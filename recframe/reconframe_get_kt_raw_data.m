function reconframe_get_kt_raw_data( rawDataDir, seriesNo, varargin )
%RECONFRAME_GET_KT_RAW_DATA  get k-t SENSE raw data using ReconFrame
%
%   Acquired data is 2D Cartesian with non-interleaved training data.
%
%   RECONFRAME_GET_KT_RAW_DATA( rawDataDir, seriesNo );
% 
%   RECONFRAME_GET_KT_RAW_DATA( ..., 'saveDataDir', '~/path/to/save/data' );
%   RECONFRAME_GET_KT_RAW_DATA( ..., 'isSaveKtData',    true );
%   RECONFRAME_GET_KT_RAW_DATA( ..., 'isSaveParamData', true );
%   RECONFRAME_GET_KT_RAW_DATA( ..., 'isSaveCsmData',   true );
%
%   RECONFRAME_GET_KT_RAW_DATA( ..., 'verbose', true );

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Setup on Beastie

%{

  run on beastie01 with pnraw01 mounted,

  e.g., 

      jfpvambp13:scp reconframe_get_kt_raw_data.m jva13@beastie01:/home/jva13/Matlab/

      jfpvambp13:ssh jva13@beastie01

      jva13@beastie01:~$ sshfs jva13@10.0.1.150:/export/pnraw/raw-ingenia ~/mnt/pnraw01-ingenia
      jva13@beastie01:~$ sshfs jva13@10.0.1.150:/export/archive-rawdata ~/mnt/pnraw01-archive

      jva13@beastie01:~$ cd ~/Matlab/
      jva13@beastie01:~$ matlab -nosplash -nodisplay -nojvm -singleCompThread
          >> reconframe_get_kt_raw_data( '/home/jva13/mnt/pnraw01-ingenia/2016_09_15/EL_100630', 25, 'saveDataDir', '/home/jva13/Matlab/data/fcmr119s25', 'verbose', true );
          >> exit
      jva13@beastie01:~$ fusermount -u ~/mnt/pnraw01-ingenia/
      jva13@beastie01:~$ fusermount -u ~/mnt/pnraw01-archive/

      jva13@beastie01:~$ exit

      jfpvambp13:realtime joshuavanamerom$ scp jva13@beastie01:/home/jva13/Matlab/data/* .

%}


%% Parse Input

p = inputParser;

default.saveDataDir     = pwd;
default.isSaveKtData    = true;
default.isSaveParamData = true;
default.isSaveCsmData   = true;
default.isVerbose       = false;

addRequired(  p, 'rawDataDir', ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
addRequired(  p, 'seriesNo', ...
    @(x) validateattributes( x, {'numeric'}, {'scalar'}, mfilename) );

addParameter( p, 'saveDataDir', default.saveDataDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );
addParameter( p, 'isSaveKtData', default.isSaveKtData, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
addParameter( p, 'isSaveParamData', default.isSaveParamData, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
addParameter( p, 'isSaveCsmData', default.isSaveCsmData, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
addParameter( p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, rawDataDir, seriesNo, varargin{:} );

saveDataDir         = p.Results.saveDataDir;
isSaveKtData        = p.Results.isSaveKtData;
isSaveParamData     = p.Results.isSaveParamData;
isSaveCsmData       = p.Results.isSaveCsmData;
isVerbose           = p.Results.verbose;


%% Add Dependencies

origPath = path;
resetPath = onCleanup( @() path( origPath ) );

addpath( '~/Matlab/MRecon-3.0.522' ), 


%% Initalise

if ( isVerbose ),
    fprintf( '\n%s()\n\n', mfilename );
end


%% Identify Raw Data Files

[ rawDataFilePath, coilSurveyFilePath, senseRefFilePath ] = id_pnraw_data( rawDataDir, seriesNo );


%% Start timer

tic,


%% Setup Save Data

if ( isVerbose ),
    fprintf( 'Setup... ' );
end

% File Names

[ ~, rawDataFileName ] = fileparts( rawDataFilePath );

get_save_file_str = @( descStr ) fullfile( saveDataDir, sprintf( '%s_%s.mat',  rawDataFileName, descStr ) );

% Files Used

filesUsed.rawData    = rawDataFilePath;
filesUsed.senseRef   = senseRefFilePath;
filesUsed.coilSurvey = coilSurveyFilePath;

% Save Directory

if ( ~exist( saveDataDir, 'dir' ) ),
    mkdir( saveDataDir ),
end


%% Anonymous Functions

% Swap dimensions 3 and 5 to match MRecon, i.e., x-y-z-c-t <-> x-y-t-c-z
swap_dim = @( data ) permute( data, [1,2,5,4,3] ); 


%% Load and Preprocess Undersampled Data

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
	fprintf( 'Loading undersampled data... ' )
end

ACQ = MRecon( rawDataFilePath );

ACQ.Parameter.Parameter2Read.typ = 1;
ACQ.Parameter.Parameter2Read.mix = 0;
ACQ.Parameter.Parameter2Read.Update;

ACQ.ReadData;
ACQ.PDACorrection;
ACQ.DcOffsetCorrection;
ACQ.RandomPhaseCorrection;
ACQ.MeasPhaseCorrection;
ACQ.SortData;


%% Derive K-Space Sampling Pattern

SMP = ACQ.Copy;
SMP.Data = single( ACQ.Data ~= 0 );


%% Load and Preprocess Training Data

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
	fprintf( 'Loading training data... ' )
end

TRN = MRecon( rawDataFilePath );

TRN.Parameter.Parameter2Read.typ = 1;
TRN.Parameter.Parameter2Read.mix = 1;
TRN.Parameter.Parameter2Read.Update;

TRN.ReadData;
TRN.PDACorrection;
TRN.DcOffsetCorrection;
TRN.RandomPhaseCorrection;
TRN.MeasPhaseCorrection;
TRN.SortData;


%% Load Sense Reference and Coil Survey Data

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
	fprintf( 'Loading sense reference and coil survey data... ' )
end


% Sense Reference

SREF = MRecon( senseRefFilePath );

SREF.Parameter.Parameter2Read.typ = 1;
SREF.Parameter.Parameter2Read.Update;

SREF.ReadData;
SREF.PDACorrection;
SREF.DcOffsetCorrection;
SREF.RandomPhaseCorrection;
SREF.MeasPhaseCorrection;
SREF.SortData;


% Coil Survey

COIL = MRecon( coilSurveyFilePath );

COIL.Parameter.Parameter2Read.typ = 1;
COIL.Parameter.Parameter2Read.Update;

COIL.ReadData;
COIL.PDACorrection;
COIL.DcOffsetCorrection;
COIL.RandomPhaseCorrection;
COIL.MeasPhaseCorrection;
COIL.SortData;


%% Load and Preprocess Noise Data

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
	fprintf( 'Loading noise data... ' )
end

NOISE = MRecon( rawDataFilePath );

NOISE.Parameter.Parameter2Read.typ = 5;
NOISE.Parameter.Parameter2Read.mix = 0;
NOISE.Parameter.Parameter2Read.Update;

NOISE.ReadData;
NOISE.PDACorrection;
NOISE.DcOffsetCorrection;
NOISE.RandomPhaseCorrection;
NOISE.MeasPhaseCorrection;
NOISE.SortData;


%% Save k-t Data

if ( isSaveKtData )

    ktdataMatFilePath = get_save_file_str( 'ktdata' );

    if ( isVerbose ),
        fprintf( '    ... %g s\n', round(toc) )
        fprintf( 'Saving data to %s ... ', ktdataMatFilePath )
    end
       
    ktAcq   = swap_dim( ACQ.Data );
    ktTrn   = swap_dim( TRN.Data );
    ktNoise = swap_dim( NOISE.Data );
        
    save( ktdataMatFilePath, 'filesUsed', 'ktAcq', 'ktTrn', 'ktNoise', '-v7.3' );
    
    clear ktAcq ktTrn ktNoise

end


%% Save Parameters


if ( isSaveParamData )

    paramMatFilePath = get_save_file_str( 'param' );

    if ( isVerbose ),
        fprintf( '    ... %g s\n', round(toc) )
        fprintf( 'Saving data to %s ... ', paramMatFilePath )
    end
    
    Param = struct( 'acq', get_recon_parameters( ACQ ), 'trn', get_recon_parameters( TRN ) );
       
    I = Param.acq.ImageInformation(1,1,2);

    dt = I.DynamicScanTime;
    if dt > 0.1,
        tr = 2 * I.EchoTime / 1000;
        nky = I.Resolution(2) / Param.acq.Scan.KtFactor;
        if ( isVerbose ),
            fprintf( '\n    Dynamic scan time (%g s) larger than expected; caluclating as num. ky lines x TR = %g.\n    ', dt, nky * tr )
        end
        dt = nky * tr;
    end
    pixdimAcq = Param.acq.Scan.AcqVoxelSize(1:2);
    
    save( paramMatFilePath, 'filesUsed', 'Param', 'dt', 'pixdimAcq', '-v7.3' );
    
    clear Param dt tr nky pixdimAcq

end


%% Coil Sensitivity Maps

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
    fprintf( 'Calculating coil sensitivity maps... ' )
end


% Create target

TGT = ACQ.Copy;
TGT.K2I;
TGT.RemoveOversampling;    % NOTE: MRsense doesn't handle oversampling correctly; oversampling removed in TGT to get correct matched dimensions.


% Create SENSE object

SENS = MRsense( SREF, TGT, COIL );     
    % NOTE: MRsense doesn't handle oversampling correctly; 
    %       oversampling removed in ACQ to get correct matched dimensions.

SENS.CalculateSensitivity = 1;
SENS.Mask = 1;
SENS.Smooth = 1;
SENS.Extrapolate = 1;
SENS.MatchTargetSize = 1;
SENS.RemoveMOversampling = 0;

SENS.Perform;


%% Save CSM Data

if ( isSaveCsmData )
      
    csm     = swap_dim( SENS.Sensitivity );             % coil sensitivity maps
    imBody  = swap_dim( SENS.ReformatedBodycoilData );  % body coil image
    imCoil  = swap_dim( SENS.ReformatedCoilData );      % array coil images
    psiSens = SENS.Psi;                                 % array coil noise covariance

    csmdataMatFilePath = get_save_file_str( 'csmdata' );
    
    if ( isVerbose ),
        fprintf( '    ... %g s\n', round(toc) )
        fprintf( 'Saving data to %s ... ', csmdataMatFilePath )
    end
       
    save( csmdataMatFilePath, 'filesUsed', 'csm', 'imBody', 'imCoil', 'psiSens', '-v7.3' );
    
    clear csm imBody imCoil psiSens
    
end


%% Add Sensitivity Maps to MRecon Objects

ACQ.Parameter.Recon.Sensitivities = SENS;
TRN.Parameter.Recon.Sensitivities = SENS;


%% Sliding Window Recon

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
    fprintf( 'Sliding window reconstruction... ' )
end

SLW = ACQ.Copy;
SLW.Data = swap_dim( kt_sliding_window( swap_dim( SLW.Data ) ) );
SLW.K2I;
SLW.RemoveOversampling;
SLW.SENSEUnfold;

xtSlw = swap_dim( SLW.Data );

xtslwMatFilePath = get_save_file_str( 'xtslw' );

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
    fprintf( 'Saving data to %s ... ', xtslwMatFilePath )
end

save( xtslwMatFilePath, 'filesUsed', 'xtSlw', '-v7.3' );

clear xtSlw


%% Cleanup

if ( isVerbose ),
    fprintf( '    ... %g s\n', round(toc) )
    fprintf( 'Done.\n\n' )
end


end  % reconframe_ktsense(...)


function S = get_recon_parameters( R )

    paramGroup = { 'Scan', 'Encoding' };

    for iP = 1:numel(paramGroup),
        fieldName = fieldnames( R.Parameter.(paramGroup{iP}));
        for iF = 1:numel(fieldName),
            S.(paramGroup{iP}).(fieldName{iF}) = R.Parameter.(paramGroup{iP}).(fieldName{iF});
        end
    end

    for iI1 = 1:size(R.Parameter.ImageInformation,1),
        for iI2 = 1:size(R.Parameter.ImageInformation,2),
            for iI3 = 1:size(R.Parameter.ImageInformation,3),
                fieldName = fieldnames( R.Parameter.ImageInformation(iI1,iI2,iI3) );
                for iF = 1:numel(fieldName),
                    S.ImageInformation(iI1,iI2,iI3).(fieldName{iF}) = R.Parameter.ImageInformation(iI1,iI2,iI3).(fieldName{iF});
                end
            end
        end
    end

end  % get_recon_parameters(...)


function [ k, ytReferenceDyn ] = kt_sliding_window( kAcq, span, method )
% KT_SLIDING_WINDOW
%
%   k = kt_sliding_window( kAcq, span, method )
%
%   k space arrays are 4-d: freq, phase, dynamic, channel,  
% 
%   span measured in frames/dynamics
% 
%   method default is nearest
% 

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Initialise

if ~exist( 'span', 'var' )
    span = size(kAcq,3);
end

if ~exist( 'method', 'var' )
    method = 'nearest';
end


%% Setup

[ nRow, nCol, nDyn, nChan ] = size( kAcq );

k = cast( zeros( size( kAcq) ), 'like', kAcq );


%% Sampling Pattern

ytSamplingPattern = squeeze( sum( sum( kAcq, 4 ), 1 ) ~= 0 );
ytReferenceDyn    = nan( size( ytSamplingPattern ) );


%% Sliding Window

for iDyn = 1:nDyn
  
  absDistToAcqLines = abs( bsxfun( @minus, ytSamplingPattern, (1:nDyn)-iDyn+1 ) );
  absDistToAcqLines( absDistToAcqLines > span | ~ytSamplingPattern ) = NaN;
  
  switch method
      
      case 'nearest'

          [ ~, iDynOfClosestAcqLine ] = min( absDistToAcqLines, [], 2 );

          for iCol = 1:nCol

            k(:,iCol,iDyn,:) = kAcq(:,iCol,iDynOfClosestAcqLine(iCol),:); 

          end

          ytReferenceDyn( :, iDyn ) = iDynOfClosestAcqLine;
  
  end
  
end


end  % kt_sliding_window(...)


