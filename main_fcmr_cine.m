function main_fcmr_cine( fcmrNo, seriesNo, safetyMargin, alpha  )
%MAIN_FCMR_CINE
%
%   MAIN_FCMR_CINE( fcmrNo, seriesNo )

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Dependencies


script_add_dependencies


%% Setup


idStr = sprintf( 'fcmr%03is%02i', fcmrNo, seriesNo );

% Number of Real-Time Frames to Reconstruct

numFrame = 96;

% k-t SENSE Recon Params 

if ~exist( 'safetyMargin', 'var' )
    safetyMargin = 27;
end

if ~exist( 'alpha', 'var' )
    alpha = 100;
end


%% Directories 


mfileDir = fileparts( mfilename( 'fullpath' ) );
resultsDir = fullfile( mfileDir, '..', 'results', idStr );


% Realtime Recon

srcDir = fullfile( mfileDir, '..', 'data', idStr );

descStr = sprintf('sm%03ia%03if%03i',safetyMargin,alpha,numFrame);
rltDirName = sprintf( '%s_%s', 'rlt', descStr );
rltDir  = fullfile( resultsDir, strcat( rltDirName, filesep ) );


% Cine Recon

if safetyMargin == 27 && alpha == 100,
  cineDir = fullfile( resultsDir, 'cine' );
else
  cineDir = fullfile( resultsDir, strcat( 'cine_', descStr ) );
end


%% Real-Time Recon


% Identify Data Sources

[ ~, fileStr ]  = fileparts( ls( fullfile( srcDir, '*_ktdata.mat' ) ) );
filePrefixStr   = fileStr(1:(end-7));

ktdataMatFilePath   = fullfile( srcDir, sprintf( '%s_ktdata.mat', filePrefixStr ) );
csmMatFilePath      = fullfile( srcDir, sprintf( '%s_csmdata.mat', filePrefixStr ) ); 
paramMatFilePath    = fullfile( srcDir, sprintf( '%s_param.mat', filePrefixStr ) );
xtslwMatFilePath    = fullfile( srcDir, sprintf( '%s_xtslw.mat', filePrefixStr ) );

maskMatFilePath     = fullfile( srcDir, 'maskfetalheart.mat' );


% Load Sliding Window Data

X = matfile( xtslwMatFilePath );
xtSlw = X.xtSlw; 


% Load K-Space Data

K = matfile( ktdataMatFilePath );
ktAcq = K.ktAcq;
ktTrn = K.ktTrn;


% Load CSM Data

C = matfile( csmMatFilePath );
csm   = C.csm;
psi   = C.psiSens;


% Load Parameters Data

P = matfile( paramMatFilePath );
Param = P.Param;
dtRlt = P.dt; 
pixdimAcq = P.pixdimAcq;



%% Close matfile Objects

clear X K C P


%% Get / Create Mask

if ( ~exist( maskMatFilePath, 'file' ) ),
    mask = uiget_roi_mask( abs( mean( xtSlw, 3 ) ), 'fetal heart' );
    save( maskMatFilePath, 'mask', '-v7.3' );
else
    M = matfile( maskMatFilePath );
    mask = M.mask;
    clear M
end


%% Drop Non-Steady State Frames

if isempty( numFrame ),
    numFrame = size(ktAcq,3);  % e.g., s = squeeze(sum(sum(bsxfun(@times,mask,abs(xtSlw)))));  figure,  plot(s),  grid on,  
end                            % or,   frameNo = find_steadystate_framenos( permute(abs(xtSlw),[1,2,4,3]), 'showFig' ); 
                               %       ktFactor = 8; 
                               %       numFrame = floor(length(frameNo)/ktFactor)*ktFactor;

ktAcq = ktAcq(:,:,(end-numFrame+1):end,:);
ktTrn = ktTrn(:,:,1:numFrame,:);
xtSlw = xtSlw(:,:,(end-numFrame+1):end);


%% Prep x-t Data

[ xtAcq, xtTrn, xtSmp, psiEst, xtBln, xtDff ] = prep_ktsense( ktAcq, ktTrn, 'removeoversampling', true );


%% Recon Real-Time

% Create output directory

mkdir( rltDir )

% Log 

diary( fullfile( rltDir, 'log.md') );
fprintf( 'recon_xtsense\n' )
fprintf( '=============\n' )
fprintf( '\n\nstart: %s\n\n\n', datestr(now) )

% Recon


[ xtRlt, safetyMargin, xfFlt, unmix, xfRlt, xfBln, xfAcq, xfTrn, xfMask, psi, csm, mask, dtRlt, alpha ] = recon_xtsense( xtDff, xtSmp, xtTrn, csm, psi, 'xtBln', xtBln, 'safetyMargin', safetyMargin, 'mask', mask, 'alpha', alpha, 'dt', dtRlt, 'makeAdaptiveFilter', true, 'verbose', true );

% Save figures
save_figs( rltDir )
close all,

% Save .mat

save( fullfile( rltDir, 'xtsense_recon' ), 'xtRlt', 'safetyMargin', 'alpha', 'xfFlt', 'xfTrn', 'mask', 'csm', 'unmix', 'psi', 'pixdimAcq', 'dtRlt', '-v7.3' )

% Save .nii

% TODO: change save_xt2nii to take input argument pixdim = [dx,dy,dz,dt]
% instead of Param and remove reading of Param from .mat file to reduce
% processing time
save_xt2nii( xtRlt, Param, fullfile( rltDir, 'realtime.nii.gz' ) );
delete( fullfile( rltDir, 'realtime.nii.gz' ) );
save_xt2nii( sum(bsxfun(@times,unmix,xtTrn),4), Param, fullfile( rltDir, 'training.nii.gz' ) );
delete( fullfile( rltDir, 'training.nii.gz' ) );
save_xt2nii( mask, Param, fullfile( rltDir, 'mask_fetalheart.nii.gz' ) );
delete( fullfile( rltDir, 'mask_fetalheart_xyt.nii.gz' ) );

% Stop Logging

fprintf( '\n\nend: %s\n\n\n', datestr(now) )
diary( 'off' );    



%% rlt2cine


[ imCine, tCine ] = rlt2cine_linear( xtRlt, dtRlt, 'mask', mask, 'pixdimAcq', pixdimAcq, 'isSaveResults', true, 'outputDir', cineDir, 'verbose', true );

if safetyMargin == 27 && alpha == 100,
  [ imCineNoMoco, tCineNoMoco ] = rlt2cine_linear( xtRlt, dtRlt, 'mask', mask, 'pixdimAcq', pixdimAcq, 'isSaveResults', true, 'doMoco', false, 'outputDir', strcat( cineDir, '_noMoco' ), 'doConvergenceTest', true, 'verbose', true );
  [ imCineNoOutrej, tCineNoOutrej ] = rlt2cine_linear( xtRlt, dtRlt, 'mask', mask, 'pixdimAcq', pixdimAcq, 'isSaveResults', true, 'doOutrej', false, 'outputDir', strcat( cineDir, '_noOutrej' ), 'doConvergenceTest', true, 'verbose', true );
end


%% create animated gifs

make_animated_gifs( srcDir, resultsDir );


end  % main_fcmr_cine(...)