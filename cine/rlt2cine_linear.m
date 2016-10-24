function [ imCine, tCine ] = rlt2cine_linear( imRlt, tRlt, varargin )
%RLT2CINE_LINEAR  reconstruct cine image sequence from real-time image sequence
%
%   [ imCine, tCine ] = RLT2CINE_LINEAR( imRlt, tRlt )
%
%       imRlt           array of complex real-time images           [x,y,t]
%       tRlt            vector of image times in seconds
%
%   Additional name-value input: 
%
%       mask            logical image mask; default all voxels included
%       tRTrigger       cardiac R-R trigger times; default estimated from
%                       imRlt
%       imRef           default imRlt; TODO: add description 
%       pixdimAcq       default [ 2.00, 2.00 ] mm; TODO: add description 
%       pixdimRcn       default [ 1.25, 1.25 ] mm; TODO: 
%       nFrameCine      default 25; TODO: add description
%       doMoco          default true; TODO: add description
%       doOutrej        default true; TODO: add description
%       isAdult         default false; TODO: add description
%       isSaveResults   default false; TODO: add description
%       outputDir       default current directory; TODO: add description
%       verbose         verbosity; default false

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Notes / WIP


% TODO: maskMoco as optional input

% TODO: doKspaceApodisation as optional innput, if used during resizing images
doKspaceApodisation = false;  


%% Parse Input


default.mask          = true(size(imRlt(:,:,1))); 
default.tRTrigger     = [];
default.imRef         = imRlt;
default.pixdimAcq     = [ 2.00, 2.00 ];
default.pixdimRcn     = [ 1.25, 1.25 ]; 
default.nFrameCine    = 25;
default.doMoco        = true;
default.doOutrej      = true;
default.doConvergenceTest = true;
default.isAdult       = false;
default.isSaveResults = false;
default.outputDir     = pwd;
default.isVerbose     = false;

p = inputParser;

addRequired(  p, 'imRlt', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );

addRequired(  p, 'tRlt', @(x) validateattributes( x, {'numeric'}, ...
        {'nonempty'}, mfilename ) );

addParameter( p, 'mask', default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(imRlt(:,:,1)),NaN]}, mfilename ) );

addParameter( p, 'tRTrigger', default.tRTrigger, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'increasing'}, mfilename ) );

addParameter( p, 'imRef', default.imRef, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );

addParameter( p, 'pixdimAcq', default.pixdimAcq, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'numel',2,'positive'}, mfilename ) );

addParameter( p, 'pixdimRcn', default.pixdimRcn, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'numel',2,'positive'}, mfilename ) );    

addParameter( p, 'nFrameCine', default.nFrameCine, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','positive'}, mfilename ) );

addParameter( p, 'doMoco', default.doMoco, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

addParameter( p, 'doOutrej', default.doOutrej, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

addParameter( p, 'doConvergenceTest', default.doConvergenceTest, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );
    
addParameter( p, 'isAdult',     default.isAdult, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );
    
addParameter( p, 'isSaveResults',     default.isSaveResults, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );
    
addParameter( p, 'outputDir',     default.outputDir, ...
        @(x) validateattributes( x, {'char'}, ...
        {'nonempty'}, mfilename ) );
    
addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

parse( p, imRlt, tRlt, varargin{:} );

mask        = p.Results.mask;
tRTrigger   = p.Results.tRTrigger;
imRef       = p.Results.imRef;
pixdimAcq   = p.Results.pixdimAcq;
pixdimRcn   = p.Results.pixdimRcn;
nFrameCine  = p.Results.nFrameCine;
doMoco      = p.Results.doMoco;
doOutrej    = p.Results.doOutrej;
doConvergenceTest = p.Results.doConvergenceTest;
isAdult     = p.Results.isAdult;
isSaveResults = p.Results.isSaveResults;
outputDir   = p.Results.outputDir;
isVerbose   = p.Results.verbose;


%% Setup


% Dependencies

% origPath  = path;
% resetPath = onCleanup( @() path( origPath ) );
% addpath( '~/Research/fcmr_cine/' )
% addpath( '~/Research/fcmr_cine/moco' )
% addpath( '~/Research/fcmr_cine/aux/nifti/' )
% addpath( genpath( '~/Research/fcmr_cine/util' ) )


% Utility (anonymous) Functions

get_mask_values = @( x, msk ) double( x( repmat( msk, [1,1,size(x,3)/size(mask,3)] ) ) );
    % param2tform     = @( tx, ty, rz ) affine2d( [cosd(rz), sind(rz), 0; -sind(rz), cosd(rz), 0; tx, ty, 1 ] );
calc_rmsd       = @( a, b ) sqrt( mean( ( a(:) - b(:) ) .^2 ) );
calc_rmsd_cine  = @( c2, c1, msk ) calc_rmsd( get_mask_values( abs(c2), mask ), get_mask_values( abs(c1), mask ) );
    % time2cphase     = @( t, rr ) 2*pi*t./rr;
    % calc_rmsd_card  = @( a, b, rr1, rr2 ) rr1 / ( 2*pi ) * sqrt( mean( angle( exp( -sqrt(-1) * ( time2cphase(a(:),rr1) - time2cphase(b(:),rr2) ) ) ) ).^2 );


% Close Open Figures

if ( isSaveResults ),
    close all,
end


%% Init


% Iterations

maxIter = 5;  % max. num. iterations


% Tolerances

    % TOL.tRrRMSD       = 0.001*60/140;  % ~0.1 pct. of R-R interval, seconds   %0.0003;                 % seconds (~0.1bpm)
    % TOL.dispRMSD      = 0.001*mean(pixdimAcq); %0.0005*mean(pixdimAcq); % 0.1 pct. of voxel size, millimetres
    
TOL.imCineRMSD    = 0.005*mean(get_mask_values(abs(imRlt),mask));  % 0.01 pct. of max. signal intensity, aritrary units


% Tests 

    % isSameFrmTime = @( tRr1, tRr2, rr1, rr2 ) ( calc_rmsd_card ( tRr1, tRr2, rr1, rr2 ) < TOL.tRrRMSD );

    % isSameDispMap = @( dispMap1, dispMap2 ) ...
    %     ( abs( calc_rmsd( get_mask_values( dispMap1, mask ), ...
    %                       get_mask_values( dispMap2, mask ) ) ) ...
    %         < TOL.dispRMSD );

isSameCineIms = @( imCine1, imCine2 ) ( calc_rmsd_cine( imCine1, imCine2, mask ) < TOL.imCineRMSD );


% Dimensions

nFrameRlt   = size(imRlt,3);
mmPerSecond = 100;  % mm per s for saving as xyt .nii files


% Timing

if ( numel( tRlt ) == 1 ),  % frame duration given instead of frame times
    dtRlt = tRlt;
    tRlt  = dtRlt * (0:1:(nFrameRlt-1));
else
    dtRlt = mean(diff(tRlt));
end


% Motion Correction Mask

rMax_mm  = 15;                              % max. dilation of ROI in mm
rMax     = rMax_mm / mean(pixdimAcq);       % max. dilation of ROI in pixels   
aMsk     = sum( mask(:) );                  % area of mask in sq. pixels
aFac     = 2;                               % area incease factor 
rInc     = (sqrt(aFac)-1)*sqrt(aMsk/pi);    % increase in radius to inrease area by aFac
rDil     = min( rMax, rInc );               % amount to dilate mask
maskMoco = bwmorph( mask, 'dilate', rDil ); % OR: maskMoco = uiget_rect_mask( abs(mean(imRlt,3)), 'motion-correction' );


% Entropy Measurements

E0 = struct();
E  = struct();


% Output Directory

if ( isSaveResults ) && ~exist( outputDir, 'dir' )
    mkdir( outputDir )
end


%% Logging


logFilePath = fullfile( outputDir, 'log.txt');

diary( logFilePath );

closeDiary = onCleanup( @() diary( 'off' ) );

fprintf( 'rlt2cine\n' )
fprintf( '========\n' )
fprintf( '\n\nstart: %s\n\n\n', datestr(now) )


%% Save Global Data


if ( isSaveResults ),
           
    pixdimAcqRltNii  = [ pixdimAcq mmPerSecond*dtRlt  ];
    pixdimRcnRltNii  = [ pixdimRcn mmPerSecond*dtRlt  ];

    save_nii( make_nii( abs(imRlt), pixdimAcqRltNii ),    fullfile( outputDir, 'rlt_xyt.nii.gz' ) );
    save_nii( make_nii( abs(imRef), pixdimAcqRltNii ),    fullfile( outputDir, 'ref_xyt.nii.gz' ) );
    save_nii( make_nii( abs(mask), pixdimAcqRltNii ),     fullfile( outputDir, 'mask_fetalheart.nii.gz' ) );
    save_nii( make_nii( abs(maskMoco), pixdimAcqRltNii ), fullfile( outputDir, 'mask_moco.nii.gz' ) );
    
end


%% Resize Images


acqDx = pixdimAcq(1);  
acqDy = pixdimAcq(2);  

rcnDx = pixdimRcn(1);
rcnDy = pixdimRcn(2);

% imRlt

xt2kt = @( xt ) fftshift( fftshift( fft2( xt ), 2 ), 1 );
kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, 2 ), 1 ) );

acqDim = size( imRlt ); 

fov = pixdimAcq .* acqDim(1:2);

rcnDim = [ acqDx/rcnDx*size(imRlt,1), acqDy/rcnDy*size(imRlt,2), nFrameRlt ];  
padDim   = round( ( rcnDim - acqDim ) / 2 );

ktRlt    = xt2kt( imRlt );  
 
if ( doKspaceApodisation ),  % TODO: FIXME
    %{
    filterEdgeLengthY = 2 * A.ktFactor;
    filterEdgeY = cos( linspace(0,pi,filterEdgeLengthY+2) )/2 + 0.5; 
    filterEdgeY = filterEdgeY( 2:(end-1) );
    filterMaskY = [ flip( filterEdgeY(:) ); ones( size(imRlt,1)-2*length(filterEdgeY), 1 ); filterEdgeY(:) ];
    filterEdgeX = cos( linspace(0,pi,round(A.dx/A.dy*length(filterEdgeY))+2) )/2 + 0.5; 
    filterEdgeX = filterEdgeX( 2:(end-1) );
    filterMaskX = [ flip( filterEdgeX(:) ); ones( size(imRlt,2)-2*length(filterEdgeX), 1 ); filterEdgeX(:) ];
    filterMask  = filterMaskY * filterMaskX';
    if ( isVerbose ),
        figure,
        imgrid(filterMask),
        colormap(parula),
        title( 'k-space apodisation' ),
    end
    ktRlt = bsxfun( @times, filterMask, ktRlt );
    %}
else
    filterMask = ones(size(ktRlt(:,:,1)));
end

ktRltQ = padarray( ktRlt, padDim );

scaleFactor = numel(ktRltQ(:,:,1))/sum(filterMask(:));

imRltQ = kt2xt( scaleFactor * ktRltQ );


% masks

x = linspace(0,fov(2),size(imRlt,2)+1); x = x(2:end)-mean(diff(x))/2;
y = linspace(0,fov(1),size(imRlt,1)+1); y = y(2:end)-mean(diff(y))/2;
xq = linspace(0,fov(2),size(imRltQ,2)+1); xq = xq(2:end)-mean(diff(xq))/2;
yq = linspace(0,fov(1),size(imRltQ,1)+1); yq = yq(2:end)-mean(diff(yq))/2;
[X,Y] = meshgrid(x,y);
[Xq,Yq] = meshgrid(xq,yq); 
maskQ = interp2(X,Y,single(mask),Xq,Yq,'linear') > 0.5;
maskMocoQ = interp2(X,Y,single(maskMoco),Xq,Yq,'linear') > 0.5;


% Save Results

if ( isSaveResults ),

    saveDir = fullfile( outputDir );
    
    save_nii( make_nii( abs(imRltQ), pixdimRcnRltNii ),    fullfile( saveDir, 'rltq_xyt.nii.gz' ) );
    save_nii( make_nii( abs(maskQ), pixdimRcnRltNii ),     fullfile( saveDir, 'maskq_fetalheart.nii.gz' ) );
	save_nii( make_nii( abs(maskMocoQ), pixdimRcnRltNii ), fullfile( saveDir, 'maskq_moco.nii.gz' ) );
    save_nii( make_nii( abs(maskQ), pixdimRcnRltNii ),     fullfile( saveDir, 'maskq_fetalheart.nii.gz' ) );
	save_nii( make_nii( abs(maskMocoQ), pixdimRcnRltNii ), fullfile( saveDir, 'maskq_moco.nii.gz' ) );

end


%% Initial Parameters


% Cardiac Synchronisation


if isempty( tRTrigger ),
    
    doCardsync = true;
    
    if ( isAdult ),
        [ P0.R.rrInterval, P0.R.tRTrigger ] = estimate_heartrate_xf( imRlt, dtRlt, 'roi', mask, 'hrRange', [45 145], 'useHarmonic', false );
    else
        [ P0.R.rrInterval, P0.R.tRTrigger ] = estimate_heartrate_xf( imRlt, dtRlt, 'roi', mask, 'useHarmonic', false );
    end
    
    P0.R.tRr = calc_cardiac_timing( tRlt, P0.R.tRTrigger );  


else
    
    P0.R.tRTrigger = tRTrigger;
    
    doCardsync = false;
    
    switch numel( tRTrigger ),        
        
        case 1,     % Calculate R-Wave Trigger Times
            
            P0.R.rrInterval = tRTrigger;
            
            nTrigger        = ceil( nFrameRlt * dtRlt / P0.R.rrInterval );
            P0.R.tRTrigger  = P0.R.rrInterval * (0:nTrigger);
            
        otherwise,  % Calculate Representative R-R Interval
            
            P0.R.tRTrigger   = tRTrigger;
            
    end
    
    % Calculate Time Since Last R-Wave Trigger and R-R Interval

    [ P0.R.tRr, ~, P0.R.rrInterval ] = calc_cardiac_timing( tRlt, P0.R.tRTrigger );   

    P0.R.rrInterval = mean( P0.R.rrInterval );
    
end


% Motion Correction

P0.T.tform      = init_tform_struct( nFrameRlt );
P0.T.dispMap    = zeros( size( imRlt ) );


% Outlier Rejection

P0.P.vox        = ones( size( imRlt ) );
P0.P.frm        = ones( nFrameRlt, 1 );
[ P0.P.imCine, P0.P.tCine ] = recon_cine( imRlt, dtRlt, P0.R.tRr, P0.R.rrInterval, P0.T.tform, maskMoco, P0.P.vox, P0.P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );


% Adjust R-Wave Triggers to Align Cardiac Cycle

if ( doCardsync ),
    
    % TODO: use different calculation of tRTrigger than taking output of
    % estimate_hr_xf and shifting to align cardiac cycle to save on
    % computations
    
    tOffset = calc_cine_timing_offset( P0.P.imCine, P0.P.tCine, mask );
    
    if ( tOffset ~= 0 ), 

        P0.R.tRTrigger = P0.R.tRTrigger - tOffset;

        while min( P0.R.tRTrigger ) > min( tRlt )
            P0.R.tRTrigger = [ min( P0.R.tRTrigger ) - P0.R.rrInterval, P0.R.tRTrigger ];
        end

        while max( P0.R.tRTrigger ) < max( tRlt )
            P0.R.tRTrigger = [ P0.R.tRTrigger, max( P0.R.tRTrigger ) + P0.R.rrInterval ];
        end

        P0.R.tRr = calc_cardiac_timing( tRlt, P0.R.tRTrigger );  

        [ P0.P.imCine, P0.P.tCine ] = recon_cine( imRlt, dtRlt, P0.R.tRr, P0.R.rrInterval, P0.T.tform, maskMoco, P0.P.vox, P0.P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );

    end

end


% Entropy

bwCrop = mask; row = find(sum(bwCrop,2)>0); col = find(sum(bwCrop,1)>0);
% E0.time  = imagemetric( P0.P.imCine(row,col,:), {'PC'} );
% E0.image = imagemetric( P0.P.imCine(row,col,:), {'Cine'} );


%% Iterative Processing   


for iIter = 1:maxIter,

    
    if ( iIter == 1 ),
        P(iIter) = P0;
    else
        P(iIter) = P(iIter-1);
    end
    
    
%% Cardiac Synchronisation


if ( doCardsync ),
           
    
    if ( isVerbose ),
        fprintf( '\n%02i a. Cardiac Synchronisation\n', iIter )
        fprintf( '-----------------------------\n\n' )
    end
    

% Apply Transformation to Real-Time Images

imRltT = transform_imseq( imRlt, P(iIter).T.tform, maskMoco, pixdimAcq );


% Estimate Cardiac Triggers

if isAdult,
    [ P(iIter).R.rrInterval, P(iIter).R.tRTrigger ] = estimate_heartrate_xf( imRltT, dtRlt, 'roi', mask, 'hrRange', [45 145], 'useHarmonic', false, 'verbose', isVerbose );
else
    [ P(iIter).R.rrInterval, P(iIter).R.tRTrigger ] = estimate_heartrate_xf( imRltT, dtRlt, 'roi', mask, 'useHarmonic', false, 'verbose', isVerbose );
end

if( isVerbose ),
    fprintf( 'Estimated R-R interval = %.2f ms  \n\n', P(iIter).R.rrInterval*1000 )
end


% Calculate Time Since Last R-Wave Trigger

P(iIter).R.tRr = calc_cardiac_timing( tRlt, P(iIter).R.tRTrigger );  


% Recon Cine

[ P(iIter).R.imCine, P(iIter).R.tCine ] = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );


% Adjust R-Wave Triggers to Align Cardiac Cycle

tOffset = calc_cine_timing_offset( P(iIter).R.imCine, P(iIter).R.tCine, mask );

if ( tOffset ~= 0 ), 
    
    P(iIter).R.tRTrigger = P(iIter).R.tRTrigger - tOffset;
    
    while min( P(iIter).R.tRTrigger ) > min( tRlt )
        P(iIter).R.tRTrigger = [ min( P(iIter).R.tRTrigger ) - P(iIter).R.rrInterval, P(iIter).R.tRTrigger ];
    end
    
    while max( P(iIter).R.tRTrigger ) < max( tRlt )
        P(iIter).R.tRTrigger = [ P(iIter).R.tRTrigger, max( P(iIter).R.tRTrigger ) + P(iIter).R.rrInterval ];
    end
    
    P(iIter).R.tRr = calc_cardiac_timing( tRlt, P(iIter).R.tRTrigger );
    
    [ P(iIter).R.imCine, P(iIter).R.tCine ] = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );

end


% Measure Entropy

% E(iIter).R.time  = imagemetric( P(iIter).R.imCine(row,col,:), {'PC'} );
% E(iIter).R.image = imagemetric( P(iIter).R.imCine(row,col,:), {'Cine'} );


% Save Results

if ( isSaveResults ),

    cardSaveDir  = sprintf( '%02ia_cardsync', iIter );
    cardSavePath = fullfile( outputDir, cardSaveDir );

    save_figs( cardSavePath ),

    close all,

    if ( isVerbose ),
        fprintf( '![](%s)  \n\n', fullfile( cardSaveDir, 'figs', 'heart_rate_estimate.png' ) )
    end
    
    recon_cine( imRltQ, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMocoQ, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimRcnRltNii, fov, outputDir );
    movefile( fullfile( outputDir, 'cineq_xyt.nii.gz' ), fullfile( outputDir, sprintf( 'cineq_xyt_%s.nii.gz', cardSaveDir ) ) )
    delete( fullfile( outputDir, 'voxprobq_xyt.nii.gz' ) ),

end


if ( ~doMoco )  
    
    doCardsync = false;  % don't need do cardiac synchronisation in subsequent iterations if input conditions don't change

    if ( ~doOutrej ),
        break,  % don't need do any more iterations if not moco or outrej being performed
    end
    
end


else

    % Recon Cine, Measure Entropy

    [ P(iIter).R.imCine, P(iIter).R.tCine ] = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
    % E(iIter).R.time  = imagemetric( P(iIter).R.imCine(row,col,:), {'PC'} );
    % E(iIter).R.image = imagemetric( P(iIter).R.imCine(row,col,:), {'Cine'} );
    

end  % if ( doCardsync ),


%% Motion-Correction


if ( doMoco ),

    
if ( isVerbose ),
    fprintf( '\n%02i b. Motion Correction\n', iIter )
    fprintf( '-----------------------\n\n' )
end


[ P(iIter).T.tform, P(iIter).T.dispMap, ~, ~, imTgt, imMvg, imFix ] = motion_correction( imRef, P(iIter).R.tRr, dtRlt, P(iIter).R.rrInterval, maskMoco, pixdimAcq, 'transformations', P(iIter).T.tform, 'voxelprob',P(iIter).P.vox, 'frameprob',P(iIter).P.frm, 'verbose', isVerbose );


% Apply Transformations to Real-Time Image Sequence

imRltT = transform_imseq( imRlt, P(iIter).T.tform, maskMoco, pixdimAcq );


% Recon Cine, Measure Entropy

P(iIter).T.imCine = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
% E(iIter).T.time  = imagemetric( P(iIter).T.imCine(row,col,:), {'PC'} );
% E(iIter).T.image = imagemetric( P(iIter).T.imCine(row,col,:), {'Cine'} );


% Save Results

if ( isSaveResults ),
    
    mocoSaveDir  = sprintf( '%02ib_moco', iIter );
    mocoSavePath = fullfile( outputDir, mocoSaveDir );
    
    save_figs( mocoSavePath ),
    
    close all,
    
    if ( isVerbose ),
        fprintf( '![](%s)  \n\n', fullfile( mocoSaveDir, 'figs', 'motion_correction.png' ) )
    end
        
    save_nii( make_nii( abs(imRltT), pixdimAcqRltNii ), fullfile( mocoSavePath, 'rlt_moco_xyt.nii.gz' ) );
    save_nii( make_nii( abs(imTgt),  pixdimAcqRltNii ), fullfile( mocoSavePath, 'tgt_xyt.nii.gz' ) );
    save_nii( make_nii( abs(imMvg),  pixdimAcqRltNii ), fullfile( mocoSavePath, 'mvg_xyt.nii.gz' ) );
    save_nii( make_nii( abs(imFix),  pixdimAcqRltNii ), fullfile( mocoSavePath, 'fix_xyt.nii.gz' ) );
    
    recon_cine( imRltQ, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMocoQ, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimRcnRltNii, fov, outputDir );
    movefile( fullfile( outputDir, 'cineq_xyt.nii.gz' ), fullfile( outputDir, sprintf( 'cineq_xyt_%s.nii.gz', mocoSaveDir ) ) )
    delete( fullfile( outputDir, 'voxprobq_xyt.nii.gz' ) ),
    
end


else  % if ( doMoco ),
    
    
    % Recon Cine, Measure Entropy
    
    P(iIter).T.imCine = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
    % E(iIter).T.time  = imagemetric( P(iIter).T.imCine(row,col,:), {'PC'} );
    % E(iIter).T.image = imagemetric( P(iIter).T.imCine(row,col,:), {'Cine'} );

    
end  % if ( doMoco ),


%% Outlier Rejection


if ( doOutrej ),

    
fprintf( '\n%02i c. Outlier Rejection\n', iIter )
fprintf( '-----------------------\n\n' )

    
[ P(iIter).P.vox, P(iIter).P.frm ] = outlier_rejection_step( imRlt, P(iIter).R.tRr, dtRlt, P(iIter).R.rrInterval, 'mask', mask, 'transformations', P(iIter).T.tform, 'transformationMask', maskMoco, 'pixdim', pixdimAcq, 'voxProb', P(iIter).P.vox, 'frmProb', P(iIter).P.frm, 'verbose', isVerbose );


% Recon Cine, Measure Entropy

P(iIter).P.imCine = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
% E(iIter).P.time  = imagemetric( P(iIter).P.imCine(row,col,:), {'PC'} );
% E(iIter).P.image = imagemetric( P(iIter).P.imCine(row,col,:), {'Cine'} );


% Save Results

if ( isSaveResults ),

    outrejSaveDir  = sprintf( '%02ic_outrej', iIter );
    outrejSavePath = fullfile( outputDir, outrejSaveDir );
    
    save_figs( fullfile( outrejSavePath ) ),
    
    close all,
    
    if ( isVerbose ),
        fprintf( '\n![](%s/figs/voxel_error_distribution.png)  \n', outrejSaveDir ),
        fprintf( '![](%s/figs/voxel_probability.png)  \n', outrejSaveDir ),
        fprintf( '![](%s/figs/frame_probability.png)  \n\n', outrejSaveDir ),
    end
    
    save_nii( make_nii( P(iIter).P.vox, pixdimAcqRltNii ), fullfile( outrejSavePath, 'voxprob_xyt.nii.gz' ) );
    
    recon_cine( imRltQ, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMocoQ, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimRcnRltNii, fov, outputDir );
    movefile( fullfile( outputDir, 'cineq_xyt.nii.gz' ), fullfile( outputDir, sprintf( 'cineq_xyt_%s.nii.gz', outrejSaveDir ) ) )
    movefile( fullfile( outputDir, 'voxprobq_xyt.nii.gz' ), fullfile( outputDir, sprintf( 'voxprobq_xyt_%s.nii.gz', outrejSaveDir ) ) )
    
end


else  % if ( doOutrej ),
    
    
    % Recon Cine, Measure Entropy

    P(iIter).P.imCine = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
    % E(iIter).P.time  = imagemetric( P(iIter).P.imCine(row,col,:), {'PC'} );
    % E(iIter).P.image = imagemetric( P(iIter).P.imCine(row,col,:), {'Cine'} );


end  % if ( doOutrej ),


%% Test for Convergence


if ( doConvergenceTest ) 
    
    if iIter == 1,
        if ( isSameCineIms( P(iIter).P.imCine,  P0.P.imCine ) ),
              break,
        end
    else
        if ( isSameCineIms( P(iIter).P.imCine,  P(iIter-1).P.imCine ) ),
              break,
        end
    end

end


end  % for iIter = 1:maxIter,


%% Save Reordered Real-Time 


if ( isSaveResults ),
    
    [ ~, indRlt2Card ] = sort( P(end).R.tRr );
       
    dtReordered = P(end).R.rrInterval / size(imRltQ,3);
    
    pixdimRcnReorderNii  = [ pixdimRcn mmPerSecond*dtReordered ];
   
    save_nii( make_nii( abs(imRltQ(:,:,indRlt2Card)), pixdimRcnReorderNii ), fullfile( saveDir, 'reorderedq_xyt.nii.gz' ) );
    
end 


%% Summarise Convergence


nIter = iIter;

cineRMSD = nan( nIter, 1 );
    % cardRMSD = nan( nIter, 1 );
    % dispRMSD = nan( nIter, 1 );

cineRMSD(1) = calc_rmsd_cine( P(1).P.imCine, P0.P.imCine, mask );
    % cardRMSD(1) = calc_rmsd_card( P(1).R.tRr, P0.R.tRr, P(1).R.rrInterval, P0.R.rrInterval );
    % dispRMSD(1) = calc_rmsd( get_mask_values( P(1).T.dispMap, mask ), get_mask_values( P0.T.dispMap, mask ) );

for iIter = 2:nIter, 
       
    cineRMSD(iIter) = calc_rmsd_cine( P(iIter).P.imCine, P(iIter-1).P.imCine, mask );
        % cardRMSD(iIter) = calc_rmsd_card( P(iIter).R.tRr, P(iIter-1).R.tRr, P(iIter).R.rrInterval, P(iIter-1).R.rrInterval );
        % dispRMSD(iIter) = calc_rmsd( get_mask_values( P(iIter).T.dispMap, mask ), get_mask_values( P(iIter-1).T.dispMap, mask ) );

end

if ( isVerbose ),
   
    figure( 'Name', 'convergence' )
    
    %{
        subplot(3,1,1),
        plot( [0,(nIter+1)], 1000*TOL.tRrRMSD*ones(2,1), 'c-', 1:nIter, 1000*cardRMSD, 'bo-', 'LineWidth', 1.5 )
        xlabel('iteration')
        ylabel('RMSD timing (ms)')
        set(gca,'XLim',[0.5,nIter+0.5])
        grid on
    
        subplot(3,1,2),
        plot( [0,(nIter+1)], TOL.dispRMSD*ones(2,1), 'c-', 1:nIter, abs(dispRMSD), 'bo-', 'LineWidth', 1.5 )
            % plot( 1:nIter, abs(dispRMSD), 'o-', 1:nIter, real(dispRMSD), 's-.', 1:nIter, imag(dispRMSD), 's-.', [0,(nIter+1)], TOL.dispRMSD*ones(2,1), 'c:', 'LineWidth', 1.5 )
            % legend('disp.','disp. x','disp. y', 'Orientation', 'horizontal' )
        xlabel('iteration')
        ylabel('RMSD displacement (mm)')
        grid on
        set(gca,'XLim',[0.5,nIter+0.5])
        
        subplot(3,1,3),
    %}
    
    plot( [0,(nIter+1)], TOL.imCineRMSD*ones(2,1), 'c-', 1:nIter, cineRMSD, 'bo-', 'LineWidth', 1.5 )
    xlabel('iteration')
    ylabel('RMSD |X| (a.u.)')
    grid on
    set(gca,'XLim',[0.5,nIter+0.5])
        
    if ( isSaveResults ),
        
        fprintf( '\nConvergence\n' )
        fprintf( '-----------\n\n' )
        
        save_figs( outputDir ),
        
        close all,
        
        fprintf( '![](figs/convergence.png)  \n' )
        
    end
    
end


%{
NOTE: if using the code below, update imCineRMSD calculation to include abs()
imCineRMSD = nan( 3*nIter, 1 );

imCineRMSD(1) = calc_rmsd( get_mask_values( P(1).R.imCine, mask ), get_mask_values( P0.P.imCine, mask ) );

for iIter = 1:nIter, 
    
    if ( iIter > 1 )
            imCineRMSD(3*(iIter-1)+1) = calc_rmsd( get_mask_values( P(iIter).R.imCine, mask ), get_mask_values( P(iIter-1).P.imCine, mask ) );
    end
    
    imCineRMSD(3*(iIter-1)+2) = calc_rmsd( get_mask_values( P(iIter).T.imCine, mask ), get_mask_values( P(iIter).R.imCine, mask ) );
    imCineRMSD(3*(iIter-1)+3) = calc_rmsd( get_mask_values( P(iIter).P.imCine, mask ), get_mask_values( P(iIter).T.imCine, mask ) );

end

if ( isVerbose ),
   
    figure( 'Name', 'cine_convergence' )
    
    plot( [0,(nIter+1)], TOL.imCineRMSD*ones(2,1), 'c-', 1:(1/3):(nIter+2/3), abs(imCineRMSD), 'bo-', 'LineWidth', 1.5 )
    xlabel('iteration')
    ylabel('RMSD |X| (a.u.)')
    set(gca,'XLim',[1-0.5,nIter+2/3+0.5])
    grid on
        
    if ( isSaveResults ),
        
        save_figs( outputDir ),
        
        close all,
        
        fprintf( '![](figs/cine_convergence.png)  \n' )
        
    end
    
end
%}


%% Summarise Entropy


%{

if ( isVerbose ),
    
    entropyTimeVsR = E0.time;
    entropyTimeVsT = E0.time;
    entropyTimeVsP = E0.time;
    
    entropyImageVsR = E0.image;
    entropyImageVsT = E0.image;
    entropyImageVsP = E0.image;
    
    for iIter = 1:nIter,
        
        entropyTimeVsR(end+1) = E(iIter).R.time;
        entropyTimeVsT(end+1) = E(iIter).T.time; 
        entropyTimeVsP(end+1) = E(iIter).P.time; 
            
        entropyImageVsR(end+1) = E(iIter).R.image;
        entropyImageVsT(end+1) = E(iIter).T.image; 
        entropyImageVsP(end+1) = E(iIter).T.image; 
        
    end
   
    figure( 'Name', 'entropy' )
    
    yyaxis left
    
    plot( 0:nIter, entropyTimeVsP, 'LineWidth', 1 ), 
    
    yyaxis right
    
    plot( 0:nIter, entropyImageVsP, 'LineWidth', 1 ), 
    
    yyaxis left
    
    hold on
    
    p1 = plot( 0:nIter, entropyTimeVsR, 'x', 'MarkerSize', 8 );
    p2 = plot( 0:nIter, entropyTimeVsT, '+', 'MarkerSize', 8 );
    p3 = plot( 0:nIter, entropyTimeVsP, 'o', 'MarkerSize', 8 );
    
    
    xlabel('iteration')
    ylabel('Time Entropy (a.u.)')
    
    set(gca,'XLim',[-0.5,maxIter+0.5],'XTick',0:maxIter,'YTick',[]), 
    
    legend('C','T','P')
    
	yyaxis right
        
    hold on
    
    p4 = plot( 0:nIter, entropyImageVsR, 'x', 'MarkerSize', 8 );
    p5 = plot( 0:nIter, entropyImageVsT, '+', 'MarkerSize', 8 );
    p6 = plot( 0:nIter, entropyImageVsP, 'o', 'MarkerSize', 8 );
  
    ylabel('Image Entropy (a.u.)')
    
    set(gca,'XLim',[-0.5,nIter+0.5],'XTick',0:nIter,'YTick',[]), 
    
    legend( [p1,p2,p3,p4,p5,p6], 'C','T','P','C','T','P')        
    
    
    if ( isSaveResults ),
        
        fprintf( '\nEntropy\n' )
        fprintf( '-------\n\n' )
        
        save_figs( outputDir ),
        
        close all,
        
        fprintf( '![](figs/entropy.png)  \n' )
        
    end
    
end

%}


%% Save Final Results


if ( isSaveResults )
    
    [ imCine, tCine ] = recon_cine( imRlt, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMoco, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimAcqRltNii, fov, '' );
    [ imCineQ, tCine ] = recon_cine( imRltQ, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMocoQ, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimRcnRltNii, fov, outputDir );
    PARAM  = P;
    PARAM0 = P0;
    save( fullfile( outputDir, 'results' ), 'imRlt', 'tRlt', 'dtRlt', 'pixdimAcq', 'imRltQ', 'pixdimRcn', 'PARAM0', 'PARAM', 'mask', 'maskQ', 'maskMoco', 'cineRMSD', 'imCine', 'imCineQ', 'tCine', 'TOL', '-v7.3' )

else
    
    [ imCine, tCine ] = recon_cine( imRltQ, dtRlt, P(iIter).R.tRr, P(iIter).R.rrInterval, P(iIter).T.tform, maskMocoQ, P(iIter).P.vox, P(iIter).P.frm, nFrameCine, pixdimRcnRltNii, fov, '' );

end


%% Finalise


if ( isVerbose ) 
    fprintf( '\n\nComplete\n--------\n\nfinish: %s\n\n\n', datestr(now) )
end

diary( 'off' )

copyfile( logFilePath, fullfile( outputDir, 'results.md' ) )


end  % rlt2cine(...)


%% Sub-Function: recon_cine
function [ imCineQ, tCine, voxProbQ ] = recon_cine( imRltQ, dtRlt, tRr, rrInterval, T, maskMocoQ, voxProb, frmProb, nFrameCine, pixdimRcnRltNii, fov, saveDir )


    % Pixel Dimensions
    
    pixdimRcn       = pixdimRcnRltNii(1:2);
    

    % Cine Frame Timing

    tCine = linspace( 0, rrInterval, nFrameCine+1 ); 
    tCine = tCine(1:(end-1));


    % Resize Voxel Probability Maps

    x = linspace(0,fov(2),size(voxProb,2)+1); x = x(2:end)-mean(diff(x))/2;
    y = linspace(0,fov(1),size(voxProb,1)+1); y = y(2:end)-mean(diff(y))/2;
    xq = linspace(0,fov(2),size(imRltQ,2)+1); xq = xq(2:end)-mean(diff(xq))/2;
    yq = linspace(0,fov(1),size(imRltQ,1)+1); yq = yq(2:end)-mean(diff(yq))/2;
    [X,Y] = meshgrid(x,y);
    [Xq,Yq] = meshgrid(xq,yq); 
    voxProbQ = zeros( size( imRltQ ) );
    for iF = 1:size(imRltQ,3),
        voxProbQ(:,:,iF) = interp2(X,Y,voxProb(:,:,iF),Xq,Yq,'linear'); 
    end
    voxProbQ(isnan(voxProbQ)) = 0;  % edge voxels probabilities will get set to NaN during interp2


    % Recon Cine

    imCineQ = imseq_kernel_smooth( tRr, transform_imseq( imRltQ, T, maskMocoQ, pixdimRcn ), tCine, 'tPeriod', rrInterval, 'kSigma', dtRlt, 'kWidth', rrInterval, 'vWeight', transform_imseq( voxProbQ, T, maskMocoQ, pixdimRcn, 'forward', 'linear' ), 'fWeight', frmProb );


    % Save as Nifti

    if ( ~isempty( saveDir ) ),

        dtCine = mean(diff(tCine));
        
        pixdimRcnCineNii = [ pixdimRcnRltNii(1:2), pixdimRcnRltNii(3)/dtRlt*dtCine ] ;
        originRcnCineNii = [ 0 0 pixdimRcnRltNii(3)*dtRlt/2 ];

        im = abs(imCineQ); 
        imMax = prctile(im(:),99.99);  % removing very large values to improve visualisation, TODO: improve, e.g., imMax = mean(im(:))+6*std(im(:)); or imMax = abs( complex( mean(real(imCineQ(:))) + 6*std(real(imCineQ(:))), mean(imag(imCineQ(:))) + 6*std(imag(imCineQ(:))) ) );
        im(im>imMax) = imMax;  

        save_nii( make_nii( im, pixdimRcnCineNii, originRcnCineNii ), fullfile( saveDir, 'cineq_xyt.nii.gz' ) );

        save_nii( make_nii( abs(voxProbQ), pixdimRcnRltNii ), fullfile( saveDir, 'voxprobq_xyt.nii.gz' ) );

    end

    
end  % recon_cine(...)


%% Sub-Function: outlier_rejection_step
function [ voxProb, frmProb, voxParam, frmParam ] = outlier_rejection_step( imSeq, tSeq, dtRlt, rrInterval, varargin )
%OUTLIER_REJECTION  
%
%       input         description
% 
%   Additional name-value input: 
%
%       mask            logical image mask; default all voxels included
%       transformations TODO
%       transformationMask  TODO
%       pixdim          TODO
%       ...             ...
%       verbose         verbosity; default false
%
%   See also: estimate_voxel_probability, estimate_frame_probability

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Parse Input

default.mask            = true(size(imSeq(:,:,1))); 
default.T               = struct([]);
default.maskMoco        = true(size(imSeq(:,:,1))); 
default.pixdim          = [ 2, 2 ];
default.voxProb         = ones(size(imSeq));
default.frmProb         = ones(1,size(imSeq,3));
default.isVerbose       = false;

p = inputParser;

addRequired(  p, 'imSeq', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );

addRequired(  p, 'tSeq', @(x) validateattributes( x, {'numeric'}, ...
        {'numel',size(imSeq,3)}, mfilename ) );

addRequired(  p, 'dtRlt', @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','positive'}, mfilename ) );

addRequired(  p, 'rrInterval', @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','positive'}, mfilename ) );

addParameter( p, 'mask', default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(imSeq(:,:,1)),NaN]}, mfilename ) );

addParameter( p, 'transformations', default.T, ...
        @(x) validateattributes( x, {'struct'}, ...
        {'numel',size(imSeq,3)}, mfilename ) );

addParameter( p, 'transformationMask', default.maskMoco, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(imSeq(:,:,1)),NaN]}, mfilename ) );

addParameter( p, 'pixdim', default.pixdim, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','size',[1,2]}, mfilename ) );

addParameter( p, 'voxProb', default.voxProb, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',size(imSeq),'>=',0,'<=',1}, mfilename ) );

addParameter( p, 'frmProb', default.frmProb, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'numel',size(imSeq,3),'>=',0,'<=',1}, mfilename ) );

addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, imSeq, tSeq, dtRlt, rrInterval, varargin{:} );

mask            = p.Results.mask;
T               = p.Results.transformations;
maskMoco        = p.Results.transformationMask;
pixdim          = p.Results.pixdim;
voxProb         = p.Results.voxProb;
frmProb         = p.Results.frmProb;
isVerbose       = p.Results.verbose;


%% Initialise

if isempty( T ),
    T        = init_tform_struct( size(imSeq,3) );
    maskMoco = mask;
end


%% Setup

pctexc = 0.1;


%% Anon Functions

nFrame      = size(imSeq,3);  % number of frames in imSeq
imSeqT      = transform_imseq( imSeq, T, maskMoco, pixdim );  % imSeq transformed

tform_pmap  = @( probMap ) transform_imseq( probMap, T, maskMoco, pixdim, 'forward', 'linear' );
calc_imref  = @( voxProb, frmProb ) transform_imseq( imseq_kernel_smooth( tSeq, imSeqT, tSeq, 'tPeriod', rrInterval, 'kSigma', dtRlt, 'kWidth', rrInterval, 'vWeight', tform_pmap( voxProb ), 'fWeight', frmProb, 'fExclude', num2cell(1:nFrame) ), T, maskMoco, pixdim, 'invert' );


%% Compute Reference Images

imRef = calc_imref( voxProb, frmProb );


%% Estimate Probabilities

[ voxProb, voxParam ] = estimate_voxel_probability( imSeq, imRef, 'mask', mask, 'pctexc', pctexc, 'verbose', isVerbose );
[ frmProb, frmParam ] = estimate_frame_probability( voxProb, 'mask', mask, 'verbose', isVerbose );


end  % outlier_rejection_step(...)


%% Sub-Function: calc_cine_timing_offset
function tOffset = calc_cine_timing_offset( imCine, tCine, mask )

nVoxMask = sum(sum(mask,2),1);

sig = squeeze( sum(sum(bsxfun(@times,abs(imCine),mask),2),1) ./ nVoxMask );

[ ~, indSig ] = max( sig );

indTgt = size( imCine, 3);

tOffset = tCine(indTgt) - tCine(indSig);

end  % calc_cine_timing_offset(...)