function make_animated_gifs( dataDir, resultsDir )

  
xf2xt = @( xf ) ifft( ifftshift( xf, 3 ), [], 3 );

dtLimit = 1/20;  % max 20 fps animated gif frame rate playable in browser per http://nullsleep.tumblr.com/post/16524517190/animated-gif-minimum-frame-delay-browser

origDir = pwd;
cdToOrig = onCleanup( @() cd( origDir ) );
cd( resultsDir )


%% Diary 

diaryFile = fullfile( resultsDir, 'summary.md'); 
if ( exist( diaryFile, 'file' ) )
    delete( diaryFile );
end
diary( diaryFile );
closeDiary = onCleanup( @() diary( 'off' ) );

fprintf( '## Real-Time\n\n' );
fprintf( '### Training Data\n\n' )
fprintf( '![](rlt_sm027a100f096/figs/trn_full.gif) ![](rlt_sm027a100f096/figs/trn_crop.gif)  \n\n' )
fprintf( '### Under-Sampled Data Sliding Window Recon\n\n' )
fprintf( '![](rlt_sm027a100f096/figs/slw_full.gif) ![](rlt_sm027a100f096/figs/slw_crop.gif)   \n\n' )
fprintf( '### _k-t_ SENSE Reconstruction\n\n' )
fprintf( '![](rlt_sm027a100f096/figs/rlt_full.gif) ![](rlt_sm027a100f096/figs/rlt_crop.gif)  \n\n' )
fprintf( '## Cine \n\n' )
fprintf( '### $X$\n\n' )
fprintf( 'L-R: $X$, $X^{(T)}$, $X^{(T)}$ reordered based on $\\theta$  \n' )
fprintf( '![](rlt_sm027a100f096/figs/rlt_crop.gif) \n' )
fprintf( '![](cine/figs/rltt_crop.gif) \n' )
fprintf( '![](cine/figs/rdr_crop.gif)  \n\n' )
fprintf( '### $Y^{[1]}$\n\n' )
fprintf( 'L-R: $Y^{[1]}$ after cardsync, cardsync+moco, cardsync+moco+outrej  \n' )
fprintf( '![](cine/figs/cine1a_crop.gif) \n' )
fprintf( '![](cine/figs/cine1b_crop.gif) \n' )
fprintf( '![](cine/figs/cine1c_crop.gif)  \n\n' )
fprintf( '### $Y^{[m]}$\n\n' )
fprintf( 'L-R: cardsync+moco, cardsync+outrej, full pipeline   \n' )
fprintf( '![](cine_noOutrej/figs/cineAB_crop.gif) \n' )
fprintf( '![](cine_noMoco/figs/cineAC_crop.gif) \n' )
fprintf( '![](cine/figs/cine_crop.gif)  \n' )
fprintf( '### _k-t_ SENSE Reconstruction (uniform regularisation)\n\n' )
fprintf( '![](rlt_sm027a001f096/figs/rlt_full.gif) ![](rlt_sm027a001f096/figs/rlt_crop.gif)  \n\n' )
fprintf( '## Uniform _k-t_ SENSE Regularisation  \n\n' )
fprintf( 'L-R: $X^{(T)}$, $X^{(T)}$ reordered based on $\\theta$  \n' )
fprintf( '![](rlt_sm027a001f096/figs/rltt_crop.gif) \n' )
fprintf( '![](rlt_sm027a001f096/figs/rdr_crop.gif)  \n\n' )
fprintf( 'L-R: $X$, $Y$  \n' )
fprintf( '![](rlt_sm027a100f096/figs/rlt_crop.gif) \n' )
fprintf( '![](cine/figs/cine1c_crop.gif)  \n\n' )


%% Make Gifs

load( fullfile( resultsDir, 'rlt_sm027a100f096', 'xtsense_recon.mat' ), 'xfTrn' );

N = dir( fullfile( dataDir, '*_xtslw.mat' ) );
load( fullfile( dataDir, N(1).name ), 'xtSlw' );

load( fullfile( resultsDir, 'cine', 'results.mat' ), 'tRlt', 'imRltQ', 'PARAM', 'maskQ', 'pixdimAcq', 'pixdimRcn', 'RESULTS' );

xtqDim = size( imRltQ );
[~,indRlt2Crd] = sort( PARAM(end).R.tRr );

dtRlt = median( diff( tRlt ) );

imSlw       = kspace_zeropad_interp( xtSlw( :, :, (end-xtqDim(3)+1):end ), xtqDim );
imTrn       = kspace_zeropad_interp( xf2xt( xfTrn( :, :, (end-xtqDim(3)+1):end ) ), xtqDim );
imRlt       = imRltQ;
imRltT      = transform_imseq( imRltQ, PARAM(end).T.tform, maskQ, pixdimRcn );
imRdr       = imRltT(:,:,indRlt2Crd);
imCine      = PARAM(end).P.imCine;
imCine1a    = PARAM(1).R.imCine;
imCine1b    = PARAM(1).T.imCine;
imCine1c    = PARAM(1).P.imCine;

% crop window indices
[ indRow, indCol ] = get_im_crop_indices( maskQ, pixdimRcn(1), pixdimRcn(1) );

% signal intensity limits
imLimit     = [ 0, prctile( reshape( abs( imRlt(indRow,indCol,:) ), 1, [] ), 98 ) ];
imLimitSlw  = [ 0, prctile( reshape( abs( imSlw(indRow,indCol,:) ), 1, [] ), 98 ) ];
imLimitTrn  = [ 0, prctile( reshape( abs( imTrn(indRow,indCol,:) ), 1, [] ), 98 ) ];

% sliding window
gifFilePath = ims2gif( abs(imSlw), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'slw_full', 't', dtRlt, 'imlimits', imLimitSlw, 'spatialScaling', 1 );
    % fprintf( '![](%s)  \n', gifFilePath ) 
gifFilePath = ims2gif( abs(imSlw(indRow,indCol,:)), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'slw_crop', 't', dtRlt, 'imlimits', imLimitSlw, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% training data
gifFilePath = ims2gif( abs(imTrn), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'trn_full', 't', dtRlt, 'imlimits', imLimitTrn, 'spatialScaling', 1 );
    % fprintf( '![](%s)  \n', gifFilePath ) 
gifFilePath = ims2gif( abs(imTrn(indRow,indCol,:)), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'trn_crop', 't', dtRlt, 'imlimits', imLimitTrn, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% real-time
gifFilePath = ims2gif( abs(imRlt), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'rlt_full', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 1 );
    % fprintf( '![](%s)  \n', gifFilePath ) 
gifFilePath = ims2gif( abs(imRlt(indRow,indCol,:)), 'filedir', fullfile( 'rlt_sm027a100f096', 'figs' ), 'filename', 'rlt_crop', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% transformed
gifFilePath = ims2gif( abs(imRltT(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'rltt_crop', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% reordered
gifFilePath = ims2gif( abs(imRdr(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'rdr_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% cine (final)
gifFilePath = ims2gif( abs(imCine(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'cine_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% cine1a
gifFilePath = ims2gif( abs(imCine1a(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'cine1a_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% cine1b
gifFilePath = ims2gif( abs(imCine1b(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'cine1b_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 

% cine1c
gifFilePath = ims2gif( abs(imCine1c(indRow,indCol,:)), 'filedir', fullfile( 'cine', 'figs' ), 'filename', 'cine1c_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    % fprintf( '![](%s)  \n', gifFilePath ) 


%% cine cardsync+moco / cine cardsync+outrej

try
    % cine cardsync+moco (final)
    SnoOutrej   = load( fullfile( resultsDir, 'cine_noOutrej', 'results.mat' ), 'PARAM', 'RESULTS' );
    imCineAB    = SnoOutrej.PARAM(end).P.imCine;   
    gifFilePath = ims2gif( abs(imCineAB(indRow,indCol,:)), 'filedir', fullfile( 'cine_noOutrej', 'figs' ), 'filename', 'cineAB_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
catch
end

try
    % cine cardsync+outrej (final)
    SnoMoco     = load( fullfile( resultsDir, 'cine_noMoco', 'results.mat' ), 'PARAM', 'RESULTS' );
    imCineAC    = SnoMoco.PARAM(end).P.imCine;
    gifFilePath = ims2gif( abs(imCineAC(indRow,indCol,:)), 'filedir', fullfile( 'cine_noMoco', 'figs' ), 'filename', 'cineAC_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
catch
end


%% rlt alpha=1

try
    
    load( fullfile( resultsDir, 'cine_rlt_sm027a001f096', 'results.mat' ), 'imRltQ', 'imCineQ' );

    dtRlt = median( diff( tRlt ) );

    imRltNoMask  = imRltQ;
    imRltNoMaskT = transform_imseq( imRltQ, PARAM(end).T.tform, maskQ, pixdimRcn );
    imRdrNoMask  = imRltNoMaskT(:,:,indRlt2Crd);
    imCineNoMask = imCineQ;
    
    % real-time 
    ims2gif( abs(imRltNoMask), 'filedir', fullfile( 'rlt_sm027a001f096', 'figs' ), 'filename', 'rlt_full', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 1 );
    ims2gif( abs(imRltNoMask(indRow,indCol,:)), 'filedir', fullfile( 'rlt_sm027a001f096', 'figs' ), 'filename', 'rlt_crop', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 2 );
    
    % transformed
    ims2gif( abs(imRltNoMaskT(indRow,indCol,:)), 'filedir', fullfile( 'cine_rlt_sm027a001f096', 'figs' ), 'filename', 'rltt_crop', 't', dtRlt, 'imlimits', imLimit, 'spatialScaling', 2 );
    
    % reordered
    ims2gif( abs(imRdrNoMask(indRow,indCol,:)), 'filedir', fullfile( 'cine_rlt_sm027a001f096', 'figs' ), 'filename', 'rdr_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
    
    % cine (final)
    ims2gif( abs(imCineNoMask(indRow,indCol,:)), 'filedir', fullfile( 'cine_rlt_sm027a001f096', 'figs' ), 'filename', 'cine_crop', 't', dtLimit, 'imlimits', imLimit, 'spatialScaling', 2 );
      
catch
end


end   % make_animated_gifs(...)