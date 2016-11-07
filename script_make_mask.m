% script_make_mask

fcmrNo = 000;

seriesNo = 00;


idStr = sprintf( 'fcmr%03is%02i', fcmrNo, seriesNo );

% Number of Real-Time Frames to Reconstruct

numFrame = 96;


srcDir = fullfile( '/Volumes/fcmr_cine/data/', idStr );

[ ~, fileStr ]  = fileparts( ls( fullfile( srcDir, '*_ktdata.mat' ) ) );
filePrefixStr   = fileStr(1:(end-7));

xtslwMatFilePath    = fullfile( srcDir, sprintf( '%s_xtslw.mat', filePrefixStr ) );

maskMatFilePath     = fullfile( srcDir, 'maskfetalheart.mat' );

% Load Sliding Window Data

X = matfile( xtslwMatFilePath );
xtSlw = X.xtSlw; 

xtSlw = xtSlw(:,:,(end-numFrame+1):end);

imagine( xtSlw );

mask = uiget_roi_mask( abs( mean( xtSlw, 3 ) ), 'fetal heart' );
save( maskMatFilePath, 'mask', '-v7.3' );

close, 

clear