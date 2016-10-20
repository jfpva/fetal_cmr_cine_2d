function save_xt2nii( xt, Param, saveFilePath )
%SAVE_XT2NII
%
%   SAVE_XT2NII( xt, Param, saveFilePath )

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Setup

origPath  = path;
resetPath = onCleanup( @() path( origPath ) );
addpath( '~/Research/fcmr_cine/aux/nifti/' ),


%% Get Acqn Parameters

I = Param.acq.ImageInformation(1,1,2);

dx = I.ACQVoxelSize(1);
dy = I.ACQVoxelSize(2);
dz = I.ACQVoxelSize(3);
dt = I.DynamicScanTime;
if dt > 0.1,
    fprintf( 'Dynamic scan time (%g s) larger than expected; caluclating as num. ky lines x TR.\n', dt )
    tr = 2 * I.EchoTime / 1000;
    nky = I.Resolution(2) / Param.acq.Scan.KtFactor;
    dt = nky * tr;
end

%% Save as Nifti Files

XYZT = make_nii( abs(permute(xt,[1,2,4,3])) ); 
XYZT.hdr.dime.pixdim(2:5) = [ dx dy dz dt ];
save_nii( XYZT, saveFilePath );

[ saveFileDir, saveFileName ] = fileparts( saveFilePath );
[ ~, saveFileName ] = fileparts( saveFileName );
saveFilePathXYT = fullfile( saveFileDir, strcat( saveFileName, '_xyt.nii.gz' ) );
XYT = make_nii( abs(xt), [ dx dy dx ] ); 
XYT.hdr.dime.pixdim(2:4) = [ dx dy dy ];
save_nii( XYT, saveFilePathXYT );


end  % save_xt2nii(...)