function U = rreg_imseq2imseq( imSeqSrc, imSeqTgt, mask, pixdim )
%RREG_IMSEQ2IMSEG  
%
%   U = RREG_IMSEQ2IMSEG( imSeqSrc, imSeqTgt, mask )
%
%   U = RREG_IMSEQ2IMSEG( imSeqSrc, imSeqTgt, mask, pixdim ) use specified
%   pixel dimensions
% 
%   See also: transform_imseq, make_imref2d

% jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Setup

if ~exist( 'pixdim', 'var' ),
    pixdim = [1,1];
end

tform2param = @(U) [ U.T(3,1), U.T(3,2), asind(U.T(1,2)) ];  % NOTE: for reference, http://uk.mathworks.com/discovery/affine-transformation.html
    % param2tform = @(tx,ty,rz) affine2d( [cosd(rz), sind(rz), 0; -sind(rz), cosd(rz), 0; tx, ty, 1 ] );

nF = size(imSeqSrc,3);


%% Frame-to-Frame Rigid Registration

U = struct();

for iF = 1:nF,
    
    U(iF).A = rreg_im2im( imSeqSrc(:,:,iF), imSeqTgt(:,:,iF), mask, pixdim );
    
    p        = tform2param( U(iF).A );
    U(iF).tx = p(1);
    U(iF).ty = p(2);
    U(iF).rz = p(3);

end


end  % rreg_imseq2imseq(...)


function A = rreg_im2im( imSrc, imTgt, bwMsk, pixdim )

bwSrc = bwMsk;

maxDisp = 15;  % mm
rMax    = floor( maxDisp / min(pixdim ) ); 
rEst    = floor(sqrt(sum(bwSrc(:))/pi));  % radius, assuming bwSrc is circular
r       = min( rMax, rEst );  
bwTgt   = bwmorph( bwSrc, 'dilate', r );

[OPT,METRIC] = imregconfig('monomodal');
OPT.GradientMagnitudeTolerance  = 1e-6;
OPT.MaximumStepLength           = 6.25e-3;
OPT.MinimumStepLength           = 1e-6;

T0 = [  1 0 0 ; ...     % null transform
        0 1 0 ; ...
        0 0 1 ];
A  = affine2d( T0 );    % initial 2d affine transform

sigma = [ 0.8, 0.6, 0.4 ]; % source image smoothing, in pixels

for iS = 1:numel(sigma),
    
    srcCpx = complex(imgaussfilt(real(imSrc),sigma(iS)),imgaussfilt(imag(imSrc),sigma(iS)));
    [ srcR, src ] = make_imref2d( abs(srcCpx), bwMsk, pixdim, bwSrc );
    
    [ tgtR, tgt ] = make_imref2d( abs(imTgt),  bwMsk, pixdim, bwTgt );
    
    A = imregtform(src,srcR,tgt,tgtR,'rigid',OPT,METRIC,'DisplayOptimization',false,'PyramidLevels',1,'InitialTransformation',A);  
        % NOTE: imregtform uses bilinear interpolation, http://uk.mathworks.com/matlabcentral/answers/73537-how-does-imregister-work

end

end  % rreg_im2im(...)