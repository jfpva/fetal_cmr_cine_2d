function imTrf = transform_imseq( imSeq, U, mask, pixdim, transformDirn, interpMethod )
%TRANSFORM_IMSEQ  
%
%   imTrf = TRANSFORM_IMSEQ( imSeq, U, mask )
%
%   imTrf = TRANSFORM_IMSEQ( imSeq, U, mask, pixdim ) use specified pixel
%   dimensions
% 
%   imTrf = TRANSFORM_IMSEQ( imSeq, U, mask, pixdim, 'invert' ) does
%   inverse transformation
%
%   imTrf = TRANSFORM_IMSEQ( imSeq, U, mask, pixdim, 'forward', interpStr ) 
%   uses interpStr interpolation during transformation

%   See also: rreg_imseq2imseq, make_imref2d

% jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Setup


if ~exist( 'pixdim', 'var' ),
    pixdim = [1,1];
end

if ~exist( 'transformDirn', 'var' ),
    transformDirn = 'forward';
end

if ~exist( 'interpMethod', 'var' ),
    interpMethod = 'cubic';
end

nF = size(imSeq,3);

R = make_imref2d( imSeq(:,:,1), mask, pixdim );

imTrf = zeros( size(imSeq), 'like', imSeq );

for iF = 1:nF,
    
    switch transformDirn,       
        case 'invert'
            A = U(iF).A.invert;
        otherwise
            A = U(iF).A;
    end

    imTrf(:,:,iF) = imwarp( imSeq(:,:,iF), R, A, interpMethod, 'OutputView', R );  
    
end

% imTrf( isnan( imTrf ) ) = 0;


end  % transform_imseq(...)