function [ R, imOut ] = make_imref2d( imIn, bwMask, pixdim, bwCrop )
%MAKE_IMREF2D  make reference 2-D image to world coordinates object
%
%   R = MAKE_IMREF2D( imIn, bwMask )
%
%   R = make_imref2d( ..., pixdim ) use specified pixel dimensions
%
%   [ R, imOut ] = make_imref2d( ..., bwCrop ) crop imIn to imOut using
%   extents of bwCrop
%
%   See also: imref2d, rreg_imseq2imseq, transform_imseq

% jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Notes

% TODO: add input parsing/validation


%% Setup

if ~exist( 'pixdim', 'var' ),
    pixdim = [1,1];
end

if ~exist( 'bwCrop', 'var' ),
    bwCrop = true( size( imIn ) );
end


%% Crop Image Sequence

iX = find( sum( bwCrop, 2) > 0 );
iY = find( sum( bwCrop, 1) > 0 );

imOut = imIn(iX,iY,:);


%% Find Centroid

props   = regionprops( bwMask, 'centroid' );
cX      = props.Centroid(2);  
cY      = props.Centroid(1);  


%% Create imref2d object

pX = pixdim(1) * ( [ iX(1)-0.5, iX(end)+0.5 ] - cX );
pY = pixdim(2) * ( [ iY(1)-0.5, iY(end)+0.5 ] - cY );

R  = imref2d( size(imOut), pY, pX );


end  % make_imref2d(...)