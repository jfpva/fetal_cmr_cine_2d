function roiMask = uiget_roi_mask( refIm, figTitle )
%UIGET_ROI_MASK  Get user-drawn ROI mask on reference image
%
%   example:
%
%       roiMask = uiget_roi_mask( refIm, figTitle )
%

if ~exist( 'figTitle', 'var' )
    figTitle = 'Draw ROI';
else
    figTitle = [ 'Draw ROI: ', figTitle ];
end

im = abs(refIm);

imLimits = [ prctile(im(:),1), prctile(im(:),99) ];

h = figure;

imshow( im, imLimits ),

title( figTitle ),

h = imfreehand;

roiMask = h.createMask;

close

end  % make_roi_mask()