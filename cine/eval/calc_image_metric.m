function metricVal = calc_image_metric( im, mask, t, rr )


% NOTE: assuming static mask and same mask for im and imRef


%% Utility Functions

get_mask_values_vec = @( im, msk ) double( abs( im( repmat( msk, [1,1,size(im,3)/size(msk,3)] ) ) ) );

get_mask_values_3d  = @( im, msk ) reshape( get_mask_values_vec( im, msk ), sum(msk(:)), 1, size(im,3) );



%% Entropy

% {

s = get_mask_values_vec( im, mask );

b = sqrt( sum( get_mask_values_vec( im, mask ).^2, 1 ) );
    % b = mean( sqrt( sum( get_mask_values_3d( im, mask ).^2, 1 ) ) );

p = bsxfun( @rdivide, s, b );

h = - mean( p(:) .* log( p(:) ) );

metricVal = h;

%}


%% Histogram Entropy

%{

metricVal = entropy( get_mask_values_vec( im, mask ) ); 

%}


%% Image Entropy, per Roy2012 

% v0
%{

metricVal = imagemetric( get_mask_values_3d( im, mask ), {'Cine'} ) / numel( get_mask_values_vec( im, mask ) );

%}


% v1
%{ 

calc_max_val_j      = @( ref, msk ) sqrt( nansum( get_mask_values_3d( ref, msk ).^2, 3 ) );

calc_norm_val       = @( im, msk, im ) bsxfun( @rdivide, get_mask_values_3d( im, mask ), calc_max_val_j( ref, mask ) );

calc_time_entropy   = @( im, msk, ref ) - sum( calc_max_val_j( ref, mask ) / sum( calc_max_val_j( ref, mask ), 1 ) .* sum( calc_norm_val( im, mask, ref ) .* log( calc_norm_val( im, mask, ref ) ), 3 ), 1 );

metricVal = calc_time_entropy( im, mask, ref );

%} 


% v2
%{

s = get_mask_values_3d( im, mask );

b = sqrt( sum( get_mask_values_3d( im, mask ).^2, 1 ) );

p = bsxfun( @rdivide, s, b );

h = sum( ( b / sum( b, 3 ) ) .* - sum( p .* log( p ), 1 ), 3 ) / numel( s );

metricVal = h;

%}


%% Entropy of Gradient

%{
dx = 1.25;
dy = 1.25;

millimetrePerRr = 1e-1*mean([dx,dy])*2*sqrt(sum(mask(:))/pi);

[g,gx,gy] = calc_imseq_gradient( im,  dx, dy, millimetrePerRr * t, millimetrePerRr * rr );
% g = sqrt( gx.^2 + gy.^2 );

calc_max_val = @( r, msk ) mean( sum( get_mask_values_3d( r, msk ), 1 ), 3 );  % calc_max_val = @( r, msk ) sqrt( sum( get_mask_values_vec( r, msk ).^2 ) );

p = g / calc_max_val( g, mask );

calc_grad_entropy = @( p, msk ) - mean( get_mask_values_vec( p, msk ) .* log2( get_mask_values_vec( p, msk ) ) );

metricVal = calc_grad_entropy ( p, mask );
%}

end  % calc_image_metric(...)





function MetricValue = imagemetric(Images,DataType)
%   MetricValue = imagemetric(Images)
%   This function takes a 3D image series (xyt), or a cropped subset of a
%   series and returns a scalar value which is the time-entropy of those
%   images
%
%   Inputs:
%   Images          - 3D series of images.  The x and y dimensions are not
%                   important (ie can be collapsed) but the 3rd must be
%                   preserved)
%
%   Outputs:
%   MetricValue     - Scaler representing the metric value


if strcmp(DataType{1},'PC')
% Check inputs
if length(size(Images)) ~= 3
    error('images is not a 3D array')
end

%% Rectify images
Images = abs(Images);

%% Calculate time-sums
PixelTotals = sqrt(nansum(Images.^2,3));
PixelTotals(PixelTotals == 0) = 1;

%% Normalize Images
NormalizedImages = Images./PixelTotals(:,:,ones(size(Images,3),1));

%% Eliminate zeros from Images
NormalizedImages(NormalizedImages == 0) = 1;

%% Compute entropy
PixelEntropies = nansum(NormalizedImages.*log(NormalizedImages),3);
if (isnan(nansum(PixelTotals(:)))==1)
1
end
if nansum(PixelTotals(:))
    PixelWeights = PixelTotals/nansum(PixelTotals(:));
else
    PixelWeights = zeros(size(PixelTotals));
end

MetricValue = -nansum(nansum(PixelEntropies.*PixelWeights,1),2);


else
    % Check inputs
if length(size(Images)) ~= 3
%     error('images is not a 3D array')
end
% Rectify images
dImages = abs(double(Images));

% Calculate Space-sums
PixelTotals = sqrt(nansum(nansum(dImages.^2,2),1));
PixelTotals(PixelTotals == 0) = 1;

% Normalize Images
NormalizedImages = dImages./PixelTotals(ones(size(dImages,1),1),ones(size(dImages,2),1),:);

% Eliminate zeros from Images
NormalizedImages(NormalizedImages == 0) = 1;

% Compute entropy
PixelEntropies = nansum(nansum(NormalizedImages.*log(NormalizedImages),2),1);
if nansum(PixelTotals(:))
    PixelWeights = PixelTotals/sum(PixelTotals(:));
else
    PixelWeights = zeros(size(PixelTotals));
end
MetricValue = -nansum(PixelEntropies.*PixelWeights,3);
end

end  % imagemetric()