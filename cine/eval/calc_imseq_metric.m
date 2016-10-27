function metricVal = calc_imseq_metric( imseq, mask )
%CALC_IMSEQ_METRIC calculate quality metric for image sequence
%
%   metricVal = CALC_IMSEQ_METRIC( imseq, mask )

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Utility Functions


get_mask_values_vec = @( im, msk ) double( abs( im( repmat( msk, [1,1,size(im,3)/size(msk,3)] ) ) ) );

get_mask_values_3d  = @( im, msk ) reshape( get_mask_values_vec( im, msk ), sum(msk(:)), 1, size(im,3) );


%% Calculate Entropy


s = get_mask_values_vec( imseq, mask );
b = sqrt( sum( get_mask_values_vec( imseq, mask ).^2, 1 ) );  
p = bsxfun( @rdivide, s, b );
h = - sum( p(:) .* log( p(:) ) ); 

metricVal = h;

    % NOTE: normalisation in usual entropy calculation would be done as 
    % b = mean( sqrt( sum( get_mask_values_3d( im, mask ).^2, 1 ) ) );
    % however, taking mean b over the total frames provided, yields an
    % entropy value that is independent of number of frames
    % BUT, this doesn't matter if all cines compared are same ROI size and
    % same number of frames, AND entropy is compared between cases
    % normalised by entropy of X[0] or mean(Y,time), so using
	% b = sqrt( sum( get_mask_values_vec( imseq, mask ).^2, 1 ) );
    
    % NOTE: 
    % h = - mean( p(:) .* log( p(:) ) );
    % instead of h = - sum( p(:) .* log( p(:) ) ); as mean 
    % gives an entropy value that is independent of number of voxels in
    % mask and number of frames, allowing for comparison between cases
    % BUT, this doesn't matter if all cines compared are same ROI size and
    % same number of frames, AND entropy is compared between cases
    % normalised by entropy of X[0] or mean(Y,time), so using
    % h = - sum( p(:) .* log( p(:) ) );

    
%% Calculate Image Entropy


%{ 
s = get_mask_values_3d( imseq, mask );

b = sqrt( sum( get_mask_values_3d( imseq, mask ).^2, 1 ) );

p = bsxfun( @rdivide, s, b );

h = sum( ( b / sum( b, 3 ) ) .* - sum( p .* log( p ), 1 ), 3 ); % / numel( s );

metricVal = h;
%}


end  % calc_imseq_metric(...)