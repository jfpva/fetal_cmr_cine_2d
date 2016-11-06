function [ indRow, indCol ] = get_im_crop_indices( mask, dR, dC )

cropWinSize = 100; % mm

pt = bwmorph( mask, 'shrink', inf );
[iR,iC] = find(pt);

rowEndPts = iR + ( cropWinSize / dR ) / 2 * [ -1, +1 ];
if rowEndPts(1) < 0, 
    rowEndPts = rowEndPts - rowEndPts(1);
end

indRow = colon( round( rowEndPts(1) ), round( rowEndPts(2) ) );

colEndPts = iC + ( cropWinSize / dC ) / 2 * [ -1, +1 ];
if colEndPts(1) < 0, 
    colEndPts = colEndPts - colEndPts(1);
end

indCol = colon( round( colEndPts(1) ), round( colEndPts(2) ) );


end  % get_im_crop_indices(...)