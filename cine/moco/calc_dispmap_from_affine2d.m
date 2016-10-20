function dispMap = calc_dispmap_from_affine2d( A, R )
%CALC_DISPMAP_FROM_AFFINE_2D  calculate displacement map for rigid transformation
%
%   dispMap = calc_dispmap( A, R ) returns displacement map, 
%   dispmap, for each voxel specified by imref2d object, R, subject to 
%   rigid tform struc, A.

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Voxel Index Vectors


xVec = linspace( R.XWorldLimits(1)/R.PixelExtentInWorldX, R.XWorldLimits(2)/R.PixelExtentInWorldX, R.ImageSize(2)+1 ); 
xVec = xVec(2:end);

yVec = linspace( R.YWorldLimits(1)/R.PixelExtentInWorldY, R.YWorldLimits(2)/R.PixelExtentInWorldY, R.ImageSize(1)+1 ); 
yVec = yVec(2:end);


%% Calculate


% Coordinate Maps

[ xMap, yMap ] = meshgrid( xVec, yVec );


% Transformed Coordinates

[ uMap, vMap ] = transformPointsForward( A, xMap, yMap );


% Displacement

txMap = uMap - xMap;
tyMap = vMap - yMap;


% Combine Displacement (x,y) as Complex (real,imag)

dispMap = complex( txMap, tyMap );


end  % calc_dispmap(...)