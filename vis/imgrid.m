function varargout = imgrid( varargin )
%IMGRID   Display a pseudocolour ('checkerboard') image plot.
%   IMGRID(IM) is a pseudocolor or 'checkerboard' plot of image IM.
%
%   IMGRID(X,Y,IM), where X and Y are vectors, makes a
%   pseudocolor plot on the grid defined by X and Y.  
%
%   IMGRID(AX,..) plots into AX instead of GCA.
%
%   H = IMGRID(...) returns a handle to the SURFACE object.
%
%   IMGRID is a special case of matgrid.
%
%   See also MATGRID, PCOLOR.

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Call matgrid

v = varargin;

switch numel(v)
    case 1
        h = matgrid( v{1} );
    case 2
        h = matgrid( v{1}, v{2} );
    case 3
        h = matgrid( v{1}, v{2}, v{3} );
    case 4
        h = matgrid( v{1}, v{2}, v{3}, v{4} );
    case 5
        h = matgrid( v{1}, v{2}, v{3}, v{4}, v{5} );
    case 6
        h = matgrid( v{1}, v{2}, v{3}, v{4}, v{5}, v{6} );
    case 7
        h = matgrid( v{1}, v{2}, v{3}, v{4}, v{5}, v{6}, v{7} );
    case 8
        h = matgrid( v{1}, v{2}, v{3}, v{4}, v{5}, v{6}, v{7}, v{8} );
    otherwise
        error( 'Number of input arguments (%i) exceeds expected.', numel(varargin) )
end


%% Set as Image Orientation

axis ij


%% Set Output

if nargout > 0,
    varargout{1} = h;
end


end  % imgrid(...)

