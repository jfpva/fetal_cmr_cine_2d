function cmaplog = colormaplog( cmap )
%COLORMAPLOG  log of colormap
%
%   cmaplog = COLORMAPLOG( cmap )
%
%   Example
%       parulalog = COLORMAPLOG( parula );
%
%   See also colormap

% jfpva (joshua.vanamerom@kcl.ac.uk)

n = size(cmap,1);
x  = linspace(1,n,n);
xq = flip( n+1 - logspace( 0, log10(n), n ) );

cmaplog = interp1( x, cmap, xq, 'linear', 'extrap' );

end  % colormaplog(...)