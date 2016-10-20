function h = matgrid( varargin )
%MATGRID   Display a pseudocolour ('checkerboard') matrix plot.
%   MATGRID(M) is a pseudocolor or 'checkerboard' plot of matrix M.
%
%   MATGRID(X,Y,M), where X and Y are vectors, makes a
%   pseudocolor plot on the grid defined by X and Y.  
%
%   MATGRID(AX,..) plots into AX instead of GCA.
%
%   H = MATGRID(...) returns a handle to the SURFACE object.
%
%   MATGRID is based on PCOLOR. 
%
%   See also IMGRID, PCOLOR.

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Parse Inputs and Check for Errors


% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 3
    error(message('MATLAB:narginchk:tooManyInputs'));
end
% do error checking before calling newplot. This argument checking should
% match the surface(x,y,z) or surface(z) argument checking.
if nargs == 2
  error(message('MATLAB:pcolor:InvalidNumberOfInputs'))
end
if isvector(args{end})
  error(message('MATLAB:pcolor:NonMatrixColorInput'));
end
if nargs == 3 && LdimMismatch(args{1:3})
  error(message('MATLAB:pcolor:InputSizeMismatch'));
end
for k = 1:nargs
  if k<nargs,
      if ~isvector(args{k})
          error( 'Grid data input must be a vector' );
      end
  end
  if ~isreal(args{k})
    error(message('MATLAB:pcolor:NonRealInputs'));
  end
end


%% Set Up Axes


cax = newplot(cax);
hold_state = ishold(cax);


%% Prep Input Arguments and Create Surface


if nargs == 1
    
    x = args{1};
    x = [ x; x(end,:) ];
    x = [ x, x(:,end) ];
    hh = surface(zeros(size(x)),x,'Parent',cax);
    [m,n] = size(x);
    lims = [ 0 n-1 0 m-1];

elseif nargs == 3

    [x,y,c] = deal(args{1:3});
    x = [ x(:); x(end)+mean(diff(x(:))) ];
    y = [ y(:); y(end)+mean(diff(y(:))) ];
    c = flip( c, 1 );
    c = [ c; c(end,:) ];
    c = [ c, c(:,end) ];
    hh = surface(x,y,zeros(size(c)),c,'Parent',cax);
    lims = [min(min(x)) max(max(x)) min(min(y)) max(max(y))];

end


%% Adjust Surface and Axis 

set(hh,'AlignVertexCenters','on');
if ~hold_state
    set(cax,'View',[0 90]);
    set(cax,'Box','on');
    if lims(2) <= lims(1)
        lims(2) = lims(1)+1;
    end
    if lims(4) <= lims(3)
        lims(4) = lims(3)+1;
    end
    axis(cax,lims);
end
if nargout == 1
    h = hh;
end

axis image off,
shading faceted,


end  % imgrid(...)


function ok = LdimMismatch(x,y,z)
[xm,xn] = size(x);
[ym,yn] = size(y);
[zm,zn] = size(z);
ok = (xm == 1 && xn ~= zn) || ...
     (xn == 1 && xm ~= zn) || ...
     (xm ~= 1 && xn ~= 1 && (xm ~= zm || xn ~= zn)) || ...
     (ym == 1 && yn ~= zm) || ...
     (yn == 1 && ym ~= zm) || ...
     (ym ~= 1 && yn ~= 1 && (ym ~= zm || yn ~= zn));
end  % LdimMismatch(...)
