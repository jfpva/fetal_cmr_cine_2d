function imOut = imseq_kernel_smooth( tSeq, imSeq, tOut, varargin )
%IMSEQ_KERNEL_SMOOTH  image sequence kernel smoothing
%
%   imOut = IMSEQ_KERNEL_SMOOTH( tSeq, imSeq, tOut );
%
%   imOut = IMSEQ_KERNEL_SMOOTH( ..., 'name', value );
%
%   Name-Value Options:
%
%       tPeriod         period, if cyclic
%       kernelMethod    kernel method; default 'sinc'
%       windowMethod    window method; default 'box'
%       kSigma          kernel 'sigma'
%       kWidth          kernel width
%       windowEdgePct   percent of window for each edge of edge+flattop window
%       vWeight         voxel weights
%       fWeight         frame weights
%       fExclude        frames to exclude
%       verbose         verbosity

% 
%   Example
%       imOut = imseq_kernel_smooth(tSeq,imSeq,tOut,'tPeriod',rr,'kSigma',dt,');

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Parse Input

% TODO: add kernel definition/method as optional input

default.tPeriod       = [];                                % period, if cyclic
default.kSigma        = [];                                % kernel 'sigma'
default.kWidth        = [];                                % kernel width
default.kernelMethod  = 'sinc';                            % sinc kernel
default.windowMethod  = 'hann';                            % tukey window
default.windowEdgePct = 15;                                % percent of window for each hann edge
default.vWeight       = ones( size( imSeq ) );             % voxel weights
default.fWeight       = ones( size( tSeq ) );              % frame weights
default.fExclude      = num2cell( zeros( size( tOut ) ) ); % frames to exclude
default.isVerbose     = false;                             % verbosity

p = inputParser;

addRequired(  p, 'tSeq', @(x) validateattributes( x, {'numeric'}, ...
        {'vector','numel',size(imSeq,3)}, mfilename ) );

addRequired(  p, 'imSeq', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );
    
addRequired(  p, 'tOut', @(x) validateattributes( x, {'numeric'}, ...
        {'vector','nonempty'}, mfilename ) );

addParameter( p, 'tPeriod', default.tPeriod, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );

addParameter( p, 'kSigma', default.kSigma, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );    
    
addParameter( p, 'kWidth', default.kWidth, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );    

addParameter( p, 'kernelMethod', default.kernelMethod, ...
        @(x) validateattributes( x, {'char'}, ...
        {'nonempty'}, mfilename ) );    

addParameter( p, 'windowMethod', default.windowMethod, ...
        @(x) validateattributes( x, {'char'}, ...
        {'nonempty'}, mfilename ) );    

addParameter( p, 'windowEdgePct', default.windowEdgePct, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','>=',0,'<=',100}, mfilename ) );    
    
addParameter( p, 'vWeight', default.vWeight, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'nonnegative','size',size(imSeq),'<=',1,'>=',0}, mfilename ) );

addParameter( p, 'fWeight', default.fWeight, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'nonnegative','vector','<=',1,'>=',0}, mfilename ) );    

addParameter( p, 'fExclude', default.fExclude, ...
        @(x) validateattributes( x, {'cell'}, ...
        {'vector','numel',numel(tOut)}, mfilename ) );    

addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, tSeq, imSeq, tOut, varargin{:} );

tPeriod       = p.Results.tPeriod;
kSigma        = p.Results.kSigma;
kWidth        = p.Results.kWidth;
kernelMethod  = p.Results.kernelMethod;
windowMethod  = p.Results.windowMethod;
windowEdgePct = p.Results.windowEdgePct;
vWeight       = p.Results.vWeight;
fWeight       = p.Results.fWeight;
fExclude      = p.Results.fExclude;
isVerbose     = p.Results.verbose;


%% Initialise

if ( isempty( tPeriod ) || tPeriod == 0 || isinf( tPeriod ) ),
    tPeriod = 0;
elseif ( ( tPeriod < max( tSeq ) ) || ( tPeriod < max( tOut ) ) ) 
        warning( 'specified input or output sequence time is longer than period' )
end

% Kernel Width
if isempty( kWidth ),
    if ( tPeriod == 0 ),
        kWidth = range( tSeq );
    else
        kWidth = tPeriod;
    end
end
if ( kWidth == 0 ),
    kWidth = eps*1e2;
end

% Kernel Sigma   % TODO: improve default kernel sigma value
if isempty( kSigma ),
    if ( tPeriod == 0 ),
        kSigma = range( tSeq ) / 5;
    else
        kSigma = tPeriod / 5;
    end
end
if ( kSigma == 0 ),
    kSigma = eps*1e2;
end


%% Setup

% Padding for Cyclic Treatment of Data

if ( tPeriod == 0 ),
	nCycPad = 0;
else
    lowerNumCycPad  = abs( floor( ( min(tOut) - kWidth/2 ) / tPeriod ) );
    upperNumCycPad  = abs( ceil( ( max(tOut) + kWidth/2 ) / tPeriod ) ) - 1;
    nCycPad         = max( lowerNumCycPad, upperNumCycPad );    % pad size
end

% Total cycles
nCycTot = 1 + 2 * nCycPad;  

% Timing
tSeqPad = reshape(bsxfun(@plus,(-nCycPad:+nCycPad)*tPeriod,repmat(tSeq(:),1,nCycTot)),1,[]);

% Frame Numbers
fSeqPad = reshape( repmat( (1:size(imSeq,3))', 1 , nCycTot ), 1, [] );


%% Define Kernel

% NOTE: may not be correctly normalised here; dealt with later

% Boxcar
rectWin = @(t,t0,dt,w) 1 * ( abs(t-t0) <= w/2 );

% Gauss, NOTE: assuming FWHM = dt, and sigma = ~ FWHM/2.355 for normal distribution
gaussKernel = @(t,t0,dt,w) ( 1 / ( sqrt(2*pi) * (dt/2.355) ) ) * exp( -( t - t0 ).^2 / ( 2 * (dt/2.355).^2 ) );

% Sinc
sincKernel = @(t,t0,dt,w) sinc((t-t0)/dt);

% Window: window and truncate using central lobe of sinc
sincWin    = @(t,t0,dt,w) rectWin(t,t0,[],w) .* sincKernel(t,t0,w/2,[]);

% Window: truncate sinc at last zero crossing
zcR =  @( t, t0, dt, w ) t0 + dt*floor((w/2)/dt);  % first zero-crossing in kernel window
zcL =  @( t, t0, dt, w ) t0 - dt*floor((w/2)/dt);  % last  zero-crossing in kernel window
truncZcWin = @( t, t0, dt, w ) rectWin( t, t0, dt, zcR(t,t0,dt,w)-zcL(t,t0,dt,w) );  

% Window: window and truncate using flat top with Hanning edges
hannWin      = @( t, t0, dt, w ) 0.5*(1+cos(2*pi*((t-t0)/(w)))) .* rectWin( t, t0, dt, w );
halfrectWinR = @( t, t0, dt, w ) rectWin(t,t0,dt,w) .* ( t-t0 >= 0 );
halfrectWinL = @( t, t0, dt, w ) rectWin(t,t0,dt,w) .* ( t-t0 <= 0 );
halfhanWinR  = @( t, t0, dt, w ) halfrectWinR(t,t0,dt,w) .* hannWin(t,t0,dt,w);
halfhanWinL  = @( t, t0, dt, w ) halfrectWinL(t,t0,dt,w) .* hannWin(t,t0,dt,w);
p = 2 * windowEdgePct/100;
hannPctWin = @( t, t0, dt, w ) halfhanWinL(t,t0-w*(1-p)/2,[],p*w) + rectWin(t,t0,dt,w*(1-p-eps)) + halfhanWinR(t,t0+w*(1-p)/2,[],p*w);

% Kernel Function
switch kernelMethod,
    case 'sinc',
        kernFn = sincKernel;
        ts = sort(tSeq(fWeight~=0)); 
        if ( max(diff([ts,ts(1)+tPeriod])) > kSigma ),
            warning( 'Large gaps in input data: results using sinc kernel may be inaccurate' ),
        end
    case 'gauss',
        kernFn = gaussKernel;
    case 'box',
        kernFn = rectWin;
    otherwise
        error( 'Invalid kernel method specified: %s', kernelMethod )
end
   
% Windowing Fuction
switch windowMethod,
    case 'sinc',
        winFn = sincWin;
    case 'box',
        winFn = rectWin;        
    case 'zerocross',
        winFn = truncZcWin;
    case 'hann', 
        winFn = hannPctWin;
    otherwise,
        error( 'Invalid window method specified: %s', windowMethod )
end        

% Kernel to Use 
kernel = @(t,t0,dt,w) winFn(t,t0,dt,w) .* kernFn(t,t0,dt,w);


%% Show Kernel

if ( isVerbose ),
    tPlot0 = linspace( min(tSeqPad(:)), max(tSeqPad(:)), 10000 );
    if ( tPeriod == 0 )
        tPlot1 = tPlot0;
    else
        tPlot1 = tPlot0/tPeriod;
    end
    figure( 'Name', sprintf( 'kernel_%s_window_%s', kernelMethod, windowMethod ) ),
    plot( tPlot1, kernel( tPlot0, tPlot0(round(numel(tPlot0)/2)), kSigma, kWidth ), 'LineWidth', 2 )
    hold on
    plot( tPlot1, kernFn( tPlot0, tPlot0(round(numel(tPlot0)/2)), kSigma, kWidth ), '--', 'LineWidth', 1.5 )
    plot( tPlot1, winFn(  tPlot0, tPlot0(round(numel(tPlot0)/2)), kSigma, kWidth ), '--', 'LineWidth', 1.5 )
    hold off
    if ( tPeriod == 0 )
        xlabel( 'time' )
    else
        xlabel( 'time (R-R interval)' )
    end
    set(gca,'yTick',0)
    grid on
    legend( 'windowed kernel', sprintf( 'kernel (%s)', kernelMethod ), sprintf( 'window (%s)', windowMethod ) )
end


%% Anon Helper Functions

normalise_time   = @(x) bsxfun( @rdivide, x, sum( x, 3 ) );
reshape_time_vec = @(x) reshape( x, 1, 1, [] );


%% Kernel 'Smoothing'

imOut = zeros( size(imSeq,1), size(imSeq,2), numel(tOut), 'like', imSeq );

% frame weight
f = reshape_time_vec( fWeight( fSeqPad ) );
    
%voxel weight
v = vWeight(:,:,fSeqPad);

for iF = 1:numel(tOut), 
    
    % kernel
    k = reshape_time_vec( kernel( tSeqPad, tOut(iF), kSigma, kWidth ) );
    
    % frames to exclude
    x = reshape_time_vec( fSeqPad ~= fExclude{iF} );
   
    % combine all weighting
    frWgt = k .* f .* x;
    vxWgt = normalise_time( bsxfun( @times, frWgt, v ) );
    
    % calculate current frame 
    imOut(:,:,iF) = sum( imSeq(:,:,fSeqPad) .* vxWgt, 3 );
    
end


%% Corrections

imOut( isnan( imOut ) ) = 0;  % NOTE: this should correct for cases where vxWgt = 0/0  
                              % TODO: verify replacing nan with zero is
                              % correct, e.g., could/should use nansum in weighted
                              % averaging insteadif imSeq has NaN values

% NOTE: correction results in errors for complex image data
%
% imOut( imOut > max(imSeq(:)) ) = max(imSeq(:));
% imOut( imOut < min(imSeq(:)) ) = min(imSeq(:));


end  % imseq_kernel_smooth(...)