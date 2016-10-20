function [ xtRlt, safetyMargin, xfFlt, unmix ] = recon_xtsense( xtAcq, xtSmp, xtTrn, csm, psi, varargin )
%[ xtRlt, xfRcn, xfBlnXY, xfDff, xfPri, xfTrn, xfMask, psi, safetyMargin, xt2xf, xf2xt, xfFlt, unmix ] = recon_xtsense( xtAcq, xtSmp, xtTrn, csm, psi, varargin )
%RECON_XTSENSE  k-t SENSE dynamic MRI reconstruction from x-t data
% 
%   xtRcn = RECON_XTSENSE( xtBln, xtDff, xtPsf, xtTrn, csm, psi )
%
%   input:
%       xtAcq:  acquired undersampled complex image-space data;     x-y-t-c
%       xtSmp:  image-space transform of k-space sampling pattern;  x-y-t-1
%       xtTrn:  training complex image-space data;                  x-y-t-c
%       csm:    coil sensitivity matrices;                          x-y-1-c
%       psi:    noise covariance matrix;                                c-c
%   optional:
%       xtBln:          image-space baseline (DC) of acquired data; x-y-1-c
%                       default zeros
%       safetyMargin:   safety margin;                                    1
%       mask:           fetal heart ROI;                            x-y-1-1 
%                       default full FOV
%       additional mask-related inputs (alpha,beta,dt,loFreq,xfMaskKernelWidth)
%
%       xToRecon:           reconstruct subset of x values; 
%                           default all x
%       makeAdaptiveFilter: create adaptive filter; 
%                           default false
%       isVerbose:          verbose output; 
%                           default false
%
%   ouput:
%       xtRcn:  reconstructed complex images;                       x-y-t-1
%
%   See also recon_ktsense, recon_xfsense

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Data Dimensions

dimX = 1;  nX = size( xtAcq, dimX );
dimY = 2;  nY = size( xtAcq, dimY );
dimT = 3;  nT = size( xtAcq, dimT );  dimF = dimT;  nF = nT;
dimC = 4;  nC = size( xtAcq, dimC );


%% Parse Input

default.xtBln        = complex( zeros( nX, nY ) );
default.safetyMargin = [];

default.mask         = true( nX, nY );
default.alpha        = 4;
default.beta         = 1;
default.dt           = 0.077;  % seconds
default.loFreq       = 0.5;    % Hz
default.xfMaskKernelWidth = 1;

default.xToRecon     = 1:nX;
default.makeAdaptiveFilter = false;

default.isVerbose    = false;

p = inputParser;

addRequired(  p, 'xtAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'xtSmp', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,nY,nT,NaN]}, mfilename ) );
    
addRequired(  p, 'xtTrn', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,NaN]}, mfilename ) );
    
addRequired(  p, 'csm', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,1,nC]}, mfilename ) );

addRequired(  p, 'psi', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );

addParameter(  p, 'xtBln', default.xtBln, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,NaN,nC]}, mfilename ) );
    
addParameter( p, 'safetyMargin', default.safetyMargin, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
   
addParameter( p, 'mask',          default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

addParameter(  p, 'alpha', default.alpha, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );

addParameter(  p, 'beta', default.beta, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','nonnegative'}, mfilename ) );

addParameter(  p, 'dt', default.dt, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','positive'}, mfilename ) );

    addParameter( p, 'loFreq', default.loFreq, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
    
addParameter( p, 'xfMaskKernelWidth', default.xfMaskKernelWidth, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive'}, mfilename ) );   

addParameter( p, 'xToRecon', default.xToRecon, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','vector','>=',0,'<=',nX}, mfilename ) );    

addParameter( p, 'makeAdaptiveFilter', default.makeAdaptiveFilter, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );
    
addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, xtAcq, xtSmp, xtTrn, csm, psi, varargin{:} );

xtBln        = p.Results.xtBln;
safetyMargin = p.Results.safetyMargin;

mask         = p.Results.mask;
alpha        = p.Results.alpha;
beta         = p.Results.beta;
loFreq       = p.Results.loFreq;
dt           = p.Results.dt;
xfMaskKernelWidth = p.Results.xfMaskKernelWidth;

xToRecon     = p.Results.xToRecon;
makeAdaptiveFilter = p.Results.makeAdaptiveFilter;

isVerbose    = p.Results.verbose;


%% Validate

% dt

if ( dt > 0.2 ),
    warning( 'Large temporal resolution (dt) specified: %g seconds', dt ),
end


% xfMaskKernelWidth

if ( numel( xfMaskKernelWidth  ) == 1 ), 
    xfMaskKernelWidth = round( xfMaskKernelWidth ) * ones(1,3);
elseif ( numel( xfMaskKernelWidth ) == 2 ),
    xfMaskKernelWidth = [ reshape( round( xfMaskKernelWidth ), 1, 2 ), 1 ];
else
    warning( 'Using xfMaskKernelWidth(1:3)' )
    xfMaskKernelWidth = reshape( round( xfMaskKernelWidth(1:3) ), 1, 3 );
end


%% Setup

if ( isVerbose ),
    fprintf( '\n%s()\n\n', mfilename );
end

% Transformation Functions

xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );
xf2xt = @( xf ) ifft( ifftshift( xf, dimT ), [], dimT );


%% Transform to x-f Space 

xfAcq = xt2xf( xtAcq );  

xfSmp = xt2xf( xtSmp );

xfBln = xtBln * nF;  % NOTE: xt2xf not required since xtBln has singular time component, so scaled to match xt2xf of xtBln with size nF in time/freq.

xfTrn = xt2xf( xtTrn );


%% Identify DC Frequency

[ ~, iFdc ] = max( abs( squeeze( sum( sum( sum( xfTrn, dimX ), dimY ), dimC ) ) ) );

fMax     = 1/(2*dt);        % Hz
df       = fMax/(nF/2);     % Hz
nFLoFreq = ceil( min( fMax, loFreq ) / df );
iLoFreq  = mod( iFdc + (-(nFLoFreq-1):(nFLoFreq)) - 1, nF ) + 1;


%% Combine Channels in Training and Baseline Data

if ( size( xfBln, dimC ) > 1 || size( xfTrn, dimC ) > 1 ),
    
    % Calculate SENSE unfolding matrix for R=1 (for coil combination), 
    % using same approach as x-f unaliasing below
    unmix = complex(zeros(nX,nY,1,nC));
    [cX,cY] = find( sum( sum( csm, dimC ), dimT ) ~= 0 );
    for iR = 1:numel(cX);       
        S = reshape( csm( cX(iR), cY(iR), : ), nC, 1 );
        unmix(cX(iR),cY(iR),1,:) = inv( S'*inv(psi)*S ) * S' * inv(psi);
    end
    
    % Create coil combination anonymous function
    combine_coils = @( x ) sum( bsxfun( @times, x, unmix ), dimC );  
    
    % Combine channels in xfBln
    if ( size( xfBln, dimC ) > 1 ),
        xfBln = combine_coils( xfBln );
    end

    % Combine channels in xfTrn
    if ( size( xfTrn, dimC ) > 1 ),
        xfTrn = combine_coils( xfTrn );
    end

end


%% Create x-f Estimate

xfPri = xfTrn;

% Remove DC if Baseline Provided

if any( xfBln(:) ),  
    xfPri(:,:,iFdc,:) = 0;
end

% Create xfMask

xfMask = ones( nX, nY, nF );

if any( ~mask ),  % only if any part of mask is not 1
  
    xfMask = beta * xfMask;

    xfMask(:,:,iLoFreq) = 1;

    xfMask( repmat( mask, [1,1,nF] ) ) = alpha;
            
    % Smooth Transition from mask to background
    
    if ( any( xfMaskKernelWidth > 1 ) ),
    
        xfKernelX = gausswin( xfMaskKernelWidth(dimX) );  % TODO: compare to Tsao2003 frequency (Hanning) filter, e.g., xfKernel1D = hann( xfMaskKernelWidth )
        xfKernelY = gausswin( xfMaskKernelWidth(dimY) ); 
        xfKernelF = gausswin( xfMaskKernelWidth(dimF) ); 
        
        xfKernel = bsxfun( @times, xfKernelX * xfKernelY.', reshape( xfKernelF, 1,1,[] ) );
        xfKernelNorm = xfKernel / sum( xfKernel( : ) );
        xfMask = convn( padarray( xfMask, floor(xfMaskKernelWidth/2), 'replicate' ), xfKernelNorm, 'valid' );
    
    end
    
end

% Apply x-f Mask to x-f Estimate

xfPri = bsxfun( @times, xfPri, xfMask );


%% Zero-Pad Baseline Data in Frequency

if ( size( xfBln, dimF ) == 1 )
    xfBln = cat( 3, zeros(nX,nY,iFdc-1,1), xfBln, zeros(nX,nY,nF-iFdc,1) );
end


%% Recon

[ xfRcn, safetyMargin, xfFlt ] = recon_xfsense( xfAcq, xfPri, xfSmp, csm, psi, 'safetyMargin', safetyMargin, 'xToRecon', xToRecon, 'makeAdaptiveFilter', makeAdaptiveFilter, 'verbose', isVerbose );

if isempty( xfFlt ),
    xfFlt = zeros(nX,nY,nF,nC);
end


%% Combine Baseline and Reconstructed 

xfRlt = xfBln + xfRcn;


%% Transform x-f to x-t

xtRlt = xf2xt( xfRlt );


%% Visualise

if ( isVerbose ),
    H = vis_xtsense( xfRlt, xfBln, xfAcq, xfTrn, xfMask, xfFlt, psi, csm, mask, dt, safetyMargin, alpha, unmix );
    fprintf('Figures:  \n\n')
    for iH = 1:length(H),
        fprintf('![](figs/%s.png)  \n',H{iH}.Name)
    end
end


end  % recon_xtsense(...)