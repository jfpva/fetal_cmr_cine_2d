function [ xtAcq, xtTrn, xtSmp, psi, xtBln, xtDff ] = prep_ktsense( ktAcq, ktTrn, varargin )
%PREP_KTSENSE  prep k-t dynamic MRI data for k-t SENSE reconstruction
%
%   [ xtAcq, xtTrn, xtSmp, psi ] = PREP_KTSENSE( ktAcq, ktTrn )
%
%   [ xtAcq, xtTrn, xtSmp, psi, xtBln, xtDff ] = PREP_KTSENSE( ... )
%
%   PREP_KTSENSE( ..., 'ktNoise', ktNoise )
%
%   PREP_KTSENSE( ..., 'removeoversampling', false )
%
%   PREP_KTSENSE( ..., 'fn_rmv_os', anon_fn_rmv_os )
%
%   PREP_KTSENSE( ..., 'isVerbose', false )
%
%   See also recon_ktsense

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Data Dimensions

dimX = 1;  nX = size( ktAcq, dimX );
dimY = 2;  nY = size( ktAcq, dimY );
dimT = 3;  nT = size( ktAcq, dimT );
dimC = 4;  nC = size( ktAcq, dimC );


%% Parse Input

default.ktNoise             = [];
default.removeoversampling  = false;
default.fn_rmv_os           = @( x ) x((round(size(x,1)/4)+1):(size(x,1)-round(size(x,1)/4)),:,:,:);
default.isVerbose           = false;

p = inputParser;

addRequired(  p, 'ktAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'ktTrn', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addParameter( p, 'ktNoise', default.ktNoise, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,NaN,NaN,nC]}, mfilename ) );
    
addParameter( p, 'removeoversampling', default.removeoversampling, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

addParameter( p, 'fn_rmv_os', default.fn_rmv_os, ...
        @(x) validateattributes( x, {'function_handle'}, ...
        {}, mfilename ) );

addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, ktAcq, ktTrn,varargin{:} );

ktNoise             = p.Results.ktNoise;
isRmvOversampling   = p.Results.removeoversampling;
fn_rmv_os           = p.Results.fn_rmv_os;
isVerbose           = p.Results.verbose;


%% Setup

if ( isVerbose ),
    fprintf( '\n%s()\n\n', mfilename );
end


%% Anon Fns

kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, 1 ), 2 ) );

phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );


%% Phase Correction

ktAcq = phase_correct( ktAcq );

ktTrn = phase_correct( ktTrn );


%% Sampling Pattern

ktSmp = single( sum( sum( ktAcq, 4 ), 1 ) ~= 0 );


%% Separate Baseline and Difference

ktBln = bsxfun( @rdivide, sum( ktAcq, 3 ), sum( ktSmp, 3 ) );

ktDff = ktAcq - bsxfun( @times, ktBln, ktSmp );


%% Transform to x-t Space

xtTrn = kt2xt( ktTrn );

xtAcq = kt2xt( ktAcq );

xtBln = kt2xt( ktBln );

xtDff = kt2xt( ktDff );

xtSmp = kt2xt( ktSmp ); 


%% Estimate Noise Covariance

if ( isempty( ktNoise ) ),
    noiseData = reshape( xtAcq([1,nX],:,:,:), [], nC ) * sqrt( nX * nY * nT );
else
    noiseData = reshape( ktNoise, [], nC );
end

psi = cov( noiseData );


%% Remove Oversampling

if ( isRmvOversampling ),
   
    xtTrn = fn_rmv_os( xtTrn );
    xtAcq = fn_rmv_os( xtAcq );
    xtBln = fn_rmv_os( xtBln );
    xtDff = fn_rmv_os( xtDff );
    
end


end   % prep_ktsense(...)