function [ xfRcn, safetyMargin, xfFlt ] = recon_xfsense( xfAcq, xfPri, psf, csm, psi, varargin )
%RECON_XFSENSE  k-t SENSE dynamic MRI reconstruction from x-f data
%
%   xfRcn = RECON_XFSENSE( xfAcq, xfPri, psf, csm, psi );
%
%   [ xfRcn, safetyMargin, xfFlt ] = RECON_XFSENSE(...);
%
%   input:
%       xfAcq:  acquired x-f data;                                  x-y-f-c
%       xfPri:  x-f prior;                                          x-y-f-1
%       psf:    point spread function of aliases in x-f; x-y-f-1 or 1-y-f-1
%       csm:    coil sensitivity matrices;                          x-y-1-c
%       psi:    noise covariance matrix;                                c-c
%   optional name-value pair:
%       safetyMargin:       safety margin, scalar; 
%                           default 2
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
%   See also recon_ktsense, recon_xtsense

% jfpva (joshua.vanamerom@kcl.ac.uk)

% based on code written by Shaihan Malik for ESMRB workshop on k-t methods


%% Data Dimensions

dimX = 1;  nX = size( xfAcq, dimX );
dimY = 2;  nY = size( xfAcq, dimY );
dimF = 3;  nF = size( xfAcq, dimF );
dimC = 4;  nC = size( xfAcq, dimC );


%% Parse Input

default.safetyMargin        = 2;
default.xToRecon            = 1:nX;
default.makeAdaptiveFilter  = false;
default.isVerbose           = false;

p = inputParser;

addRequired(  p, 'xfAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nF,nC]}, mfilename ) );

addRequired(  p, 'xfPri', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nF,1]}, mfilename ) );

addRequired(  p, 'psf', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,nY,nF]}, mfilename ) );
    
addRequired(  p, 'csm', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,1,nC]}, mfilename ) );

addRequired(  p, 'psi', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );

addParameter( p, 'safetyMargin', default.safetyMargin, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );

addParameter( p, 'xToRecon', default.xToRecon, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','vector','>=',0,'<=',nX}, mfilename ) );    

addParameter( p, 'makeAdaptiveFilter', default.makeAdaptiveFilter, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );
    
addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p,  xfAcq, xfPri, psf, csm, psi, varargin{:} );

safetyMargin        = p.Results.safetyMargin;
xToRecon            = p.Results.xToRecon;
makeAdaptiveFilter  = p.Results.makeAdaptiveFilter;
isVerbose           = p.Results.verbose;


%% Validate Input

% safetyMargin

if isempty( safetyMargin ),
    safetyMargin = default.safetyMargin;
end

% xToRecon

if isempty( xToRecon ),
   xToRecon = default.xToRecon; 
end

% makeAdaptiveFilter

xfFlt = []; 
if ( makeAdaptiveFilter ),
    fprintf( 'Constructing adaptive filter, reconstruction may take ~5x longer.\n' ),
end


%% Setup

if ( isVerbose ),
    fprintf( '\n%s()\n\n', mfilename );
end

% Suppress warnings for inv()

origWarningState = warning;
warning( 'off', 'MATLAB:nearlySingularMatrix' );
resetWarnings = onCleanup( @() warning( origWarningState ) );


%% Generate Aliasing Pattern and Calculate Acceleration Factor

% Point spread function

xfPsf   = abs( psf ) > 0.7 * max( abs( psf(:) ) );
nA      = sum( xfPsf(:) );   % number of aliases

% x-f Aliasing Pattern

[sX,sY,sF] = ind2sub( size(xfPsf), find(xfPsf) );

xShift = mod( sX - sX(1) + 1, nX ) - 1;
yShift = mod( sY - sY(1) + 1, nY ) - 1;
fShift = mod( sF - sF(1) + 1, nF ) - 1;

% Alias signs

aliasSign = zeros(nA,1);
for iA = 1:nA
    aliasSign(iA) = sign( real( psf(sX(iA),sY(iA),sF(iA)) ) );
end

%% Calculate x-f Prior Scaling

priScaleFactor = safetyMargin;  % alternatively, could specify regStrength, with safetyMargin = 1 / ( regStrength^2 ); 


%% Initialise Output

xfRcn = complex(zeros(nX,nY,nF));    

if ( makeAdaptiveFilter ),
    xfFlt      = complex(zeros(nX,nY,nF,nC));        % adaptive filter
    xfFltAlias = complex(zeros(nA,nC,nX,nY/nA,nF));  % adaptive filter in aliased space
end

%% x-f Unaliasing

%{

Each unique point in x-f (x,y,f) = ( iX=1:nX, iY=1:(nY/acFactor), iF=1:nF ), 
has aliases at locations, iA, (x,y,f) = ( aX, aY, aF ).

Build unaliased x-f images, xfRcn, from point-by-point reconstructions, 
rhoRcn, i.e., for x-y-f point (iX,iY,iF), rhoRcn = xfRcn(aX,aY,aF). 

For input data,
    xfAcq   acquired, undersampled data
    xfPri   training data prior

Assemble the following at each unique point in x-f space,
    rhoAls  xfAcq(iX,iY,iF,:)
    rhoEst  xfPri(aX,aY,aF) 
    M       abs( diag( rhoEst ) ) * scaleFactor
    theta   M.^2
    S       csm(aX,aY,1,:)

And solve,

    rhoRcn = theta*S' * inv( S*theta*S' + psi ) * rhoAls;

%}


for iY = 1:(nY/nA),  % TODO: could replace iX-iY loops with iR loop and change to parfor
    
    for iX = xToRecon,
        
        % Get indices of aliased locations for current iX-iY
        aX = 1 + mod( iX-xShift(:)-1, nX );
        aY = 1 + mod( iY-yShift(:)-1, nY );
        
        % Get coil sensitivity data for current iX-iY
        S = zeros( nC, nA );
        for iC = 1:nC,  
            for iA = 1:nA,
                S(iC,iA) = csm( aX(iA), aY(iA), 1, iC );
            end
        end
        
        % Loop over frequency
        for iF = 1:nF,
            
            % Get indices of aliased locations for current iX-iY-iF
            aF = 1 + mod( iF-fShift(:)-1, nF ); 
            
            % Convert to linear indices of aliased aX-aY-aF points 
            aL = aX + (aY-1)*nX + (aF-1)*(nX*nY);  % NOTE: aL = sub2ind([nX,nY,nC],aX,aY,aF) is slower
            
            % rhoAls, scaled
            rhoAcq = nA * reshape( xfAcq(iX,iY,iF,:), nC, 1 );
           
            % rhoEst, M
            rhoPri = xfPri( aL );
            M      = priScaleFactor * abs( diag( rhoPri ) ); 

            % Create adaptive filter for current alias problem
            F0 = M.^2*S' * inv( S*M.^2*S' + psi );  
            
            % Adjust sign for each alias
            F = bsxfun( @times, aliasSign, F0 );

            % Solve
            rhoRcn = F * rhoAcq; 
            
            % Insert in reconstructed x-f array
            xfRcn(aL) = rhoRcn(:);
            
            % Insert in adaptive filter array
            if ( makeAdaptiveFilter ),
                xfFltAlias(:,:,iX,iY,iF) = F;  % NOTE: fastest to put into 5D array and rearrange later
            end
                
        end
        
    end
    
end


%% Rearrange Adaptive Filter

%{ 

fastest to put adaptive filter values into 5D array during reconstruction 
and rearrange later;

vectorised indexing isn't actually faster, 
    e.g., 
      aLxyfc = reshape(bsxfun(@plus,aL,((1:nC)-1)*(nX*nY*nF)),[],1);   
    xfFlt(aLxyfc) = F(:);
    
nor is looping over aliases,
    e.g., 
      for iC = 1:nC,
          for iA = 1:nA,
              xfFlt(aX(iA),aY(iA),aF(iA),iC) = F(iA,iC);  
          end
      end

%}


if ( makeAdaptiveFilter ),
    [~,oY] = sort(aY);
    for iA = 1:nA,
        iY = (1:(nY/nA))+(nY/nA)*(oY(iA)-1);
        xfFlt(:,iY,:,:) = circshift( permute( xfFltAlias(iA,:,:,:,:), [3,4,5,2,1] ), -fShift(iA), dimF );
    end
end


end  % recon_xfsense(...)