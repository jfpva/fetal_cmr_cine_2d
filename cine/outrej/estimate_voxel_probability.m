function [ voxProb, P ] = estimate_voxel_probability( img, est, varargin )
%ESTIMATE_VOXEL_PROBABILITY  estimate voxel-wise posterior probabilities
%for voxels in image sequence
%
%   [ voxProb, PARAM ] = ESTIMATE_VOXEL_PROBABILITY( img, est ), computes 
%   complex error between voxels in input image sequence, img, and matched 
%   reference data, est, and returns voxel-wise posterior probabilities, 
%   voxProb, based on maximum likelihood estimate of inlier-outlier mixture
%   model parameters, PARAM.
%
%   In-/Out-lier Class Definintions
%       Error:         e_jk         = Img_jk - Est_jk
%       Inlier Class:  P_{in}(e)    = Gaussian_{s1,s2}(e)
%       Outlier Class: P_{out}(e)   = Uniform_{m}(e)
%       Likelihood:    P(error|s1,s2,c) = c*P_{in}(e) + (1-c)*P_{out}(e)
%
%
%       img             array of images                             [x,y,t]
%       est             array of estimated image priors             [x,y,t]
% 
%   Additional name-value input: 
%
%       mask            logical image mask; default all voxels included
%       pctexc          percentile of data to exclude from parameter 
%                       initalisation (i.e., remove < pctexc/2 th percentile 
%                       and > 100-pctexc/2 th precentile; default zero
%       verbose         verbosity; default false
%
%   See also: estimate_frame_probability

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Parse Input

default.mask          = true(size(img)); 
default.pctexc        = 0;
default.isVerbose     = false;

p = inputParser;

addRequired(  p, 'img', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );

addRequired(  p, 'est', @(x) validateattributes( x, {'numeric'}, ...
        {'size',size(img)}, mfilename ) );
    
addParameter( p, 'mask', default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(img(:,:,1)),NaN]}, mfilename ) );

addParameter( p, 'pctexc', default.pctexc, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'scalar','>=',0,'<',100}, mfilename ) );
    
addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, img, est, varargin{:} );

mask        = p.Results.mask;
pctexc      = p.Results.pctexc;
isVerbose   = p.Results.verbose;


%% Init

if isreal( img ),
    warning( 'expected complex images' ),
    img = complex( img, img );
end

if isreal( est ),
    warning( 'expected complex images' ),
    est = complex( est, est );
end


%% Setup

% Anonymous function to extract values in masked regions

get_mask_values = @( x, msk ) double( x( repmat( msk, [1,1,size(x,3)/size(msk,3)] ) ) );


% Error maps

err  = img - est;  


% Error vector

e    = get_mask_values( err, mask );


%% In-/Out-lier Estimation

% In-/Out-lier Class Definintions
%   Inlier Class:  P_{in}(e)    = G_{s1,s2}(e)
%   Outlier Class: P_{out}(e)   = U_{m}(e)
%   Likelihood:    P(e|s1,s2,c) = c*P_{in}(e) + (1-c)*P_{out}(e)


% Prepare Error Data for Parameter Initialisation
% separate real and imaginary, exclude extreme values

pctLo   = pctexc / 2;
pctHi   = 100 - pctLo;

eReal = real( e );  
eReal = eReal( ( eReal >= prctile(eReal,pctLo) ) & ( eReal <= prctile(eReal,pctHi ) ) );
eImag = imag( e );  
eImag = eImag( ( eImag >= prctile(eImag,pctLo) ) & ( eImag <= prctile(eImag,pctHi ) ) );


% Inlier 

u1 = 0;             % mean of inlier class real pdf
s1 = std(eReal);    % standard deviation of inlier class real pdf
u2 = 0;             % mean of inlier class imag pdf
s2 = std(eImag);    % standard deviation of inlier class imag pdf


% Outlier

r1 = max(abs(eReal));       % maximum (radius) of real error values
r2 = max(abs(eImag));       % maximum (radius) of imag error values
m  = 1 / ( pi * r1 * r2 );  % density of uniform outlier pdf
r = @(x) (r1*r2) ./ sqrt( r1^2*sin(angle(x)).^2 + r2^2*cos(angle(x)).^2 ); % radius at angle(x) 


% Mix

c  = 0.9;


% PDF Definitions

vPdfIn       = @(x,s1,s2,c) reshape( ( c ) * mvnpdf([real(x(:)),imag(x(:))],[u1,u2],[s1,s2].^2), size(x) );  
vPdfOut      = @(x,c)       (1-c) * m * ones( size( x ) );
vPdfMix      = @(x,s1,s2,c) ( vPdfIn(x,s1,s2,c) + vPdfOut(x,c) );

vPdfMask     = @(x)         ( abs(x) <= r(x) );
vPdfOutMask  = @(x,c)       ( vPdfOut(x,c) .* vPdfMask(x) );
vPdfMixMask  = @(x,s1,s2,c) ( vPdfIn(x,s1,s2,c) + vPdfOutMask(x,c) );
vWeight      = @(x,s1,s2,c) ( vPdfIn(x,s1,s2,c) ./ vPdfMix(x,s1,s2,c) ) .* vPdfMask(x); 


% Log Likelihood Definitions

negLogLhd    = @(p,x)      -sum( log( vPdfMix(x,p(1),p(2),p(3)) ) );         % negative log-likelihood
conNegLogLld = @(p,x,lb,ub) max( negLogLhd(p,x), inf*(any(p<lb|p>ub)-0.5) ); % constrained negative log-likelihood


% PDF Parameters

paramStart = [  s1,  s2,   c ]; 
paramLower = [   0,   0, 0.0 ];  
paramUpper = [ inf, inf, 1.0 ];  


% Search Options

opts = optimset( statset('mlecustom'), 'GradObj','off', 'FunValCheck','off');

% if ( isVerbose ),
%     fprintf( 'Iter.      s1       s2       c\n' ),
%     outputFn = @(p,optimValues,state) fprintf('%-5i %8.4f %8.4f %8.4f\n',optimValues.iteration,p(1),p(2),p(3)) < 0;
%     opts = optimset(opts, 'Display', 'on', 'OutputFcn', outputFn);
% end


% Maximum Likelihood Estimation

paramHat = fminsearch( @(p) conNegLogLld(p,e,paramLower,paramUpper), paramStart, opts );

s1 = paramHat(1);
s2 = paramHat(2);
c  = paramHat(3);


%% Results 

if ( isVerbose ),

    
w = vWeight(e,s1,s2,c);
nVoxIn  = sum(sum(sum(w>0.5)));
nVoxOut = sum(sum(sum(w<=0.5)));
nVoxTot = numel(w);
pVoxIn  = 100 * nVoxIn  / nVoxTot;
pVoxOut = 100 * nVoxOut / nVoxTot;
   
fprintf( '\n\t\testimate_voxel_probability\n' ),
fprintf( '\t\t\tinlier class:  c * mvnpdf([real(e),imag(e)],[0,0],[s1,s2]), s1 = %6.4f, s2 = %6.4f, %7i/%-7i voxels (%.2f%%)\n', s1, s2, nVoxIn, nVoxTot, pVoxIn )
fprintf( '\t\t\toutlier class: (1-c) * m,                                    m = %6.4f,  c = %6.4f, %7i/%-7i voxels (%.2f%%)\n', m, c, nVoxOut, nVoxTot, pVoxOut )

% use histcount to determine bin spacing

[ ~, eBinEdgeReal0 ] = histcounts( real(e), 'BinMethod', 'auto', 'Normalization', 'pdf' );   
[ ~, eBinEdgeImag0 ] = histcounts( imag(e), 'BinMethod', 'auto', 'Normalization', 'pdf' );   


% adjust bin edges to be evenly distributed around zero

eBinEdgeMax         = max(abs([eBinEdgeReal0,eBinEdgeImag0]));
eBinWidthRescaleF   = 1.7;
eBinWidthReal       = eBinWidthRescaleF * mean(diff(eBinEdgeReal0));
eBinWidthImag       = eBinWidthRescaleF * mean(diff(eBinEdgeImag0));
eBinEdgeFn          = @( w, m ) [ -flip((0:w:m)+w/2), (0:w:m)+w/2 ];
eBinEdgeReal        = eBinEdgeFn( eBinWidthReal, eBinEdgeMax );  
eBinEdgeImag        = eBinEdgeFn( eBinWidthImag, eBinEdgeMax );


% calculate bin centres

eBinCentreReal  = eBinEdgeReal(2:end) - eBinWidthReal/2;   
eBinCentreImag  = eBinEdgeImag(2:end) - eBinWidthImag/2;   
[eGridReal,eGridImag] = meshgrid( eBinCentreReal, eBinCentreImag ); 
eGrid           = complex( eGridReal, eGridImag );


% 2d histogram

eBinCount       = permute( hist3( [real(e),imag(e)], { eBinCentreReal, eBinCentreImag } ), [2,1]);  
eBinArea        = mean(diff(eBinCentreReal)) * mean(diff(eBinCentreImag));
eBinDensity     = ( eBinCount / sum(eBinCount(:)) ) / eBinArea;

pInBinDensity   = vPdfIn(  eGrid, s1, s2, c );
pOutBinDensity  = vPdfOutMask( eGrid, c );
pMixBinDensity  = vPdfMixMask( eGrid, s1, s2, c );
wBinDensity     = vWeight( eGrid, s1, s2, c );


% Visualise Error Distribution and Maximum Likelihood Estimate

show_hist2d = @( v ) imshow(v,[],'XData',eBinCentreReal,'YData',eBinCentreImag,'InitialMagnification',500,'Colormap',parula);

figure( 'Name', 'voxel_error_distribution', 'Position', [279 1 724 705] ),

vOffset = 0.01;

subplot(4,4,[1,2,5,6]),
hAx1 = gca;
show_hist2d( eBinDensity ) 
axis on xy,
xlabel('\Re(error)'),
ylabel('\Im(error)'),
hCb = colorbar('Location','east','Color',0.85*[1,1,1]);
hCb.Label.String = 'density';
hCb.Position(3) = hCb.Position(3)/2;
hCb.Position(1) = hCb.Position(1) + hCb.Position(3);
title('Error Distribution')
hAx1.YTick = hAx1.XTick;

subplot(4,4,[3,4]),
hAx2 = gca;
hB = bar( eBinCentreImag, sum(eBinDensity,2) );
hB.FaceColor = [ 0.5 0.5 1.0 ];
ylabel('density'),
title('Error Distribution'),
hAx2.YAxisLocation = 'right';
hAx2.XLim = hAx1.YLim;

subplot(4,4,[7,8])
hAx3 = gca;
hold on
plot(eBinCentreImag,sum(pMixBinDensity,2),'LineWidth',3)
plot(eBinCentreImag,sum(pInBinDensity,2),'g--','LineWidth',2)
plot(eBinCentreImag,sum(pOutBinDensity,2),'r--','LineWidth',2)
hold off
box on
xlabel('\Im(error)'),
ylabel('density'),
legend('mix','in','out')
title('Maximum Likelihood Estimate'),
hAx3.YAxisLocation = 'right';
hAx3.XLim = hAx1.YLim;

subplot(4,4,[9,10])
hAx4 = gca;
hB = bar( eBinCentreReal, sum(eBinDensity,1) );
hB.FaceColor = [ 0.5 0.5 1.0 ];
ylabel('density'),
title('Error Distribution'),
hAx4.Position(2) = hAx4.Position(2)-vOffset;
hAx4.XLim = hAx1.XLim;

subplot(4,4,[13,14])
hAx5 = gca;
hold on
plot(eBinCentreReal,sum(pMixBinDensity,1),'LineWidth',3)
plot(eBinCentreReal,sum(pInBinDensity,1),'g--','LineWidth',2)
plot(eBinCentreReal,sum(pOutBinDensity,1),'r--','LineWidth',2)
hold off
box on
xlabel('\Re(error)'),
ylabel('density'),
legend('mix','in','out')
title('Maximum Likelihood Estimate'),
hAx5.Position(2) = hAx5.Position(2)-vOffset;
hAx5.XLim = hAx1.XLim;

subplot(4,4,[11,12,15,16])
hAx6 = gca;
show_hist2d( vPdfMix(eGrid,s1,s2,c) ),
axis on xy,
xlabel('\Re(error)'),
ylabel('\Im(error)'),
hAx6.YAxisLocation = 'right';
hAx6.Position(2) = hAx6.Position(2)-vOffset;
hCb = colorbar('Location','east','Color',0.85*[1,1,1]);
hCb.Label.String = 'density';
hCb.Position(3) = hCb.Position(3)/2;
hCb.Position(1) = hCb.Position(1) + hCb.Position(3);
title('Maximum Likelihood Estimate')
hAx6.YTick = hAx6.XTick;

dLim1D = max( [ hAx2.YLim; hAx3.YLim; hAx4.YLim; hAx5.YLim ] );
hAx2.YLim = dLim1D;
hAx3.YLim = dLim1D;
hAx4.YLim = dLim1D;
hAx5.YLim = dLim1D;

dLim2D = max( [ hAx1.CLim; hAx6.CLim ] );
hAx1.CLim = dLim2D;
hAx6.CLim = dLim2D;


% Visualise Posterior Probability

hFig2 = figure( 'Name', 'voxel_probability' );

show_hist2d( wBinDensity ),
hAx = gca;
axis on xy,
xlabel('\Re(error)'),
ylabel('\Im(error)'),
set(gca,'CLim',[0,1])
title('Posterior Probability')
colormap(pink),
hFig2.Position = [364 162 554 384];
hAx.YTick = hAx.XTick;
hCb = colorbar;
hCb.Label.String = 'probability';

end


%% Output

% Voxel Weights

voxProb = vWeight(err,s1,s2,c);

% Parameters

P.u1 = u1;  % PARAM.inlier.real.mean = u1;
P.s1 = s1;  % PARAM.inlier.real.stdv = s1;
P.u2 = u2;	% PARAM.inlier.imag.mean = u2;
P.s2 = s2;  % PARAM.inlier.imag.stdv = s2;
P.m  = m;   % PARAM.outlier.density  = m;
P.c  = c;   % PARAM.mix.c            = c; 


end  % estimate_voxel_probability(...)
