function [ frmProb, P ] = estimate_frame_probability( voxProb, varargin )
%ESTIMATE_FRAME_PROBABILITY  estimate frame-wise posterior probabilities
%given voxel-wise posterior probabilities in image sequence
%
%   [ frmProb, PARAM ] = ESTIMATE_FRAME_PROBABILITY( voxProb ), returns
%   frame-wise posterior probabilities for frames with voxel probabilities,
%   voxProb, based on maximum likelihood estimate of inlier-outlier mixture
%   model parameters, PARAM.
%
%   In-/Out-lier Class Definintions
%       Potential:     f_k          = sqrt( sum( voxProb_k.^2 / N_k ) )
%       Inlier Class:  P_{in}(f)    = Rician_{1-u,s}(1-f)
%       Outlier Class: P_{out}(f)   = Uniform_{m}(f)
%       Likelihood:    P(f|u,s,c)   = c*P_{in}(f) + (1-c)*P_{out}(f)
%
% 
%       voxProb         array of voxel-wise probabilities           [x,y,t]
% 
%   Additional name-value input: 
%
%       mask            logical image mask; default all voxels included
%       verbose         verbosity; default false
%
%   See also: estimate_voxel_probability

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Parse Input

default.mask          = true(size(voxProb)); 
default.isVerbose     = false;

p = inputParser;

addRequired(  p, 'voxProb', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );

addParameter( p, 'mask', default.mask, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(voxProb(:,:,1)),NaN]}, mfilename ) );

addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, voxProb, varargin{:} );

mask        = p.Results.mask;
isVerbose   = p.Results.verbose;


%% Init

% Mask
if size( mask, 3 ) < size( voxProb, 3 ),
    mask = repmat( mask, [1,1,size(voxProb,3)/size(mask,3)] );
end

% Frame Potential 
% f_k = sqrt( sum( voxProb_k.^2 / N_k ) ), for frame k with N_k voxels in mask
f = double(squeeze(sqrt(sum(sum(bsxfun(@times,mask,voxProb.^2)))./sum(sum(mask)))));


%% Identify Frame Outliers

% inlier class:  rician

fPdfIn  = @( x, u, s, c ) c * pdf('Rician', 1-x, 1-u, s );   


% outlier class: uniform

m = 1/range(f);

v = linspace(prctile(f,5),prctile(f,95),numel(f)*10);

fPdfOut = @( x, u, s, c ) (1-c) * m * ( ones( size( x ) ) - ((x>v(find(fPdfIn( v, u, s, 1 )==max(fPdfIn( v, u, s, 1 )),1,'last'))).*(max(fPdfIn( v, u, s, 1 ))-fPdfIn( x, u, s, 1 ))/max(fPdfIn( v, u, s, 1 ) ) ) );


% mixture

fPdfMix = @( x, u, s, c ) fPdfIn( x, u, s, c ) + fPdfOut( x, u, s, c );
fLogPdfMix = @( x, u, s, c ) log( fPdfIn( x, u, s, c ) + fPdfOut( x, u, s, c ) );


% frame weight (posterior probability)

fWeight = @( x, u, s, c ) fPdfIn( x, u, s, c ) ./ fPdfMix( x, u, s, c ); 


% parameter initial estimates and bounds

fC0  = 0.9;
fIn  = f(f>prctile(f,100*(1-fC0)));
paramStart = [ mean(fIn) std(fIn) fC0 ];
paramLower = [ 0   0     0   ];
paramUpper = [ 1   inf   1   ];


% maximum likelihood estimation

opt = optimset( 'MaxFunEvals', 800, 'MaxIter', 400 );

fParamEst = mle( f, 'logpdf', fLogPdfMix, 'start', paramStart, 'lowerbound', paramLower, 'upperbound', paramUpper, 'Options', opt );

u = fParamEst(1);
s = fParamEst(2);
c = fParamEst(3);


% frame-wise probability

frmProb = fWeight(f,u,s,c);


% in-/out-lier counts and percentages

nFrmIn  = sum(sum(sum(frmProb>0.5)));
nFrmOut = sum(sum(sum(frmProb<=0.5)));
nFrmTot = numel(frmProb);
pFrmIn  = 100 * nFrmIn  / nFrmTot;
pFrmOut = 100 * nFrmOut / nFrmTot;


% Verbose Output

if ( isVerbose ),

[ fBinDensity, fBinEdges ] = histcounts( f, 'BinMethod', 'auto', 'Normalization', 'pdf' );
if ( numel(fBinDensity) < sqrt(numel(f)) ),
    [ fBinDensity, fBinEdges ] = histcounts( f, 'BinMethod', 'fd', 'Normalization', 'pdf' );
end
fBinCentre = fBinEdges(2:end) - mean(diff(fBinEdges))/2;

figure('Name','frame_probability')
yyaxis left
hB = bar( fBinCentre, fBinDensity );
hB.FaceColor = [.5 .5 1 ];
xlabel( 'frame potential' )
x = linspace(min(f),max(f),10000);  
hold on
plot(x,fPdfMix(x,u,s,c),'LineWidth',2)
plot(x,fPdfIn(x,u,s,c),'g--','LineWidth',1)
plot(x,fPdfOut(x,u,s,c),'r--','LineWidth',1)
grid on
ylabel('density')
yyaxis right
plot(x,fWeight(x,u,s,c),'LineWidth',2)
ylabel('probability')
title( 'Frame-Wise Outlier Rejection' )
        
fprintf( '\n\t\testimate_frame_probability\n' ),
fprintf( '\t\t\tinlier class:  c * pdf(''Rician'', 1-x, 1-u, s ), u = %6.4f, s = %6.4f, %3i/%-3i frames (%.2f%%)\n', u, s, nFrmIn, nFrmTot, pFrmIn )
fprintf( '\t\t\toutlier class: (1-c) * m,                       m = %6.4f, c = %6.4f, %3i/%-3i frames (%.2f%%)\n', m, c, nFrmOut, nFrmTot, pFrmOut )

end


% Parameter Output

P = struct( 'u', u, 's', s, 'm', m, 'c', c );


%% OR... Identify Frame Outliers Using Mixture of Gaussians 
% inlier class:  gaussian
% outlier class: gaussian

% fPdfIn  = @( x, u, s, c ) c * normpdf( x, u, s );  
% fPdfIn  = @( x, u, s, c ) c * pdf('Rician', 1-x, 1-u, s );   
% fPdfOut = @( x, u, s, c ) (1-c) * normpdf( x, u, s );
% 
% fPdfMix = @( x, u1, s1, u2, s2, c ) fPdfIn( x, u1, s1, c ) + fPdfOut( x, u2, s2, c );
% fLogPdfMix = @( x, u1, s1, u2, s2, c ) log( fPdfIn( x, u1, s1, c ) + fPdfOut( x, u2, s2, c ) );
% 
% fWeight = @( x, u1, s1, u2, s2, c ) fPdfIn( x, u1, s1, c ) ./ fPdfMix( x, u1, s1, u2, s2, c );
% 
% % param         u1  s1    u2  s2     c
% fC0  = 0.9;
% fIn  = f(f>prctile(f,100*(1-fC0)));
% fOut = f(f<prctile(f,100*(1-fC0)));
% paramStart = [ mean(fIn) std(fIn) mean(fOut) std(fOut) fC0 ];
% paramLower = [ 0   0     0   0   0   ];
% paramUpper = [ 1   inf   1   inf 1   ];
% 
% opt = optimset( 'MaxFunEvals', 800, 'MaxIter', 400 );
% [fParamEst,fParamCI] = mle( f, 'logpdf', fLogPdfMix, 'start', paramStart, 'lowerbound', paramLower, 'upperbound', paramUpper, 'Options', opt );
% 
% u1 = fParamEst(1);
% s1 = fParamEst(2);
% u2 = fParamEst(3);
% s2 = fParamEst(4);
% c = fParamEst(5);
% 
% frmProb = fWeight(f,u1,s1,u2,s2,c);
% 
% nFrmOut = sum(sum(sum(frmProb<0.5)));
% nFrmTot = numel(frmProb);
% 
% 
% % Verbose Output
% 
% if ( isVerbose ),
%     
% [ fBinDensity, fBinEdges ] = histcounts( f, 'BinMethod', 'auto', 'Normalization', 'pdf' );
% if ( numel(fBinDensity) < sqrt(numel(f)) ),
%     [ fBinDensity, fBinEdges ] = histcounts( f, 'BinMethod', 'fd', 'Normalization', 'pdf' );
% end
% fBinCentre = fBinEdges(2:end) - mean(diff(fBinEdges))/2;
% 
% figure('Name','frame_probabilities')
% yyaxis left
% hB = bar( fBinCentre, fBinDensity );
% hB.FaceColor = [.5 .5 1 ];
% xlabel( 'frame potential' )
% x = linspace(min(f),max(f),10000);  % x = linspace(min(f),max(f),10000);
% hold on
% plot(x,fPdfMix(x,u1,s1,u2,s2,c),'LineWidth',2)
% plot(x,fPdfIn(x,u1,s1,c),'g--','LineWidth',1)
% plot(x,fPdfOut(x,u2,s2,c),'r--','LineWidth',1)
% grid on
% ylabel('density')
% yyaxis right
% plot(x,fWeight(x,u1,s1,u2,s2,c),'LineWidth',2)
% ylabel('probability')
% title( 'Frame-Wise Outlier Rejection' )
%         
% fprintf( 'frame outliers: u1 = %6.4f, s1 = %6.4f, u2 = %6.4f, s2 = %6.4f, c = %6.4f, pct outlier frames %.2f (%i/%i)\n', u1, s1, u2, s2, c, 100*nFrmOut/nFrmTot, nFrmOut, nFrmTot ) % fprintf( 'frame outliers: u1 = %7.4f (%7.4f,%7.4f), s1 = %7.4f (%7.4f,%7.4f), u2 = %7.4f (%7.4f,%7.4f), s2 = %7.4f (%7.4f,%7.4f), c = %7.4f (%7.4f,%7.4f), outlier frames (p_{k}<0.5) = %i/%i (%.2f%%)\n', u1, fParamCI(1,1), fParamCI(2,1), s1, fParamCI(1,2), fParamCI(2,2), u2, fParamCI(1,3), fParamCI(2,3), s2, fParamCI(1,4), fParamCI(2,4), c, fParamCI(1,5), fParamCI(2,5), nFrmOut, nFrmTot, 100*nFrmOut/nFrmTot )
% 
% end
% 
% 
% % Output
% 
% P = struct( 'u1', u1, 's1', s1, 'u2', u2, 's2', s2, 'c', c );


end  % estimate_frame_probability(...)