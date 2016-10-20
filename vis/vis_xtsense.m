function H = vis_xtsense( xfRcn, xfBln, xfDff, xfTrn, xfMask, xfFlt, psi, csm, mask, dt, safetyMargin, alpha, unmixMap )

%% Initialise

H = cell(0);

%% Dimensions

dimX = 1;  nX = size( xfDff, dimX );
dimY = 2;  nY = size( xfDff, dimY );
dimT = 3;  nT = size( xfDff, dimT );  
dimF = dimT; nF = nT;
dimC = 4;  nC = size( xfDff, dimC );

%% Pt

pt = bwmorph( mask, 'shrink', inf );
[iX,iY] = find(pt);
[~,iC]=max(abs(squeeze(csm(iX,iY,:))));


%% Colormaps

addpath( '~/Documents/MATLAB/Toolboxes/phasemap_v1.0/phasemap/' )

cMapSize = 1024; 

% Magnitude Colour Map (Log-Scale)
cMapMag = colormaplog( parula( cMapSize ) );
    %     cMapRef = parula( cMapSize );
    %     cMapMag = interp1( linspace(1,cMapSize,cMapSize), cMapRef, flip( cMapSize+1 - logspace(0, log10(cMapSize), cMapSize+1) ) );
    %     cMapMag = cMapMag(2:end,:);

% Angle Colour Map

cMapPha = phasemap( cMapSize );


%% Image Display

% Montage Multi-Channel Images

chshow = @( imCh ) montage( imCh, 'DisplayRange', [] );
chxfshow = @( imCh, iX ) montage( permute(imCh(iX,:,:,:),[2,3,1,4]), 'DisplayRange', [] );
%     chshow_range = @( imCh, imLim ) montage( imCh, 'DisplayRange', imLim );
%     chshow_prctile = @( imCh, pLo, pHi ) montage( imCh, 'DisplayRange', [prctile(imCh(:),pLo),prctile(imCh(:),pHi)] );

% Display x-f Space

xfRef = xfTrn; 
xfRef(:,:,(nF/2+1)+(-5:5)) = 0;
xfLim = [0,max(abs(xfRef(:)))];  
clear xfRef;
xfLimMask = [0,max(4,max(xfMask(:)))];
xfshow = @( xf, ix, xflim ) imshow(squeeze(abs(xf(ix,:,:))),xflim,'InitialMagnification',500,'Colormap',cMapMag);
xfcshow = @( xf, ix, xflim, ic ) xfshow( xf(:,:,:,ic), ix, xflim );

% x-f Threshold Limits and Maps

xfPriLim      = xfLim;
xfRssLim      = xfLim;
xfCutLim      = [0,50];
xfCutThresh   = 1/sqrt(2);
xfCutVal      = linspace(xfCutLim(1),xfCutLim(2),cMapSize);

% cMapCut = parula( cMapSize ); 
% cMapCut = flip( colormaplog( flip( parula( cMapSize ), 1 ) ), 1 ); 
cMapCut = colormaplog( parula( cMapSize ) ); 
cMapCut(xfCutVal<xfCutThresh,:) = repmat(permute(rgb2gray(permute(cMapCut(xfCutVal<xfCutThresh,:),[1,3,2])),[1,3,2]),[1,3]);

%% csm

imCsm = permute( csm, [2,1,3,4] );

% hFigCsmAbs
H{end+1} = figure( 'Name', 'csm_abs' );
chshow(abs(imCsm)), colormap(cMapMag), title('|csm|'), colorbar

% hFigCsmAng 
H{end+1} = figure( 'Name', 'csm_ang' );
chshow(angle(imCsm)), colormap(cMapPha), title('\anglecsm'), colorbar

%% psi

% hFigPsiAbs 
H{end+1}  = figure( 'Name', 'psi_abs' );
imgrid(abs(psi)),colorbar,colormap(cMapMag),title('|psi|')

% hFigPsiAng 
H{end+1}  = figure( 'Name', 'psi_ang' );
imgrid(angle(psi)),colorbar,colormap(cMapPha),title('\anglepsi')

%% x-f threshold

% coil combination

rss     = @( x ) sqrt(sum(real(x).^2+imag(x).^2,4));  
unmix   = @( x ) sum( bsxfun( @times, x, unmixMap ), 4 );

combine_coils = @( x ) unmix( x );

% calculate regularisation map

calc_regmap = @(sm0,smRoi,msk) (1/sm0^2)*(~msk)+(1/(smRoi))^2*(msk);

% x-f prior/support (including baseline)

xfPri = abs( xfTrn ).^2;

% noise variance vector

psiVec = reshape( diag( psi ), 1, 1, 1, [] );
psiMap = abs( unmix( psiVec ) );
psi1   = mean(psiMap(psiMap~=0));

% map of S \Theta S^H in k-t SENSE recon

pri = bsxfun( @times, csm, bsxfun( @times, xfPri , conj(csm) ) );

% map of \Lambda \Psi in k-t SENSE recon

reg0 = bsxfun( @times, calc_regmap(safetyMargin,safetyMargin,mask), psiVec ); 
reg1 = bsxfun( @times, calc_regmap(safetyMargin,safetyMargin*alpha,mask), psiVec ); 

% filter

f = xfPri ./ bsxfun( @plus, xfPri, calc_regmap(safetyMargin,safetyMargin,mask) .* psiMap );

% figures

% hFigXfPriAbs 
H{end+1} = figure( 'Name', 'xfPri_abs' );
xfshow(abs(xfPri),iX,xfPriLim), colormap(cMapMag), title(sprintf('x-f Prior\nx=%i',iX)), hCb = colorbar; hCb.Label.String = '|\Theta|'; xlabel('f'), ylabel('y'), 

% hFigPriAbs 
H{end+1} = figure( 'Name', 'pri_abs_separate_channels' );
chxfshow(abs(pri),iX), set(gca,'CLim',xfPriLim), colormap(cMapMag), title(strcat('coil-by-coil x-f support',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = '|S\ThetaS^H|'; xlabel('f'), ylabel('y'), 

% hFigRegAbs 
H{end+1} = figure( 'Name', 'reg_separate_channels' );
chxfshow(repmat(reg1,[1,1,nF,1]),iX), set(gca,'CLim',xfPriLim), colormap(cMapMag), title(strcat('spatially-varying regularisation term',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = '|\Lambda\Psi|'; xlabel('f'), ylabel('y'), 

% hFigxfCutUniReg 
H{end+1} = figure( 'Name', 'xfCut_uniformReg_separate_channels' );
chxfshow(abs(bsxfun(@rdivide,pri,reg0)),iX), set(gca,'CLim',xfCutLim), colormap(cMapCut), title(strcat('| x-f support / uniform reg. term |',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = '|S\ThetaS^H|\oslash|\lambda\Psi|'; hCb.Ticks = sort( [ xfCutThresh hCb.Ticks ] ); xlabel('f'), ylabel('y'), 

% hFigxfCutVarReg 
H{end+1} = figure( 'Name', 'xfCut_variableReg_separate_channels' );
chxfshow(abs(bsxfun(@rdivide,pri,reg1)),iX), set(gca,'CLim',xfCutLim), colormap(cMapCut), title(strcat('| x-f support / varying reg. term |',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = '|S\ThetaS^H|\oslash|\Lambda\Psi|'; hCb.Ticks = sort( [ xfCutThresh hCb.Ticks ] ); xlabel('f'), ylabel('y'), 

% hFigPriRss 
% H{end+1} = figure( 'Name', 'pri_rss' );
% xfshow(abs(combine_coils(pri)),iX,xfRssLim), colormap(cMapMag), title(strcat('RSS( x-f support )',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = 'RSS(S\ThetaS^H)'; xlabel('f'), ylabel('y')

% hFigRegRss 
% H{end+1} = figure( 'Name', 'reg_rss' );
% xfshow(abs(repmat(combine_coils(reg1),[1,1,nF,1])),iX,xfRssLim), colormap(cMapMag), title(strcat('RSS( spatially-varying regularisation term)',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = 'RSS(\Lambda\Psi)'; xlabel('f'), ylabel('y')

% hFigXfCutUniReg 
H{end+1} = figure( 'Name', 'xfCut_uniformReg' );
xfshow(combine_coils(bsxfun(@rdivide,pri,reg0)),iX,xfCutLim), colormap(cMapCut), title(strcat('x-f support / uniform reg. term' ,sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = 'S\ThetaS^H)\oslash(\lambda\Psi)'; hCb.Ticks = sort( [ xfCutThresh hCb.Ticks ] ); xlabel('f'), ylabel('y')

% hFigXfCutVarReg 
H{end+1} = figure( 'Name', 'xfCut_variableReg' );
xfshow(combine_coils(bsxfun(@rdivide,pri,reg1)),iX,xfCutLim), colormap(cMapCut), title(strcat('x-f support / varying reg. term',sprintf('\nx=%i',iX))), hCb = colorbar; hCb.Label.String = 'RSS(S\ThetaS^H)\oslashRSS(\Lambda\Psi)'; hCb.Ticks = sort( [ xfCutThresh hCb.Ticks ] ); xlabel('f'), ylabel('y')

%% adaptive filter

% hFigCsmAbsXf
H{end+1} = figure( 'Name', 'csm_abs_xf' );
chxfshow(abs(repmat(csm,[1,1,nF,1])),iX), colormap(cMapMag), title(strcat('|csm|',sprintf(' x=%i',iX))), colorbar, xlabel('f'), ylabel('y')

% hFigCsmAngXf
H{end+1} = figure( 'Name', 'csm_ang_xf' );
chxfshow(angle(repmat(csm,[1,1,nF,1])),iX), colormap(cMapPha), title(strcat('\anglecsm',sprintf(' x=%i',iX))), colorbar, xlabel('f'), ylabel('y')

% hFigFltAbs
H{end+1} = figure( 'Name', 'flt_abs' );
chxfshow(abs(xfFlt),iX), colormap(cMapMag), title(strcat('|xfFlt|',sprintf('x=%i',iX))), colorbar, xlabel('f'), ylabel('y')

% hFigFltAng
H{end+1} = figure( 'Name', 'flt_ang' );
chxfshow(angle(xfFlt),iX), colormap(cMapPha), title(strcat('\anglexfFlt',sprintf(' x=%i',iX))), colorbar, xlabel('f'), ylabel('y')

%% x-f training

% hFigXfTrnAbs
H{end+1} = figure( 'Name', 'xfTrn_abs' );
xfshow( xfTrn, iX, xfLim ), colorbar, title(strcat('|xfTrn|',sprintf('_{x=%i}',iX))), xlabel('f'), ylabel('y')

% hFigXfTrnAng
H{end+1} = figure( 'Name', 'xfTrn_ang' );
xfshow( angle(xfTrn), iX, pi*[-1,+1] ), colormap(cMapPha), colorbar, title(strcat('\anglexfTrn',sprintf('_{x=%i}',iX))), xlabel('f'), ylabel('y')

%% x-f mask

% hFigXfMask
H{end+1} = figure( 'Name', 'xfMask' );
xfshow( xfMask, iX, xfLimMask ), colorbar, title(strcat('xfMask',sprintf('_{x=%i}',iX))), xlabel('f'), ylabel('y')

%% x-f reconstructed

% hFigXfRcnAbs
H{end+1} = figure( 'Name', 'xfRcn_abs' );
xfshow( xfRcn, iX, xfLim ), colorbar, title(strcat('|xfRcn|',sprintf(' x=%i',iX))), xlabel('f'), ylabel('y')

% hFigxfRcnAng
H{end+1} = figure( 'Name', 'xfRcn_ang' );
xfshow( angle(xfRcn), iX, pi*[-1,+1] ), colormap(cMapPha), colorbar, title(strcat('\anglexfRcn',sprintf('_{x=%i}',iX))), xlabel('f'), ylabel('y')


%% x-f difference

% hFigXfDffAbs
H{end+1} = figure( 'Name', 'xfDff_abs' );
xfcshow( xfDff, iX, xfLim, iC ), colorbar, title(strcat('|xfDff|',sprintf('_{x=%i,c=%i}',iX,iC))), xlabel('f'), ylabel('y')

% hFigXfDffAbs
H{end+1} = figure( 'Name', 'xfDff_ang' );
xfcshow( angle(xfDff), iX, pi*[-1,+1], iC ), colormap(cMapPha), colorbar, title(strcat('|xfDff|',sprintf('_{x=%i,c=%i}',iX,iC))), xlabel('f'), ylabel('y')


%% Spatial Location of High-Frequency Content

[ ~, iFdc ] = max( abs( squeeze( sum( sum( sum( xfTrn, dimX ), dimY ), dimC ) ) ) );

fMax     = 1/(2*dt);        % Hz
df       = fMax/(nF/2);     % Hz
loFreq   = 1.8;             % Hz
nFLoFreq = ceil( loFreq / df );
iLoFreq  = mod( iFdc + (-(nFLoFreq-1):(nFLoFreq)) - 1, nF ) + 1;

xfHiF = xfTrn; 
xfHiF(:,:,iLoFreq) = 0;

imBln  = permute( abs(xfBln(:,:,iFdc)), [2,1] );
imHiF  = permute( max(abs(xfHiF),[],3), [2,1] );

% hFigXfBln
H{end+1} = figure( 'Name', 'xf_bln' );
imshow( imBln, [], 'InitialMagnification', 500 )
title('x-f baseline')
colorbar,
line([iX,iX],[1,size(imBln,2)],'Color',[0.5,0.5,0.5],'LineWidth',2)

% hFigXfHiF
H{end+1} = figure( 'Name', 'xy_hifreq' );
imshow( imHiF, xfLim, 'Colormap', cMapMag, 'InitialMagnification', 500 ), 
title( sprintf( 'max x-f signal in training data (f>%gHz)', loFreq ) )
colorbar
xlabel('x')
ylabel('y')

% Frequency of max signal in x-f...
%
%     fKernel = reshape(gausswin(5),[],1); fKernel = fKernel/sum(fKernel);
%     xfRef = xfTrn; 
%     xfRef(:,:,iFdc) = 0; 
%     iHiFreq = [(1:40),nF-(39:-1:0)]; 
%     xfRef(:,:,iHiFreq) = 0;
%     xfRef = convn( xfRef, fKernel, 'same' ); 
%     xfRef(:,:,iFdc) = 0;
% 
%     fMaxBW = abs(xfRef) == repmat( max(abs(xfRef),[],3), [1,1,nF] ); fMaxBW(abs(xfRef)==0) = false;
% 
%     for iX = 1:nX, for iY = 1:nY, f = find(fMaxBW(iX,iY,:)); if isempty(f), f=NaN; end, imMaxF(iX,iY) = f; end, end,
% 
%     hFigFMax = figure;
% 
%     imshow(df*abs(permute(imMaxF,[2,1])-iFdc),[0,df*nF/2],'Colormap',cMapM,'InitialMagnification',500),colorbar

%% Figure 2

H{end+1} = figure( 'Name', 'fig2a_xfBln' );
imshow( imBln, [], 'InitialMagnification', 500 )
line([iX,iX],[1,size(imBln,2)],'Color',[0.5,0.5,0.5],'LineWidth',2)

H{end+1} = figure( 'Name', 'fig2b_xfTrn' );
xfshow( xfTrn, iX, xfLim ), 

H{end+1} = figure( 'Name', 'fig2c_xfCut0' );
xfshow(bsxfun(@rdivide,combine_coils(pri),combine_coils(reg0)),iX,xfCutLim), colormap(cMapCut), 

H{end+1} = figure( 'Name', 'fig2d_xfCutRoi' );
xfshow(bsxfun(@rdivide,combine_coils(pri),combine_coils(reg1)),iX,xfCutLim), colormap(cMapCut), 

end