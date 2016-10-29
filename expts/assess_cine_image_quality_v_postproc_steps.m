function [ hFig, R, D ] = assess_cine_image_quality_v_postproc_steps( resultsDir )
%assess_cine_image_quality_v_postproc_steps
%
%   [ hFig, R, D ] = assess_cine_image_quality_v_postproc_steps( resultsDir );


dataSrcStr = 'mat';
initCardPhaseStr = 'zero';


%% Dependencies


script_add_dependencies


%% Identify Data

D = dir( resultsDir );

isInc = cat(1,D.isdir);  % include all directories

% loop to identify hidden directories
ind = find(isInc);
for iD = 1:numel(ind), 
   % on OSX, hidden directories start with a dot
   isInc(ind(iD)) = ~strcmp(D(ind(iD)).name(1),'.');
   if isInc(iD) && ispc
       % check for hidden Windows directories - only works on Windows
       [~,stats] = fileattrib(fullfile(resultsDir,D(iD).name));
       if stats.hidden
          isInc(iD) = false;
       end
   end
end

D = D(isInc);

idStrCell = { D.name };


%% Load Results


indFailed = false( 1, numel(idStrCell) );

for iD = 1:numel(idStrCell),
    
    try

        idStr = idStrCell{iD};

        D(iD).idStr = idStr;
        
        switch dataSrcStr
            
            case 'mat'
                               
                M = matfile( fullfile( resultsDir, idStr, 'cine', 'results.mat' ) );
                
                P = M.PARAM;
                
                D(iD).mask = M.maskQ;
                
                D(iD).rrInterval = P(end).R.rrInterval;
                
                D(iD).tRlt = M.tRlt;
                D(iD).imRlt = M.imRltQ;
                [ ~, indRlt2Card ] = sort( P(end).R.tRr );
                D(iD).tRdr = P(end).R.tRr(indRlt2Card);
                D(iD).imRdr = D(iD).imRlt(:,:,indRlt2Card);
                D(iD).tCine = M.tCine;
                
                D(iD).dtRlt = M.dtRlt;
                                
                n = size(D(iD).imRlt,3);
                rng('default');
                
                switch initCardPhaseStr
                    case 'zero',  % init with card phase = 0 for all
                        tRr = zeros(1,n);
                    case 'randi', % random cardiac phase 
                        tRr = (randi(n,1,n)-1)/n * D(iD).rrInterval;  
                    case 'randperm',  % randperm cardiac phase
                        tRr = randperm(n,n-1)/n * D(iD).rrInterval; 
                    case 'linear',  % linear cardiac phase
                        tRr = linspace(0,D(iD).rrInterval,n+1); tRr = tRr(1:n);   
                    otherwise 
                        warning( 'initial cardiac phase type (initCardPhaseStr = ''%s'') not an option\n', initCardPhaseStr )
                        indFailed(iD) = true;
                end
                D(iD).imCine0 = imseq_kernel_smooth( tRr, D(iD).imRlt, D(iD).tCine, 'tPeriod', D(iD).rrInterval, 'kSigma', D(iD).dtRlt, 'kWidth', D(iD).rrInterval );

                D(iD).imCineCardSync = P(1).R.imCine;
                
                MnoMoco = matfile( fullfile( resultsDir, idStr, 'cine_noMoco', 'results.mat' ) );
                PnoMoco = MnoMoco.PARAM;
                D(iD).imCineNoMoco = PnoMoco(end).P.imCine;
                
                MnoOutrej = matfile( fullfile( resultsDir, idStr, 'cine_noOutrej', 'results.mat' ) );
                PnoOutrej = MnoOutrej.PARAM;
                D(iD).imCineNoOutrej = PnoOutrej(end).P.imCine;
                
                D(iD).imCine = P(end).P.imCine;
                
            case 'nii'
                 
                 N = load_untouch_nii( fullfile( resultsDir, idStr, 'cine', 'maskq_fetalheart.nii.gz' ) );
                 D(iD).mask = logical( N.img );
                 
                 M = matfile( fullfile( resultsDir, idStr, 'cine', 'results.mat' ) );
                 P = M.PARAM;
                 
                 D(iD).rrInterval = P(end).R.rrInterval;
                 
                 pixdimAcq = M.pixdimAcq;
                 pixdimRcn = M.pixdimRcn;
                 xt2kt = @( xt ) fftshift( fftshift( fft2( xt ), 2 ), 1 );
                 kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, 2 ), 1 ) );
                 imRlt = M.imRlt;
                 acqDx = pixdimAcq(1);
                 acqDy = pixdimAcq(2);
                 rcnDx = pixdimRcn(1);
                 rcnDy = pixdimRcn(2);
                 acqDim = size( imRlt );
                 % fov = pixdimAcq .* acqDim(1:2);
                 rcnDim = [ acqDx/rcnDx*size(imRlt,1), acqDy/rcnDy*size(imRlt,2), acqDim(3) ];
                 padDim   = round( ( rcnDim - acqDim ) / 2 );
                 ktRlt    = xt2kt( imRlt );
                 filterMask = ones(size(ktRlt(:,:,1)));
                 ktRltQ = padarray( ktRlt, padDim );
                 scaleFactor = numel(ktRltQ(:,:,1))/sum(filterMask(:));
                 imRltQ = kt2xt( scaleFactor * ktRltQ );
                 D(iD).tRlt = M.tRlt;
                 D(iD).imRlt = imRltQ;
                 [ ~, indRlt2Card ] = sort( P(end).R.tRr );
                 D(iD).tRdr = P(end).R.tRr(indRlt2Card);
                 D(iD).imRdr = D(iD).imRlt(:,:,indRlt2Card);
                 D(iD).tCine = M.tCine;
                 % imCineMean = imseq_kernel_smooth( zeros(size(D(iD).tRlt)), imRltQ, 0 );
                 % D(iD).imCine0 = repmat( abs(imCineMean), 1, 1, numel(D(iD).tCine) );
                 dtRlt = mean(diff(D(iD).tRlt));
                 % n = size(imRltQ,3);
                 % rng('default');
                 % tRr = (randi(n,1,n)-1)/n * D(iD).rrInterval;
                 % tRr = (randperm(size(imRltQ,3),size(imRltQ,3))-1)/size(imRltQ,3) * D(iD).rrInterval;
                 % tRr = linspace(0,D(iD).rrInterval,n+1); tRr = tRr(1:n);
                 % D(iD).imCine0 = imseq_kernel_smooth( tRr, imRltQ, D(iD).tCine, 'tPeriod', D(iD).rrInterval, 'kSigma', dtRlt, 'kWidth', D(iD).rrInterval );
                 D(iD).imCine0 = imseq_kernel_smooth( zeros(size(D(iD).tRlt)), imRltQ, D(iD).tCine, 'tPeriod', D(iD).rrInterval, 'kSigma', dtRlt, 'kWidth', D(iD).rrInterval );
                 
                 N = load_untouch_nii( fullfile( resultsDir, idStr, 'cine', 'cineq_xyt_01a_cardsync.nii.gz' ) );
                 D(iD).tCine = M.tCine;
                 D(iD).imCineCardSync = N.img;
                 
                 N = load_untouch_nii( fullfile( resultsDir, idStr, 'cine_noMoco', 'cineq_xyt.nii.gz' ) );
                 D(iD).imCineNoMoco = N.img;
                 
                 N = load_untouch_nii( fullfile( resultsDir, idStr, 'cine_noOutrej', 'cineq_xyt.nii.gz' ) );
                 D(iD).imCineNoOutrej = N.img;
                 
                 N = load_untouch_nii( fullfile( resultsDir, idStr, 'cine', 'cineq_xyt.nii.gz' ) );
                 D(iD).imCine = N.img;
                
            otherwise
                
                warning( 'data source type (%s) not an option\n', dataSrcStr )
                indFailed(iD) = true;
                
        end
        
    catch ME
        
        fprintf( 'couldn''t load data in %s\n', fullfile( resultsDir, D(iD).name ) )
        
        disp(ME),
        
        indFailed(iD) = true;

    end
    
end

D = D(~indFailed);

idStrCell = { D.idStr };


%% Calc Entropy

E = [];

for iD = 1:numel(D),

    E(iD,1) = calc_imseq_metric( D(iD).imRlt, D(iD).mask );

    E(iD,2) = calc_imseq_metric( D(iD).imRdr, D(iD).mask );
    
    E(iD,3) = calc_imseq_metric( D(iD).imCine0, D(iD).mask );
	
    E(iD,4) = calc_imseq_metric( D(iD).imCineCardSync, D(iD).mask );
   
	E(iD,5) = calc_imseq_metric( D(iD).imCineNoMoco, D(iD).mask );
    
	E(iD,6) = calc_imseq_metric( D(iD).imCineNoOutrej, D(iD).mask );
	
    E(iD,7) = calc_imseq_metric( D(iD).imCine, D(iD).mask );

end


%% Relative Entropy

R = bsxfun(@rdivide,E,E(:,3));


%% Plot

hFig{1} = figure('Name','entropy_plot');
hAx{1} = gca;

plot(R(:,4:end)','o:','MarkerSize',8,'LineWidth',2,'MarkerSize',12),

ylabel('relative entropy' )  % ylabel('H(X^{[m]}) / H(X^{[m]})')
legend(idStrCell,'Location','SW','FontSize',9)

set(gca,'XLim',[0.17,4.83],'XTick',[1,2,3,4],'XTickLabel',{'cardsync','+outrej','+moco','+outrej+moco'})
set(gca,'YTick',[])
set(gca,'FontSize',16)


%% Boxplot

hFig{2} = figure('Name','entropy_boxplot');
hAx{2} = gca;

labels = {'cardsync','+outrej','+moco','+outrej+moco'};
h = boxplot(R(:,4:end),'Labels',labels,'Widths',0.33);
for ih = 1:size(h,1),
    set(h(ih,:),'LineWidth',2,'LineStyle','-','MarkerSize',12,'Color','k');
end

ylabel('relative entropy' )  % ylabel('H(X^{[m]}) / H(X^{[m]})')

set(gca,'XLim',[0.17,4.83],'XTick',[1,2,3,4],'XTickLabel',{'cardsync','+outrej','+moco','+outrej+moco'})
set(gca,'YTick',[])
set(gca,'FontSize',16)
set(gca,'TickLength',[0,0]);



%% Match YLim

set( hAx{1}, 'YLim', get(hAx{2},'YLim') )


end  % assess_cine_image_quality_v_postproc_steps(...)
