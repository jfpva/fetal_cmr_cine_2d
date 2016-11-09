function generate_correction_v_iteration_summary_per_case( resultsDir, fcmrStr )


%% Identify Results to Include

D = dir( resultsDir );

isInc = [ D.isdir ];  % include all directories

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
   if isInc(ind(iD))
       if numel( D(ind(iD)).name ) < numel( fcmrStr ),
           isInc(ind(iD)) = false;
       else
           isInc(ind(iD)) = strcmp( D(ind(iD)).name(1:numel(fcmrStr)), fcmrStr );
       end
   end
end

D = D(isInc);


%% Get Corrections

get_mask_values_vec = @( x, msk ) double( x( repmat( msk, [1,1,size(x,3)/size(msk,3)] ) ) );
get_mask_values_3d  = @( im, msk ) reshape( get_mask_values_vec( im, msk ), sum(msk(:)), 1, size(im,3) );

isFailed = false( size(D) );

fprintf( '\nLoading data...\n' )

for iD = 1:numel(D), 
    
    clear PARAM RESULTS mask 
    
    try
        
    fprintf( '%s\n', D(iD).name )

    load( fullfile( resultsDir, D(iD).name, 'cine', 'results.mat' ), 'PARAM', 'mask', 'RESULTS' );

    D(iD).rrInterval    = RESULTS.rrInterval(end);
    D(iD).rrRange       = range( [ RESULTS.rrInterval ] ) ;
    D(iD).disp          = squeeze( mean( get_mask_values_3d( abs( PARAM(end).T.dispMap ), mask ), 1 ) );
    D(iD).maxDisp       = max( D(iD).disp );
    D(iD).meanDisp      = RESULTS.meanDisp(end);
    D(iD).pctTotOutlier = RESULTS.pctTotOutlier(end);
    D(iD).meanVoxProb   = squeeze( mean(  get_mask_values_3d( PARAM(end).P.vox, mask ), 1 ) );
    D(iD).frmProb       = PARAM(end).P.frm;
    D(iD).meanProb      = squeeze( mean(  get_mask_values_3d( bsxfun( @times, PARAM(end).P.vox, reshape( PARAM(end).P.frm, 1, 1, [] ) ), mask ), 1 ) );
    
    catch ME
        
        isFailed(iD) = true;
        
        warning( 'Error processing %s:', D(iD).name )
        
        disp( ME ),
        
    end
    
end

indFailed = find ( isFailed );

if ~isempty( indFailed )
    fprintf( '\n\n\nFailed to process:\n' )
    for iD = indFailed,
        fprintf( '   %s\n', D(iD).name )
    end
    D = D(~isFailed);
end


%% Create Markdown Summary

logFilePath = fullfile( resultsDir, sprintf( 'summary_correction_v_iteration_%s.md', fcmrStr ) );

if exist( logFilePath, 'file' )
    delete( logFilePath ),
end

diary( logFilePath );

closeDiary = onCleanup( @() diary( 'off' ) );

fprintf( 'summary corrections\n' )
fprintf( '===================\n\n\n' )

fprintf( '## final corrections\n\n' )
fprintf( 'ID         | HR (bpm) | mean/max disp. (mm) | %% p < 0.5\n' )
fprintf( '---------- | -------- | ------------------- | ----------\n' )
for iD = 1:numel(D)
    fprintf( '%s | %7.1f  |     %3.1f / %3.1f     | %7.1f\n', D(iD).name, 60000/(1000*D(iD).rrInterval), D(iD).meanDisp, D(iD).maxDisp, D(iD).pctTotOutlier )
end
fprintf( '\n\n' )

fprintf( '## final corrections (sorted by mean displacement)\n\n' )
fprintf( 'ID         | HR (bpm) | mean disp. (mm) | %% p < 0.5\n' )
fprintf( '---------- | -------- | --------------- | ----------\n' )
[~,indDispOrder] = sort( [ D.meanDisp ] );
for iD = flip(indDispOrder)
    fprintf( '%s | %7.1f  | %9.1f       | %7.1f\n', D(iD).name, 60000/(1000*D(iD).rrInterval), D(iD).meanDisp, D(iD).pctTotOutlier )
end
fprintf( '\n\n' )

fprintf( '## final corrections (sorted by percent outliers)\n\n' )
[~,indOutliersOrder] = sort( [ D.pctTotOutlier ] );
fprintf( 'ID         | HR (bpm) | mean disp. (mm) | %% p < 0.5\n' )
fprintf( '---------- | -------- | --------------- | ----------\n' )
for iD = flip(indOutliersOrder)
    fprintf( '%s | %7.1f  | %9.1f       | %7.1f\n', D(iD).name, 60000/(1000*D(iD).rrInterval), D(iD).meanDisp, D(iD).pctTotOutlier )
end
fprintf( '\n\n' )

for iD = 1:numel(D)
    fprintf( '## %s\n\n', D(iD).name )
    fprintf( '![](%s/cine/figs/correction_v_iteration.png)  \n', D(iD).name )
    fprintf( '![](%s/cine/03a_cardsync/figs/heart_rate_estimate.png)  \n', D(iD).name )
end


end  % generate_correction_v_iteration_summary(...)