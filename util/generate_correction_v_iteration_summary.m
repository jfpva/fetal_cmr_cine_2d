function generate_correction_v_iteration_summary( resultsDir )


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
       [~,stats] = fileattrib(fullfile(dataDir,D(iD).name));
       if stats.hidden
          isInc(iD) = false;
       end
   end
   if isInc(ind(iD))
       isInc(ind(iD)) = strcmp( D(ind(iD)).name(1:4), 'fcmr' );
   end
end

D = D(isInc);


%% Create Markdown Summary

logFilePath = fullfile( resultsDir, 'summary_correction_v_iteration.md');

diary( logFilePath );

closeDiary = onCleanup( @() diary( 'off' ) );

fprintf( 'summary correction v iteration_summary\n' )
fprintf( '======================================\n' )

for iD = 1:numel(D);
    fprintf( '## %s\n\n', D(iD).name )
    fprintf( '![](%s/cine/figs/correction_v_iteration.png)  \n\n\n', D(iD).name )
end


end  % generate_correction_v_iteration_summary(...)