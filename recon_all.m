% recon_all( dataDir )

%% Get list data directories to process

if ~exist( 'dataDir', 'var' ),
    dataDir = fullfile( fileparts( mfilename ), '..', 'data' );
    if ~exist( dataDir, 'dir' )
        error( 'data directory doesn''t exist, %s\n', dataDir )
    end
end

D = dir( dataDir );

isInc = [ D.isdir ];  % include all directories

% loop to identify hidden directories
ind = find(isInc);
for iD = ind, 
   % on OSX, hidden directories start with a dot
   isInc(iD) = ~strcmp(D(iD).name(1),'.');
   if isInc(iD) && ispc
       % check for hidden Windows directories - only works on Windows
       [~,stats] = fileattrib(fullfile(dataDir,D(iD).name));
       if stats.hidden
          isInc(iD) = false;
       end
   end
   if isInc(iD),
       isInc(iD) = strcmp( D(iD).name(1:4), 'fcmr' );
   end
end

D = D(isInc);


%% Reconstruct

isFailed = false( 1, numel(D) );

for iD = 1:numel(D), 
    
    fcmrNo   = 0;
    seriesNo = 0;
    
    try 
    
    fcmrNo   = str2double( D(iD).name(5:7) );
    seriesNo = str2double( D(iD).name(9:10) );
    
    main_fcmr_cine( fcmrNo, seriesNo );
    
    catch ME
        
        isFailed(iD) = true;
        
        warning( 'Caught error reconstructing %s:', D(iD).name )
        
        disp( ME ),
        
    end
    
end

fprintf( '\n\n\nFailed to reconstruct:\n' )
indFailed = find ( isFailed );
for iF = indFailed,
    fprintf( '   %s\n', D(iF).name )
end

% end