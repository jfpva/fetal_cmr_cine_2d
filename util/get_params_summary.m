% get_params_summary( dataDir )

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


%% Get Params


isFailed = false( size(D) );

fprintf( '\nLoadin data...\n' )

for iD = 1:numel(D), 
    
    clear dt pixdimAcq
    
    try
        
    fprintf( '%s\n', D(iD).name )

    M = dir( fullfile( dataDir, D(iD).name, '*_param.mat' ) );

    load( fullfile( dataDir, D(iD).name, M(1).name ), 'dt', 'pixdimAcq' );

    D(iD).dx = pixdimAcq(1,1);
    D(iD).dy = pixdimAcq(1,2);
    D(iD).dt = dt;
    
    catch ME
        
        isFailed(iD) = true;
        
        warning( 'Caught error loading %s:', D(iD).name )
        
        disp( ME ),
        
    end
    
end

indFailed = find ( isFailed );

if ~isempty( indFailed )
    fprintf( '\n\n\nFailed to load:\n' )
    for iD = indFailed,
        fprintf( '   %s\n', D(iD).name )
    end
    D = D(~isFailed);
end


%% Summarise


fprintf( '\n\nID      dx  dy  dt  \n' )

for iD = 1:numel(D),
    fprintf( '%s   %.1f   %.1f   %.3f\n', D(iD).name, D(iD).dx, D(iD).dy, D(iD).dt )
end
