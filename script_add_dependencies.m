% script_add_dependencies

% NOTE: this file must be located project base directory 
% (e.g., fcmr_cine/src) to perform as expected

resetPath = onCleanup( sprintf( '@() path( %s )', path ) );

[StackTrace,indWorkspace] = dbstack;
dependenciesDir = fileparts( which( StackTrace(indWorkspace).file ) );

addpath( genpath( dependenciesDir ) );

excludeDir = dir( fullfile( dependenciesDir, '_*' ) );
for iD = 1:numel(excludeDir),
    if ( excludeDir(iD).isdir ),
        rmpath( genpath( fullfile( dependenciesDir, excludeDir(iD).name ) )  )
    end
end

clearvars -except resetPath
