function save_figs( workingDir )

figDir = fullfile( workingDir, 'figs' );
matfigDir = fullfile( figDir, 'matfigs' );

if ~exist( figDir, 'dir' )
    mkdir( figDir )
end

if ~exist( matfigDir, 'dir' )
    mkdir( matfigDir );
end

h = get(0,'children');

offset = 0;
d = dir( fullfile( matfigDir, '*.fig' ) );
if ~isempty(d),
    s = strsplit(d(end).name,'.');
    offset = str2num(s{1});
end
offset = offset + 1;

for i=1:length(h)
    figName = h(i).Name;
    if isempty( figName )
        figName = sprintf( '%03i', offset );
        offset = offset + 1;
    end    
    saveas( h(i), fullfile( matfigDir, figName ), 'fig' );
    saveas( h(i), fullfile( figDir, figName ), 'png' );
end
