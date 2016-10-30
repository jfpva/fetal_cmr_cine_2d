% script_run_multiple_reconframe_get_kt_raw_data

generate_idstrcell;

isFailed = false( size( idStrCell ) );
for iS = 1:numel(idStrCell)
    idStr = idStrCell{iS};
    fcmrNo = str2double(idStr(5:7));
    seriesNo = str2double(idStr(9:10));
    try
        run_reconframe_get_kt_raw_data(fcmrNo,seriesNo);
    catch ME
        isFailed(iS) = true;
        warning( 'failed to recon %s', idStr )
        disp( ME ),
    end
end

indFailed = find( isFailed );
for iF = indFailed,
    fprintf( 'failed to recon %s\n', idStrCell{iF} )
end