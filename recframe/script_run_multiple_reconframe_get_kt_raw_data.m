% script_run_multiple_reconframe_get_kt_raw_data

generate_idstrcell;

idStrCell = { ... 
'fcmr053s26', ...
'fcmr063s28', ...
'fcmr081s21', ...
'fcmr100s34' };

for iS = 1:numel(idStrCell)
idStr = idStrCell{iS};
fcmrNo = str2double(idStr(5:7));
seriesNo = str2double(idStr(9:10));
run_reconframe_get_kt_raw_data(fcmrNo,seriesNo);
end