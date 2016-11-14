function run_reconframe_get_kt_raw_data( fcmrNo, seriesNo ),
%RUN_RECONFRAME_GET_KT_RAW_DATA
%
%	RUN_RECONFRAME_GET_KT_RAW_DATA( fcmrNo, seriesNo ),

%% Identify Raw Data Dir

switch fcmrNo,
    
    case  53,  % s26
        rawDataDirDateStr   = '2015_08_04';
        rawDataDirName      = 'KH_80512';
    case 63,  % s28
        rawDataDirDateStr   = '2015_10_27';
        rawDataDirName      = 'BI_97800';
    case  64,  % s25, s26, s37, s38
        rawDataDirDateStr   = '2015_11_05';
        rawDataDirName      = 'LO_100002';
    case  65,  % s29
        rawDataDirDateStr   = '2015_11_06';
        rawDataDirName      = 'KO_100900';
    case  66,  % s26
        rawDataDirDateStr   = '2015_11_19';
        rawDataDirName      = 'TW_103100';
    case  67,  % s27, s28
        rawDataDirDateStr   = '2015_12_10';
        rawDataDirName      = 'EL_108201';
    case  68,  % s24, s26, s28, s30
        rawDataDirDateStr   = '2016_01_05';
        rawDataDirName      = 'JA_112304';
    case  69,  % s34, s35
        rawDataDirDateStr   = '2016_01_07';
        rawDataDirName      = 'LE_112801';
    case  70,  % s28, s29
        rawDataDirDateStr   = '2016_01_07';
        rawDataDirName      = 'AL_112802';
    case  71,  % s29
        rawDataDirDateStr   = '2016_01_14';
        rawDataDirName      = 'RI_114600';
    case  72,  % s49, s54, s55
        rawDataDirDateStr   = '2016_01_21';
        rawDataDirName      = 'BE_116600';
    case  73,  % s26
        rawDataDirDateStr   = '2016_01_29';
        rawDataDirName      = 'BO_119900';
    case  74,  % s29, s30 
        rawDataDirDateStr   = '2016_01_29';
        rawDataDirName      = 'IG_119901';
    case  75,  % s28
        rawDataDirDateStr   = '2016_02_08';
        rawDataDirName      = 'BR_122033';
    case  76,  % s29
        rawDataDirDateStr   = '2016_02_11';
        rawDataDirName      = 'GO_123100';
    case  77,  % s24 
        rawDataDirDateStr   = '2016_02_11';
        rawDataDirName      = 'CR_123101';
    case  78,  % s23, s28, s30
        rawDataDirDateStr   = '2016_02_18';
        rawDataDirName      = 'CA_125000';
    case  79,  % s29
        rawDataDirDateStr   = '2016_02_26';
        rawDataDirName      = 'HI_126401';
    case  81,  % s21
        rawDataDirDateStr   = '2016_03_15';
        rawDataDirName      = 'GR_130400';
    case  82,  % s25, s28, s30, 36
        rawDataDirDateStr   = '2016_03_17';
        rawDataDirName      = 'RO_131111';
    case  83,  % s26
        rawDataDirDateStr   = '2016_03_17';
        rawDataDirName      = 'DA_131112';
    case  84,  % s25, s27
        rawDataDirDateStr   = '2015_12_11';
        rawDataDirName      = 'MA_108800';
    case  85,  % s26, s29
        rawDataDirDateStr   = '2015_12_14';
        rawDataDirName      = 'HA_109212';
    case  86,  % s26, s27, s28
        rawDataDirDateStr   = '2015_12_18';
        rawDataDirName      = 'CA_111003';
    case  87,  % s26, s27
        rawDataDirDateStr   = '2015_12_18';
        rawDataDirName      = 'SP_111004';
    case  88,  % 28
        rawDataDirDateStr   = '2015_12_21';
        rawDataDirName      = 'KH_111036';
    case  89,  % s33, s37
        rawDataDirDateStr   = '2016_04_07';
        rawDataDirName      = 'HO_700';
    case  90,  % s32
        rawDataDirDateStr   = '2016_04_19';
        rawDataDirName      = 'BI_2901';
    case  92,  % s28
        rawDataDirDateStr   = '2016_04_28';
        rawDataDirName      = 'MC_1689';
    case  93,  % s20, s22
        rawDataDirDateStr   = '2016_05_19';
        rawDataDirName      = 'WA_14771';
    case  94,  % s28, s31
        rawDataDirDateStr   = '2016_05_19';
        rawDataDirName      = 'HA_15009';
    case  95,  % s31, s33, s35
        rawDataDirDateStr   = '2016_05_26';
        rawDataDirName      = 'JA_20069';
    case  96,  % s27
        rawDataDirDateStr   = '2016_05_27';
        rawDataDirName      = 'BR_21609';
    case  97,  % s30
        rawDataDirDateStr   = '2016_06_10';
        rawDataDirName      = 'DE_29333';
    case  98,  % s32
        rawDataDirDateStr   = '2016_06_13';
        rawDataDirName      = 'HO_30022';
    case  99, 
        rawDataDirDateStr   = '2016_06_16';
        rawDataDirName      = 'MA_32402';
    case 100,  % s34
        rawDataDirDateStr   = '2016_06_20';
        rawDataDirName      = 'TE_34437';
    case 101,  % s23, s24, s25
        rawDataDirDateStr   = '2016_06_27';
        rawDataDirName      = 'BU_39097';
    case 102,  % s24, s27
        rawDataDirDateStr   = '2016_06_30';
        rawDataDirName      = 'BU_42322';
    case 103,  % s22, s24, s53
        rawDataDirDateStr   = '2016_07_05';
        rawDataDirName      = 'FR_46424';
    case 105,  % s23
        rawDataDirDateStr   = '2016_07_28';
        rawDataDirName      = 'WE_66966';
    case 106,  % s27, s31
        rawDataDirDateStr   = '2016_08_11';
        rawDataDirName      = 'RA_77501';
    case 107,  % s22 s24
        rawDataDirDateStr   = '2016_08_11';
        rawDataDirName      = 'RA_77887';
    case 108,  % s31
        rawDataDirDateStr   = '2016_08_12';
        rawDataDirName      = 'NI_79199';
    case 109,  % s25
        rawDataDirDateStr   = '2016_08_18';
        rawDataDirName      = 'LY_82088';
    case 111,  % s45
        rawDataDirDateStr   = '2016_08_25';
        rawDataDirName      = 'BE_86437';
    case 112,  % s29
        rawDataDirDateStr   = '2016_08_25';
        rawDataDirName      = 'TH_86739';
    case 113,  % s23
        rawDataDirDateStr   = '2016_08_26';
        rawDataDirName      = 'TH_88070';
    case 114,  % s26
        rawDataDirDateStr   = '2016_08_26';
        rawDataDirName      = 'UN_87798';
    case 116,  % s21, s23
        rawDataDirDateStr   = '2016_09_08';
        rawDataDirName      = 'WH_94316';
    case 117,  % s29
        rawDataDirDateStr   = '2016_09_09';
        rawDataDirName      = 'BA_95633';
    case 119,  % s25, s27
        rawDataDirDateStr   = '2016_09_15';
        rawDataDirName      = 'EL_100630';   
    case 120,  % s28
        rawDataDirDateStr   = '2016_09_22';
        rawDataDirName      = 'AL_107918';   
    case 122,  % s36, s40
        rawDataDirDateStr   = '2016_10_13';
        rawDataDirName      = 'SN_125152';   
    case 123,  % s22, s24
        rawDataDirDateStr   = '2016_10_17';
        rawDataDirName      = 'BA_128047'; 
    case 124,  % s37-48
        rawDataDirDateStr   = '2016_10_20';
        rawDataDirName      = 'WH_131485'; 
    case 126,  % s22-30
        rawDataDirDateStr   = '2016_10_21';
        rawDataDirName      = 'ST_133099'; 
    case 127,  % s22-30, 32-41
        rawDataDirDateStr   = '2016_10_28';
        rawDataDirName      = 'MA_137139'; 
    case 130,  % s21-30, 33-40
        rawDataDirDateStr   = '2016_11_04';
        rawDataDirName      = 'FI_141988'; 
    case 131,  % s19-28, 30-38
        rawDataDirDateStr   = '2016_11_07';
        rawDataDirName      = 'PI_142913'; 
%     case 000,  % s00
%         rawDataDirDateStr   = '0000_00_00';
%         rawDataDirName      = 'AA_000000'; 
    otherwise
        error( 'Raw data directory undefined for fcmr%03i', fcmrNo )
end



rawDataDir = fullfile( strcat( filesep, 'home'), 'jva13', 'mnt', 'pnraw01-archive', 'archive-ingenia', rawDataDirDateStr, rawDataDirName );
fprintf( 'Trying raw data directory %s ', rawDataDir ),
if ~exist( rawDataDir, 'dir' )    
    fprintf( ' ...cannot be found.\n' ),
    rawDataDir = fullfile( strcat( filesep, 'home'), 'jva13', 'mnt', 'pnraw01-ingenia', rawDataDirDateStr, rawDataDirName );
    fprintf( 'Trying raw data directory %s.', rawDataDir ),
    if ~exist( rawDataDir, 'dir' )    
        fprintf( ' ...cannot be found.\n' ),
        warning( 'Raw data directory %s cannot be found.', rawDataDir ),
        error( 'No raw data directory found.' )
    else
        fprintf( ' ...complete.\n' ),
    end
else
    fprintf( ' ...complete.\n' ),
end



%% Define Save Data Dir


saveDataDir = fullfile( strcat( filesep, 'scratch'), 'jva13', 'kt_raw_data', sprintf( 'fcmr%03is%02i', fcmrNo, seriesNo ) );


%% Get Data


reconframe_get_kt_raw_data( rawDataDir, seriesNo, 'saveDataDir', saveDataDir, 'verbose', true );  



end  % run_reconframe_get_ktdata(...)