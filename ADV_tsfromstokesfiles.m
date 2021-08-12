function [tsFMT, tsvrai]  = ADV_tsfromstokesfiles(PSD, TI, TF, restime,typedate)
%%% Makes a time series of the Stokes drift files needed 
%%% for the simulation (use TI, TF and the files names)

%%% MIT License, copyright (c) 2021 Maxime Mouyen, https://github.com/maximeOSO/laps#license

% the full time series of SD as online
time_first_StDr = [1992 01 01 0 0 0];
restimeStDr = 1/8; % new velocity field every restimeStDr days
StDrfull = datenum(time_first_StDr):restimeStDr:(now-100); 
StDrfull_mth0 = datevec(StDrfull);
StDrfull_mth = unique([StDrfull_mth0(:,1)*100 + StDrfull_mth0(:,2)]);

% the user StDr time series fmt_SD_201801_1440x720xtime3h.nc
Ufn0 = dir([PSD,'fmt*.nc']);
Ufn1 = char(Ufn0.name);
UStDruser = str2num(Ufn1(:,8:11))*100 + str2num(Ufn1(:,12:13));

% the needed time series
tneed3h = datenum(TI):restime:datenum(TF);
tneed3hFMT = datevec(tneed3h);
tneedMONTH = tneed3hFMT(:,1:2);
tneedSORT = unique(tneedMONTH,'stable','rows');
tsvrai = tneed3h;
tfiles = tneedSORT(:,1)*100 + tneedSORT(:,2);
tsFMT = UStDruser;

% Check
m = 0;
disp(['Checking if all needed STOKES DRIFT files are in ', PSD, ' ...'])
for KK= 1:length(tfiles)
    if isempty(find(UStDruser == tfiles(KK)))
        disp (['! File fmt_SD_',num2str(tfiles(KK)),'_1440x720xtime3h.nc is missing'])
        m = m+1;
    end
end

if m>0
    disp('************** ')
    disp(['ERROR: ', num2str(m),' velocity files are missing in ', PSD])
    disp('Get them here --> ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST/GLOBAL/')
    disp('************** ')
    return
elseif m ==0 && strcmpi(typedate,'ADV') == 1
    disp('OK')
end

