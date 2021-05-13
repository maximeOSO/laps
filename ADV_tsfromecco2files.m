function [tsFMT, tsvrai]  = ADV_tsfromecco2files(PECCO2, TI, TF, typedate)
%%% Fait une serie temporelle des fichiers necessaires a la simulation a
%%% partir des TI et TF de simul et des noms de fichier.

% the full time series of ECCO2 as online
time_first_ecco2 = [1992 01 02 0 0 0];
restimeecco2 = 3; % new velocity field every restimeecco2 days
ecco2full = datenum(time_first_ecco2):restimeecco2:(now-40);

% the user U ecco2 time series
Ufn0 = dir([PECCO2,'UVEL*']);
Ufn1 = char(Ufn0.name);
Uecco2user = datenum([str2num(Ufn1(:,18:21)) str2num(Ufn1(:,22:23)) ...
    str2num(Ufn1(:,24:25)) zeros(size(Ufn1,1),1) ...
    zeros(size(Ufn1,1),1) zeros(size(Ufn1,1),1)]);

% the user V ecco2 time series
Vfn0 = dir([PECCO2,'VVEL*']);
Vfn1 = char(Vfn0.name);
Vecco2user = datenum([str2num(Vfn1(:,18:21)) str2num(Vfn1(:,22:23)) ...
    str2num(Vfn1(:,24:25)) zeros(size(Vfn1,1),1) ...
    zeros(size(Vfn1,1),1) zeros(size(Vfn1,1),1)]);

% the user W ecco2 time series
Wfn0 = dir([PECCO2,'WVEL*']);
Wfn1 = char(Wfn0.name);
Wecco2user = datenum([str2num(Wfn1(:,18:21)) str2num(Wfn1(:,22:23)) ...
    str2num(Wfn1(:,24:25)) zeros(size(Wfn1,1),1) ...
    zeros(size(Wfn1,1),1) zeros(size(Wfn1,1),1)]);

% the needed time series
%tneed = datenum(TI):restimeecco2:datenum(TF)  !!! FAUX, ca peut zapper un
%temps de calcul si le delta TI/TF n'est pas divisble par restimeecco2
tneed = datenum(TI):datenum(TF);
tsvrai = tneed;

[~,udtECCO2] = histc(tneed,ecco2full); % the files needed
udtECCO2 = unique(udtECCO2); % elimine les multiples (3 dates peuvent avoir le meme indice)

if strcmpi(typedate,'CHK') == 1
    % Check
    %     datesn = [];
    %     disp('ECCO2 U, V and W components are needed on dates:')
    %     for kk = 1:length(udtECCO2)
    %     datesn = [datesn, datestr(ecco2full(udtECCO2(kk)),'yyyymmdd'), ' '];
    %     end
    %     disp(datesn)
    % disp(['Checking if all needed ECCO2 files are in ', PECCO2, ' ...'])
    m = 0;
    for KK= 1:length(udtECCO2)
        if isempty(find(Uecco2user == ecco2full(udtECCO2(KK))))
            disp (['! File UVEL.1440x720x50.',datestr(ecco2full(udtECCO2(KK)),'yyyymmdd'),'.nc is missing'])
            m = m+1;
        end
        if isempty(find(Vecco2user == ecco2full(udtECCO2(KK))))
            disp (['! File VVEL.1440x720x50.',datestr(ecco2full(udtECCO2(KK)),'yyyymmdd'),'.nc is missing'])
            m = m+1;
        end
        if isempty(find(Wecco2user == ecco2full(udtECCO2(KK))))
            disp (['! File WVEL.1440x720x50.',datestr(ecco2full(udtECCO2(KK)),'yyyymmdd'),'.nc is missing'])
            m = m+1;
        end
    end
    
    if m>0
        disp('************** ')
        disp(['ERROR: ', num2str(m),' velocity files are missing in ', PECCO2])
        disp('Get them here --> ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/')
        disp('************** ')
        return
    elseif m ==0 && strcmpi(typedate,'CHK') == 1
        disp('ECCO2 files OK')
    end
    
elseif strcmpi(typedate,'ADV') == 1
    % If all files exist -> final series, can be any of U, V, W.
    % The time series needed from ECCO2
    tsFMT = [ecco2full(udtECCO2) ecco2full(udtECCO2(end)+1)];
    % NOTE : the last value in tsFMT is NOT needed for computation but it has
    % to be given as the upper bound of the time, because one ECCO2 file is
    % valid from the day in its name to 3 days after. If we don't do that then
    % the code will skip the remaining 3 days.
elseif strcmpi(typedate,'INJ') == 1
    tsFMT = ecco2full(udtECCO2);
end

