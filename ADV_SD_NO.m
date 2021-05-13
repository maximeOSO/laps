function [OUTPUTFgd, OUTPUTFtk] = ADV_SD_NO(configfile, PECCO2, PSD, TI, TSI, TF, ...
    INPUTF, OUTPUTP, MODE, GRZ, RHOP, TRK, AGES, PPSINK, PVSINK)
%%% Simulation du transport de particules dans l'ocean avec modele ECCO2 et
%%% settling velocity fonction de la taille de particule.

%%% INPUT:
%%% PECCO2: path to ECCO2 files
%%% PSD: path to Stokes Drift files (*** NOT USED HERE ***)
%%% TI: Advection and Injection START
%%% TSI: Injection STOP
%%% TF: Advection STOP
%%% INPUTF: path-to-file of particle input: X, Y, Z(def = 0)
%%% OUTPUTP: path-to-folder for the output files
%%% SV: settling velocity (Stokes - particle sinks?): YES/NO
%%% GRZ: Grain Size (only used if SV = YES)
%%% RHOP: Particle density (only used if SV = YES)
%%% TRK: Track the particles all along their travel : YES/NO
%%% RESINJ: Injection time step, in days (so 1hr = 1/24)

%%% OUTPUT:
%%% - OUTPUTFgd: path to matrix formatted data:
%%%             _ LON0 : mesh lon
%%%             _ LAT0 : mesh lat
%%%             _ dep : depth limits
%%%             _ COUNT3D: 3D matrix (LON, LAT, DEP), number of particles
%% Geo var
dt = 0.25/24; % 15min in day units
dts = dt * 86400; % dt in sec
rad_earth = 6371000; % earth radius m
%% Temps d'avection
ti = datenum(TI); % date de debut de simulation (datenum = jours)
tsf = datenum(TSI); % date de fin d'injection
tf = datenum(TF); % date de fin de simulation
resfile = 3; % resolution temporelle des champs de vitesse (jours)
[tseccoDN, ~] = ADV_tsfromecco2files(PECCO2, TI, TF,'adv');
tsecco_n = str2num(datestr(tseccoDN,'yyyymmdd'));
%% read input particles
fidpart = fopen(INPUTF);
inputparticles = textscan(fidpart, '%f %f %f %s');
fclose(fidpart);
x0 = inputparticles{1};
y0 = inputparticles{2};
z0 = inputparticles{3};
label_orig = inputparticles{4};
%% initalize particle tracking
TRACK = [];
trki = [];
trkf = [];
%% time series vs files name
[tsinjFMT, tsinjVRAI] = ADV_tsfromecco2files(PECCO2, TI, TSI,'inj');
tsecco_sel = str2num(datestr(tsinjFMT,'yyyymmdd'));
%% Init param avec un fichier ECCO2 (n'importe lequel)
ecco2varU = [PECCO2,'UVEL.1440x720x50.',datestr(tsinjFMT(1),'yyyymmdd'),'.nc'] ;
%% general
lat0 = double(ncread(ecco2varU,'LATITUDE_T'));
lon0 = double(ncread(ecco2varU,'LONGITUDE_T'));
dep = double(ncread(ecco2varU,'DEPTH_T'));
resll = lon0(2)-lon0(1);
edgeslat0 = min(lat0)-resll/2:resll:max(lat0)+resll/2;
edgeslon0 = min(lon0)-resll/2:resll:max(lon0)+resll/2;
resdep = diff(dep);
edgesdep = [dep(1)-resdep(1)/2;dep(2:end)-resdep/2;dep(end)+resdep(end)/2];
[LAT0, LON0] = meshgrid(lat0,lon0);
% Settling velocity
if strcmpi(MODE,'SED') == 1
    g = 9.8;
    rhop = RHOP;
    rhof = 1030;
    Dp = GRZ;
    mu = 1.4e-3 ; % Pa.s
    vsi = g*(rhop-rhof)*Dp^2/(18*mu);
else
    vsi = 0;
end
v0 = x0*0+vsi; % reinitialise les vitesses de chute a chaque injection
%% Init
COUNT = zeros(size(LAT0,1),size(LAT0,2));
COUNT3D = zeros(size(LAT0,1),size(LAT0,2),length(dep));
%% Initialisation
% la boucle d'injection est dans le while
dtstep = (datenum(TSI) - tsinjVRAI(1))/dt;
t0name = tsecco_sel(1) ; % format date dans nom de fichier
t0 = tsinjFMT(1);
kk = 1;
t = t0 ;
trk0 = [x0*0+t0 x0 y0 z0];
trki = [trki; trk0];
tname = t0name;
tnameold = tname;
a = 0 ;
u_can_move = 1 ; % initialisation

% Prepare sous-region a importer
mrg = 6; % marge en degres
% cas particulier latitude: si on depasse +/-90
if min(y0)-mrg <-90
    lati = -89.8750; latf = max(y0)+mrg;
    chglo = 1; % il faut changer les bornes des longitudes car on passe de l'autre cote du globe
elseif max(y0)+mrg>90
    lati = min(y0)-mrg; latf = 89.8750;
    chglo = 1; % il faut changer les bornes des longitudes car on passe de l'autre cote du globe
else
    lati = min(y0)-mrg; latf = max(y0)+mrg;
    chglo = 0;
end
% cas particulier longitude: si on depasse -0/+360
if min(x0)-mrg < 0 || max(x0)+mrg > 360 || chglo == 1
    loni = 0.125; lonf = 359.875;
else
    loni = min(x0)-mrg; lonf = max(x0)+mrg;
end
[~,ulati] = histc(lati,edgeslat0); [~,ulatf] = histc(latf,edgeslat0);
[~,uloni] = histc(loni,edgeslon0); [~,ulonf] = histc(lonf,edgeslon0);
[~,ud] = histc(z0,edgesdep);
startlon = uloni;
startlat = ulati;
cntlat = ulatf - ulati;
cntlon = ulonf - uloni;
mrgdep = 50;  % max to avoid error l268: idx = sub2ind(size(U), ulo, ula, ud');
if max(ud)<length(dep)-mrgdep
    cntdep = max(ud)+mrgdep;
else
    cntdep = length(dep);
end
lat = double(ncread(ecco2varU,'LATITUDE_T',startlat,cntlat));
lon = double(ncread(ecco2varU,'LONGITUDE_T',startlon,cntlon));
edgeslat = min(lat)-resll/2:resll:max(lat)+resll/2;
edgeslon = min(lon)-resll/2:resll:max(lon)+resll/2;
[LAT, LON] = meshgrid(lat,lon);
srt = [startlon startlat 1 1];
cnt = [cntlon cntlat cntdep 1]; %length(dep)
% champ de vitesses initial
U = double(ncread([PECCO2,'UVEL.1440x720x50.',num2str(tname),'.nc'],'UVEL',srt,cnt));
V = double(ncread([PECCO2,'VVEL.1440x720x50.',num2str(tname),'.nc'],'VVEL',srt,cnt));
W = double(ncread([PECCO2,'WVEL.1440x720x50.',num2str(tname),'.nc'],'WVEL',srt,cnt));
W(W<-1000) = 0;
% initialize les coordonnees des particules
x = [];
y = [];
z = [];
vs = [];
age = [];
sink = [];
timetrack = [];
t_prev = t0; % == ti normalement
%% Main advection
disp(' ')
disp('Advection')
% Advection principale, tourne tant que:
% - on n'a pas atteint la date finale de simul
% - des particules peuvent etre deplacees
while t < tf - dt && isempty(u_can_move) == 0
    % age of the particle
    age = age + dts/86400; % age is in days
    % If microplastic debris
    if strcmpi(MODE,'MPD') == 1
        n_mpd_can_sink = find(age>AGES & sink==0); % older than age and should not be sinking already
        if ~isempty(n_mpd_can_sink)
            % randolmy choose PSINK % of these MDP and set them to sink = 1
            % disp('Proba sink')
            chgs = randi(length(n_mpd_can_sink), [round(length(n_mpd_can_sink)*PPSINK) 1]);
            sink(chgs) = 1;
            vs(chgs) = PVSINK; % set these partciles to sink at VSINK velocity
            % disp([num2str(length(chgs)),'   -   ',num2str(length(age)),'   -   ',num2str(length(vs))])
        end
    end
    % tracking 
    if TRK > 0
        TRK_day = TRK/24;
        if mod(a,TRK_day/dt) == 0
            disp(['--- RECORD TRACK on ',datestr(t,'dd-mmm-yyyy HH'),' ---'])
            if a == 0
                xto = x0;
                idp = (1:length(x0))';
                TRACK = [TRACK; x0*0+t x0 y0 z0 idp];
            else
                if length(xto) == length(x)  % le nombre de particule est reste le meme sur ce dt donc pas dinjection, donc pas de changement d'indice
                    idp;
                    TRACK = [TRACK; x*0+t x y z idp];
                    xto = x;
                else   % le nombre de particules a change, donc injection, donc nouveaux indices
                    % There are n_newpart = length(x) - length(xto) new particles
                    n_newpart = length(x) - length(xto); % amount of new particles
                    idp = [idp; (1:n_newpart)' + max(idp)];
                    TRACK = [TRACK; x*0+t x y z idp];
                    xto = x;
                end
            end
        end
    end
    % disp(['npart=', num2str(length(x)),' - xto=',num2str(length(xto)), ' - idp range = ', num2str(idp(1)), ' -> ', num2str(idp(end))])
    % injection every day
    if t<tsf && (t==t0 || floor(t)-floor(t_prev)==1) % every day
        disp(['Injection: ' datestr(t,'dd/mm/yyyy'),' done on ',datestr(now,'dd/mm/yyyy HH:MM')])
        x = [x; x0];
        y = [y; y0];
        z = [z; z0];
        vs = [vs; v0];
        age = [age; x0*0]; % the new batch is aged 0
        sink = [sink; x0*0]; % the new batch is not sinking
        timetrack = [timetrack; x0*0 + t];  % ces t sont les temps d'INJECTION
    end
    kk = kk + 1;
    if tnameold ~= tname
        % Prepare sous-region a importer
        % cas partic latitude
        if min(y)-mrg <-90
            lati = -89.8750; latf = max(y)+mrg;
            chglo = 1; % il faut changer les bornes des longitudes car on passe de l'autre cote du globe
        elseif max(y)+mrg>90
            lati = min(y)-mrg; latf = 89.8750;
            chglo = 1; % il faut changer les bornes des longitudes car on passe de l'autre cote du globe
        else
            lati = min(y)-mrg; latf = max(y)+mrg;
            chglo = 0;
        end
        % cas partic longitude
        if min(x)-mrg < 0 || max(x)+mrg > 360 || chglo == 1
            loni = 0.125; lonf = 359.875;
        else
            loni = min(x)-mrg; lonf = max(x)+mrg;
        end
        [~,ulati] = histc(lati,edgeslat0); [~,ulatf] = histc(latf,edgeslat0);
        [~,uloni] = histc(loni,edgeslon0); [~,ulonf] = histc(lonf,edgeslon0);
        [~,ud] = histc(z,edgesdep);
        startlon = uloni;
        startlat = ulati;
        cntlat = ulatf - ulati;
        cntlon = ulonf - uloni;
        if max(ud)<length(dep)-mrgdep
            cntdep = max(ud)+mrgdep;
        else
            cntdep = length(dep);
        end
        lat = double(ncread(ecco2varU,'LATITUDE_T',startlat,cntlat));
        lon = double(ncread(ecco2varU,'LONGITUDE_T',startlon,cntlon));
        edgeslat = min(lat)-resll/2:resll:max(lat)+resll/2;
        edgeslon = min(lon)-resll/2:resll:max(lon)+resll/2;
        [LAT, LON] = meshgrid(lat,lon);
        srt = [startlon startlat 1 1];
        cnt = [cntlon cntlat cntdep 1];
        U = double(ncread([PECCO2,'UVEL.1440x720x50.',num2str(tname),'.nc'],'UVEL',srt,cnt));
        V = double(ncread([PECCO2,'VVEL.1440x720x50.',num2str(tname),'.nc'],'VVEL',srt,cnt));
        W = double(ncread([PECCO2,'WVEL.1440x720x50.',num2str(tname),'.nc'],'WVEL',srt,cnt));
        W(W<-1000) = 0;
    end
    a = a + 1 ;
    % trouve l'indice de profondeur
    [~,ud] = histc(z,edgesdep);
    % trouve les indices de position de la
    % particule dans le mesh afin d'en tirer les
    % vitesses u et v associees
    [Nla,lala] = histc(y,edgeslat);
    ula = lala' ;
    [Nlo,lolo] = histc(x,edgeslon);
    ulo = lolo';
    % cas ou y est dans l'intervalle 89.8750 - 90
    iju0_avant = find(y>89.750);
    iju0_apres = find(y<-89.75);
    ula(iju0_avant) = length(edgeslat)-1;
    ula(iju0_apres) = 1;
    % cas ou x est dans l'intervalle 359.875 - 0.125
    ku0 = find(ulo == 0);
    iku0_avant = find(x>359.750);
    iku0_apres = find(x<0.125);
    ulo(iku0_avant) = max(ulo);
    ulo(iku0_apres) = 1;
    % Position dans la matrice LON/LAT/DEPTH/TIME(=1)
    % U, V et W ont le meme format, meme indice pour les 3
    idx0 = sub2ind(size(U), ulo, ula, 1+ula*0); % indice sur la premiere slice de la matrice, pour LAT --- 4 eme dim 1+ula*0
    idx = sub2ind(size(U), ulo, ula, ud');
    ue = U(idx);
    ve = V(idx);
    we = W(idx); % ajoute pas vs ici pour check if vitesse ECO2 nulle
    vel_mag = sqrt(ue.^2 + ve.^2 + we.^2);
    % Check conditions d'arret
    u_not_move = find(vel_mag' ==0); % cas ou U=0 (echouage, touche cotes)
    u_can_move = find(vel_mag ~=0); % U!=0 --> bouge encore en horiz
    % update coord
    x = x + ue'*dts./(pi/180*rad_earth*sind(90-LAT(idx0)'));
    y = y + ve'*dts./(pi/180*rad_earth*sind(90-LAT(idx0)'));
    z(u_can_move) = z(u_can_move) + vs(u_can_move)*dts + we(u_can_move)'*dts; % z-axis positive downward
    z(u_not_move) = z(u_not_move); % profondeur change pas
    z(z<0) = 0.001; % particules au dessus de l'eau (cause vitesse verticale) sont remises a 0
    % passe la limite -90/90 en latitude
    x(y>90 | y<-90) = x(y>90 | y<-90)+180; % commence par changer les longitudes avant de modifier y (sinon on perd l'info ><90)
    y(y>90)=180-y(y>90);
    y(y<-90)=-180-y(y<-90);
    % passe la limite 0/360 en longitude
    x(x<0)=360+(x(x<0));
    x(x>=360)=x(x>=360)-360;
    tnameold = tname ;
    t_prev = t; % memory this t
    t = t + dt;  % t au format datenum Matlab (jours)
    [~,udt] = histc(t,tseccoDN);
    tname = tsecco_n(udt);
end
%% Wrap-up particule final positions
% COUNT
[xsortc,sxc] = sort(x);
[ysortc,syc] = sort(y);
% trouve les indices dans le mesh COUNT
Nlac = histc(y,edgeslat0);
ula0c = find(Nlac>0);
ula1c = [ula0c Nlac(ula0c)];
ulac_ok = [];
for uula = 1:size(ula1c,1)
    test_ulac = zeros(1,ula1c(uula,2))+ula1c(uula,1);
    ulac_ok = [ulac_ok test_ulac];
end
ulac_ok(syc) = ulac_ok; % remet les indices dans l'ordre initial des particules (pas de tri croissant !!!)
ulac = ulac_ok;
Nloc = histc(x,edgeslon0);
ulo0c = find(Nloc>0);
ulo1c = [ulo0c Nloc(ulo0c)];
uloc_ok = [];
for uulo = 1:size(ulo1c,1)
    test_uloc = zeros(1,ulo1c(uulo,2))+ulo1c(uulo,1);
    uloc_ok = [uloc_ok test_uloc];
end
uloc_ok(sxc) = uloc_ok; % remet les indices dans l'ordre initial des particules (pas de tri croissant !!!)
uloc = uloc_ok;
idxc = sub2ind(size(COUNT), uloc, ulac);
% il y a des doublons dans idxc car des particules peuvent
% se retrouver dans la meme zone --> doit boucler pour
% toutes les compter
for k = 1:length(idxc)
    COUNT(idxc(k)) = COUNT(idxc(k)) + 1;
end
% COUNT3D
[zsortc,szc] = sort(z);
Ndep = histc(z,edgesdep);
ude0c = find(Ndep>0);
ude1c = [ude0c Ndep(ude0c)];
udec_ok = [];
for uude = 1:size(ude1c,1)
    test_udec = zeros(1,ude1c(uude,2))+ude1c(uude,1);
    udec_ok = [udec_ok test_udec];
end
udec_ok(szc) = udec_ok; % get the good index order, not sorted ascend
udec = udec_ok;
idxc3D = sub2ind(size(COUNT3D), uloc, ulac, udec);
% il y a des doublons dans idxc car des particules peuvent
% se retrouver dans la meme zone --> doit boucler pour
% toutes les compter
for k = 1:length(idxc3D)
    COUNT3D(idxc3D(k)) = COUNT3D(idxc3D(k)) + 1;
end

%% save final positions
OUTPUTFgd = [OUTPUTP, ...
    'MAT_' datestr(ti,'yyyymmdd'),'_',datestr(TSI,'yyyymmdd'),...
    '_',datestr(tf,'yyyymmdd'),'_',MODE, '_', configfile,'.mat'];
save(OUTPUTFgd,'LON0', 'LAT0','dep', 'COUNT3D')
disp(['Final positions: ', OUTPUTFgd])

%% save tracking
if TRK > 0
    timestr = datestr(dateshift(datetime(datestr(TRACK(:,1))), 'start', 'minute', 'nearest'), 'yyyymmdd HHMMSS');
    n_step = length(TRACK)/length(x0);
    L = {};
    for k = 1:n_step
        L = [L; label_orig];
    end
    OUTPUTFtk = [OUTPUTP, 'TRK_' datestr(ti,'yyyymmdd'),'_',datestr(TSI,'yyyymmdd'),...
        '_',datestr(tf,'yyyymmdd'),'_',MODE, '_', configfile];
    fid = fopen(OUTPUTFtk,'w');
    for k = 1:length(TRACK)
        fprintf(fid,'%s %.6f %.6f %.3f %s%06d \n',timestr(k,:),TRACK(k,2:4), L{k}, TRACK(k,5));
    end
    fclose(fid);
    disp(['Tracking: ', OUTPUTFtk])
else
    OUTPUTFtk = [OUTPUTP,'null.tmp'];
end