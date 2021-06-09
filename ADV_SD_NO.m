function [OUTPUTFgd, OUTPUTFtk] = ADV_SD_NO(configfile, PECCO2, PSD, TI, TSI, TF, ...
    INPUTF, OUTPUTP, MODE, GRZ, RHOP, TRK, AGES, PPSINK, PVSINK)
%%% LAPS without Stokes drift velocities

%%% INPUT: cf config file

%%% OUTPUT:
%%% - OUTPUTFgd: path to matrix formatted data:
%%%             _ LON0 : mesh lon
%%%             _ LAT0 : mesh lat
%%%             _ dep : depth limits
%%%             _ COUNT3D: 3D matrix (LON, LAT, DEP), number of particles

%% Parameters
dt = 0.25/24; % 15 min in day units. The particle position is updatede every dt
dts = dt * 86400; % dt in sec
rad_earth = 6371000; % earth radius m
%% Advection times
ti = datenum(TI); % simulation/injection start (datenum = days)
tsf = datenum(TSI); % injection stops
tf = datenum(TF); % simulation stops
resfile = 3; % Time resolution of ECCO2 files (days)
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
v0 = x0*0+vsi; % re-initialize settling velocity at each injection
%% Initialisation
% (injection loop is in the while loop below)
COUNT = zeros(size(LAT0,1),size(LAT0,2));
COUNT3D = zeros(size(LAT0,1),size(LAT0,2),length(dep));
dtstep = (datenum(TSI) - tsinjVRAI(1))/dt;
t0name = tsecco_sel(1) ; % date format in the file name
t0 = tsinjFMT(1);
kk = 1;
t = t0 ;
trk0 = [x0*0+t0 x0 y0 z0];
trki = [trki; trk0];
tname = t0name;
tnameold = tname;
a = 0 ;
u_can_move = 1 ;

% Prepare sub-region to load: this allows to not load the entire netcdf
% file but only the area were particles are, plus a safety margin around
mrg = 6; % margin in degrees
% special case latitude: if greater/lower than +/-90
if min(y0)-mrg <-90
    lati = -89.8750; latf = max(y0)+mrg;
    chglo = 1; % must change longitudes boundaries because particle moves on the other side of the earth
elseif max(y0)+mrg>90
    lati = min(y0)-mrg; latf = 89.8750;
    chglo = 1; % must change longitudes boundaries because particle moves on the other side of the earth
else
    lati = min(y0)-mrg; latf = max(y0)+mrg;
    chglo = 0;
end
% special case longitude: if greater/lower than 0/360
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
mrgdep = 50;  % max to avoid error line 257: idx = sub2ind(size(U), ulo, ula, ud');
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
cnt = [cntlon cntlat cntdep 1]; % length(dep)
% Initial velocity field
U = double(ncread([PECCO2,'UVEL.1440x720x50.',num2str(tname),'.nc'],'UVEL',srt,cnt));
V = double(ncread([PECCO2,'VVEL.1440x720x50.',num2str(tname),'.nc'],'VVEL',srt,cnt));
W = double(ncread([PECCO2,'WVEL.1440x720x50.',num2str(tname),'.nc'],'WVEL',srt,cnt));
W(W<-1000) = 0;
% initialize particules coordinates
x = [];
y = [];
z = [];
vs = [];
age = [];
sink = [];
timetrack = [];
t_prev = t0;
%% Main advection
disp(' ')
disp('Advection')
% Main advection, runs while:
% - the time of "simulation stops" as not been reached
% - particles can still move (not beached nor settled on the seafloor)
while t < tf - dt && isempty(u_can_move) == 0
    % age of the particle
    age = age + dts/86400; % age is in days
    % If microplastic debris
    if strcmpi(MODE,'MPD') == 1
        n_mpd_can_sink = find(age>AGES & sink==0); % older than age and should not be sinking already
        if ~isempty(n_mpd_can_sink)
            % randolmy choose PSINK percentage of these MPD and set them to sink = 1
            chgs = randi(length(n_mpd_can_sink), [round(length(n_mpd_can_sink)*PPSINK) 1]);
            sink(chgs) = 1;
            vs(chgs) = PVSINK; % set these partciles to sink at VSINK velocity
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
                if length(xto) == length(x)  % the amount of particles remains the same over that dt (no injection) so ID does not change
                    idp;
                    TRACK = [TRACK; x*0+t x y z idp];
                    xto = x;
                else   % there are new particles (injection), hence new IDs 
                    % There are n_newpart = length(x) - length(xto)  new particles
                    n_newpart = length(x) - length(xto); % amount of new particles
                    idp = [idp; (1:n_newpart)' + max(idp)];
                    TRACK = [TRACK; x*0+t x y z idp];
                    xto = x;
                end
            end
        end
    end
    % !!!!!!!!!!!!!!!!!!!!
    % injection every day
    % !!!!!!!!!!!!!!!!!!!!
    if t<tsf && (t==t0 || floor(t)-floor(t_prev)==1) % every day
        disp(['Injection: ' datestr(t,'dd/mm/yyyy'),' done on ',datestr(now,'dd/mm/yyyy HH:MM')])
        x = [x; x0];
        y = [y; y0];
        z = [z; z0];
        vs = [vs; v0];
        age = [age; x0*0]; % the new batch is aged 0
        sink = [sink; x0*0]; % the new batch is not sinking
        timetrack = [timetrack; x0*0 + t];  % these are INJECTION times
    end
    kk = kk + 1;
    if tnameold ~= tname
        % Prepare sub-region to load:
        % special case latitude
        if min(y)-mrg <-90
            lati = -89.8750; latf = max(y)+mrg;
            chglo = 1; % must change longitudes boundaries because particle moves on the other side of the earth
        elseif max(y)+mrg>90
            lati = min(y)-mrg; latf = 89.8750;
            chglo = 1; % must change longitudes boundaries because particle moves on the other side of the earth
        else
            lati = min(y)-mrg; latf = max(y)+mrg;
            chglo = 0;
        end
        % special case longitude
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
    % find the depth's index
    [~,ud] = histc(z,edgesdep);
    % Get the indices of the particle in the mesh in order to know the
    % velocity of the current (U and V components) where the particle is. 
    [Nla,lala] = histc(y,edgeslat);
    ula = lala' ;
    [Nlo,lolo] = histc(x,edgeslon);
    ulo = lolo';
    % case interval 89.8750 - 90
    iju0_avant = find(y>89.750);
    iju0_apres = find(y<-89.75);
    ula(iju0_avant) = length(edgeslat)-1;
    ula(iju0_apres) = 1; % if outside -0/+360
    % case interval 359.875 - 0.125
    ku0 = find(ulo == 0);
    iku0_avant = find(x>359.750);
    iku0_apres = find(x<0.125);
    ulo(iku0_avant) = max(ulo);
    ulo(iku0_apres) = 1;
    % Position in the matrix LON/LAT/DEPTH/TIME(=1)
    % U, V and W have the same format, hence same index for all three
    idx0 = sub2ind(size(U), ulo, ula, 1+ula*0); % index first slice of the matrix, for LAT --- 4st dim 1+ula*0
    %%
    idx = sub2ind(size(U), ulo, ula, ud');
    ue = U(idx);
    ve = V(idx);
    we = W(idx); % Do not add vs yet to check if velocity ECCO2 == 0
    vel_mag = sqrt(ue.^2 + ve.^2 + we.^2);
    % Check stop condition
    u_not_move = find(vel_mag' ==0 | z>edgesdep(end)-0.1); % velocity == 0: particle beached or down on the seafloor
    u_can_move = find(vel_mag ~=0); % velocity != 0: particle still in the current and moving
    % update coordinates
    x = x + ue'*dts./(pi/180*rad_earth*sind(90-LAT(idx0)'));
    y = y + ve'*dts./(pi/180*rad_earth*sind(90-LAT(idx0)'));
    z(u_can_move) = z(u_can_move) + vs(u_can_move)*dts + we(u_can_move)'*dts; % z-axis positive downward
    z(u_not_move) = z(u_not_move); % depth does not change for beached or settled particles
    z(z<0) = 0.001; % if particle goes above the sea surface (just safety check in case W fields at the surface are <0), then they are put back in the water
    z(z>edgesdep(end)) = edgesdep(end)-0.05; % if particle goes below the maximal depth defined by ECCO2, move it back to the seafloor to avoid histc error (no error with histcounts though)
    % special case -90/90 latitude
    x(y>90 | y<-90) = x(y>90 | y<-90)+180; % first change longitudes and then the latitude (or loose the info ><90)
    y(y>90)=180-y(y>90);
    y(y<-90)=-180-y(y<-90);
    % special case 0/360 longitude
    x(x<0)=360+(x(x<0));
    x(x>=360)=x(x>=360)-360;
    tnameold = tname ;
    t_prev = t; % memory this t
    t = t + dt;  % t with format datenum Matlab (days)
    [~,udt] = histc(t,tseccoDN);
    tname = tsecco_n(udt);
end
%% Wrap-up particule final positions
% COUNT
[xsortc,sxc] = sort(x);
[ysortc,syc] = sort(y);
% find indices in the mesh COUNT
Nlac = histc(y,edgeslat0);
ula0c = find(Nlac>0);
ula1c = [ula0c Nlac(ula0c)];
ulac_ok = [];
for uula = 1:size(ula1c,1)
    test_ulac = zeros(1,ula1c(uula,2))+ula1c(uula,1);
    ulac_ok = [ulac_ok test_ulac];
end
ulac_ok(syc) = ulac_ok; % keep particle order (no ascending order !!!)
ulac = ulac_ok;
Nloc = histc(x,edgeslon0);
ulo0c = find(Nloc>0);
ulo1c = [ulo0c Nloc(ulo0c)];
uloc_ok = [];
for uulo = 1:size(ulo1c,1)
    test_uloc = zeros(1,ulo1c(uulo,2))+ulo1c(uulo,1);
    uloc_ok = [uloc_ok test_uloc];
end
uloc_ok(sxc) = uloc_ok; % keep particle indices order (no ascending order !!!)
uloc = uloc_ok;
idxc = sub2ind(size(COUNT), uloc, ulac);
% there are "doubles" in idxc because particles can end up in the same
% voxel. Must loop in order to count all particles
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
% there are "doubles" in idxc because particles can end up in the same
% voxel. Must loop in order to count all particles
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
