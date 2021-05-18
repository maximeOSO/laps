function ADV_main(configfile)
%%% Main function of the LAPS program

%%% INPUT: see configfile

%%% OUTPUT:
%%% - adv_grd: full path to the matrix formatted data:
%%%             _ LON0 : mesh lon
%%%             _ LAT0 : mesh lat
%%%             _ dep : depth limits
%%%             _ COUNT3D: 3D matrix (LON, LAT, DEP), number of particles
%%% - adv_trk: full path to the table with initial time and position of the
%%% advected particles and of their final XYZ position


disp(" ************************************** ")
disp(" *               LAPS                 * ")
disp(" ************************************** ")


%% Check and Read user's config file
disp(['configuration file: ', configfile])
[R, PECCO2, TI, TF, TSI, INPUTF,...
    OUTPUTP, MODE, GRZ, RHOP, TRK, PSD,...
    STOKESDRIFT, AGES, PPSINK, PVSINK] = ADV_checkinput(configfile);
if R == 1
    error('Program stopped due to 1 error in inputs')
elseif R > 1
    error(['Program stopped due to ',num2str(R),' errors in inputs'])
else
    disp('Inputs OK')
end

copyfile(configfile, OUTPUTP)

%% Check files ECCO2
typedate = 'CHK';
ADV_tsfromecco2files(PECCO2, TI, TF, typedate);

%% Main call
tic
switch STOKESDRIFT
    case 0
        [adv_grd, adv_trkHR] = ADV_SD_NO(configfile, PECCO2, PSD, TI, TSI, TF, ...
            INPUTF, OUTPUTP, MODE, GRZ, RHOP, TRK, AGES, PPSINK, PVSINK);
    case 1
        [adv_grd, adv_trkHR] = ADV_SD_YES(configfile, PECCO2, PSD, TI, TSI, ...
            TF, INPUTF, OUTPUTP, MODE, GRZ, RHOP, TRK, AGES, PPSINK, PVSINK);
end
toc


mkfig = 1;
if mkfig == 1
    disp('Making figure...')
    set(figure,'Position',[0 0 1200 800])
    set(gcf,'PaperPositionMode','auto')
    % load input
    fidpart = fopen(INPUTF);
    inputparticles = textscan(fidpart, '%f %f %f %s');
    fclose(fidpart);
    x0 = inputparticles{1};
    y0 = inputparticles{2};
    z0 = inputparticles{3};
    label_orig = inputparticles{4};
    % load MAT
    load(adv_grd)
    sC = sum(COUNT3D,3);
    sC(sC == 0) = nan;
    % make map
    worldmap world
    p1 = pcolorm(LAT0, LON0, sC);
    hold on
    partexistid = find(~isnan(sC));
    load coastlines
    [latcells, loncells] = polysplit(coastlat, coastlon);
    p2 = plotm(coastlat, coastlon,'k');
    p3 = plotm(LAT0(partexistid),LON0(partexistid), '.');
    cb = colorbar;  
    cb.Label.String = 'Number of particles';
    if TRK >= 0
        % load TRK
        fidpart = fopen(adv_trkHR);
        trk = textscan(fidpart, '%d %d %f %f %f %s');
        fclose(fidpart);
        dd = trk{1};
        tt = trk{2};
        xtrk = trk{3};
        ytrk = trk{4};
        ztrk = trk{5};
        ltrk = trk{6};
        p4 = plotm(ytrk,xtrk,'k.');
    end
    p5 = plotm(y0,x0,'ro','MarkerFaceColor','r');
    legend([p5 p4 p1], {'initial position', 'track', 'final position'})
end
end
