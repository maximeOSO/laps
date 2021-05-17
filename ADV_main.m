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
    % STOKESDRIFT, AGES, PPSINK, PVSINK] = ADV_checkinput('ADV_config.txt');
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


mkfig = 0;
if mkfig == 1
disp('Making figure...')
    set(figure,'Position',[0 0 1200 1200])
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
    pcolorm(LAT0, LON0, sC)
    hold on
    partexistid = find(~isnan(sC));   
    load coastlines
    [latcells, loncells] = polysplit(coastlat, coastlon);
    plotm(coastlat, coastlon,'k')
    plotm( LAT0(partexistid),LON0(partexistid), '.')
    colorbar
    %patchm(coastlat, coastlon,[0.5 0.5 0.5])
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
        plotm(ytrk,xtrk,'k.')       
    end
    plotm(y0,x0,'ro','MarkerFaceColor','r')
end









mkfig = 0;
if mkfig == 1
    
    
    %% reload and process results
    load(adv_grd)
    % load(adv_trkHR)
    sC = sum(COUNT3D,3);
    sC(sC == 0) = nan;
    
    % ix = find(isnan(sC)==0);
    % tosave = [LON0(ix) LAT0(ix) sC(ix)];
    % save('/home/maxime/Desktop/densite.txt', 'tosave', '-ASCII')
    %
    % fid = fopen('/home/maxime/Desktop/densite2.txt', 'w');
    % for k = 1:length(ix)
    %     fprintf(fid,'%f \t %f \t %f \n', tosave(k,:));
    % end
    % fclose(fid);
    
    % % si tu veux aussi les prof
    % ix2 = find(COUNT3D>0);
    % [i1, i2, i3] = ind2sub(size(COUNT3D), ix2);
    
    %% format output to text files
    disp('Convert *.mat to *.txt files ...')
    ttres = 24; %% modif: time resolution of the track files in counts of dt
    % tsrec = 1/dt * ttres/24; %% modif
    if TRK>0
        for ii = 1:size(adv_trkHR,1)
            load(adv_trkHR(ii,:))
            fileName = [adv_trkHR(ii,:),'_FMT.txt'];
            fid = fopen(fileName, 'w');
            header = ['particle-ID \t date-time \t',...
                'Longitude \t Latitude \t Depth'];
            fprintf(fid,[header,'\n']);
            npart = size(TRACK,1);
            ntime = size(TRACK,3);
            % INIT: nombre de particule * nombre de temps
            ct = 1;
            tosave = zeros(npart*round(ntime/tsrec),5); %% modif
            for kk = 1:npart
                t = squeeze(TRACK(kk,1,1:tsrec:end)); %% modif
                lon = squeeze(TRACK(kk,2,1:tsrec:end)); %% modif
                lat = squeeze(TRACK(kk,3,1:tsrec:end)); %% modif
                z = squeeze(TRACK(kk,4,1:tsrec:end)); %% modif
                ID = ii*1e6 + kk + 0*t;
                % final format : ID year month day hour lon lat z,
                tosave(ct:ct+length(t)-1,:) = [ID t lon lat z]; %% modif
                ct = ct + length(t); %% modif
            end
            dateFMT0_str = datestr(tosave(:,2),'yyyy mm dd HH MM SS');
            dateFMT0 = str2num(dateFMT0_str);
            tosaveFMT = [tosave(:,1) dateFMT0 tosave(:,3:5)];
            fprintf(fid,['%d \t %d %02d %02d %02d %02d %02d'...
                '\t %.6f \t %.6f \t %.3f \n'],tosaveFMT');
            fclose(fid);
        end
    end
    
    %% Map input check
    % [LATmap ,LONmap ,velmap] = input_map_visu(PECCO2, TI, INPUTF);
    
    %% simple visu
    disp('Making figure...')
    set(figure,'Position',[0 0 1200 1200])
    set(gcf,'PaperPositionMode','auto')
    fidpart = fopen(INPUTF);
    inputparticles = textscan(fidpart, '%f %f %f %s');
    fclose(fidpart);
    x0 = inputparticles{1};
    y0 = inputparticles{2};
    z0 = inputparticles{3};
    label_orig = inputparticles{4};
    worldmap world
    %pcolorm(LATmap,LONmap,velmap)
    pcolorm(LAT0,LON0,sC)
    hold on
    load coastlines
    [latcells, loncells] = polysplit(coastlat, coastlon);
    plotm(coastlat, coastlon,'k')
    colorbar
    %patchm(coastlat, coastlon,[0.5 0.5 0.5])
    if TRK > 0
        for ii = 1:size(adv_trkHR,1)
            load(adv_trkHR(ii,:))
            for kk = 1:size(TRACK,1)
                plotm(squeeze(TRACK(kk,3,:)),squeeze(TRACK(kk,2,:)))
                hold on
            end
        end
    end
    plotm(y0,x0,'ro','MarkerFaceColor','r')
    % dl = 5;
    % setm(gca,'MapLonLimit',[min(x0)-dl max(x0)+dl],'MapLatLimit',[min(y0)-dl max(y0)+dl])
end
end
