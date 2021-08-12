clear
close all
clc

% TASK = 1; % prepare 4 inputs
% TASK = 2; % edit netcdf
% TASK = 3; % check synth files
% TASK = 4; % check LAPS results

TASK = 4;

switch TASK
    case 1 % create synthetic velcoty fields from ECCO2 and WWIII
        eU1 = 'v_orig/UVEL.1440x720x50.20180220.nc';

        %% general
        lat = double(ncread(eU1,'LATITUDE_T'));
        lon = double(ncread(eU1,'LONGITUDE_T'));
        dep = double(ncread(eU1,'DEPTH_T'));
        
        uvel1 = double(ncread(eU1,'UVEL'));
        uvel1(uvel1==0) = nan;
        idnan = find(isnan(uvel1));
        
        uvel1(   1:360,  :, :) = 1;
        uvel1( 361:720,  :, :) = -1;
        uvel1( 721:1080, :, :) = 0.5;
        uvel1(1081:1440, :, :) = 0;
        
        uvel1(idnan) = nan;

        imagesc(lon, lat, uvel1(:,:,1)')
        colorbar
        caxis([-2 2])
        hold on
        plot(70, -20, 'k+')
        plot(100, -20, 'k+')
        plot(220, -20, 'k+')
        plot(340, -20, 'k+')
        
    case 2 % edit netcdf - ncwrite
        eU1 = 'v_synth/UVEL.1440x720x50.20180217.nc';
        eV1 = 'v_synth/VVEL.1440x720x50.20180217.nc';
        eW1 = 'v_synth/WVEL.1440x720x50.20180217.nc';
        
        eU2 = 'v_synth/UVEL.1440x720x50.20180220.nc';
        eV2 = 'v_synth/VVEL.1440x720x50.20180220.nc';
        eW2 = 'v_synth/WVEL.1440x720x50.20180220.nc';
        
        sd = 'v_synth/fmt_SD_201802_1440x720xtime3h.nc';
        
        %% general
        lat = double(ncread(eU1,'LATITUDE_T'));
        lon = double(ncread(eU1,'LONGITUDE_T'));
        dep = double(ncread(eU1,'DEPTH_T'));
        
        uvel1 = double(ncread(eU1,'UVEL'));
        vvel1 = double(ncread(eV1,'VVEL'));
        wvel1 = double(ncread(eW1,'WVEL'));
        
        uvel2 = double(ncread(eU2,'UVEL'));
        vvel2 = double(ncread(eV2,'VVEL'));
        wvel2 = double(ncread(eW2,'WVEL'));
        
        uuss = double(ncread(sd,'UUSS'));
        vuss = double(ncread(sd,'VUSS'));
        
        %% modif 1
        uvel1(   1:360,  :, :) = 1;
        uvel1( 361:720,  :, :) = -1;
        uvel1( 721:1080, :, :) = 0;
        uvel1(1081:1440, :, :) = 0;
        
        vvel1(   1:360,  :, :) = 0;
        vvel1( 361:720,  :, :) = 0;
        vvel1( 721:1080, :, :) = -1;
        vvel1(1081:1440, :, :) = 1;
        
        wvel1( :, :, 1:10) = 0;
        wvel1( :, :, 11:end) = 1;
               
        ncwrite(eU1, 'UVEL', uvel1)
        ncwrite(eV1, 'VVEL', vvel1)
        ncwrite(eW1, 'WVEL', wvel1)
        
        %% modif 2
        uvel2(   1:360,  :, :) = 0;
        uvel2( 361:720,  :, :) = 0;
        uvel2( 721:1080, :, :) = 1;
        uvel2(1081:1440, :, :) = -1;
        
        vvel2(   1:360,  :, :) = -1;
        vvel2( 361:720,  :, :) = 1;
        vvel2( 721:1080, :, :) = 0;
        vvel2(1081:1440, :, :) = 0;
        
        wvel2( :, :, 1:10) = 0;
        wvel2( :, :, 11:end) = 1;
        
        ncwrite(eU2, 'UVEL', uvel2)
        ncwrite(eV2, 'VVEL', vvel2)
        ncwrite(eW2, 'WVEL', wvel2)
        
        %% SD
        uuss(   1:360,  :, :) = 0;
        uuss( 361:720,  :, :) = 0;
        uuss( 721:1080, :, :) = 0.5;
        uuss(1081:1440, :, :) = 0.5;
        
        vuss(   1:360,  :, :) = 0.5;
        vuss( 361:720,  :, :) = 0.5;
        vuss( 721:1080, :, :) = 0;
        vuss(1081:1440, :, :) = 0;
        
        ncwrite(sd, 'UUSS', uuss)
        ncwrite(sd, 'VUSS', vuss)        
        
        
    case 3 % check synth netcdf
        eU1 = 'v_synth/UVEL.1440x720x50.20180217.nc';
        eV1 = 'v_synth/VVEL.1440x720x50.20180217.nc';
        eW1 = 'v_synth/WVEL.1440x720x50.20180217.nc';
        
        eU2 = 'v_synth/UVEL.1440x720x50.20180220.nc';
        eV2 = 'v_synth/VVEL.1440x720x50.20180220.nc';
        eW2 = 'v_synth/WVEL.1440x720x50.20180220.nc';
        
        sd = 'v_synth/fmt_SD_201802_1440x720xtime3h.nc';
        
        
        %% general
        lat = double(ncread(eU1,'LATITUDE_T'));
        lon = double(ncread(eU1,'LONGITUDE_T'));
        dep = double(ncread(eU1,'DEPTH_T'));
        
        uvel1 = double(ncread(eU1,'UVEL'));
        vvel1 = double(ncread(eV1,'VVEL'));
        
        uvel2 = double(ncread(eU2,'UVEL'));
        vvel2 = double(ncread(eV2,'VVEL'));
        
        uuss = double(ncread(sd,'UUSS'));
        vuss = double(ncread(sd,'VUSS'));
        
        subplot(3,2,1)
        imagesc(lon, lat, uvel1(:,:,1)')
        colorbar
        subplot(3,2,2)
        imagesc(lon, lat, vvel1(:,:,1)')
        colorbar
        subplot(3,2,3)
        imagesc(lon, lat, uvel2(:,:,1)')
        colorbar
        subplot(3,2,4)
        imagesc(lon, lat, vvel2(:,:,1)')
        colorbar
        
        subplot(3,2,5)
        imagesc(lon, lat, uuss(:,:,1)')
        colorbar
        subplot(3,2,6)
        imagesc(lon, lat, vuss(:,:,1)')
        colorbar
        
    case 4 % check LAPS results no SD
        adv_trkHR = 'res_reference/TRK_20180217_20180217_20180222_MPD_TMP_cfg_verif.txt'; 
        INPUTF = 'inp.txt';
        
        % load input
        fidpart = fopen(INPUTF);
        inputparticles = textscan(fidpart, '%f %f %f %s');
        fclose(fidpart);
        x0 = inputparticles{1};
        y0 = inputparticles{2};
        z0 = inputparticles{3};
        label_orig = inputparticles{4};

        fidpart = fopen(adv_trkHR);
        trk = textscan(fidpart, '%d %d %f %f %f %s');
        fclose(fidpart);
        dd = trk{1};
        tt = trk{2};
        xtrk = trk{3};
        ytrk = trk{4};
        ztrk = trk{5};
        ltrk = trk{6};
        
        disp('Making figure...')
        set(figure,'Position',[0 0 1200 400])
        set(gcf,'PaperPositionMode','auto')
        mrg = 3;
        uvel1 = [1 -1 0 0]; % see uvel1 in TASK 2
        vvel1 = [0 0 -1 1]; % see vvel1 in TASK 2
        uvel2 = [0 0 1 -1]; % see uvel2 in TASK 2
        vvel2 = [-1 1 0 0]; % see vvel2 in TASK 2
        offset = 4;
        sc = 1;
        for kk = 1:4
            subplot(1,5,kk)
            hold on
            p5 = plot(x0,y0,'ro','MarkerFaceColor','r','DisplayName','initial position');
            p4 = plot(xtrk,ytrk,'k.','DisplayName','track');
            
            headWidth = 8;
            headLength = 8;
            ah1 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
            set(ah1,'parent',gca);
            set(ah1,'position',[x0(kk),y0(kk)+offset,uvel1(kk),vvel1(kk)]);
            
            ah2 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
            set(ah2,'parent',gca);
            set(ah2,'position',[x0(kk)+uvel1(kk),y0(kk)+vvel1(kk)+offset,uvel2(kk),vvel2(kk)]);
            
            axis equal
            xlim([x0(kk)-mrg, x0(kk)+mrg])
            ylim([y0(kk)-mrg, y0(kk)+mrg+offset])
            
            legend('Initial position', 'trajectory')
            
            npart = 5;
            id_part = kk:npart:length(xtrk);
            disp(cell2mat(trk{6}(id_part(1))))
            
            distx = (max(xtrk(id_part)) - min(xtrk(id_part)))*pi/180*6371000*sind(90-mean(ytrk(id_part))); % geo deg to m
            disty = (max(ytrk(id_part)) - min(ytrk(id_part)))*pi/180*6371000*sind(90-mean(ytrk(id_part))); % geo deg to m
            
            distx1_theor = abs(uvel1(kk))*3600*(24*3)/1000; % 3 days travel in sec x velocity [m/s]
            disty1_theor = abs(vvel1(kk))*3600*(24*3)/1000;
            
            distx2_theor = abs(uvel2(kk))*3600*(24*2-1)/1000; % 2 days travel in sec x velocity [m/s]
            disty2_theor = abs(vvel2(kk))*3600*(24*2-1)/1000;
            
            disp(['X Velocity = 1 m/s. Simulated distance = ', num2str(round(distx/1000),'%d'),' km, theoretical distance = ', num2str(round(max(distx1_theor, distx2_theor)),'%d'),' km'])
            disp(['Y Velocity = 1 m/s. Simulated distance = ', num2str(round(disty/1000),'%d'),' km, theoretical distance = ', num2str(round(max(disty1_theor, disty2_theor)),'%d'),' km'])
            
            if kk == 4
                disp(['The simulated distance implies conversion from degrees to ' ...
                'km which is function of latitude.'])
                disp(['We use only the average of the '... '
                'particle latitude, hence a introducing a slight bias in the '...
                'distance computation, which should yet remain small since particles ' ...
                'do not travel far for in this sverification'])
            end
        end
        subplot(1,5,5)
        kk = 5;
        hold on
        p5 = plot(x0,y0,'ro','MarkerFaceColor','r','DisplayName','initial position');
        p4 = plot(xtrk,ytrk,'k.','DisplayName','track');
        plot([90, 90],[-13, -23],'k--')
        
        ah1 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah1,'parent',gca);
        set(ah1,'position',[x0(kk)-2,y0(kk)+offset,uvel1(1),vvel1(1)]);
        
        ah2 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah2,'parent',gca);
        set(ah2,'position',[x0(kk)+uvel1(1)-2,y0(kk)+vvel1(1)+offset,uvel2(1),vvel2(1)]);
        
        ah3 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah3,'parent',gca);
        set(ah3,'position',[x0(kk)+1,y0(kk)+offset,uvel1(2),vvel1(2)]);
        
        ah4 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah4,'parent',gca);
        set(ah4,'position',[x0(kk)+uvel1(2)+1,y0(kk)+vvel1(2)+offset,uvel2(2),vvel2(2)]);
        
        axis equal
        xlim([x0(kk)-mrg, x0(kk)+mrg])
        ylim([y0(kk)-mrg, y0(kk)+mrg+offset])
        legend('Initial position', 'trajectory')
               
        case 5 % check LAPS results with SD
        adv_trkHR = 'res_reference/TRK_20180217_20180217_20180222_MPD_TMP_cfg_verif_SD.txt'; 
        INPUTF = 'inp.txt';
        
        % load input
        fidpart = fopen(INPUTF);
        inputparticles = textscan(fidpart, '%f %f %f %s');
        fclose(fidpart);
        x0 = inputparticles{1};
        y0 = inputparticles{2};
        z0 = inputparticles{3};
        label_orig = inputparticles{4};

        fidpart = fopen(adv_trkHR);
        trk = textscan(fidpart, '%d %d %f %f %f %s');
        fclose(fidpart);
        dd = trk{1};
        tt = trk{2};
        xtrk = trk{3};
        ytrk = trk{4};
        ztrk = trk{5};
        ltrk = trk{6};
        
        disp('Making figure...')
        set(figure,'Position',[0 0 1200 400])
        set(gcf,'PaperPositionMode','auto')
        mrg = 4;
        uvel1 = [1 -1 0 0]*2 + [0 0 0.5 0.5]*2; % see uvel1 + uuss in TASK 2
        vvel1 = [0 0 -1 1]*2 + [0.5 0.5 0 0]*2; % see vvel1 + vvss in TASK 2
        uvel2 = [0 0 1 -1]*2 + [0 0 0.5 0.5]*2; % see uvel2 + uuss in TASK 2
        vvel2 = [-1 1 0 0]*2 + [0.5 0.5 0 0]*2; % see vvel2 + vvss in TASK 2
        offset = 3;
        sc = 1;
        for kk = 1:4
            subplot(1,4,kk)
            hold on
            p5 = plot(x0,y0,'ro','MarkerFaceColor','r','DisplayName','initial position');
            p4 = plot(xtrk,ytrk,'k.','DisplayName','track');
            
            headWidth = 8;
            headLength = 8;
            ah1 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
            set(ah1,'parent',gca);
            set(ah1,'position',[x0(kk),y0(kk)+offset,uvel1(kk),vvel1(kk)]);
            
            ah2 = annotation('arrow', 'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
            set(ah2,'parent',gca);
            set(ah2,'position',[x0(kk)+uvel1(kk),y0(kk)+vvel1(kk)+offset,uvel2(kk),vvel2(kk)]);
            
            axis equal
            xlim([x0(kk)-mrg, x0(kk)+mrg])
            ylim([y0(kk)-mrg, y0(kk)+mrg+offset])
            
            legend('Initial position', 'trajectory','location','south')
        end
        
end