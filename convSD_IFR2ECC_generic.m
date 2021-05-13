clear
close all
format long
clc

%%% Read netcdf file
% This function transform the netcdf files of the Stokes drift velocity
% fields into netcdf files at the same format and spatial resolution of the
% ECCO2 files.
% it requires the gridfit package that must be downloaded here:
% https://se.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit

dateECCO2 = '20180103'  % can be any date in the name of your ECCO2 files
path2ecco2 = '/your/path/to/ECCO2/'
path2stokesdrift = '/your/path/to/Sokesdrift/'
path2gridfitdir = '/your/path/to/gridfitdir'
year_to_convert = [2015 2018 2019]  % example [2010 2011 2014], the code will run over all months (1:12). If you only need some months, then adjust the jj variable (second for-loop below)
addpath(path2gridfitdir) 

%% ECCO2
ecU = [path2ecco2, 'UVEL.1440x720x50.',dateECCO2,'.nc'];
ecV = [path2ecco2, 'VVEL.1440x720x50.',dateECCO2,'.nc'];
ecW = [path2ecco2, 'WVEL.1440x720x50.',dateECCO2,'.nc'];
elat = double(ncread(ecU,'LATITUDE_T'));
elon = double(ncread(ecU,'LONGITUDE_T'));
eums = double(ncread(ecU,'UVEL'));
evms = double(ncread(ecV,'VVEL'));
ewms = double(ncread(ecW,'WVEL'));
vcont = 9.9e22;
eams = sqrt(eums.^2 + evms.^2 + ewms.^2);
eams1 = eams(:,:,1,1);
eams1(eams1>vcont) = NaN;
eresll = elon(2)-elon(1);
edgeselat = min(elat)-eresll/2:eresll:max(elat)+eresll/2;
edgeselon = min(elon)-eresll/2:eresll:max(elon)+eresll/2;
[ELAT, ELON] = meshgrid(elat,elon);

%% SD
for ii = year_to_convert
    for jj = 1:12
        sd = [path2stokesdrift, 'WW3-GLOB-30M_',num2str(ii*100+jj),'_uss.nc'];
        slat = double(ncread(sd,'latitude'));
        slon = double(ncread(sd,'longitude'));
        sums = double(ncread(sd,'uuss'));
        svms = double(ncread(sd,'vuss'));
        sresll = slon(2)-slon(1);
        edgesslat = min(slat)-sresll/2:sresll:max(slat)+sresll/2;
        edgesslon = min(slon)-sresll/2:sresll:max(slon)+sresll/2;
        [SLAT, SLON] = meshgrid(slat,slon);
        %% SD -180:180 --> 0:360
        % Attention at 360, make an overlap on that border (double the last line) or 
        % the spatial interpolation will return NAN there
        ui = find(slon(slon<0),1,'first');
        uf = find(slon(slon<0),1,'last');
        % UUSS
        sumsNEG = cat(1,sums(ui:uf,:,:),sums(uf,:,:)); % double the last line to make an overlap at 360
        sumsPOS = sums(uf+1:end,:,:);
        sums360 = cat(1,sumsPOS,sumsNEG);
        % VUSS
        svmsNEG = cat(1,svms(ui:uf,:,:),svms(uf,:,:)); % double the last line to make an overlap at 360
        svmsPOS = svms(uf+1:end,:,:);
        svms360 = cat(1,svmsPOS,svmsNEG);
        % new coord
        SLON360 = [SLON(uf+1:end,:);SLON(ui:uf,:)+360;SLON(uf,:)*0+360];
        SLAT360 = [SLAT;SLAT(end,:)];
        %% Reformat SD at ECCO2's standard
        load extrapol_StokesDriftMap; % island1 & ue_fill
        SD_U = zeros(length(elon),length(elat),size(sums360,3));
        SD_V = zeros(length(elon),length(elat),size(sums360,3));
        for kk = 1:size(sums360,3)
            disp([num2str(ii),'-',num2str(jj,'%02d'),'-',num2str(kk,'%03d'),'/',num2str(size(sums360,3))])
            sums1 = squeeze(sums360(:,:,kk));
            sums1R = interp2(SLAT360,SLON360,sums1,ELAT,ELON);
            svms1 = squeeze(svms360(:,:,kk));
            svms1R = interp2(SLAT360,SLON360,svms1,ELAT,ELON);
            %% Coast needs extrapolation
            % extrapolate everywhere
            all_extrap_U = gridfit(ELON,ELAT,sums1R,elon,elat);
            all_extrap_V = gridfit(ELON,ELAT,svms1R,elon,elat);
            % reformat at sums1R
            all_extrap_U = all_extrap_U';
            all_extrap_V = all_extrap_V';
            % keep extrapolated values only on ue_fill
            sel_extrap_U = all_extrap_U(island1(ue_fill));
            sums1R(island1(ue_fill)) = sel_extrap_U;
            sel_extrap_V = all_extrap_V(island1(ue_fill));
            svms1R(island1(ue_fill)) = sel_extrap_V;
            % save
            SD_U(:,:,kk) = sums1R;
            SD_V(:,:,kk) = svms1R;
        end

        sf = 1e-4;
        fv = -32767;
        disp([path2stokesdrift, '/fmt_SD_',num2str(ii*100+jj),'_1440x720xtime3h.nc'])
        ncid = netcdf.create([path2stokesdrift, '/fmt_SD_',num2str(ii*100+jj),'_1440x720xtime3h.nc'],'NETCDF4');
        londimid = netcdf.defDim(ncid,'lon',size(SD_U,1));
        latdimid = netcdf.defDim(ncid,'lat',size(SD_V,2));
        timdimid = netcdf.defDim(ncid,'t3h',size(SD_V,3));
        var_uuss = netcdf.defVar(ncid,'UUSS','NC_SHORT',[londimid latdimid timdimid]);
        var_vuss = netcdf.defVar(ncid,'VUSS','NC_SHORT',[londimid latdimid timdimid]);
        netcdf.defVarDeflate(ncid,var_uuss,true,true,5);
        netcdf.defVarDeflate(ncid,var_vuss,true,true,5);
        netcdf.endDef(ncid);
        Utosave = SD_U;
        Vtosave = SD_V;
        Utosave(isnan(Utosave)) = fv;
        Vtosave(isnan(Vtosave)) = fv;
        Utosave = SD_U/sf;
        Vtosave = SD_V/sf;
        netcdf.putVar(ncid,var_uuss,Utosave);
        netcdf.putVar(ncid,var_vuss,Vtosave);
        % Re-enter define mode.
        netcdf.reDef(ncid);
        % Create an attribute associated with the variable.
        netcdf.putAtt(ncid,var_uuss,'scale_factor',sf);
        netcdf.putAtt(ncid,var_uuss,'_Fillvalue',fv);
        netcdf.putAtt(ncid,var_uuss,'valid_min',-9900);
        netcdf.putAtt(ncid,var_uuss,'valid_max',9900);
        netcdf.putAtt(ncid,var_uuss,'add_offset',0.0);
        
        netcdf.putAtt(ncid,var_vuss,'scale_factor',sf);
        netcdf.putAtt(ncid,var_vuss,'_Fillvalue',fv);
        netcdf.putAtt(ncid,var_vuss,'valid_min',-9900);
        netcdf.putAtt(ncid,var_vuss,'valid_max',9900);
        netcdf.putAtt(ncid,var_vuss,'add_offset',0.0);
        netcdf.close(ncid);
        
    end
end
