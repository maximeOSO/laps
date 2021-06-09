function [R, PECCO2, TI, TF, TSI, INPUTF,...
    OUTPUTP, MODE, GRZ, RHOP, TRK, PSD,...
    STOKESDRIFT, AGES, PPSINK, PVSINK] = ADV_checkinput(configfile)

fid = fopen(configfile);
a = textscan(fid, '%s %s');
fclose(fid);

keyvar = a{1};
valvar = a{2};

id_pecco2 = find(strcmp(keyvar,'PECCO2') == 1);
id_ti = find(strcmpi(keyvar,'TI') == 1);
id_tf = find(strcmpi(keyvar,'TF') == 1);
id_tsi = find(strcmpi(keyvar,'TSI') == 1);
id_inpf = find(strcmpi(keyvar,'INPUTF') == 1);
id_outp = find(strcmpi(keyvar,'OUTPUTP') == 1);
id_mode = find(strcmpi(keyvar,'MODE') == 1);
id_grz = find(strcmpi(keyvar,'GRZ') == 1);
id_rhop = find(strcmpi(keyvar,'RHOP') == 1);
id_trk = find(strcmpi(keyvar,'TRK') == 1);
id_psd = find(strcmpi(keyvar,'PSD') == 1);
id_sd = find(strcmpi(keyvar,'STOKESDRIFT') == 1);
id_ages = find(strcmpi(keyvar,'AGES') == 1);
id_ppsink = find(strcmpi(keyvar,'PPSINK') == 1);
id_pvsink = find(strcmpi(keyvar,'PVSINK') == 1);

% Check if config file is complete
if isempty(id_pecco2) == 1; disp('Path to ecco2 files "PECCO2" is missing in the config file'); return; end
if isempty(id_ti) == 1; disp('Start time "TI" is missing in the config file'); return; end
if isempty(id_tsi) == 1; disp('Time when injection stops "TSI" is missing in the config file'); return; end
if isempty(id_tf) == 1; disp('Time when simulation stops "TF" is missing in the config file'); return; end
if isempty(id_inpf) == 1; disp('Path to input file "INPUTF" is missing in the config file'); return; end
if isempty(id_outp) == 1; disp('Path for output file "OUTPUTP" is missing in the config file'); return; end
if isempty(id_mode) == 1; disp('Indication whether particles are sediment (SED) or microplastic debris (MPD) is missing in the config file'); return; end
if isempty(id_trk) == 1; disp('Indication whether particle tracking should be recorded or not "TRK" is missing in the config file (set TRK=0 for no or TRK=an integer number of hours, the time setp of the tracking record, for yes)'); return; end
if isempty(id_sd) == 1; disp('Indication whether Stokes drift should be considered or not "SD" [0/1] is missing in the config file'); return; end

% also check the path , must end with / or \
PECCO2 = valvar{id_pecco2};
TI = datevec(datestr(valvar{id_ti},'yyyy-mm-dd HH:MM:SS'));
TF = datevec(datestr(valvar{id_tf},'yyyy-mm-dd HH:MM:SS'));
TSI = datevec(datestr(valvar{id_tsi},'yyyy-mm-dd'));
INPUTF = valvar{id_inpf};
OUTPUTP = valvar{id_outp};
MODE = valvar{id_mode};
if strcmpi(MODE,'SED')
    disp('Sediment mode, particle sinks / Stokes settling velocity')
    if isempty(id_grz) == 1; disp('Particles can sink --> needs to set the particle size "GRZ" [m] in the config file'); return; end
    if isempty(id_rhop) == 1; disp('Particles can sink --> needs to set the particle density "RHOP" in the config file'); return; end
    GRZ = str2num(valvar{id_grz});
    RHOP = str2num(valvar{id_rhop});
    AGES = 0.; % not used
    PPSINK = 0.; % not used
    PVSINK = 0.; % not used
else
    disp('Plastic mode, particle sinks following user-defined settling velocity')
    if isempty(id_ages) == 1; disp('Particles can sink --> needs to set the ages "AGES" [days] in the config file'); return; end
    if isempty(id_ppsink) == 1; disp('Particles can sink --> needs to set the particle sinking probability "PPSINK" in the config file'); return; end
    if isempty(id_pvsink) == 1; disp('Particles can sink --> needs to set the particle setling velocity "PVSINK" [m.s-1] in the config file'); return; end   
    GRZ = 0; % not used
    RHOP = 0; % not used
    AGES = str2num(valvar{id_ages});
    PPSINK = str2num(valvar{id_ppsink});
    PVSINK = str2num(valvar{id_pvsink});
end
TRK = str2num(valvar{id_trk});
if TRK > 0
    disp(['Particle tracking: every ',num2str(TRK),' hours'])
else
    disp('No particle tracking')
end
STOKESDRIFT = str2num(valvar{id_sd});
if STOKESDRIFT == 1
    disp('With Stokes drift')
    if isempty(id_psd) == 1; disp('Stokes drift considered --> needs to set the path to Stokes drift velocity files "PSD" in the config file'); return; end
    PSD = valvar{id_psd};
elseif STOKESDRIFT == 0
    disp('No Stokes drift')
    PSD = 'NULL';
else
    disp('ERROR: STOKESDRIFT must be 0 or 1');
end

%% Check input
R = 0;
if ischar(PECCO2) == 0
    disp('ERROR: invalid path to ECCO2'); R = R + 1;
%    return
end
if ischar(PECCO2) == 1
    ECCO2FolderInfo = dir(PECCO2);
    if length(ECCO2FolderInfo) == 2 % only . and ..
        disp(['ERROR: no files in ', PECCO2]); R = R + 1;
        disp('Get them here --> ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/')
        %        return
    elseif isempty(ECCO2FolderInfo)
        disp(['ERROR: no directory ', PECCO2]); R = R + 1;
        %        return
    end
end
if ischar(PSD) == 0
    disp('ERROR: invalid path to Stokes drift data'); R = R + 1;
%    return
end
if ischar(INPUTF) == 0
    disp('ERROR: invalid path to particule input file'); R = R + 1;
%    return
end
if ischar(OUTPUTP) == 0
    disp('ERROR: invalid path to results'); R = R + 1;
%    return
end
if ischar(MODE) == 0 || strcmpi(MODE,'SED') == 0 && strcmpi(MODE,'MPD') == 0
    disp('ERROR: MODE must be "SED" (sediment) or "MPD" (microplastic debris) '); R = R + 1;
%    return
end
if TRK <0    
    disp('ERROR: TRK must be positive, in hours, can be < 1, eg 0.25 = 15 min tracking (0 if no tracking) '); R = R + 1;
%    return
end
if size(TI,1) ~= 1 || size(TI,2) ~= 6
    disp('ERROR: TI must be a matrix [yyyy mm dd HH MM SS]'); R = R + 1;
%    return
end
if size(TF,1) ~= 1 || size(TF,2) ~= 6
    disp('ERROR: TF must be a matrix [yyyy mm dd HH MM SS]'); R = R + 1;
%    return
end
if size(TSI,1) ~= 1 || size(TSI,2) ~= 6
    disp('ERROR: TSI must be a matrix [yyyy mm dd HH MM SS]'); R = R + 1;
%    return
end
if datenum(TI) > datenum(TF)
    disp('ERROR: TI must be lower than TF'); R = R + 1;
%    return
end
if datenum(TI) > datenum(TSI)
    disp('ERROR: TI must be lower than TSI'); R = R + 1;
%    return
end
if datenum(TSI) > datenum(TF)
    disp('ERROR: TSI must be lower than TF'); R = R + 1;
%    return
end
if isfloat(RHOP) == 0
    disp('ERROR: RHOP must be a float [kg/m3]'); R = R + 1;
%    return
end
if isfloat(GRZ) == 0
    disp('ERROR: GRZ must be a float [meter]'); R = R + 1;
%    return
end
if isfloat(AGES) == 0
    disp('ERROR: AGES must be a float [day unit]'); R = R + 1;
%    return
end
if isfloat(PPSINK) == 0
    disp('ERROR: PPSINK must be a float [0<=PPSINK<=1]'); R = R + 1;
%    return
end
if isfloat(PVSINK) == 0
    disp('ERROR: PVSINK must be a float [m.s-1]'); R = R + 1;
%    return
end
end