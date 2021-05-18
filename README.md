# Laps
LAPS stands for Lagrangian Advection of Particles at Sea. LAPS simulates the trajectories of particlesat sea and is programmed with Matlab. It runs on any operating system. Use paths ``` /path/to/folder/ ``` for Linux/Mac and ``` C:\path\to\folder\ ``` for Windows. In it's most basic form, the only script to use is ```ADV_main('your_config_file.txt')```.

1. Get this Github repository.
2. Set it as your working directory in Matlab.
3. Edit your configuration file, one line per keyword (see keywords meanings and the example configuration file).
4. Download the sea velocities you need (see next setion).
5. Run ```ADV_main('your_config_file.txt')``` with Matlab.

### Sea current velocities and Stokes drift
LAPS reads the sea velocities provided by ECCO2 available through https://urs.earthdata.nasa.gov/ , you need to set a free account in order to download such data. You must at least have ECCO2 data. ECCO2 files have a 3-days resolution and there is one file per velocity component (East/North/Depth). Suppose you want to simulate some particle advection that occured between 01 May 2016 and 03 May 2016, then you need to get 3 files: 
```
- UVEL.1440x720x50.20160501.nc
- VVEL.1440x720x50.20160501.nc
- WVEL.1440x720x50.20160501.nc
```
Optionnaly, if you intend to account for the Stokes drift effects (the influence of the wind on the current velocities, mostly at the sea surface), then you must also get the Stokes drift velocities, specfically the WaveWatch 3 products availble on Ifremer's ftp. These files have a 3-hours resolution but contains all components (North/East) at monthly intervals. There is no depth component but LAPS takes care of that. So for that same simulation between 01 May 2016 and 03 May 2016, you need:
```
WW3-GLOB-30M_201605_uss.nc
```
Where do I get the files? At the time thie README is written:
- ECCOS2: https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/UVEL.nc/UVEL.1440x720x50.20160501.nc  (edit accordingly for VVEL and WVEL)
- Stokes drift: ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST/GLOBAL/2016_ECMWF/uss/WW3-GLOB-30M_201605_uss.nc  (edit accordingly for the year)

The Stokes drift velocity files must be converted to match the format of the ECCO2 files (it saves computation time). Use the code ```convSD_IFR2ECC_generic.m```, which must be edited to set your paths.

### Example of what should be written in a configuration file
```
TI 2018-01-03
TSI 2018-01-15
TF 2018-01-31
PECCO2                  /path/to/ECCO2/           # or C:\path\to\ECCO2\          for Windows 
INPUTF                  /path/to/input.xyz        # or C:\path\to\input.xyz       for Windows
OUTPUTP                 /path/to/result/folder/   # or C:\path\to\result\folder\  for Windows
PSD                     /path/to/Stokesdrift/     # or C:\path\to\Stokesdrift\    for Windows
TRK                     24
STOKESDRIFT             1
MODE                    MPD
#sedimentmode_SED (autoset to 0 if MDP)
GRZ                     1e-6
RHOP                    2500
#plasticmode_MicroPlasticDebris_MPD (autoset to 0 if SED)
PVSINK                  0.016
AGES                    10
PPSINK                  0.5
```

### Description of each keyword
Keyword | Description
------------ | -------------
TI		|Time when advection and injection start: YYYY-MM-DD
TSI		|Time when injection stops: YYYY-MM-DD
TF 		|Time when advection stops: YYYY-MM-DD
PECCO2		|path to ECCO2 files
INPUTF 		|path to the file of particle input (an ascii file with 4 columns X, Y, Z(def = 0), ID: X = lon, Y = lat [degree.decimal], Z = depth [m], ID = a 3-letters identification
OUTPUTP 	|path to result files
PSD 		|path to Stokes drift files
TRK 		|Record the path of the particles (not only their final position): 0 for no tracking or an integer meaning the tracking time step in hours
STOKESDRIFT 	|Activate Stokes drift: 0 or 1
MODE            |Advection of sediment particles (=SED) or micro plastic debris (=MPD)

Depending on the mode chosen, only the relevant parameters are used

If MODE is SED:
Keyword | Description
------------ | -------------
GRZ             |Particle size [m]
RHOP            |Particle density in [kg.m-3]

If MODE is MPD:
Keyword | Description
------------ | -------------
PVSINK          |The settling velocity of MDP (complicated to know, default value: 0.016 m.s-1 according to a bunch of papers: https://iopscience.iop.org/article/10.1088/1748-9326/aa8e8b/meta (mainly) https://www.sciencedirect.com/science/article/pii/S0269749116300264?via%3Dihub
AGES            |The age in [days] after which MDP might start sinking
PPSINK          |The probability [0=<p<=1] that the MDP older than AGES actually start sinking. (Old enough MDP are chossen randomly). The concept of MDP starting to setle with some propbility after a certain age threshold is inspired from https://www.sciencedirect.com/science/article/pii/S0025326X18301000
