# Laps

1. Put all the Matlab files (.m, .mat) in the same directory.
2. Set this directory as your working directory in Matlab.
3. Edit your configuration file, one line per keyword (see keywords meanings and the example configuration file).
4. Run ```ADV_main('your_config_file.txt')``` with Matlab
## Description of each keyword

Keyword | Description
------------ | -------------
PECCO2		|path to ECCO2 files
TI		|Time when advection and injection start: YYYY-MM-DD
TSI		|Time when injection stops: YYYY-MM-DD
TF 		|Time when advection stops: YYYY-MM-DD
INPUTF 		|path to the file of particle input (an ascii file with 4 columns X, Y, Z(def = 0), ID: X = lon, Y = lat [degree.decimal], Z = depth [m], ID = a 3-letters identification
OUTPUTP 	|path to result files
TRK 		|Record the path of the particles (not only their final position): 0 for no tracking or an integer meaning the tracking time step in hours
PSD 		|path to Stokes drift files
STOKESDRIFT 	|Activate Stokes drift (wave action <-- wind): 0 or 1
MODE            |Advection of sediment particles (=SED) or micro plastic debris (=MPD)

Depending on the mode chosen, only the relevant parameters are used
### If MOD is SED:
Keyword | Description
------------ | -------------
GRZ             |Particle size [m]
RHOP            |Particle density in [kg.m-3]
### If MOD is MPD:
Keyword | Description
------------ | -------------
PVSINK          |The settling velocity of MDP (complicated to know, default value: 0.016 m.s-1 according to a bunch of papers: https://iopscience.iop.org/article/10.1088/1748-9326/aa8e8b/meta (mainly) https://www.sciencedirect.com/science/article/pii/S0269749116300264?via%3Dihub
AGES            |The age in [days] after which MDP might start sinking
PPSINK          |The probability [0=<p<=1] that the MDP older than AGES actually start sinking. (Old enough MDP are chossen randomly). The concept of MDP starting to setle with some propbility after a certain age threshold is inspired from https://www.sciencedirect.com/science/article/pii/S0025326X18301000


## Example:
```
TI 2018-01-03
TSI 2018-01-15
TF 2018-01-31
PECCO2                  /home/maxime/Work/ECCO2/
INPUTF                  /home/maxime/Work/advection/input/test_label.xyz
OUTPUTP                 /home/maxime/Work/advection/res/
PSD                     /home/maxime/Work/Stokesdrift/
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
