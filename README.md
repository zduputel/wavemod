# Wavemod
WaveMod calculates synthetic waveforms using various tools. Three main classes are currently available:
- `WaveDB` deals with synthetic waveforms computed from Green's function database.
- `WaveInt` is a wrapper of Bob Hermann's Wavenumber Integration code.
- `FastSweep` is a simple 2D eikonal solver used for rupture front propagation on simple faults.

## Some instructions
- The WaveMod module must be located in a directory that is incuded in the PYTHONPATH environment.
- To use `WaveDB`, Green's function databases must be organized following the *W-phase format*. We use SAC files stored in sub-directories as follows:
   `H????.?/??/GF.????.SY.LH?.SAC`
where 
   - `H????.?` indicates depth (in km)
   - `/??/` indicates the moment tensor component (PP,RP,RR,RT,TP,TT)
   - `????.SY` is the epicentral distance x 10 (in km or deg)
   - `LH?` is the component (either LHL or LHZ)  
- SAC files are handled using the sacpy module available at https://github.com/eost/sacpy
- To use WaveInt, the package "Computer programs in seismology" must be installed and path to associated binaries must be included in the PATH environment variable : http://www.eas.slu.edu/eqc/eqccps.html
