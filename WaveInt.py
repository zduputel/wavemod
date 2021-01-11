'''
A class that deals with Bob Hermann's wavenumber integration code

Written by Z. Duputel, September 2013
'''

# Personals
from .utils import rm
from .sacpy import sac

# Externals
import sys
import os
import numpy as np
from subprocess import call
from glob       import glob
from copy       import deepcopy


class WaveInt(object):
    '''
    A class that deals with Bob Hermann's wavenumber integration code
    '''

    def __init__(self,model_file,npts,delta,T0=0.0,Vred=0.0,Xs=None,
                 stat=None,Xr=None,dist_file='dfile'):
        '''
        Args:
             model_file:  Earth model filename
             npts:        Number of samples 
             delta:       Sampling period
             T0:          time of the first sample if Vred=0.0 (optional, default=0.0) 
             Vred:        Reduction velocity (optional, default=0.0)
             Xs:          Source [X,Y,Z] coordinate in km (optional, default=None)
             stat:        Station names (optional, default=None)
             Xr:          Receivers XY coordinates in km (optional, default=None)
             dist_file:   Name of dist filename (optional, default='dfile')
        '''

        self.name = 'Waveform Integration'

        # Check stat and Xr if both of them are specified
        if stat!=None and Xr!=None:
            assert len(stat)==len(Xr), 'Xr and stat must have same length'

        # Convert Xr
        if Xr!=None and not isinstance(Xr,np.ndarray):
            Xr = np.array(Xr,dtype='float64')
            
        # Assign attributes from input parameters
        self.model_file = os.path.abspath(model_file)
        self.dist_file  = dist_file
        self.npts       = npts
        self.delta      = delta
        self.T0         = T0
        self.Vred       = Vred
        self.Xs         = Xs
        self.stat       = stat
        self.Xr         = Xr
        
        # Output synthetic seismograms
        self.synth    = None

        # Hard-wired assignements
        self.dist       = None

        # All done
        return
    
    def readStatXY(self,station_file):
        '''
        Read station file and populate the Xr attribute (station coordinates)
        Args:
              station_file: station filename including station coordinates
        file format:
        STNAME  X_COORD Y_COORD
        '''
        
        # Assert if station file exists
        assert os.path.exists(station_file), 'Cannot read %s (no such file)'%(station_file)
        
        # Read the file and fill-up Xr
        self.stat = []; self.Xr = []
        for l in open(station_file):
            if l.strip()[0]=='#':
                continue
            items = l.strip().split()
            self.stat.append(items[0].strip())
            self.Xr.append([float(items[1]),float(items[2])])
        self.Xr = np.array(self.Xr,dtype='float64')
        
        # All done
        return

    def setXr(self,stat_names,X,Y):
        '''
        Set X and Y receiver values
        Args:
            *stat_names: list of station names
            * X: list or array of x coordinates
            * Y: list or array of y coordinates
        '''

        # Check if X and Y are lists or arrays
        if type(X)!=list and type(X)!=np.ndarray:
            X = [X]
        if type(Y)!=list and type(Y)!=np.ndarray:
            Y = [Y]
        if type(stat_names)!=list:
            stat_names = [stat_names]
        
        # Check length of stat_names, X and Y
        assert len(X)==len(Y), 'X and Y should have same length'
        assert len(stat_names)==len(X), 'stat_names and X should have same length'

        # Assign X and Y to Xr
        self.stat = deepcopy(stat_names)
        self.Xr = np.array([[x,y] for x,y in zip(X,Y)],dtype='float64')

        # All done
        return

    def checkXs(self):
        '''
        Check and convert source coordinates to the correct format
        '''
        
        # Check if it has been correctly assigned
        assert self.Xs!=None,   'Xs is not assigned correctly (must be a (3,) ndarray)'
        assert len(self.Xs)==3, 'Xs is not assigned correctly (must be a (3,) ndarray)'

        # Convert to ndarray
        if not isinstance(self.Xs,np.ndarray):            
            self.Xs = np.array(self.Xs,dtype='float64')
        
        # All done
        return

    def calcDist(self,station_file=None):
        '''
        Calculate distances and populate the dist attribute
        Args:
             station_file:  file including station coordinate (optional, default=None)
        '''
        
        # Read station file
        if station_file!=None:
            self.readStatXY(station_file)

        # Check/Convert source coordinates
        self.checkXs()

        # Assert if Xr is correctly assigned
        assert isinstance(self.Xr,np.ndarray), 'Xr is not assigned correctly (must be a (N,2) ndarray)'
        
        # Compute distance and create the dist array attribute
        dX = self.Xr-self.Xs[:2]
        self.dist = np.sqrt((dX*dX).sum(axis=1))
        
        # All done
        return
    
    def writeDistFile(self):
        ''' 
        Create distance file from the dist array attribute
        '''
        
        # Assert if dist is correct
        assert isinstance(self.dist,np.ndarray), 'The array dist is not correct'
        
        # Find unique elements of the dist array attribute
        unique_dist = np.unique(self.dist)
        
        # Write the file and close
        f = open(self.dist_file,'wt')
        format = '%10.5f %6.3f %6d %7.3f %6.3f\n'
        for d in unique_dist:
            f.write(format%(d,self.delta,self.npts,self.T0,self.Vred))
        f.close()
        
        # All done
        return

    def writeModelFile(self,Vp,Vs,Rho,H):
        '''
        Create model file from input Vp, Vs and thickness (H)
        '''
    
        # Assert if input parameters have the correct size
        assert len(Vp)==len(Vs)==len(Rho)==len(H), 'Vp,Vs,Rho,H must have same dimensions'

        # Write the file and close
        f = open(self.model_file,'wt')
        f.write('%d 1000.\n'%(len(h)))
        for vp,vs,rho,h in zip(Vp,Vs,Rho,H):
            f.write('%7.4f %7.4f %7.4f %7.4f\n'%(rho,vp,vs,h))
        f.close()
        
        # All done
        return

    def calcKernel(self,ofd=sys.stdout,efd=sys.stderr):
        '''
        Calculate Green's functions in the frequency domain
        Args:
             ofd: stream for standard output (optional, default=sys.stdout)
             efd: stream for standard error  (optional, default=sys.stdout)
        '''
        
        # Check/convert source coordinates
        self.checkXs()

        # Assert if the distance file exists and if Xs is correct
        assert os.path.exists(self.dist_file), 'Cannot read %s (no such file)'%(self.dist_file)
        
        # Define commands
        hprep_cmd = 'hprep96 -M %s -d %s -HS %f -HR 0.0'%(self.model_file,self.dist_file,self.Xs[2])
        hspec_cmd = 'hspec96'
        
        # Calculate the GF spectrums
        call(hprep_cmd, shell=True,stdout=ofd,stderr=efd)
        call(hspec_cmd, shell=True,stdout=ofd,stderr=efd)
        
        # All done
        return

    def writeRfile(self,rfile_name,STF):
        '''
        Write STF to rfile_name
        '''
        f  = open(rfile_name,'wt')
        ns = len(STF)
        f.write('%5d%10.3f\n'%(ns,self.delta))
        for i in xrange(len(STF)):
            f.write('%15.7e'%STF[i])				
        f.write('\n')
        f.close()
        
        # All done
        return
        
    def synthKernelSDR(self,out_type,strike,dip,rake,M0,stf_type,duration=None,rfile=None,
                      ofd=sys.stdout,efd=sys.stderr):
        '''
        Calculate synthetics from pre-calculated Green's functions for a given Strike, Dip, Rake, M0 and a given STF
        Args:
             out_type   output type ('D'=displacement, 'V'=velocity, 'A'=acceleration)
             strike:    Stike angle in degrees
             dip:       Dip angle in degrees
             rake:      Rake angle in degrees
             M0:        Moment in dyne-cm
             stf_type:  STF type: 'triangle', 'boxcar', 'dirac', 'rfile'=user supplied STF in rfile
             duration:  STF duration (optional, used if stf_type is not 'dirac' or 'rfile')
             rfile:     name of file including a user supplied normalized pulse
             ofd:       stream for standard output (optional, default=sys.stdout)
             efd:       stream for standard error  (optional, default=sys.stdout)        
        '''
        
        # Assert if input parameters are correct
        if stf_type=='rfile':
            assert rfile!=None, 'rfile filename should be specified (%s)'%(stf_type)
            assert os.path.exists(rfile), 'Cannot read %s (no such file)'%(rfile)
        if stf_type!='dirac' and stf_type!='rfile':
            assert duration!=None, 'STF duration missing (%s)'%(stf_type)
        
        # Assert if Xr, Xs and dist are correctly assigned
        assert isinstance(self.Xr,np.ndarray), 'Xr is not assigned correctly (must be a (N,2) ndarray)'
        assert isinstance(self.Xs,np.ndarray), 'Xs is not assigned correctly (must be a (3,) ndarray)'
        assert isinstance(self.dist,np.ndarray), 'The array dist is not correct'
        assert len(self.Xr)==len(self.dist), 'dist and Xr must have the same length (%d vs %d)'%(len(self.Xr),len(self.dist))

        # Assign self.stat if not specified
        if self.stat==None:
            self.stat = []
            for i in range(len(self.Xr)):
                self.stat.append('STA_%d'%(i+1))

        # Find unique elements of the dist array attribute
        unique_dist = np.unique(self.dist)
        
        # From the half-duration, check if the source is a dirac
        if duration!=None:
            half_dur   = duration/2.   # Half-duration
            if half_dur < self.delta: 
                efd.write('Warning: Using half-duration<sampling rate: will use dirac STF\n')
                stf_type = 'dirac'
                half_dur_factor = None
            else:
                half_dur_factor = int(np.round(half_dur/self.delta))
        
        # Initialize the hpulse command
        hpulse_f96 = 'hpulse_file96' # hpulse output file96
        if stf_type=='dirac':      # Dirac
            hpulse_cmd = 'hpulse96 -%c -i > %s'%(out_type,hpulse_f96)
        elif stf_type=='triangle': # Triangle
            hpulse_cmd = 'hpulse96 -%c -t -l %d -Z > %s'%(out_type,half_dur_factor,hpulse_f96)
        elif stf_type=='rfile':    # User supplied STF
            hpulse_cmd = 'hpulse96 -%c -F %s > %s'%(out_type,rfile,hpulse_f96)
        elif stf_type=='boxcar':     # Boxcar STF
            rfile = 'boxcar'
            hpulse_cmd = 'hpulse96 -%c -F %s > %s'%(out_type,rfile,hpulse_f96)
            # build a rfile with a boxcar
            ns  = int(round(duration/self.delta))
            STF = np.ones(ns)/duration
            self.writeRfile(rfile,STF)
        else:
            raise TypeError('stf_type')
        
        # STF convolution
        call(hpulse_cmd,shell=True,stdout=ofd,stderr=efd)
        
        # Compute azimuth and backazimuth
        dX  = self.Xr-self.Xs[:2]
        Az  = np.arctan2( dX[:,0], dX[:,1])*180./np.pi
        BAz = np.arctan2(-dX[:,0],-dX[:,1])*180./np.pi
        
        # Compute data from mechanism and azimuth and convert to SAC
        fsel_f96  = 'fsel_file96'
        fmech_f96 = 'fmech_file96'
        fsel_cmd_format  = 'fsel96 -NS %d < %s > %s'
        fmech_cmd_format = 'fmech96 -S %.4f -D %.4f -R %.4f -M0 %.5e -A %.4f -B %.4f < %s > %s'
        for j in xrange(unique_dist.size):
            for k in xrange(self.dist.size):
                if self.dist[k]==unique_dist[j]:
                    # Select Green's functions
                    fsel_cmd  = fsel_cmd_format%(j+1,hpulse_f96,fsel_f96)
                    call(fsel_cmd, shell=True,stdout=ofd,stderr=efd)
                    # Moment tensor product
                    fmech_cmd = fmech_cmd_format%(strike,dip,rake,M0,Az[k],BAz[k],fsel_f96,fmech_f96)
                    call(fmech_cmd, shell=True,stdout=ofd,stderr=efd)
                    # Conversion to sac
                    f96tosac_cmd = 'f96tosac '+fmech_f96 # Conversion to Sac
                    call(f96tosac_cmd, shell=True,stdout=ofd,stderr=efd)
                    # Cleanup
                    os.rename('B00101Z00.sac','%s_Z.SAC'%(self.stat[k]))
                    os.rename('B00102N00.sac','%s_N.SAC'%(self.stat[k]))
                    os.rename('B00103E00.sac','%s_E.SAC'%(self.stat[k]))
                    rm(glob('*.sac'))
                    rm([fsel_f96,fmech_f96])
        rm(hpulse_f96)

        # All done
        return

    def synthSDR(self,out_type,strike,dip,rake,M0,stf_type,duration=None,rfile=None,
                 calc_dist=True,ofd=sys.stdout,efd=sys.stderr):
        '''
        Calculate Synthetics from scratch for a given Strike, Dip, Rake, M0 and a given STF
        Args:
             out_type   output type ('D'=displacement, 'V'=velocity, 'A'=acceleration)
             strike:    Stike angle in degrees
             dip:       Dip angle in degrees
             rake:      Rake angle in degrees
             M0:        Moment in dyne-cm
             stf_type:  STF type: 'triangle', 'boxcar', 'dirac', 'rfile'=user supplied STF in rfile
             duration:  STF duration (optional, used if stf_type is not 'dirac' or 'rfile')
             rfile:     name of file including a user supplied normalized pulse
             calc_dist: if True, calculate epicentral distances (optional, default=False)
             ofd:       stream for standard output (optional, default=sys.stdout)
             efd:       stream for standard error  (optional, default=sys.stdout)        
        '''
        
        # Calculate distances
        if calc_dist:
            self.calcDist()
        
        # Write Distance File
        self.writeDistFile()
        
        # Calculate Green's functions in the frequency domain
        self.calcKernel(ofd=ofd,efd=efd)

        # Calculate synthetics from pre-calculated kernels
        self.synthKernelSDR(out_type,strike,dip,rake,M0,stf_type,duration,rfile,ofd=ofd,efd=efd)
        
        self.synth = {}
        for stat in self.stat:
            self.synth[stat]={}
            for c in 'ZNE':
                sacfile = sac()
                sacfile.read('%s_%c.SAC'%(stat,c))
                # Conversion from cm -> m
                sacfile.depvar *= 1.0e-2
                self.synth[stat][c] = deepcopy(sacfile)

        # All done
        return
        
        
