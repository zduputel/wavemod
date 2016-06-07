'''
A class that deals with Kikuchi-Kanamori Teleseismic data and Green's functions 

Written by Z. Duputel, May 2016
'''

# Personals
from .utils import rm
from .sacpy import sac

# Externals
import sys
import os
import pyproj as pp
import numpy as np
from subprocess import call
from glob       import glob
from copy       import deepcopy

def nextpow2(i):
    n = 1
    while n < i:
        n *= 2
    return n

class WaveKK(object):
    '''
    A class that deals with Kikuchi-Kanamori Teleseismic data and Green's functions 
    '''

    def __init__(self,evid=None,delta=0.5,T0=-10.0,TL=100.,dtl=1.,scale=1.0e-18):
        '''
        Args:
             delta:       Sampling period
             T0:          time of the first sample 
             TL:          duration of each trace
             dtl:         triangle width
        '''

        self.name = 'Kikuchi Kanamori'

        # Assign attributes from input parameters
        self.evid       = evid
        self.delta      = delta
        self.T0         = T0
        self.ixa        = int(T0/delta) # safety delay
        self.TL         = TL
        self.npts       = int((TL)/delta)
        self.dtl        = dtl
        self.nchan      = None
        self.scale      = scale

        # Data
        self.data  = None
        self.chans = None

        # Green's function database
        self.GF   = {}        
        self.GFdb = None
        self.GFparams = None
        
        # Earth model
        self.vmodelname = None
        self.vp  = None
        self.vs  = None
        self.rho = None
        self.dep = None
        
        # Hard-wired assignements
        self.dist       = None

        # All done
        return
    
        
    def wimom3(self,fname,H0,strike,dip,rake):
        '''
        Writes i_mom3 file
        Args:
            * fname: imom3 filename
            * H0:    depth
            * strike,dip,rake: strike dip rake angles
        '''
        # Write i_mom3 file
        f = open(fname,'wt')
        if self.evid is not None:
            f.write(self.evid+'\n')
        else:
            f.write('None\n')
        f.write('fort.1\n')
        f.write('%.4f %.4f %.4f %.4f %.4f %.4f %.4f\n'%(self.T0,self.TL,self.delta,H0,strike,dip,rake))
        f.write('1 1 1.0 1 1 1.0 3.0 %.4f %.4f %.4f 1 0.5\n'%(-self.T0,self.dtl,self.dtl))
        f.write('%d\n'%(self.nchan))
        c = 1
        for i in range(self.nchan):
            f.write('1 ')
            c += 1
            if c ==4:
                f.write('\n')
                c = 1
        if c!=1:
            f.write('\n')        
        f.write('%d\n'%(nextpow2(self.npts)))
        f.close()

        # All done
        return

    def readfort2(self,ifile='fort.2'):
        '''
        Read kk model file (fort.2 formated file)
        Args:
            * ifile: input file
        '''
        f = open(ifile,'rt')
        self.vmodelname = f.readline().strip()
        items = f.readline().strip().split()
        self.tqp = float(items[0])
        self.tqs = float(items[1])
        N = int(items[2])
        self.vp  = []
        self.vs  = []
        self.rho = []
        self.dep = []
        h = 0
        for i in range(N):
            self.vp.append(float(items[-4]))
            self.vs.append(float(items[-3]))
            self.rho.append(float(items[-2]))
            h += float(items[-1])
            self.dep.append(h)
            items = f.readline().strip().split()            
        f.close()
        # All done
        return
        
    
    def readfort1(self,ifile='fort.1',GFfile=False):
        '''
        Read waveforms from a kk fort.1 formated file
        Args:
            * ifile: input file
            * GFfile: if True, ifile includes Green's functions
        '''
        # Initialize dictionaries and lists
        W = {}
        Wlist = []

        # Main loop
        f = open(ifile,'rt')
        l = f.readline()
        if self.evid==None:
            self.evid = l.strip().split()[0]
        while True:
            l = f.readline()
            if not l:
                break
            items = l.strip().split()

            # Parse file
            stat   = items[0]
            if GFfile:
                strike = np.round(float(items[1]),1)
                dip    = np.round(float(items[2]),1)
                rake   = np.round(float(items[3]),1)
                depth  = float(items[4])
            items = f.readline().strip().split()
            az   = float(items[0])
            az2  = float(items[1])
            dist = float(items[2])
            p    = float(items[3])
            g    = float(items[4])
            ix   = int(items[5])+self.ixa 
            assert ix >= 0, 'ix-ixa must be larger than 0 (ix=%d)'%(ix)
            items = f.readline().strip().split()
            ib = int(items[1])
            assert ib>=1 and ib<=4, 'ib must be >=1 and <=4'
            f.readline()
            f.readline()
            items = f.readline().strip().split()
            ym = float(items[0])
            N  = int(items[1])
            dt = float(items[2])
            ie   = ix + self.npts
            assert dt == self.delta, 'delta in %s (%f vs %f)'%(ifile,dt,self.delta)
            assert ie <= N, 'Waveform too short in %s (%d vs %d)'%(ifile,ie,N)
            wave = np.array([],dtype='float64')
            while len(wave)<N:
                items = list(map(float,f.readline().strip().split()))
                wave = np.append(wave,items)
            wave *= ym
            assert len(wave)==N, 'Incorrect waveform length for %s in %s'%(stat,ifile)
            if GFfile==False:
                wave = wave[ix:ie].copy()
            else:
                wave = wave[ix:].copy()
            # Set waveform id
            if ib == 1:
                Wtype = 'P'
            elif ib == 2:
                Wtype = 'SV'
            elif ib == 3:
                Wtype = 'SH'
            else:
                Wtype = 'PP'
            Wid = stat+Wtype
            Wsac = sac()
            items = Wid.strip().split('.')
            Wsac.knetwk = items[0]
            Wsac.kstnm  = items[1]
            if len(items)==4:
                Wsac.khole  = items[2]
            elif len(items)==3:
                Wsac.khole  = ''
            else:
                sys.stderr.write('Incorrect waveform id\n')
                sys.exit(1)
            Wsac.kcmpnm = Wtype
            Wsac.b      = self.T0
            Wsac.delta  = dt
            Wsac.npts   = len(wave)
            Wsac.az     = az
            Wsac.baz    = az2
            Wsac.gcarc  = dist
            Wsac.user[0] = p
            Wsac.depvar = wave.copy()

            # Fill up dictionary
            assert Wid not in Wlist, 'Multiple entries for %s'%(Wid)
            Wlist.append(Wid)
            W[Wid] = Wsac.copy()
        f.close()

        # Set channel list (if not a GF file)
        if not GFfile:
            self.nchan = len(Wlist)
            self.chans = deepcopy(Wlist)
        elif GFfile:
          assert len(Wlist)==self.nchan, 'Incorrect number of channels'  
            
        # All done
        return W

    def readData(self,ifile):
        '''
        Read data from kk fort.1 formated file
        '''
        self.data = self.readfort1(ifile)
    
        # All done
        return

    def computeGFdb(self,Hs,Strikes,Dips,Rakes):
        '''
        Compute GFs database
        Args:
            * Hs,Strikes,Dips,Rakes: lists of depth, strikes, dips, rakes
        '''

        # Check things and cleanup
        assert self.nchan is not None, 'Must read data first (self.readData)'
        
        # Main loop
        cmd = 'compute_green' # Patched version of mom3_large_v4 (by Z. Duputel)
        self.GFdb = []
        GFparams = []
        for ho,s,d,r in zip(Hs,Strikes,Dips,Rakes):

            # Create i_mom3 file
            imom3 = 'i_mom3_h'
            self.wimom3(imom3,ho,s,d,r)

            # Run compute_green
            if os.path.exists('fort.69'):
                os.remove('fort.69')
            ofile = 'o_GFs.txt'
            ifd = open(imom3,'rt')
            ofd = open('o_compute_green','wt')
            call(cmd,shell=True,stdin=ifd,stdout=ofd,stderr=sys.stderr)
            ofd.close()
            ifd.close()

            # Read results
            W = self.readfort1('fort.69',GFfile=True)
            os.remove('fort.69')

            # Append to list
            self.GFdb.append(deepcopy(W))
            GFparams.append([ho,s,d,r])
        self.GFparams = np.array(GFparams)
        
        # All Done
        return

    def computeGF(self,lon0,lat0,H0,lons,lats,Hs,Strikes,Dips,Rakes,ellps='WGS84',causal=True):
        '''
        Compute Green's functions from GF database
        Args:
            * lon0,lat0,H0: Hypocenter coordinates
            * lons,lats,Hs,Strikes,Dips,Rakes: Coordinates/Orientation of patches
            * ellps: Reference ellipsoid for distance/azimuth calculation
        '''
        # Check that everything is ready
        assert self.data is not None, 'Must read data first (use self.readData)'
        assert self.GFdb is not None, 'GF database not ready (use computeGFdb)'
        assert self.vmodelname is not None, 'Velocity model not loaded (use readfort2)'
        

        # Set up stuff to compute distances and angles
        geod = pp.Geod(ellps=ellps)
        rad  = np.pi/180.

        # Get Vp/Vs at the hypocenter
        for ll in range(len(self.vp)):
            if (H0 - self.dep[ll]) < 0:
                break

        # Create velocity dictionary (i.e., Vp or Vs for each waveform)
        V = {}
        for dkey in self.data:
            if self.data[dkey].kcmpnm[0] == 'P':
                V[dkey] = self.vp[ll]
            elif self.data[dkey].kcmpnm[0] == 'S':
                V[dkey] = self.vs[ll]
            else:
                sys.stderr.write('Error: incorrect channel name format for %s in self.data\n'%(dkey))
                sys.exit(1)        

        # Main loop
        self.GF = []        
        for lon,lat,h,s,d,r in zip(lons,lats,Hs,Strikes,Dips,Rakes):

            # Compute patch-to-patch azimuth and distance
            az_s,baz_s,dist_s = geod.inv(lon0,lat0,lon,lat,radians=False)
            dist_s *= 1.0e-3
            
            # Find appropriate GF in the GF database
            DH = np.abs(self.GFparams[:,0]-h)
            DS = np.abs(self.GFparams[:,1]-s)
            DD = np.abs(self.GFparams[:,2]-d) 
            DR = np.abs(self.GFparams[:,3]-r)           
            i = np.where((DH.min() == DH) & (DS.min() == DS) & (DD.min() == DD) & (DR.min() == DR))[0]
            assert len(i)==1, 'Lookup table issue (%d occurences) h=%f s=%f d=%f r=%f'%(len(i),h,s,d,r)
            assert DH.min() <= DH.max() * 1.0e-4, 'Lookup table issue (depth residual: %f %f)'%(DH.min())
            assert DS.min() <= DS.max() * 1.0e-4, 'Lookup table issue (strike residual: %f %f)'%(DS.min(),DS.max())
            assert DD.min() <= DD.max() * 1.0e-4, 'Lookup table issue (dip residual: %f)'%(DD.min())
            assert DR.min() <= DR.max() * 1.0e-4, 'Lookup table issue (rake residual: %f)'%(DR.min())
            
            GFdb = self.GFdb[i[0]]
            self.GF.append({})
            # Loop over stations
            for wid in GFdb:
                self.GF[-1][wid] = {}
                # Extract relevant information
                p  = GFdb[wid].user[0]
                az = GFdb[wid].az
                dh = h - H0
                #dh = 0.
                v  = V[wid]
                # Time-shift
                ts1 = dist_s * p * np.cos((az-az_s)*rad)
                ts2 = dh * np.sqrt(1./(v*v) - p*p)
                #its  = int(-(self.T0 + ts1 + ts2)/self.delta)
                if causal:
                    tg = -(ts1 + ts2) + np.arange(GFdb[wid].npts)*GFdb[wid].delta
                else:
                    tg = -(ts1 + ts2) + np.arange(GFdb[wid].npts)*GFdb[wid].delta - self.dtl
                td = self.T0 + np.arange(self.data[wid].npts)*self.data[wid].delta
                # Time-shifted GF
                self.GF[-1][wid] = self.data[wid].copy()
                self.GF[-1][wid].depvar *= 0.
                self.GF[-1][wid].user[0] = p
                self.GF[-1][wid].depvar = np.interp(td,tg,GFdb[wid].depvar*self.scale)
                #for j in range(self.data[wid].npts):
                #    i1 = j - its
                #    if (i1 >= 0) and (i1 < self.data[wid].npts):                        
                #         self.GF[-1][wid].depvar[j] = GFdb[wid].depvar[i1]
                            
                            
