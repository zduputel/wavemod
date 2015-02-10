'''
A class that deals with Green's function database

Written by Z. Duputel, May 2014
'''


# Personals
import sacpy 


# Externals
import numpy as np
import os.path as path
from glob import glob


class WaveDB(object):
    '''
    A class that deals with Green's function database
    
    The database is organized as W-phase GFs databases. 
    We use SAC files stored in sub-directories as follows:
       H????.?/??/GF.????.SY.LH?.SAC
    where 
       - H????.?: depth
       - ??: moment tensor component (PP,RP,RR,RT,TP,TT)
       - ????: epicentral distance x 10 (in km or deg)
       - LH?: component (either LHL or LHZ)    

    '''

    def __init__(self,GF_path,scale):
        '''
        Args:
             GF_path: path to GF database
             scale:   scalar factor for GFs
        '''

        self.name = 'GFs database'

        # Assign GF_path and scale
        assert path.exists(GF_path), '%s not found, no such directory'%(GF_path)
        self.GF_path = path.abspath(GF_path)
        self.scale   = scale

        # Set depths
        self.setDepths()

        # All done
        return

    
    def setDepths(self):
        '''
        List available source depths in the GF database
        '''
        
        # List all depths
        d = []
        for h_dir in glob(path.join(self.GF_path,'H*/')):
            items = h_dir.strip('/').split('/')
            d.append(float(items[-1][1:]))

        # Check length of d
        assert len(d) != 0, 'No H* directory found in %s'%(self.GF_path)
            
        # Convert to array and store it
        self.depths = np.array(d)
        
        # All done
        return
        

    def getMT(self,M0,strike_deg,dip_deg,rake_deg):
        '''
        Set moment tensor components from M0, strike, dip, rake 
        '''
        
        # Convert to rad
        deg2rad = np.pi/180.
        strike = strike_deg * deg2rad
        dip    = dip_deg    * deg2rad
        rake   = rake_deg   * deg2rad
        
        # Assign Sines/Cosines
        sr  = np.sin(rake)
        cr  = np.cos(rake)
        sd  = np.sin(dip)
        cd  = np.cos(dip)
        sd2 = np.sin(2*dip)
        cd2 = np.cos(2*dip)
        ss  = np.sin(strike)
        cs  = np.cos(strike)
        ss2 = np.sin(2*strike)
        cs2 = np.cos(2*strike)

        # Init moment tensor 
        MT = {}
        MT['RR'] =  M0 * sd2*sr                       # Mrr
        MT['TT'] = -M0 * (sd*cr*ss2 + sd2*sr*ss*ss)   # Mtt
        MT['PP'] =  M0 * (sd*cr*ss2 - sd2*sr*cs*cs)   # Mpp
        MT['RT'] = -M0 * (cd*cr*cs  + cd2*sr*ss)      # Mrt
        MT['RP'] =  M0 * (cd*cr*ss  - cd2*sr*cs)      # Mrp
        MT['TP'] = -M0 * (sd*cr*cs2 + 0.5*sd2*sr*ss2) # Mtp
        
        # All done
        return MT

    
    def rotMT(self,az_deg,MT):
        '''
        Rotate moment tensor MT according to az
        '''

        # Compute cosines/sines
        az = az_deg*np.pi/180.
        co  = np.cos(az)
        si  = np.sin(az)
        co2 = np.cos(2*az)
        si2 = np.sin(2*az)

        # Rotate MT
        MT_rot = {}
        MT_rot['RR'] = self.scale*MT['RR']
        MT_rot['TT'] = self.scale*(MT['TT']*co*co + MT['PP']*si*si - MT['TP']*si2)
        MT_rot['PP'] = self.scale*(MT['PP']*co*co + MT['TT']*si*si + MT['TP']*si2)
        MT_rot['RT'] = self.scale*(MT['RT']*co    - MT['RP']*si)
        MT_rot['RP'] = self.scale*(MT['RT']*si    + MT['RP']*co)
        MT_rot['TP'] = self.scale*((MT['TT'] - MT['PP'])*si2/2. + MT['TP']*co2)
        
        # All done
        return MT_rot

    
    def rotTraces(self,L_sac,T_sac,baz,cmpaz):
        '''
        Rotate L_sac and T_sac according to alpha
        '''
        
        # Calculate sines/cosines
        alpha_deg = baz-cmpaz
        alpha = alpha_deg * np.pi/180.
        co = np.cos(alpha)
        si = np.sin(alpha)
        
        # Rotate
        rot_sac = L_sac.copy()
        rot_sac.depvar = -co*L_sac.depvar -si*T_sac.depvar
        rot_sac.cmpaz = cmpaz

        # All done
        return rot_sac

    
    def bestDepth(self,depth):
        '''
        Find the best source depth in the database
        '''

        # Check depth
        assert self.depths != None, 'depths attribute must be assigned'
        assert len(self.depths) != 0, 'depths attribute must be assigned'

        # Find best depth
        dd = np.abs(depth-self.depths)
        i = np.where(dd==dd.min())[0][0]
        best_depth = self.depths[i]
        
        # All done
        return best_depth

    
    def bestDist(self,H_path, dist):
        '''
        Find the best source-station distance in the GF database 
        Args:
            * H_path: directory in which distances are listed
            * dist:   source-station distance
        '''
        
        # List all distances
        dists = []        
        for GF_file in glob(path.join(H_path,'TP','*LHT.SAC')):
            items = path.basename(GF_file).strip().split('.')
            dists.append(float(items[1])/10.)
        dists = np.array(dists)
            
        # Check dist
        assert dists.min()<dist, 'Too short source-station distance %d'%(dist)
        assert dists.max()>dist, 'Too large source-station distance %d'%(dist)
            
        # Find the best distance
        dd = np.abs(dists-dist)
        i  = np.where(dd == dd.min())[0][0]
        best_dist = dists[i]

        # All done
        return best_dist
        
        
    def sum_up(self,best_depth,dist,MT_rot):
        ''' 
        Sum up Green's functions 
        Args:
            * best_depth: best depth in the database (in km)
            * dist:       distance (in km or deg, depending on the d)
            * MT_rot:     rotated moment tensor
        Outputs:
            * Z_sac: sac object for component Z
            * L_sac: sac object for component L
            * T_sac: sac object for component T
        '''
        
        # Init sac
        Z_sac = sacpy.sac()
        L_sac = sacpy.sac()
        T_sac = sacpy.sac()
        GF_Z_sac = sacpy.sac()
        GF_L_sac = sacpy.sac()
        GF_T_sac = sacpy.sac()

        # Get prefix
        H_path   = path.join(self.GF_path,'H%05.1f'%(best_depth))
        
        # Get best distance
        best_dist = self.bestDist(H_path, dist)
        
        # GF_filenames
        Z_file = 'GF.%04.0f.SY.LHZ.SAC'%(best_dist*10.)
        L_file = 'GF.%04.0f.SY.LHL.SAC'%(best_dist*10.)
        T_file = 'GF.%04.0f.SY.LHT.SAC'%(best_dist*10.)

        # Sum up Green's functions for vertical/longitudinal components
        for MT_cmp in ['RR','TT','PP','RT']:
            # Read GFs
            GF_Z_file = path.join(H_path,MT_cmp,Z_file)
            GF_L_file = path.join(H_path,MT_cmp,L_file)
            GF_Z_sac.rsac(GF_Z_file)
            GF_L_sac.rsac(GF_L_file)
            # Set Output sacs
            if MT_cmp=='RR':
                Z_sac = GF_Z_sac.copy()
                L_sac = GF_L_sac.copy()
                Z_sac.depvar *= MT_rot[MT_cmp]
                L_sac.depvar *= MT_rot[MT_cmp]
            else:
                assert Z_sac.npts == GF_Z_sac.npts, 'Inconsistent number of samples'
                assert L_sac.npts == GF_L_sac.npts, 'Inconsistent number of samples'
                Z_sac.depvar += GF_Z_sac.depvar * MT_rot[MT_cmp]
                L_sac.depvar += GF_L_sac.depvar * MT_rot[MT_cmp]

        # Sum up Green's functions for transverse components
        for MT_cmp in ['RP','TP']:
            # Read GFs
            GF_T_file = path.join(H_path,MT_cmp,T_file)
            GF_T_sac.rsac(GF_T_file)
            # Set Output sacs
            if MT_cmp=='RP':
                T_sac = GF_T_sac.copy()
                T_sac.depvar *= MT_rot[MT_cmp]
            else:
                assert T_sac.npts == GF_T_sac.npts, 'Inconsistent number of samples: %d %d'%(T_sac.npts,GF_T_sac.npts)
                T_sac.depvar += GF_T_sac.depvar * MT_rot[MT_cmp]
        
        # All done
        return Z_sac, L_sac, T_sac

            
    def synth(self,depth,az,dist,MT):
        '''
        Compute synthetic waveforms
        Args:
            * depth: source depth
            * az:    station azimuth
            * dist:  distance
            * MT:    moment tensor
        '''

        # Rotate the moment tensor
        MT_rot = self.rotMT(az,MT)
        
        # Get the best depth
        best_depth = self.bestDepth(depth)

        # Sum up Green's functions
        Z_sac,L_sac,T_sac = self.sum_up(best_depth,dist,MT_rot)
        Z_sac.az = az
        L_sac.az = az
        T_sac.az = az

        #ituple = (MT['TT'],MT['PP'],MT['RR'],-MT['TP'],MT['RT'],-MT['RP'],az)
        #print best_depth,az,dist,'-XX %e -YY %e -ZZ %e -XY %e -XZ %e -YZ %e -A %.3f'%ituple
        
        # All done
        return Z_sac,L_sac,T_sac

    
    def synthSDR(self,depth,az,dist,M0,strike,dip,rake):
        '''
        Compute synthetic waveforms from strike, dip, rake and M0
        Args:
           * depth:  source depth
           * az:     station azimuth
           * dist:   distance
           * M0:     seismic moment in dyne-cm
           * strike: strike angle in degrees
           * dip:    dip angle in degrees
           * rake:   rake angle in degrees
        '''
        
        # Get the moment tensor from strike, dip and rake
        MT = self.getMT(M0,strike,dip,rake)
        #print MT['RR'],MT['TT'],MT['PP'],MT['RT'],MT['RP'],MT['TP']
        #exit(1)
        
        # Get depth, az, dist and MT from
        Z_sac,L_sac,T_sac = self.synth(depth,az,dist,MT)
        
        # All done
        return Z_sac,L_sac,T_sac

