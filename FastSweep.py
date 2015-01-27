'''
A set of classes that uses the Fast sweeping method to solve the eikonal equation in 2D:
    Zhao, Hongkai (2004), A FAST SWEEPING METHOD FOR EIKONAL EQUATIONS, ...
    MATHEMATICS OF COMPUTATION, 74(250), 603-627, S 0025-5718(04)01678-3

Written by Z. Duputel, March 2014
'''

# Externals
import copy
import numpy as np
import sys

class Hypocenter(object):
    '''
    A class that defines hypocenter coordinates
    '''

    def __init__(self,x=None,y=None):
        '''
        Args:
            * Hypocenter x coordinates
            * Hypocenter y coordinates
        '''

        self.setHypo(x,y)
        
        # All done
        return


    def setHypo(self,x,y):
        '''
        Set hypocenter coordinates
        '''
        # Check input parameters
        assert type(x)==float, 'hypo x coordinate must be float'
        assert type(y)==float, 'hypo y coordinate must be float'
        
        # Set coordinates
        self.x = x
        self.y = y


    def copy(self):
        '''
        Returns a copy of the Hypocenter
        '''
        
        return copy.deepcopy(self)





class Grid(object):
    '''
    A class that defines the grid for fast sweeping
    '''
    
    def __init__(self,x=None,y=None,vr=None,name=''):
        '''
        Args:
            * x  coordinates in the grid
            * y  coordinates in the grid
            * vr rupture velocity 
            * h  grid spacing
        '''
        
        self.name = name
        
        # Init Grid
        if x==None or y==None or vr==None:
            self.x  = x.copy()
            self.y  = y.copy()
            self.vr = vr.copy()
            self.h  = None
            self.nx = None
            self.ny = None
        else:
            self.setGrid(x,y,vr)

        # All done
        return


    def setGrid(self,x,y,vr):
        '''
        Setting up the grid
        Args:
            * x coordinates in the grid
            * y coordinates in the grid
            * vr rupture velocity 
        '''     
        print('Setting up the {} grid'.format(self.name))
        
        # Check input parameters
        assert type(x) == np.ndarray, 'x must be numpy.ndarray'
        assert type(y) == np.ndarray, 'y must be numpy.ndarray'
        assert type(vr)== np.ndarray, 'vr must be numpy.ndarray'

        # Init Grid
        self.x    = x.copy()
        self.y    = y.copy()
        self.vr   = vr.copy()

        # Define grid shape and spacing
        self.ny,self.nx = vr.shape
        self.shape = (self.ny,self.nx)
        self.h = self.x[0,1]-self.x[0,0]

        # Checking grid attributes
        print('Checking {} grid'.format(self.name))
        self.checkGrid()

        # All done
        return


    def checkGrid(self):
        '''
        Check the grid attributes
        '''

        # Check sizes
        assert self.nx == self.x.shape[1],  'Incorrect attribute for x size'
        assert self.ny == self.y.shape[0],  'Incorrect attribute for y size'
        assert self.x.shape==self.y.shape,  'Shapes of x and y must be consistent'
        assert self.x.shape==self.vr.shape, 'Shape of vr must be consistent'
        
        # Check spacing
        hr = round(self.h,2)
        for i in range(self.ny): 
            for j in range(self.nx-1):
                dx = round(self.x[i,j+1]-self.x[i,j],4)
                assert hr == dx, 'x grid spacing must be regular'
        for i in range(self.ny-1):
            for j in range(self.nx):
                dy = round(self.y[i+1,j]-self.y[i,j],4)
                assert hr == dy, 'y grid spacing must be regular'
        

    def copy(self):
        '''
        Returns a copy of the Grid
        '''
        
        return copy.deepcopy(self)



class FastSweep(object):

    def __init__(self,grid=None,hypo=None):
        '''
        Args:
            grid: grid object (optional)
            hypo: hypocenter object (optional)
        '''

        self.name = 'Fast sweeping'

        # Init attributes
        self.grid = None
        self.hypo = None
        self.pad_grid = None
        self.sr = None

        # Check grid and hypo
        if grid != None:
            assert grid==Grid, 'grid must be a Grid object'
            self.grid == grid.copy()
        if hypo != None:
            assert hypo==Hypocenter, 'grid must be a Hypocenter object'
            self.hypo == hypo.copy()

        # All done
        return

    
    def setGrid(self,x_vec,y_vec,vr_mat):
        '''
        set Grid
        Args:
            * x_vec: vector of x coordinates
            * y_vec: vector of y coordinates
            * vr_max: matrix of rupture velocities at (x,y) 
        '''
    
        x_mat,y_mat = np.meshgrid(x_vec,y_vec)
        self.grid = Grid(x_mat,y_mat,vr_mat,'original')
        
        # All done
        return


    def setGridFromFault(self,fault,grid_space=1.):
        '''
        set Fast sweeping Grid from kinematic fault object 
        (assuming a planar fault and rectangular patches)
        Args:
            * fault: kinematic fault object
            * grid_space: spacing for the grid used to solve eikonal
              (if == None, will use 1 grid point at the center of each patch)
        '''

        ## Check fault strike dip rake
        #assert fault.f_strike != None, 'Fault strike must be assigned'
        #assert fault.f_dip    != None, 'Fault dip must be assigned'

        # Loop over each patch
        Np = len(fault.patch)
        g_dip  = []; g_strike  = []
        g_dipc = []; g_strikec = []
        g_vr     = []        
        for p in range(Np):

            # Get patch location and geometry
            p_x,p_y,p_z,p_W,p_L,p_strike,p_dip = fault.getpatchgeometry(p,center=True)
            if p==0:
                patch_size = p_W

            ## Check that the fault is planar
            #assert np.round(p_strike,2)==np.round(fault.f_strike,2), 'Fault must be planar' 
            #assert np.round(p_dip,2)   ==np.round(fault.f_dip,2)   , 'Fault must be planar' 
            if p==0:
                width  = np.round(p_W,2)
                length = np.round(p_L,2)
            assert np.round(p_W,2)==width,  'Patch width  must be homogeneous over the fault'
            assert np.round(p_L,2)==length, 'Patch length must be homogeneous over the fault'

            # get coordinates along fault
            g_dip_c, g_strike_c = fault.getHypoToCenter(p,True)
            g_dip_e = [np.round(g_dip_c-p_W/2.,2),np.round(g_dip_c+p_W/2.,2)] 
            g_strike_e = [np.round(g_strike_c-p_L/2.,2),np.round(g_strike_c+p_L/2.,2)]
            g_dip.append(g_dip_e)
            g_strike.append(g_strike_e)
            g_dipc.append(g_dip_c)
            g_strikec.append(g_strike_c)
            if grid_space == None:
                assert p_W==patch_size, 'Patch size must be uniform'
                assert p_L==patch_size, 'Patch size must be uniform'
        g_dip    = np.array(g_dip)
        g_strike = np.array(g_strike)
        g_dipc    = np.array(g_dipc)
        g_strikec = np.array(g_strikec)
        
        # Dip is x, Strike is y
        if grid_space != None:
            x = np.arange(g_dip.min()+grid_space/2.,g_dip.max(),grid_space)
            y = np.arange(g_strike.min()+grid_space/2.,g_strike.max(),grid_space)        
        else:
            x = np.arange(g_dipc.min(),g_dipc.max()+patch_size,patch_size)
            y = np.arange(g_strikec.min(),g_strikec.max()+patch_size,patch_size)                    
        x = np.round(x,2)
        y = np.round(y,2)

        # Assign vr to grid
        vr_mat = np.zeros((y.size,x.size),dtype='float64')                        
        for p in range(Np):
            i = np.where((y>=g_strike[p][0]) & (y<=g_strike[p][1]))[0]
            j = np.where((x>=g_dip[p][0])    & (x<=g_dip[p][1]))[0]
            vr_mat[i.min():i.max()+1,j.min():j.max()+1] = fault.vr[p]

        # Check vr matrix
        assert not (vr_mat==0.).any(), 'incorrect vr assigment'

        # Set grid
        self.setGrid(x,y,vr_mat)

        # Set Hypo
        self.setHypo(0.,0.)

        # All done
        return
        
        
    def setHypo(self,x,y):
        '''
        Set Hypocenter
        '''

        self.hypo = Hypocenter(x,y)

        # All done 
        return
        

    def calcHypoDist(self,pad=True):
        '''
        Computes distance between hypocenter and grid nodes
        Args:
            * pad: if True, will compute distances for the padding grid
        '''
        
        # Check hypo
        assert self.hypo != None, 'Hypocenter coordinates must be assigned'
        
        # Coordinates with respect to hypocenter
        if pad:
            assert self.pad_grid != None, 'Padding grid must be assigned for pad=True'
            X = self.pad_grid.x - self.hypo.x
            Y = self.pad_grid.y - self.hypo.y
        else:
            assert self.grid != None, 'Grid must be assigned'
            X = self.grid.x - self.hypo.x
            Y = self.grid.y - self.hypo.y
        
        # Distances
        self.hypo_dist = np.sqrt(X*X+Y*Y)

        # All done
        return


    def gridPadding(self):
        '''
        Defines a mirror padding grid from the original grid. 
        Padding velocities are equal to those on the edges of the grid
        '''
        
        # Allocate padding grid attributes
        pad_grid_size = (self.grid.ny+2,self.grid.nx+2)
        X  = np.zeros(pad_grid_size,dtype='float64')
        Y  = np.zeros(pad_grid_size,dtype='float64')
        Vr = np.zeros(pad_grid_size,dtype='float64')
        
        # Copy grid into the padding grid
        X[1:-1,1:-1]  = self.grid.x.copy()
        Y[1:-1,1:-1]  = self.grid.y.copy()
        Vr[1:-1,1:-1] = self.grid.vr.copy()
        
        # Padding grid edges for X coordinates
        X[1:-1, 0] = self.grid.x[:,0]  - self.grid.h
        X[1:-1,-1] = self.grid.x[:,-1] + self.grid.h
        X[0 ,:]    = X[ 1,:]
        X[-1,:]    = X[-2,:]
        # Padding grid edges for Y coordinates
        Y[0 ,1:-1] = self.grid.y[ 0,:] - self.grid.h
        Y[-1,1:-1] = self.grid.y[-1,:] + self.grid.h
        Y[:, 0]    = Y[:, 1]
        Y[:,-1]    = Y[:,-2]
        # Padding grid edges for rupture velocities
        Vr[ 0,:]   = Vr[1 ,:]
        Vr[-1,:]   = Vr[-2,:]
        Vr[:, 0]   = Vr[:, 1]
        Vr[:,-1]   = Vr[:,-2]
        
        # Instanciate padding grid
        self.pad_grid = Grid(X,Y,Vr,'padding')
        
        # All done
        return
        

    def initT0(self):
        '''
        Initialize T0 values
        '''
        
        # Check sr
        assert self.sr != None, 'Slowness matrix must be assigned' 

        # Initialize to some large value
        self.t0 = np.ones(self.pad_grid.shape)*1.0e6
        
        # Get indexes of the 4 smallest distances in the Grid
        hypo_dist = self.hypo_dist.flatten()
        hypo_dist.sort()
        i,j = np.where(self.hypo_dist<=hypo_dist[3])
        i = i[:4]
        j = j[:4]        

        # For these grid points, set t0 to hypo_dist/vr
        self.t0[i,j] = self.hypo_dist[i,j]*self.sr[i,j]
        
        # All done 
        return


    def eq_solve(self,a,b,f,h):
        '''
        This solves [(x-a)^+]^2 + [(x-b)^+]^2 = f^2 * h^2
                      | z, if z>0
        where (z)^+ = | 
                      | 0, if z<=0

        The (unique) solution is given by
                 | min(a,b) + f*h,                        if |a-b|>= f*h
        x_sol =  |
                 |0.5 * [a+b+sqrt(2*f^2*h^2 - (a-b)^2)],  if |a-b| < f*h
        '''    
        
        if np.abs(a-b) >= f*h:
            min_ab = a
            if a>=b:
                min_ab=b
            x_sol = min_ab + f*h
        else:
            amb = a-b
            x_sol = a+b+np.sqrt(2*f*f*h*h-amb*amb)
            x_sol /= 2.;
    
        # All done
        return x_sol


    def upwind(self,i,j):
        '''
        Perform upwind difference
        Args:
            * i: grid index (row)
            * j: grid index (column)
        '''

        # Get indexes of neighbors
        i1 = i-1; i2 = i+1
        j1 = j-1; j2 = j+1
        if i1 < 0:
            i1 = 0
        if i2 > self.pad_grid.ny-1:
            i2 = self.pad_grid.ny-1
        if j1 < 0:
            j1 = 0
        if j2 > self.pad_grid.nx-1:
            j2 = self.pad_grid.nx-1
        
        # get min value of neighbors
        t0_xmin = self.t0[i1,j]
        t0_ymin = self.t0[i ,j1]
        if self.t0[i1,j] >= self.t0[i2,j]:
            t0_xmin = self.t0[i2,j]
        if self.t0[i,j1] >= self.t0[i,j2]:
            t0_ymin = self.t0[i,j2]

        # Solve the equation locally
        t0_new = self.eq_solve(t0_xmin,t0_ymin,self.sr[i,j],self.pad_grid.h)
        
        # Update the t0 grid
        if self.t0[i,j]>t0_new:
            self.t0[i,j] = t0_new
        
        # All done
        return

    
    def iterGaussSeidel(self,sweep_order_x,sweep_order_y):
        '''
        Gauss-Seidel iterations with alternating sweeping orderings
        '''

        # set X order
        if sweep_order_x>0:
            ixs = np.arange(self.pad_grid.nx)
        else:
            ixs = np.arange(self.pad_grid.nx-1,-1,-1)
        # set Y order
        if sweep_order_y>0:
            iys = np.arange(self.pad_grid.ny)
        else:
            iys = np.arange(self.pad_grid.ny-1,-1,-1)
        
        # Loop over grid points
        for i in iys:
            for j in ixs:
                self.upwind(i,j)

        # All done
        return


    def getT0(self,point_x,point_y):
        '''
        Get T0 value in point_x and point_y
        Args:
            * point_x,point_y (1D array or list)
        '''
        
        # Check length of point_x and point_y
        assert type(point_x)==type(point_y), 'point_x and point_y must have same type'
        assert len(point_x)==len(point_y), 'point_x and point_y must have same length'

        # Find t0 for each point
        point_t0 = []
        for x,y in zip(point_x,point_y):
            dx = np.abs(self.pad_grid.x-x)
            dy = np.abs(self.pad_grid.y-y)
            #dist = np.sqrt(dx*dx+dy*dy)
            #i,j=np.where(dist==dist.min())
            i,j = np.where( (dx==dx.min()) & (dy==dy.min()) )
            point_t0.append(self.t0[i[0],j[0]])
        if type(point_x)==np.ndarray:
            point_t0 = np.array(point_t0)

        # Check length of output array
        assert len(point_t0)==len(point_x), 'Some points are missing'

        # All done
        return point_t0

    def interpolateT0(self,point_x,point_y):
        '''
        Get T0 value in point_x and point_y
        Args:
            * point_x,point_y (1D array or list)
        '''
        
        # Check length of point_x and point_y
        assert type(point_x)==type(point_y), 'point_x and point_y must have same type'
        assert len(point_x)==len(point_y), 'point_x and point_y must have same length'

        # Find t0 for each point
        dgx = self.pad_grid.x[0,1]-self.pad_grid.x[0,0]
        dgy = self.pad_grid.y[1,0]-self.pad_grid.y[0,0]
        assert np.round(dgx,2) == np.round(dgy,2), 'Fast sweeping mesh must be uniform (%f vs %f)'%(dgx,dgy)
        point_t0 = []
        for x,y in zip(point_x,point_y):
            dx  = x-self.pad_grid.x
            dy  = y-self.pad_grid.y
            dxa = np.abs(dx)
            dya = np.abs(dy)            
            i,j = np.where( (dxa==dxa.min()) & (dya==dya.min()))
            i = i[0]
            j = j[0]
            if dx[i,j]<0:
                j -= 1
            if dy[i,j]<0:
                i -= 1
            xr = dx[i,j]/dgx
            yr = dy[i,j]/dgy
            f11 = self.t0[i  ,j  ]
            f21 = self.t0[i  ,j+1]
            f12 = self.t0[i+1,j  ]
            f22 = self.t0[i+1,j+1]
            point_t0.append(f11*(1.0-xr)*(1.0-yr) + f21*xr*(1.0-yr) + f12*(1.0-xr)*yr + f22*xr*yr)
            
        if type(point_x)==np.ndarray:
            point_t0 = np.array(point_t0)

        # Check length of output array
        assert len(point_t0)==len(point_x), 'Some points are missing'

        # All done
        return point_t0



    def getT0FromFault(self,fault,g_x,g_y,g_z):
        '''
        Get T0 value in point_x and point_y (only working for planar faults)
        Args:
            * fault: fault object
            * g_x,g_y,g_z: coordinates of points on the fault (UTM)
        '''

        # Check strike/dip
        assert fault.f_strike != None, 'Fault strike must be assigned'
        assert fault.f_dip    != None, 'Fault dip must be assigned'

        # Check length of point_x and point_y
        assert type(g_x)==type(g_y)==type(g_z), 'x, y and z must have same type'
        if type(g_x)==list or type(g_x)==np.ndarray:
            assert len(g_x)==len(g_y)==len(g_z),    'x, y and z must have same length'
            scalar = False
        else:
            g_x = [g_x]
            g_y = [g_y]
            g_z = [g_z]
            scalar = True

        # Find t0 for each point
        point_x = []
        point_y = []
        for x,y,z in zip(g_x,g_y,g_z):
            x -= fault.hypo_x
            y -= fault.hypo_y
            z -= fault.hypo_z
            # Get fault coordinates
            point_x.append(z / np.sin(fault.f_dip))
            point_y.append(x * np.sin(fault.f_strike) + y * np.cos(fault.f_strike))
        
        # Get T0s
        point_t0 = self.getT0(point_x,point_y)

        # All done
        if scalar:
            return point_t0[0]
        else:
            return point_t0


    def printT0(self,caption=None,pad=False,t0_format='%7.2f'):
        '''
        Display T0 values
        Args:
            * caption: print caption before plotting the matrix
            * pad: if pad==True, plot all pad_grid
            * t0_format: print format 
        '''

        if caption!=None:
            print(caption)
        
        if pad:
            bx = 0
            by = 0
            ex = self.pad_grid.nx
            ey = self.pad_grid.ny
        else:
            bx = 1
            by = 1
            ex = self.pad_grid.nx-1
            ey = self.pad_grid.ny-1
            
        for i in range(by,ey):
            for j in range(bx,ex):
                print(t0_format%(self.t0[i,j])),
            print('')        
            
        # All done
        return


    def fastSweep(self,num_iter=4,verbose=False):
        '''
        Calculate rupture times using the fast sweeping method 
        Args:
            * num_iter: number of iteration
            * verbose: verbose mode (True or False)
        '''
        
        # Define a padding grid
        self.gridPadding()

        # Calculate hypocentral distances
        self.calcHypoDist(pad=True)
        
        # Defines rupture slowness vector
        self.sr = 1./self.pad_grid.vr
        
        # Initialize T0
        print('Initialize T0')
        self.initT0()
        if verbose:
            self.printT0('initialized T0: ',pad=True)
        
        # Main loop
        print('Fast Sweeping')
        for iter in range(num_iter):
            if verbose:
                print('iteration: {}'.format(iter+1))
            for order_x,order_y in zip([+1,+1,-1,-1],[+1,-1,-1,+1]):  
                self.iterGaussSeidel(order_x,order_y)
                if verbose:
                    self.printT0('order {} {}'.format(order_x,order_y),pad=True)
        # All done
        return
        
        
    def copy(self):
        '''
        Returns a copy of the FastSweep solver
        '''
        
        return copy.deepcopy(self)
