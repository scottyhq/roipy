"""
functions for adding and characterizing noise in interferograms
"""
import numpy as np


# ======================================================
# ADDING "NOISE"
# ======================================================
def add_noise(data, type='white', min=-0.1, max=0.1, mean=0.0, std=0.05, Lc=100,
              cortype='powerlaw', topotype='linear', dem='radar_32rlks.unw',
              toposlope=-1e-4):
    """Generate various type of noise for an interferogram
    mean & std in meters. toposlope --> is the slope from a plot of phase vs.
    elevation. -1e-4 means the deformation drops from 0 to -3cm as elevation
    increases from 0 to 1km. Note, min max values scaled by random numbers
    generated between 0 & 1, therefore given in meters"""        
    
    # Uniform velocities
    if type == 'uniform':
        noise = min + (max-min)*np.random.random_sample(data.shape)

    # Gaussian white noise
    if type == 'white':
        noise = np.random.normal(mean,std,data.shape)
    
    # Correlated Noise (from Lohman & Simons 2005)
    if type == 'correlated':
        covd = make_covariance(data.shape,type=cortype,Lc=Lc,std=std)
        noise = make_correlated_noise(covd)
    
    # Signal Correlated with Topography
    if type == 'topo':
        #NOTE: should only be applied as a constant over small regions?
        #NOTE: there is probably a more sophisticated way of going about this
        elev = rp.tools.load_half(dem)
        b = np.median(elev)
        noise = b + (toposlope * elev)
    
    return noise


def make_covariance(shape, type='powerlaw', Lc=100, std=0.01):
    """ NOTE: this takes a long time. Lc in pixels. Modified from Rowena
    Lohman's scripts"""
    var = std**2
    nr, nc = shape
    X,Y = np.indices(shape) #like meshgrid
    X = X.flatten()
    Y = Y.flatten()
   
    # FULL COVARIANCE matrix.. it's a beast
    covd = var * np.eye(X.size)
    for i in range(X.size):
        for j in range(i+1, Y.size):
            distance = np.hypot(X[i] - X[j], Y[i] - Y[j])
            #C[r,c] = exp(-distance/Lc)) #power-law covariance matrix
            covd[i,j] = var * (10**(-distance/Lc))
            covd[j,i] = covd[i,j] #symmetric
    
    return covd
    
    
    
def make_correlated_noise(covd, mean=0.0, std=0.1):
    """ given uncorrelated noise and a covariance model, make an array of
    correlated noise. 0.1 = 10cm """
    u, v = np.linalg.eig(covd) #NOTE: how to interpret real & imaginary parts, slow for big arrays!
    noise = np.random.normal(mean, std, X.shape[0]) # Uncorrelated Noise
    noise_cor = np.dot(np.dot(v,np.sqrt(np.diag(u))),noise)
    noise_cor = np.real(noise_cor)
    #reshape to image dimensions
    #noise = noise.reshape(shape)
    noise_cor = noise_cor.reshape(shape)
    
    return noise_cor



def add_ramp(Interferogram, type='quadratic',std=0.02):
    """ generate a long-wavelength surface in the interferogram
    type = 'linear','dc','quadratic'"""
    #add_ramp('quad')
    x,y = np.indices(ig.Shape)
    x = x.flatten()
    y = y.flatten()
    
    # 6-parameter 2nd order quadratic
    if type == 'quadratic':
        G = np.array([np.ones(x.size), x, y, x**2, x*y, y**2])
        # generate random quadratic surface coefficients
        #m = np.array([1e-2, 1e-5, 1e-5, 1e-7, 1e-7, 1e-7]) * np.random.random(6)
        # Rowena's automated approach
        randos = np.random.random(6) / 6 #divide by 6 necessary?
        m = np.array([std,
                      std/np.std(x),
                      std/np.std(y),
                      std/np.std(x**2),
                      std/np.std(x*y),
                      std/np.std(y**2)])
        m = randos * m
    
    # 3-parameter linear ramp
    if type == 'linear':
        G = np.array([np.ones(ig.size), x, y])
        randos = np.random.random(3) / 3 #divide by 6 necessary?
        m = np.array([std,
                      std/np.std(x),
                      std/np.std(y),
                      std/np.std(x**2),
                      std/np.std(x*y),
                      std/np.std(y**2)])
        m = randos * m
        
    # dc-shift
    if type == 'dc':
        G = np.ones(ig.size)
        m = std*np.random.random(1)
        
    ramp = np.dot(G.T,m) 

    return surface



def add_interseismic(Interferogram):
    """ For long-track interferograms subtract estimated interseismic signal
    See Fournier et al 2011
    """
    print('Not implemented yet')
    #figure3 shows interseismic ramps on the order of 3mm



# ======================================================
# CHARACTERIZING "NOISE"
# ======================================================
def calc_covariance(Interferogram, array=None):
    """ return the FULL covariance matrix for an interferogram """
    self = Interferogram
    if array == None:
        array = rp.tools.load_half(self,2) #load phase by default
        
    
    
    
    
    
    