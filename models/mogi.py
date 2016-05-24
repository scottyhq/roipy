"""
Functions for forward volcano-geodesy analytic models

Author: Scott Henderson
Date: 8/31/2012

TODO:
-benchmark codes against paper results
-test for georeferenced coordinate grids?
-add function to convert to los
-add sphinx docstrings
"""
import numpy as np
import matplotlib.pyplot as plt
import roipy.tools


# =====================
# Inverse Models
# =====================
def invert_fullres(X,Y,look,head,xcen,ycen,depth,dV,nu):
    """
    Adjust arguments of calc_mogi to work with scipy.omptimize.curvefit convention
    Assumes UTM input for X and Y
    """
    Xshift = X - xcen
    Yshift = Y - ycen

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(Xshift,Yshift) # surface angle and radial distance
    R = np.hypot(depth,rho) # radial distance from source

    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3
    uz = C * depth / R**3
    ux, uy = pol2cart(th, ur)
    dataVec = np.dstack([ux, uy, uz]) #shape = (1107, 890, 3)

    # Get LOS transform
    cart2los = get_cart2los(look,head)
    los = np.sum(dataVec * cart2los, axis=2)

    return los.ravel() #flattened arrays required by scipy.optimize

def invert_resample(easting,northing,cart2los,xcen,ycen,depth,dV,nu=0.25):
    """
    Adjust arguments of mogi.forward() to work with scipy.omptimize.curvefit convention
    Assumes UTM input for X and Y
    """
    Xshift = easting - xcen
    Yshift = northing - ycen

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(Xshift,Yshift) # surface angle and radial distance
    R = np.hypot(depth,rho) # radial distance from source

    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3
    uz = C * depth / R**3
    ux, uy = pol2cart(th, ur)
    dataVec = np.dstack([ux, uy, uz]) #shape = (1107, 890, 3)

    los = np.sum(dataVec * cart2los, axis=2)

    return los.ravel() #flattened arrays required by scipy.optimize




# =====================
# Forward Models
# =====================

def forward(x,y,xoff=0,yoff=0,d=3e3,dV=1e6,nu=0.25,output='cyl'):
    """
    Calculates surface deformation based on point source

    References: Mogi 1958, Segall 2010 p.203

    Args:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)

    Kwargs:
    -----------------
    xoff: y-offset of point source epicenter (m)
    yoff: y-offset of point source epicenter (m)
    d: depth to point (m)
    dV: change in volume (m^3)
    nu: poisson's ratio for medium
    output: 'cart' (cartesian), 'cyl' (cylindrical)

    Returns:
    -------
    (ux, uy, uz) if output='cart'
    (ur,uz) if output='cyl'


    Examples:
    --------
    forward(0,0) #gives uz_max for default source parameters
    """
    # Center coordinate grid on point source
    x = x - xoff
    y = y - yoff

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(x,y) # surface angle and radial distance
    R = np.hypot(d,rho) # radial distance from source

    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3    # horizontal displacement, m
    uz = C * d / R**3   # vertical displacement, m

    # Convert surface cylindrical to cartesian
    if output == 'cyl':
        return ur, uz
    else:
        ux, uy = pol2cart(th, ur)
        return ux, uy, uz


def forward_dp(x,y,xoff=0,yoff=0,d=3e3,a=500,dP=100e6,mu=4e9,nu=0.25,output='cyl'):
    """
    dP instead of dV, NOTE: dV = pi * dP * a**3 / mu
    981747.7 ~ 1e6
    """
    dV = np.pi * dP * a**3 / mu
    if output == 'cyl':
        return forward(x,y,xoff,yoff,d,dV,nu,output='cyl')
    else:
        return forward(x,y,xoff,yoff,d,dV,nu)



def calc_linmax(x,y,tn,xoff=0,yoff=0,d=3e3,a=500.0,dP=100e6,mu=4e9,nu=0.25,output='cyl'):
    """ Solution for spherical source in a Maxwell Solid viscoelastic halfspace
    based on Bonafede & Ferrari 2009 (Equation 14).

    Simplified equations for z=0 free surface, for instantaneous dP at t=0.

    Note that displacements are in m, but inputs are normalized time.

    Required arguments:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)
    t: normalized time (t/tau)

    Keyword arguments:
    -----------------
    xoff: y-offset of point source epicenter (m)
    yoff: y-offset of point source epicenter (m)
    d: depth to point (m)
    dP: change in pressure at t=0 (Pa)
    K: short-term elastic incompressibility () NOTE=(5/3)mu in poisson approximation
    mu: short-term elastic rigidity (Pa)
    eta: effectice viscosity (Pa*s)
    output: 'cart' (cartesian), 'cyl' (cylindrical)

    """
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff

    # convert to surface cylindrical coordinates
    th, r = cart2pol(x,y)
    R = np.hypot(d,r)

    # Common variables
    K = (2.0*mu*(1 + nu)) / (3*(1 - 2*nu))
    A = 1 + 2.0*(mu/K)
    B = (2.0*mu**2) / (K*(3*K + mu))
    alpha = (3.0*K + mu) / (3*K)
    #tau = nu / mu #maxwell relaxation time
    #tau_a = alpha * tau
    #Rstar = R #NOTE, this is only true for solution on free surface

    term = 1 + A - (B*np.exp(-tn/alpha)) + (2*tn)

    C = (dP * a**3) / (4*mu)
    ur = C * (r/(R**3)) * term
    uz = C * (d/(R**3)) * term

    print('uz_max = {:.4f}'.format(uz.max()))
    print('ur_max = {:.4f}'.format(ur.max()))

    # Convert surface cylindrical to cartesian
    if output == 'cart':
        ux, uy = pol2cart(th, ur)
        return ux, uy, uz
    elif output == 'cyl':
        return ur, uz


def calc_linmax_dPt(tn,dVdt,xoff=0,yoff=0,d=3e3,a=500.0,dP=100e6,mu=4e9,
                         nu=0.25,output='cyl'):
    """ Instead of constant pressure, have pressure determined by a constant
    supply rate of magma

    From Bonafede 2009 Equation 16

    NOTE: Only Uzmax b/c full solution not given in Segaall
    """
    K = (2.0*mu*(1 + nu)) / (3*(1 - 2*nu))
    tau = nu / mu
    #tauB = ((3*K + 4*mu) / (3*K)) * tau
    tauA = ((3*K + mu) / (3*K)) * tau

    C = (dVdt*tau) / (2*pi*d**2)
    term1 = t/tau
    term2 = (mu / (3*K)) * (1 - np.exp(-t/tauA))
    uzmax = C * (term1 - term2)
    return uzmax


def calc_genmax(x,y,t,xoff=0,yoff=0,d=4e3,dP=100e6,a=700,nu=0.25,G=30e9,
                    mu1=0.5,eta=2e16,output='cyl',**kwargs):
    """ Solution for spherical source in a generalized maxwell viscoelastic
    halfspace based on Del Negro et al 2009.

    Required arguments:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)
    t: time (s)

    Keyword arguments:
    -----------------
    xoff: y-offset of point source epicenter (m)
    yoff: y-offset of point source epicenter (m)
    d: depth to point (m)
    dV: change in volume (m^3)
    K: bulk modulus (constant b/c incompressible)
    E: Young's moduls
    G: total shear modulus (Gpa)
    mu0: fractional shear modulus (spring part)
    mu1: fractional shear modulus (dashpot part)
    eta: viscosity (Pa s)
    output: 'cart' (cartesian), 'cyl' (cylindrical)

    """
    #WARNING: mu0 != 0
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff

    # convert to surface cylindrical coordinates
    #th, r = cart2pol(x,y)
    r = np.hypot(x,y) #surface radial distance
    R = np.hypot(d,r) #radial distance from source center

    # Calculate displacements
    #E = 2.0 * G * (1+nu)
    #K = E / (3.0* (1 - 2*nu)) #bulk modulus = (2/3)*E if poisson solid
    K = (2.0*G*(1+nu)) / (3*(1-(2*nu)))
    mu0 = 1.0 - mu1
    alpha = (3.0*K) + G #recurring terms
    beta = (3.0*K) + (G*mu0)

    # Maxwell times
    try:
        tau0 = eta / (G*mu1)
    except:
        tau0 = np.inf
    tau1 = (alpha / beta) * tau0
    tau2 = tau0 / mu0

    #print('relaxation times:\nT0={}\nT1={}\nT2={}'.format(tau0,tau1,tau2))

    term1 = ((3.0*K + 4*G*mu0) / (mu0*beta))
    term2 = ((3.0 * G**2 * np.exp(-t/tau1))*(1-mu0)) / (beta*alpha)
    term3 = ((1.0/mu0) - 1) * np.exp(-t/tau2)

    A = (1.0/(2*G)) * (term1 - term2 - term3)
    C = (dP * a**3) / R**3
    ur = C * A * r
    uz = C * A * d

    return ur, uz



def calc_mctigue(x,y,xoff=0,yoff=0,d=3e3,dP=10e6,a=1500.0,nu=0.25,mu=4e9,
                 terms=1, output='cyl'):
    """
    3d displacement field from dislocation point source (McTigue, 1987)
    Caution: analysis done in non-dimensional units!
    see also Segall Ch7 p207
    Same as forward, except change in pressure and chamber radius are specified

    Keyword arguments:
    ------------------
    xoff: y-offset of point source epicenter (m)
    yoff: y-offset of point source epicenter (m)
    d: depth to point (m)
    rad: chamber radius (m)
    dV: change in volume (m^3)
    dP: change in pressure (Pa)
    nu: poisson's ratio for medium
    mu: shear modulus for medium (Pa)
    order: highest order term to include (up to 2)
    output: 'cart' (cartesian), 'cyl' (cylindrical)

    Set terms=1 to reduce to Mogi Solution
    NOTE: eps**6 term only significant if eps > 0.5
    """
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff

    # dimensionless scaling term
    scale = dP * d / mu
    eps = a / d #NOTE eps = 0.5 # mctigue fig3

    # convert to surface cylindrical coordinates
    th, r = cart2pol(x,y)
    r = r / d #dimensionless radial distance

    # 1st order mctigue is essentially Mogi solution
    uz = eps**3 * ((1-nu) * (1 / np.hypot(r,1)**3))
    ur = eps**3 * ((1-nu) * (r / np.hypot(r,1)**3))

    # 2nd order term
    if terms==2:
        print('adding eps**6 term')
        A = ((1 - nu) * (1 + nu)) / (2 * (7 - 5*nu))
        B = (15 * (2 - nu) * (1 - nu)) / (4 * (7 - 5*nu))
        uz2 =  -eps**6 * ((A * (1 / np.hypot(r,1)**3)) - (B * (1 / np.hypot(r,1)**5)))
        ur2 =  -eps**6 * ((A * (r / np.hypot(r,1)**3)) - (B * (r / np.hypot(r,1)**5)))
        uz += uz2
        ur += ur2

    # Convert back to dimensional variables
    uz = uz * scale
    ur = ur * scale

    # Convert surface cylindrical to cartesian
    if output == 'cart':
        ux, uy = pol2cart(th, ur)
        return ux, uy, uz
    elif output == 'cyl':
        return ur, uz


def calc_viscoshell(x,y,t,xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,dP=100e6,
                             mu=30e9,nu=0.25,eta=2e16,output='cyl'):
    """ Spherical Source surronded by a viscoelastic shell in an elastic halfspace
    Derivation of equations 7.105 in Segall Ch.7 p245
    NOTE: good approximation if duration of intrusion << relation time tR

    Required arguments:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)
    t: time (s)

    Keyword arguments:
    -----------------
    Same as forward_dp() plus:
    a: inner chamber radius
    b: extent of viscoelastic region
    eta: viscosity
    """
    # characteristic relaxation time
    tR = (3*eta*(1-nu)*b**3) / (mu*(1+nu)*a**3)

    #avoid ZeroDivisionError
    if tR == 0:
        scale = 1
    else:
        scale = (np.exp(-t/tR) + ((b/a)**3)*(1.0 - np.exp(-t/tR)))

    #rho = np.hypot(x,y)
    rho = np.hypot(x,y) / d #dimensionless!
    uz = (((1-nu)*dP*a**3) / (mu*d**2)) * scale * (1 + rho**2)**(-3/2.0)
    ur = rho * uz

    #print('uz_max = {:.4f}'.format(uz.max()))
    #print('ur_max = {:.4f}'.format(ur.max()))

    return ur, uz


def calc_viscoshell_dPt(x,y,t,P0,tS,xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,
                             mu=30e9,nu=0.25,eta=2e16,output='cyl'):
    """
    Viscoelastic shell with a exponentially decaying pressure source
    from Segall 2010

    P0 = initial pressure at t0
    tS = relaxation time of pressure source
    NOTE: tS & tR should not equal zero, or else division by zero
    NOTE: eq. 7.113 has an error, when tS=tR, should have t/tR in the formula
    """
    # viscoelastic relaxation time
    tR = (3*eta*(1-nu)*b**3) / (mu*(1+nu)*a**3)
    # pressure source relaxation time = ts

    # Pressure history
    P = P0 * (1 - np.exp(-t/tS))

    # Lambda coefficient
    Lambda = ((tS/tR)*(b/a)**3 - 1) / ((tS/tR) - 1)

    # Scale term
    if tR == tS:
        scale = ( (1-(b/a)**3) * (t/tR) * np.exp(-t/tR) +
                 ((b/a)**3) * (1.0 - np.exp(-t/tR)))
    else:
        scale = (Lambda*np.exp(-t/tR) -
                 Lambda*np.exp(-t/tS) +
                 ((b/a)**3)*(1.0 - np.exp(-t/tR)))

    #rho = np.hypot(x,y)
    rho = np.hypot(x,y) / d #Dimensionless radius!

    uz = (((1-nu)*P0*a**3) / (mu*d**2)) * (1 + rho**2)**(-3/2.0) * scale
    ur = rho * uz

    print('uz_max = {:.4f}'.format(uz.max()))
    print('ur_max = {:.4f}'.format(ur.max()))

    return ur, uz, P




# ===================
# Utility Functions
# ===================
def cart2pol(x1,x2):
    #theta = np.arctan(x2/x1)
    theta = np.arctan2(x2,x1) #sign matters -SH
    r = np.hypot(x2,x1)
    return theta, r


def pol2cart(theta,r):
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    return x1,x2

def shift_utm(X,Y,xcen,ycen):
    ''' Avoid large numbers in UTM grid by creating local (0,0) origin '''
    x0 = X.min()
    y0 = Y.min()

    X = X - x0
    Y = Y - y0
    xcen = xcen - x0
    ycen = ycen - y0

    return X,Y,xcen,ycen


def get_cart2los(look,head):
    ''' coefficients for projecting cartesian displacements into LOS vector '''
    incidence = np.deg2rad(look)
    heading = np.deg2rad(head)

    EW2los = np.sin(heading) * np.sin(incidence)
    NS2los = np.cos(heading) * np.sin(incidence)
    Z2los = -np.cos(incidence)

    cart2los = np.dstack([EW2los, NS2los, Z2los])

    return cart2los
