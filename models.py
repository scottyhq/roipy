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
    

# ===================
# Forward Models
# =====================
def calc_mogi(x,y,xoff=0,yoff=0,d=3e3,dV=1e6,nu=0.25,output='cyl'):
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
    calc_mogi(0,0) #gives uz_max for default source parameters
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


def calc_mogi_dp(x,y,xoff=0,yoff=0,d=3e3,a=500,dP=100e6,mu=4e9,nu=0.25,output='cyl'):
    """
    dP instead of dV, NOTE: dV = pi * dP * a**3 / mu
    981747.7 ~ 1e6
    """
    
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff
    
    # convert to surface cylindrical coordinates
    th, rho = cart2pol(x,y) 
    #rho = np.hypot(x,y) 
    R = np.hypot(d,rho) #radial distance from source center 
    
    # Mogi displacement calculation
    C = (1-nu) * dP * a**3 / mu
    ur = C * rho / R**3    
    uz = C * d / R**3
    
    print 'uz_max = {:.4f}'.format(np.abs(uz).max()) 
    print 'ur_max = {:.4f}'.format(np.abs(ur).max())
    
    # Convert surface cylindrical to cartesian
    if output == 'cart':
        ux, uy = pol2cart(th, ur)
        return ux, uy, uz
    elif output == 'cyl':
        return ur, uz



def calc_mogi_linmax(x,y,tn,xoff=0,yoff=0,d=3e3,a=500.0,dP=100e6,mu=4e9,nu=0.25,output='cyl'):
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
    
    print 'uz_max = {:.4f}'.format(uz.max())
    print 'ur_max = {:.4f}'.format(ur.max())
    
    # Convert surface cylindrical to cartesian
    if output == 'cart':
        ux, uy = pol2cart(th, ur)
        return ux, uy, uz
    elif output == 'cyl':
        return ur, uz


def calc_mogi_linmax_dPt(tn,dVdt,xoff=0,yoff=0,d=3e3,a=500.0,dP=100e6,mu=4e9,
                         nu=0.25,output='cyl'):
    """ Instead of constant pressure, have pressure determined by a constant
    supply rate of magma
    
    From Bonafede 2009 Equation 16
    
    NOTE: can only get Uzmax b/c full solution not given
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


def calc_mogi_genmax(x,y,t,xoff=0,yoff=0,d=4e3,dP=100e6,a=700,nu=0.25,G=30e9,
                    mu1=0.5,eta=2e16,output='cyl'):
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
    
    #print 'relaxation times:\nT0={}\nT1={}\nT2={}'.format(tau0,tau1,tau2)
    
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
    Same as calc_mogi, except change in pressure and chamber radius are specified
    
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
        print 'adding eps**6 term'
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


def calc_mogi_viscoshell(x,y,t,xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,dP=100e6,
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
    Same as calc_mogi_dp() plus:
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
    
    print 'uz_max = {:.4f}'.format(uz.max()) #sub-millimeter precision at best!
    print 'ur_max = {:.4f}'.format(ur.max())
    
    return ur, uz
    
    
def calc_mogi_viscoshell_dPt(x,y,t,P0,tS,xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,
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
    
    print 'uz_max = {:.4f}'.format(uz.max()) #sub-millimeter precision at best!
    print 'ur_max = {:.4f}'.format(ur.max())
    
    return ur, uz, P


def spheroid(a,b,c,matrl,phi,theta,P):
    """
    Function called by calc_yang
    """
    
    #Unpack material parameters
    lamda, mu, nu = matrl
    
    # Model expressions
    ac    = (a - c) / (a + c)
    L1    = np.log(ac)
    #L1 = np.log((a-c)/(a+c))
    iia   = 2/a/c**2 + L1/c**3
    iiaa  = 2/3/a**3/c**2 + 2/a/c**4 + L1/c**5
    
    coef1 = -2*np.pi*a*b**2
    Ia    = coef1*iia
    Iaa   = coef1*iiaa
    
    u = 8*np.pi*(1-nu)
    Q = 3/u
    R = (1-2*nu)/u
    
    a11 = 2*R*(Ia-4*np.pi)
    a12 = -2*R*(Ia+4*np.pi)
    a21 = Q*a**2*Iaa + R*Ia - 1
    a22 = -(Q*a**2*Iaa + Ia*(2*R-Q))
    
    coef2 = 3*lamda+2*mu
    w     = 1/(a11*a22-a12*a21)
    e11   = (3*a22-a12)*P*w/coef2
    e22   = (a11-3*a21)*P*w/coef2
    Pdila = 2*mu*(e11-e22)
    Pstar = lamda*e11 + 2*(lamda+mu)*e22
    a1    = -2*b**2*Pdila
    # !PL version had (1-nu) in 2nd term!
    b1    = 3*b**2*Pdila/c**2 + 2*(1-2*nu)*Pstar  
    
    sph = np.array([a,b,c,phi,theta,Pstar,Pdila,a1,b1,P])
    
    return sph


def yang(x,y,z0,sph,xi,matrl,e_theta,coeffs):
    """
    Translated from matlab script:
    %function [u1,u2,u3]=yang(sph,xi,z0,x,y,matrl,e_theta,coeffs)
    % Calculate the double force (star) and dilatation (dila) displacements U
    % for a SPHEROIDAL pressure source in an elastic halfspace 
    % (based on Yang et al., vol 93, JGR, 4249-4257, 1988) with arbitrary plunge (theta)
    % of the long axis of the spheroid (theta = 90, prolate theta = 0, oblate).
    % Evaluate at for xi.
    %
    % Inputs: theta: dip angle of source
    %             P: pressure change in magma chamber
    %            a: semimajor axis of spheroid
    %             b: semiminor axis of spheriod
    %            xi: evaluate integrals at +- c
    %            z0: depth of source
    %             x: x location of point on surface
    %             y: y location of point on surface
    % Output: rd: calculated range displacement
    % NOTE: the x, y locations assume a source at origin
    % ALSO: the units need to be in mks units so input x, y, and z0
    %       in km will be changed into meters
    % NOTE: In the expressions of Yang et al. the spheroid dips parallel to the y axis
    %       at x=0. We will assume (initially) that it is also centered at y=0.
    """
    # extract spheroid information
    a,b,c,phi,theta,Pstar,Pdila,a1,b1,P = sph

    sinth, costh = e_theta
    
    # Poisson's ratio, Young's modulus, and the Lame coeffiecents mu and lamda
    nu    = matrl[2]
    nu4   = coeffs[1]
    nu2   = 1 - 2*nu
    nu1   = 1 - nu
    coeff = a * b**2 / c**3 * coeffs[0]
    
    # Introduce new coordinates and parameters (Yang et al., 1988, page 4251):
    xi2   = xi * costh
    xi3   = xi * sinth
    y0    = 0.0
    z     = 0.0
    x1    = x
    x2    = y - y0
    x3    = z - z0
    xbar3 = z + z0
    y1    = x1
    y2    = x2 - xi2
    y3    = x3 - xi3
    ybar3 = xbar3 + xi3
    r2    = x2*sinth - x3*costh
    q2    = x2*sinth + xbar3*costh
    r3    = x2*costh + x3*sinth
    q3    = -x2*costh + xbar3*sinth
    rbar3 = r3 - xi
    qbar3 = q3 + xi
    R1    = np.sqrt(y1**2 + y2**2 + y3**2)
    R2    = np.sqrt(y1**2 + y2**2 + ybar3**2)
    
    # C0 = y0*costh+z0*sinth  % check this!!!!!!!!!!!!!
    C0 = z0 / sinth
    
    betatop    = (costh * q2 + (1 + sinth)*(R2 + qbar3))
    betabottom = costh * y1
    
    #atnbeta = atan2(betatop,betabottom)
    nz          = (np.abs(betabottom) != 0)
    atnbeta     = np.pi/2 * sign(betatop)
    atnbeta[nz] = np.atan(betatop[nz] / betabottom[nz])
    
    # Set up other parameters for dipping spheroid (Yang et al., 1988, page 4252):
    # precalculate some repeatedly used natural logs:
    Rr  = R1 + rbar3
    Rq  = R2 + qbar3
    Ry  = R2 + ybar3
    lRr = np.log(Rr)
    lRq = np.log(Rq)
    lRy = np.log(Ry)
    
    A1star    =  a1 / (R1*Rr) + b1*(lRr + (r3 + xi)/Rr)
    Abar1star = -a1 / (R2*Rq) - b1*(lRq + (q3 - xi)/Rq)
    A1        = xi/R1 + lRr
    Abar1     = xi/R2 - lRq
    A2        = R1 - r3*lRr
    Abar2     = R2 - q3*lRq
    A3        = xi*rbar3/R1 + R1
    Abar3     = xi*qbar3/R2 - R2
    
    B      = xi*(xi + C0)/R2 - Abar2 - C0*lRq
    Bstar  = a1/R1 + 2*b1*A2 + coeffs[1]*(a1/R2 + 2*b1*Abar2)
    F1     = 0.0
    F1star = 0.0
    F2     = 0.0
    F2star = 0.0
    
    # Skip if displacement calculated at surface (z = 0)
    if z != 0:
        F1 = (-2*sinth*z* (xi*(xi+C0)/R2**3 +
                          (R2+xi+C0)/(R2*(Rq)) +
                          4*(1-nu)*(R2+xi)/(R2*(Rq))
                          )
             )
      
        F1star = (2*z*(costh*q2*(a1*(2*Rq)/(R2**3*(Rq)**2) - b1*(R2 + 2*xi)/(R2*(Rq)**2)) +
                       sinth*(a1/R2**3 -2*b1*(R2 + xi)/(R2* (Rq)))
                       )
                 )
      
        F2 = -2*sinth*z*(xi*(xi+C0)*qbar3/R2**3 + C0/R2 + (5-4*nu)*Abar1)
      
        F2star = 2*z*(a1*ybar3/R2**3 - 2*b1*(sinth*Abar1 + costh*q2*(R2+xi)/(R2*Rq)))
    
    # Calculate little f's
    ff1 = (xi*y1/Ry +
           3/(costh)**2*(y1*lRy*sinth -y1*lRq + 2*q2*atnbeta) +
           2*y1*lRq -
           4*xbar3*atnbeta/costh
          )
    
    ff2 = (xi*y2/Ry +
           3/(costh)**2*(q2*lRq*sinth - q2*lRy + 2*y1*atnbeta*sinth + costh*(R2-ybar3)) -
           2*costh*Abar2 +
           2/costh*(xbar3*lRy - q3*lRq)
          )
    
    ff3 = ((q2*lRq - q2*lRy*sinth + 2*y1*atnbeta)/costh +
            2*sinth*Abar2 + q3*lRy - xi
          )
    
    # Assemble into x, y, z displacements (1,2,3):
    u1 = coeff*(A1star + nu4*Abar1star + F1star)*y1
    
    u2 = coeff*(sinth*(A1star*r2+(nu4*Abar1star+F1star)*q2) +
                costh*(Bstar-F2star) + 2*sinth*costh*z*Abar1star)
    
    u3 = coeff*(-costh*(Abar1star*r2+(nu4*Abar1star-F1star)*q2) +
                sinth*(Bstar+F2star) + 2*(costh)**2*z*Abar1star)
    
    u1 = u1 + 2*coeff*Pdila*((A1 + nu4*Abar1 + F1)*y1 - coeffs[2]*ff1)
    
    u2 = u2 + 2*coeff*Pdila*(sinth*(A1*r2+(nu4*Abar1+F1)*q2) -
                             coeffs[2]*ff2 + 4*nu1*costh*(A2+Abar2) +
                             costh*(A3 - nu4*Abar3 - F2))
    
    u3 = u3 + 2*coeff*Pdila*(costh*(-A1*r2 + (nu4*Abar1 + F1)*q2) + coeffs[2]*ff3 +
                             4*nu1*sinth*(A2+Abar2) +
                             sinth*(A3 + nu4*Abar3 + F2 - 2*nu4*B))
    return u1,u2,u3



def calc_yang(x,y,xoff=0,yoff=0,P=10e6,z0=3e3,theta=90,strike=0,a=1e3,b=5e2,nu=0.25):
    """
    based on Yang et al., vol 93, JGR, 4249-4257, 1988
    Keyword arguments:
    -----------------
    x: x location of point on surface
    y: y location of point on surface
    theta: dip angle of source
    P: pressure change in magma chamber (used inPstar_dila.m) 
    a: semimajor axis of spheroid
    b: semiminor axis of spheriod
    z0: depth of source
    """
    
    mu = 1 #normalized shear modulus
    matrl = np.array([2*nu/(1-2*nu)*mu,
                      mu,
                      nu])
    
    # Store some commonly used coefficients
    coeffs = np.array([1/(16*mu*(1-nu)),
                        3-4*nu,
                        4*(1-nu)*(1-2*nu)])
    
    # Convert to radians & store geometric variables
    #strike  =  -strike* pi/180 + pi/2
    strike =  -strike * np.pi/180
    theta  =  theta * np.pi/180
    coss   =  np.cos(strike)
    sins   =  np.sin(strike)
    
    e_theta = np.array([np.sin(theta),
                        np.cos(theta)])
    rotx = x * coss + y * sins
    roty = -x * sins + y*coss
    
    # Evaluate integrals at xi  =  +c to -c:
    c = np.sqrt(a**2 - b**2)
    if c==0:
        print 'ERROR a must not equal b! use calc_mctigue instead'
    
    # Get spheroid parameters
    phi = 0
    sph = spheroid(a,b,c,matrl,phi,theta,P)
    
    # Run Yang calculations
    xi    = c
    Up1,Up2,Up3 = yang(sph,xi,z0,rotx,roty,matrl,e_theta,coeffs)
    xi = -xi
    Um1,Um2,Um3 = yang(sph,xi,z0,rotx,roty,matrl,e_theta,coeffs)
    
    uxj = -Up1 + Um1
    uyj = -Up2 + Um2
    uz  = Up3 - Um3
    
    ux = -uyj*sins + uxj*coss
    uy = uxj*sins + uyj*coss
    
    
    return ux,uy,uz


def calc_okada():
    """
    Crack dislocation
    TODO
    """
    return ux,uy,uz

    
def calc_fialko():
    """
    Penny-shaped crack
    """
    return ux,uy,uz




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

def cart2los(ux,uy,uz,s1,s2,s3):
    #TODO
    los = 0
    return los

