# Inversion tests!

def mogi_invert_insar_justXY( (X,Y), xcen, ycen, depth, dV, nu, look, head):
    """
    Adjust arguments of calc_mogi to work with scipy.omptimize.curvefit convention
    Assumes UTM input for X and Y
    Automatically scales inputs to avoid very large and very small number computations
    """
    # Scale distances and change dtype to avoid numerical problems
    # in this case, lengthscale = 1 km and shift UTM grid to start at 0,0
    X,Y,xcen,ycen = shift_utm(X,Y,xcen,ycen)
    
    dscale = 1e3 #set as kwarg?
    X = X.astype('float64') / dscale
    Y = Y.astype('float64') / dscale
    
    xcen = xcen / dscale
    ycen = ycen / dscale
    depth = depth / dscale
    dV = dV / (dscale**3) 
    
    # Center coordinate grid on point source
    Xshift = X - xcen
    Yshift = Y - ycen
    
    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(Xshift,Yshift) 
    R = np.hypot(depth,rho) 
    
    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3   
    uz = C * depth / R**3  
    ux, uy = pol2cart(th, ur)
    dataVec = np.dstack([ux, uy, uz]) * dscale #shape = (1107, 890, 3)
    
    cart2los = get_cart2los(look,head)
    los = -np.sum(dataVec * cart2los, axis=2)

    return los.ravel() 



def mogi_invert_insar_scale( (X,Y,look,head), xcen, ycen, depth, dV, nu, dscale=1e3):
    """
    Adjust arguments of calc_mogi to work with scipy.omptimize.curvefit convention
    Assumes UTM input for X and Y
    Automatically scales inputs to avoid very large and very small number computations
    """
    # Scale distances and change dtype to avoid numerical problems
    # in this case, lengthscale = 1 km and shift UTM grid to start at 0,0
    X,Y,xcen,ycen = shift_utm(X,Y,xcen,ycen)
    
    #dscale = 1e3 #set as kwarg
    X = X.astype('float64') / dscale
    Y = Y.astype('float64') / dscale
    
    xcen = xcen / dscale
    ycen = ycen / dscale
    depth = depth / dscale
    dV = dV / (dscale**3) 
    
    # Center coordinate grid on point source (previously corner=(0,0))
    Xshift = X - xcen
    Yshift = Y - ycen
    
    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(Xshift,Yshift) 
    R = np.hypot(depth,rho) 
    
    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3   
    uz = C * depth / R**3  
    ux, uy = pol2cart(th, ur)
    dataVec = np.dstack([ux, uy, uz]) * dscale #shape = (1107, 890, 3)
    
    cart2los = get_cart2los(look,head)
    los = -np.sum(dataVec * cart2los, axis=2)

    return los.ravel() 



def mogi_invert_insar_test( (X,Y,look,head),xcen,ycen,depth,dV,nu):
    """
    Adjust arguments of calc_mogi to work with scipy.omptimize.curvefit convention
    Assumes UTM input for X and Y
    """
    # Scale distances and change dtype to avoid numerical problems
    # in this case, lengthscale=1 km
    
    #print np.nanmin(look)
    #print np.nanmin(head)
    # Center coordinate grid on point source
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
    #dataVec = np.array([ux,uy,uz]) #shape = (3, 1107, 890)
    dataVec = np.dstack([ux, uy, uz]) #shape = (1107, 890, 3)
    
    # Get LOS transform
    #look = np.deg2rad(look)
    #head = np.deg2rad(head)
    #EW2los = np.sin(head) * np.sin(look)
    #NS2los = np.cos(head) * np.sin(look) 
    #Z2los = -np.cos(look)
    #cart2los = np.dstack([EW2los, NS2los, Z2los])
    cart2los = get_cart2los(look,head)
    
    # Examine how individual components are effected by viewing geometry
    #ux_los = ux * EW2los
    #uy_los = uy * NS2los
    #uz_los = uz * Z2los
    
    # Project each component into LOS
    # NOTE: follows convention of uplift positive, subsidece negative
    los = -np.sum(dataVec * cart2los, axis=2)

    #return los
    return los.ravel() #flattened arrays required by scipy.optimize


# =====================
# MISC
# ===================
def calc_visco_shell_uzmax(x,y,tn,xoff=0,yoff=0,depth=3e3,dP=100e6,a=500,mu=4e9,nu=0.25,
                    output='cyl'):
    """ Max vertical uplift from Bonafede 2009 eq#29 , based on simplification of
    Dragoni & Magnanensi 1989 *mu1=mu2=mu, K1=K2=(5/3)mu --> poisson solid
    
    Required arguments:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)
    tn: normalized time (t/tau_gamma)
    
    Keyword arguments:
    -----------------
    xoff: y-offset of point source epicenter (m)
    yoff: y-offset of point source epicenter (m)
    depth: depth to point (m)
    R2: outer radius (m)
    R1: inner radius (m)
    dP: change in pressure at t=0 (Pa)
    K1: short-term elastic incompressibility () NOTE=(5/3)mu in poisson approximation
    mu1: short-term elastic rigidity (Pa)
    eta: effectice viscosity (Pa*s)
    mu2: half-space shear modulus (Pa)
    K2: half-space bulk modulus
    output: 'cart' (cartesian), 'cyl' (cylindrical)
    
    """
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff
    d = depth
    
    # convert to surface cylindrical coordinates
    th, r = cart2pol(x,y)
    R = np.hypot(d,r)
    
    # max vertical displacement (normalized time)
    #tau = eta / mu
    #gamma = (9.0/5)*(R2/R1)**3
    C = (3*dP*R2**3) / (4*mu*d**2)
    term = 1 - (1-(R1/R2)**3)*np.exp(-tn)
    uzmax = C * term    
    
    return uzmax


def calc_visco_vdot(x,y,t,xoff=0,yoff=0,depth=3e3,vdot=1.94,a=500,mu=4e9,nu=0.25,
                    output='cyl'):
    """ Max vertical uplift from Bonafede 2009 eq#, based on constant volumetric
    magma supply rate vdot (m^3/s). time t given in seconds(s)
    Vo = initial volume
    """
    # center coordinate grid on point source
    x = x - xoff
    y = y - yoff
    d = depth
    
    # convert to surface cylindrical coordinates
    th, r = cart2pol(x,y)
    R = np.hypot(d,r)
    
    # max vertical displacement (normalized time)
    V0 = (4/3)*np.pi*r**3
    tau = eta / mu
    beta = (3*K + 4*mu) / (3*K)   
    C1 = (4*vdot*mu*tau) / (3*V0)
    dp = C1 * (1 - np.exp(-t/(tau*beta)))
    C2 = (vdot*tau) / (2*np.pi*d**2)
    uzmax = C2 * ( (t/tau) - (mu/(3*K))*(1 - exp(-t/(tau*alpha))))
    
    return dp, uzmax



# NOTE: Seems to be a small error with these yang routines,,, not sure where:
# So copying routines that do work to roipy

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
    atnbeta     = np.pi/2 * np.sign(betatop)
    atnbeta[nz] = np.arctan(betatop[nz] / betabottom[nz])
    
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



def calc_yang(x,y,xoff,yoff,d,P,a,b,strike,dip,nu):
    """
    pressure source based on Yang et al., vol 93, JGR, 4249-4257, 1988
    
    Keyword arguments:
    -----------------
    xoff:   x location of point on surface [km]
    yoff:   y location of point on surface [km]
    d:      depth to source center [km]
    P:      pressure change in magma chamber [mu*10^(-5) Pa]    
    a:      semimajor axis of spheroid [km]
    b:      semiminor axis of spheriod [km]
    strike: strike angle of source [degrees] / phi
    dip:    dip angle of source [degrees] / theta
    
    Output:
    ----------------
    ux,uy,uz
    """
    xn = x - xoff
    yn = y - yoff
    
    mu = 1 #normalized shear modulus
    #1st Lame Paramter, 2nd Lame parameter (Shear Modulus), 
    matrl = np.array([2*nu/(1-2*nu)*mu, mu, nu])
    
    # Store some commonly used coefficients
    coeffs = np.array([1/(16*mu*(1-nu)),
                        3-4*nu,
                        4*(1-nu)*(1-2*nu)])
    
    # Convert to radians & store geometric variables
    strike =  np.deg2rad(strike)
    dip  =  np.deg2rad(dip)
    coss   =  np.cos(strike)
    sins   =  np.sin(strike)
    e_theta = np.array([np.sin(dip),
                        np.cos(dip)])
    
    rotx = xn * coss + yn * sins
    roty = yn * coss - xn * sins
    
    c = np.sqrt(a**2 - b**2)
    if c==0:
        print 'ERROR: a must not equal b! use calc_mctigue instead'
        return
    
    # Get spheroid parameters (a,b,c,matrl,phi,theta,P)
    sph = spheroid(a,b,c,matrl,strike,dip,P)
    
    # Run Yang calculations (x,y,z0,sph,xi,matrl,e_theta,coeffs)
    xi    = c
    Up1,Up2,Up3 = yang(rotx,roty,d,sph,xi,matrl,e_theta,coeffs)
    xi = -xi
    Um1,Um2,Um3 = yang(rotx,roty,d,sph,xi,matrl,e_theta,coeffs)
    
    uxj = -Up1 + Um1
    uyj = -Up2 + Um2
    
    ux = -uyj*sins + uxj*coss
    uy = uxj*sins + uyj*coss
    uz  = Up3 - Um3
    
    
    return ux,uy,uz

def yang_benchmark():
    '''
    Reproduce demo from matlab codes on Yuri Fialko's website
    # http://sioviz.ucsd.edu/~fialko/software.html
    
    NOTE: input units in km and normalized pressure.
    '''
    #x,y,  xoff,yoff,d,P,a,b,strike,dip,nu
    params=[20.0, 30.0, 15.0, 10.0, 12.0, 4.0, 30.0, 40.0, 0.25]
    
    # Intitialize computation grid
    x0=0
    y0=0
    dx=0.5
    dy=0.5
    ndat=100
    mdat=100
    xx = np.arange(x0, ndat*dx, dx)
    yy = np.arange(y0, mdat*dy, dy)
    x,y = np.meshgrid(xx,yy)
    
    ux,uy,uz = m.calc_yang(x,y,*params)
    
    # Subsample for quiver plot
    dhx=np.round(ndat / 10)
    dhy=np.round(mdat / 10)
    xsub=x[:-1:dhx,:-1:dhy]
    ysub=y[:-1:dhx,:-1:dhy]
    uxsub = ux[0:-1:dhx,0:-1:dhy]     
    uysub = uy[0:-1:dhx,0:-1:dhy]

    plt.figure()
    plt.imshow(uz, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.quiver(xsub,ysub,uxsub,uysub)
    plt.title('Yang Uz & Ur Surface Deformation')
    plt.xlabel('EW Distance [km]')
    plt.ylabel('NS Distance [km')
    
    plt.show()
    
    
    





