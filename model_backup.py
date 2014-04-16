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