"""
Reproduce plots in Segall Ch7 or plots from original papers
"""
import numpy as np
import matplotlib.pyplot as plt
import roipy.models as m

plt.rcParams['figure.figsize'] = (11,8.5)
plt.rcParams['font.size'] = 14
plt.style.use('seaborn-white')

# =====================
# Benchmarks
# =====================
def mogi_viscoshell(xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,dP=100e6,mu=30e9,
                    nu=0.25,eta=2e16):
    """
    Mogi surrounded by viscoelastic shell in an elastic half-space
    (Segall figure 7.38)
    """
    # Set-up domain 5km, 100m pixels
    x = np.linspace(0,15e3,1e2)
    y = np.zeros_like(x)

    norm = ((1-nu)*dP*a**3) / (mu*d**2)
    tR = (3*eta*(1-nu)*b**3) / (mu*(1+nu)*a**3)
    times = np.array([0.0,1.0,3.0,10.0]) #Normalized!

    fig = plt.figure()
    ax = fig.add_subplot(211)
    uzmax = np.zeros(5)
    for i,tn in enumerate(times):
        #print(tn)
        t = tn * tR
        ur,uz = m.mogi.calc_viscoshell(x,y,t,xoff,yoff,d,a,b,dP,mu,nu,eta)
        plt.plot(x/d, uz/norm, label=str(tn))
        #print(uz.max())
        uzmax[i+1] = uz.max()
    #for segall reproduction
    plt.xlim(0,3)
    plt.ylim(0,2)
    plt.ylabel('Normalized Uz')
    plt.xlabel('Normalized distance (r/d)')
    plt.legend(title=r'$t/\tau$')
    plt.grid(True)

    ax = fig.add_subplot(212)
    #print(uzmax)
    plt.plot([0,0,1,3,10],uzmax/norm,'b.-')
    #to reproduce segall figure with b/a=1.2
    plt.xlim(-0.25,3)
    plt.ylim(0,2)
    plt.xlabel('Normalized time (t/tR)')
    plt.ylabel('Normalized Uz_max')
    plt.grid(True)
    plt.suptitle('Viscoelastic Shell Benchmark: Segall Figure 7.38')


def mogi_viscoshell_dPt(P0=1.0,xoff=0,yoff=0,d=4e3,a=1000.0,b=1200.0,mu=30e9,
                     nu=0.25,eta=2e16):
    """
    Mogi surrounded by viscoelastic shell in an elastic half-space with an
    exponentially decaying pressure source
    (Segall figure 7.39)
    """
    timesN = np.array([0,0.05,0.2,0.4,0.6,0.8,1,2,3,4,5])
    tR = (3*eta*(1-nu)*b**3) / (mu*(1+nu)*a**3)
    tS = tR * np.array([0.1,1.0,3.0])
    times = tR * timesN

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    norm = ((1-nu)*P0*a**3) / (mu*d**2)
    for ts in tS:
        U = np.zeros_like(times)
        P = np.zeros_like(times)
        for i,t in enumerate(times):
            ur,uz,p = m.mogi.calc_viscoshell_dPt(0,0,t,P0,ts)
            U[i] = uz
            P[i] = p
        ax1.plot(times/tR, P/P0, lw=2, label=str(t/tR))
        ax2.plot(times/tR, U/norm, lw=2)

    ax1.set_ylabel('Normalized Pressure')
    ax1.set_xlabel(r'Normalized time ($t/\tau$)')
    ax1.legend(['0.1','1','3'], title=r'$t_s/\tau$',loc='lower right')
    ax1.set_ylim(0,1.1)
    ax1.grid(True)

    #ax2.set_xlim(0,5)
    ax2.set_ylim(0,2)
    ax2.set_xlabel(r'Normalized time ($t/\tau$)')
    ax2.set_ylabel('Normalized Uz_max')
    ax2.grid(True)
    plt.suptitle('Viscoelastic Shell Benchmark: Segall Figure 7.39')



def mogi_linmax(d=3e3,a=500.0,dP=100e6,mu=4e9,nu=0.25):
    """
    Mogi in linear maxwell viscoelastic halfspace
    (Bonafede 2009 Figure 3)

    NOTE: I think the 'parameters used' are off in the
    """
    # Poisson solid
    K = (5/3.0)*mu

    # Set-up domain 5km, 100m pixels
    x = np.linspace(0,5e3,1e2)
    y = np.linspace(0,5e3,1e2)
    #X,Y = np.meshgrid(x,y)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    tn = np.array([1e-5,0.5,1.0,2.0])
    colors=['b','g','r','k','y']
    urmax=[]
    uzmax=[]
    for color,time in zip(colors,tn):
        #ur,uz = m.mogi.calc_linmax(X,Y,time) #map view
        #ur,uz = m.mogi.forward_dp(x,y,0,0,d,a,dP,mu,nu) #elastic solution
        ur,uz = m.mogi.calc_linmax(x,y,time,0,0,d,a,dP,mu,nu)
        urmax.append(ur.max())
        uzmax.append(uz.max())
        ax.plot(x, uz, c=color, lw=2, label=str(time))
        ax.plot(x, ur, c=color, ls='--',lw=2)

    plt.legend()
    plt.grid(True)
    plt.title('Bonafede et. al. 2009 Fig. 3a')
    plt.xlabel('radial distance (m)')
    plt.ylabel('Uz, Ur (m)')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(tn,uzmax,'b.-',label='Uz_max')
    ax.plot(tn,urmax, 'b.--', label='Ur_max')
    plt.xlabel('t/tau')
    plt.ylabel('displacement (m)')
    plt.title('Bonafede et. al. 2009 Fig. 3b')
    plt.legend()


def mogi_genmax():
    """
    Mogi in general maxwell viscoelastic halfspace
    (Del Negro 2009 Fig 2)
    """
    # Set parameters (p.301)
    params = dict(xoff = 0,
                  yoff = 0,
                  d = 4e3, #m
                  dP = 100e6, #Pa
                  a = 700, #m
                  nu = 0.25,
                  G = 30e9, #Pa
                  eta = 2e16, #Pa*s
                  mu1 = 0.5)

    # Set-up domain 0 x 10km with 100m pixels
    x = np.linspace(0,10e3,1e2)
    y = np.zeros_like(x)

    # time intervals in seconds
    times = np.array([0,1,3,5,7,9,11,13,15]) * 3600 * 24

    # plot solutions in a grid
    fig = plt.figure()
    for i,t in enumerate(times):
        ax = fig.add_subplot(3,3,i+1)
        dr,dz = m.mogi.calc_genmax(x,y,t)
        ax.plot(x/1000, dz,'-', lw=2, label=str(t))
        ax.plot(x/1000, dr,'--', lw=2, label=str(t))
        plt.tick_params(labelbottom=0, labeltop=0, labelleft=0, labelright=0)
        plt.ylim(0,0.1)
        if i==6:
            plt.tick_params(labelbottom=1, labeltop=0, labelleft=1, labelright=0)
            plt.ylabel('Uz (m)')
            plt.xlabel('radial distance (km)')
        plt.grid(True)
        plt.title('time={} days'.format(t/(3600*24)))

    plt.suptitle('Del Negro et. al. 2009 fig 2')
    plt.show()


def mogi():
    """
    Mogi Source in an elastic halfspace
    (Segall Figure 7.5)
    """
    # Set parameters
    params = dict(xoff = 0,
                yoff = 0,
                d = 3e3, #m
                dV = 1e6, #m^3
                nu = 0.25)
    depth = params['d']

    # 10km x 10km with 100m pixels
    x = np.linspace(-15e3,15e3,1e2)
    y = np.linspace(-15e3,15e3,1e2)
    X,Y = np.meshgrid(x,y)

    # Run mogi model with delta volume input
    dr,dz = m.mogi.forward(X,Y,**params)
    # Run mogi model with delta pressure input
    #dr,dz = m.mogi.forward_dp(X,Y,xoff=0,yoff=0,
    #                     depth=3e3,
    #                     dP=10e6,
    #                     a=500,
    #                     nu=0.25,
    #                     mu=4e9,
    #                     output='cyl')

    # Normalize results
    z = dz[50, 50:] / dz.max()
    r = dr[50, 50:] / dz.max()
    x = x[50:] / depth

    # Reproduce the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, z,'b-', lw=3, label='dz')
    ax.plot(x, r,'b--', lw=3, label='dr')
    plt.legend()
    plt.grid(True)
    plt.title('Mogi Displacements')
    plt.xlabel('normalized distance (r/d)')
    plt.ylabel('normalized displacement (dxi / dz.max)')
    plt.show()


def mctigue():
    """
    Spherical Source in an elastic halfspace
    (Segall Figure 7.6B)
    """
    # Set parameters
    params = dict(xoff = 0,
                  yoff = 0,
                  d = 3e3,
                  dP = 10e6,
                  a = 1500.0,
                  nu = 0.25,
                  mu = 4e9)
    depth = params['d']

    # 10km x 10km with 100m pixels
    x = np.linspace(-15e3,15e3,1e2)
    y = np.linspace(-15e3,15e3,1e2)
    X,Y = np.meshgrid(x,y)

    # Run McTigue for first term solution (matches mogi)
    #(x,y,xoff=0,yoff=0,depth=3e3,dP=10e6,a=1500.0,nu=0.25,mu=4e9,terms=1, output='cyl'):
    dr1,dz1 = m.mogi.calc_mctigue(X,Y,terms=1,output='cyl', **params)
    # Two-term solution
    dr2,dz2 = m.mogi.calc_mctigue(X,Y,terms=2,output='cyl', **params)

    # Normalize results
    z1 = dz1[50, 50:] / dz1.max()
    r1 = dr1[50, 50:] / dz1.max()
    z2 = dz2[50, 50:] / dz1.max()
    r2 = dr2[50, 50:] / dz1.max()
    x = x[50:] / depth #normalized distance

    # Reproduce the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, z1,'b-', lw=3, label=r'$order \epsilon^3$')
    ax.plot(x, r1,'b--', lw=3)
    ax.plot(x, z2,'g-', lw=3, label=r'$order \epsilon^6$')
    ax.plot(x, r2,'g--', lw=3)
    plt.legend()
    plt.grid(True)
    plt.title('Mctigue profiles')
    plt.xlabel('normalized distance (r/d)')
    plt.ylabel('normalized displacement (u*(p*d/mu))')
    plt.show()


def okada():
    """
    From script converted from matlab
    """
    # Set parameters
    params = dict(xoff=0, yoff=0,
            depth=5e3, length=1e3, width=1e3,
            slip=0.0, opening=10.0,
            strike=0.0, dip=0.0, rake=0.0,
            nu=0.25)

    # Make grid NOTE: odd number so that center is symmetrical
    n = 201
    x = np.linspace(-25e3,25e3,n)
    y = np.linspace(-25e3,25e3,n)
    X,Y = np.meshgrid(x,y)

    #ux,uy,uz = m.okada.calc_okada(**params)
    ux,uy,uz = m.okada.forward(X,Y,**params)

    # resample grid for quiver plot
    nx = ny = 10

    plt.figure()
    im = plt.imshow(uz, extent=[-25, 25, -25, 25])
    plt.quiver(X[::nx, ::ny]/1e3, Y[::nx, ::ny]/1e3,
               ux[::nx, ::ny], uy[::nx, ::ny],
                units='x', color='w')
    plt.title('Okada profiles')
    plt.xlabel('Easting [km]')
    plt.ylabel('Northing [km]')
    cb = plt.colorbar(im)
    cb.set_label('Vertical Deformation [m]')

    # make sure profile looks OK
    plt.figure()
    mid = int(n/2)
    plt.plot(x/1e3,uz[:,mid], label='vertical')
    plt.plot(x/1e3,uy[:,mid], label='ux')
    plt.plot(x/1e3,ux[mid,:], 'ro', label='uy', markevery=3)
    plt.axhline(color='k')
    plt.axvline(color='k')
    #plt.axhline(0.5, color='k', linestyle='dashed')
    plt.ylabel(' Deformation [m]')
    plt.legend()

    return ux, uy, uz

'''
def okada_fialko():
    """
    Unfortunately Segall Ch3 figures don't list all parameters,
    can instead use dMODELS matlab code (2013) or presentation from T Wright

    U, x, y, nu, delta, d, length, W, fault_type, strike, tp)
    """
    # Set parameters
    params = dict(U = -10.0, # opening [m] NOTE: negative
                  #xoff = -1e3, # offset
                  #yoff = 0, # offset
                  nu = 0.25, # poisson ratio
                  delta = 0.001, # dip angle (set very close to zero to avoid numerical issue)
                  d = 5e3,  # depth to bottom [m]
                  length = 2e3, # [m]
                  W = 2e3, #width [m]
                  fault_type = 3, #opening [bool]
                  strike = 0.0, #[degrees]
                  tp = 0.0 #topo vector [m]
                  )

    # Make grid NOTE: odd number so that center is symmetrical
    n = 201
    x = np.linspace(-25e3,25e3,n)
    y = np.linspace(-25e3,25e3,n)
    X,Y = np.meshgrid(x,y)
    params['x'] = X -1e3
    params['y'] = Y

    ux,uy,uz = m.okada.calc_okada(**params)
    #ux,uy,uz = m.okada.calc_okada_sill(**params)

    # resample grid for quiver plot
    nx = ny = 10

    plt.figure()
    im = plt.imshow(uz, extent=[-25, 25, -25, 25])
    plt.quiver(X[::nx, ::ny]/1e3, Y[::nx, ::ny]/1e3,
               ux[::nx, ::ny], uy[::nx, ::ny],
                units='x', color='w')
    plt.title('Okada profiles')
    plt.xlabel('Easting [km]')
    plt.ylabel('Northing [km]')
    cb = plt.colorbar(im)
    cb.set_label('Vertical Deformation [m]')

    # make sure profile looks OK
    plt.figure()
    mid = int(n/2)
    plt.plot(x/1e3,uz[:,mid], label='vertical')
    plt.plot(x/1e3,uy[:,mid], label='ux')
    plt.plot(x/1e3,ux[mid,:], label='uy')
    plt.axhline(color='k')
    plt.axvline(color='k')
    #plt.axhline(0.5, color='k', linestyle='dashed')
    plt.ylabel(' Deformation [m]')
    plt.legend()

    return ux, uy, uz
'''


def fialko():
    """
    TODO
    """
    print('work in progress')


# ==================
# PLOTTING FUNCTIONS
# ==================

def plot_image(u):
    """ imshow array"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(u)
    cb = plt.colorbar(im)
    cb.set_label('meters')
    plt.title('Synthetic')
    plt.xlabel('distance (m)')
    plt.ylabel('distance (m)')
    plt.show()


def plot_pcolor(x,y,u):
    """ imshow array"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.pcolor(x,y,u)
    cb = plt.colorbar(im)
    cb.set_label('meters')
    plt.title('Synthetic')
    plt.xlabel('distance (m)')
    plt.ylabel('distance (m)')
    plt.show()

def mogi_profiles(dz,dr, x, depth=3e3):
    #full symmetric
    #z = dz[50]
    #r = dr[50]
    #half
    z = dz[50, 50:]
    r = dr[50, 50:]
    x = x[50:]
    #normalized a la segal
    #z = dz[50, 50:] / dz.max()
    #r = dr[50, 50:] / dz.max()
    #x = x[50:] / depth

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, z,'b-', lw=2, label='dz')
    ax.plot(x, r,'g-', lw=2, label='dr')
    plt.legend()
    plt.grid(True)
    plt.title('Mogi profiles')
    #plt.xlabel('distance (r/d)')
    #plt.ylabel('displacement (% dzmax)')
    plt.xlabel('normalized distance (m)')
    plt.ylabel('normalized displacement (m)')
    plt.show()

def mctigue_profiles(dz1,dr1,dz2,dr2, x, depth=3e3):
    #duplicate figure 7.6 in segal ch7
    # NOTE: what is correct normalization?
    z1 = dz1[50, 50:] / dz1.max()
    r1 = dr1[50, 50:] / dz1.max()
    z2 = dz2[50, 50:] / dz1.max()
    r2 = dr2[50, 50:] / dz1.max()
    x = x[50:] / depth #normalized distance

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, z1,'b-', lw=2, label=r'$order \epsilon^3$')
    ax.plot(x, r1,'b-', lw=2)
    ax.plot(x, z2,'g-', lw=2, label=r'$order \epsilon^6$')
    ax.plot(x, r2,'g-', lw=2)
    plt.legend()
    plt.grid(True)
    plt.title('Mctigue profiles')
    #plt.xlabel('distance (r/d)')
    #plt.ylabel('displacement (% dzmax)')
    plt.xlabel('normalized distance (r/d)')
    plt.ylabel('normalized displacement (u*(p*d/mu))')
    plt.show()




if __name__ == '__main__':
    print('Runnning all benchmarks')
    #mogi()
    #mctigue()
    #mogi_linmax()
    #mogi_genmax()
    #mogi_viscoshell()
    ux, uy, uz = okada()
    print('Done!')
