# Set Okada Parameters here to run in model 
xcen = 0            # [m] NS offset of source from center of grid   
ycen = 0            # [m] EW offset of source from center of grid
U = 1.0             # [m] U is slip
d = 1e-3            # [m] depth (positive down)
nu = 0.27           # [unitless] Poisson ratio
delta = 89.99        # [degrees] delta is dip angle, 90.0 exactly might cause numerical issues?
strike = 90.0        # [degrees] counter clockwise from north
length = 70e3        # [m] # len,W are the fault length and width, resp.
width = 30e3        # [m]
fault_type = 1        # fault_type is 1 2 3 for strike, dip, and opening