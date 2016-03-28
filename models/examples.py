# Set Okada Parameters here to run in model 
xcen = 0            # [m] NS offset of source from center of grid   
def okada(example='strike slip'):
	'''
	examples = ['strike slip', 'thrust', normal', 'finite sill', 'point sill', 'dyke']
	
	ycen = # [m] EW offset of source from center of grid
	U = # [m] U is slip
	d = # [m] depth (positive down)
	nu = # [unitless] Poisson ratio
	delta = # [degrees] delta is dip angle, 90.0 exactly might cause numerical issues?
	strike = # [degrees] counter clockwise from north
	length = # [m] # len,W are the fault length and width, resp.
	width = # [m]
	fault_type = # fault_type is 1 2 3 for strike, dip, and opening
	
	Usage:
	params = rp.models.examples.okada('point sill')
	
	'''
	if example == 'strike slip':
		print('1 m of left-lateral slip on NS-striking vertical fault (70km down to 14km) that ruptures surface')
		xcen=0
		ycen=0
		U = 1.0 
		d = 1e-3 
		nu = 0.27 
		delta = 90.0 
		strike = 90.0 
		length = 70e3 
		width = 15e3
		fault_type = 1 
		
	if example == 'thrust':
		print('1m of slip, 70x30km fault, 30-dipping to W with top at 1km below surface')
		xcen = 0           
		ycen = 0            
		U = 1.0            
		d = 1e3            
		nu = 0.25           
		delta = 30.0      
		strike = 180.0        
		length = 70e3        
		width = 30e3       
		fault_type = 2 
		
	
	if example == 'finite sill':
		print('1m of opening, 50x50km sill at 10km depth') 
		xcen = 0             
		ycen = 0           
		U = -1.0      
		d = 10e3          
		nu = 0.25        
		delta = 0.0       
		strike = 0.0        
		length = 50e3        
		width = 50e3        
		fault_type = 3       
	
	
	if example == 'point sill':
		print('Sill dimensions (3x3km) much less than depth (20km)')
		xcen = 0             
		ycen = 0           
		U = -1.0      
		d = 20e3          
		nu = 0.25        
		delta = 0.0       
		strike = 0.0        
		length = 3e3        
		width = 3e3        
		fault_type = 3       
	
	
	if example == 'normal fault':
		print('1m of slip, 70x30km fault, 60-dipping to E with top at 1km below surface')
		xcen = 0           
		ycen = 0            
		U = -1.0            
		d = 1e3            
		nu = 0.25           
		delta = 60.0      
		strike = 0.0        
		length = 70e3        
		width = 30e3       
		fault_type = 2 
	
	
	if example == 'dyke':
		print('1m of opening on vertical NS-striking dyke w/ top at 5km depth') 
		xcen = 0           
		ycen = 0         
		U = -1.0       
		d = 5e3           
		nu = 0.25         
		delta = 90.0     
		strike = 0.0        
		length = 70e3      
		width = 30e3        
		fault_type = 3        
		
	
	return xcen,ycen,U,d,nu,delta,strike,length,width,fault_type
