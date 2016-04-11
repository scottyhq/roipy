"""
Useful plots to use in ipython sessions
"""

from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as pltdate
from matplotlib.widgets import Slider
from matplotlib.offsetbox import AnchoredText

try:
	from mpl_toolkits.basemap import Basemap, cm
except:
	print('Matplotlb Basemap required for map plots')

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle, Polygon
from matplotlib.colors import Normalize
#for lasso selection - note moved to stand-alone plotting utility
#from matplotlib.widgets import Lasso
#from matplotlib.nxutils import points_inside_poly #deprecated w/ matplotlib 1.3.0
#from matplotlib.colors import colorConverter
#from matplotlib.collections import RegularPolyCollection

#try:
from mpl_toolkits.axes_grid1 import ImageGrid
#except:
#    from mpl_toolkits.axes_grid import ImageGrid #NOTE: for old matplotlib
from mpl_toolkits.axes_grid.inset_locator import inset_axes

#from osgeo import gdal, ogr
import gdal
import ogr
import os
import numpy as np
import numpy.ma as ma
import scipy.stats #nice linear regression command
import scipy.interpolate
import sys

import datetime
from roipy import tools #NOTE: should keep this separate? plot can't import tools
#import roipy.tools
import roipy.data



def pig(Interferogram, array=2, title='', units='radar', cunits='phase', ax=None):
    """
    bring up a pseudocolor image of any array associated with an interferogram
    possible= phs,amp,look,head,cor. Array=nparray,1,2"
    units = 'radar','ground'
    cunits = 'phase', 'cm'
    """
    self = Interferogram
    
    if cunits == 'cm':
        convert = True
    else:
        convert = False
        
    data = tools.load_half(self,array, convert2cm=convert)
    
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    im = ax.imshow(data,cmap=plt.cm.jet)
    cb = plt.colorbar(im)
    
    if units == 'km':
        maxx = self.Width * self.col2km
        maxy = self.Length * self.row2km
        im.set_extent((0,maxx,0,maxy))
        # maintain pixel aspect ratio
        #ax.set_aspect(self.Length/self.Width)
        #ax.autoscale_view(scalex=False, scaley=False)
        
    plt.xlabel('range ({})'.format(units))
    plt.ylabel('azimuth ({})'.format(units))
    plt.title(title)
    cb.set_label(cunits)
    plt.show()
    
    return cb


def pcolor(array, title='', ig=None, units='pixels', cmap='jet', ax=None):
	""" display a colormap of an interferogram. If units='km', use RANGE_PIXEL_SIZE
	and AZIMUTH_PIXEL_SIZE in the .rsc to scale the display"""
	#print title, units, ax
	if ax==None:
		fig = plt.figure()
	
	#http://matplotlib.org/examples/color/colormaps_reference.html
	# for insar, best use RdBu_r, bwr, seismic (each differ just in terms of extreme saturation val 
	im = plt.imshow(array, cmap=plt.get_cmap(cmap)) 
	cb = plt.colorbar()
	if title:
		plt.title(title)
	if units == 'km':
		maxx = ig.Width * ig.col2km
		maxy = ig.Length * ig.row2km
		im.set_extent((0,maxx,0,maxy))
	
	plt.xlabel('range ({})'.format(units))
	plt.ylabel('azimuth ({})'.format(units))
	plt.title(title)
	#cb.set_label('dLOS (cm)') #not necessarily... can also plot dems etc
	plt.show()  



def pcolor_ma(dataFile, maskFile=None):
    """load separate data and mask and plot"""
    data = np.load(dataFile)
    mask = np.load(maskFile)             
    maskArray = ma.array(data, mask=mask)
    fig = plt.figure()
    im = plt.imshow(maskArray,origin='lower',cmap=plt.cm.jet)
    cb = plt.colorbar()
    

def plot_components(ux,uy,uz,ulos,n=10,cmode='each', cloc='right'):
    """ 'n' is number of pixels per quiver arrow """
    fig = plt.figure(figsize=(11,8))
    fig.suptitle('Components of Deformation', fontsize=14, fontweight='bold')
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (2, 2),
                direction="row",
                axes_pad = 0.5,
                add_all=True,
                label_mode = '1', #'all', 'L', '1'
                share_all = True,
                cbar_location=cloc, #top,right
                cbar_mode=cmode, #each,single,None
                cbar_size=0.1,#"7%",
                cbar_pad=0.0#,"1%",
                )
    # needed for 'quiver' plot
    Z = np.zeros_like(ux)
    R, C = np.indices(ux.shape)
    #ur = np.hypot(ux,uy)
    #uy = fix_uy(ux,uy,ur,uz)
    #datas = [uz,ur,ux,uy]
    
    datas = [ulos,uz,ux,uy]
    titles = ['LOS', 'Vertical', 'East-West', 'North-South' ]
    master = Normalize(-1,1)
    #for ax,data,title,clim in zip(grid,datas,titles,clims):
    for ax,data,title in zip(grid,datas,titles):   
        #ax = grid[i]
        im = ax.imshow(data, cmap=plt.cm.jet)
        #im = ax.imshow(data, cmap=plt.cm.jet, vmin=clim[0],vmax=clim[1])
        ax.grid(True)
        if cmode == 'each':
            ax.cax.colorbar(im)
        else:
            im.set_norm(master)
        #cmin = data[np.isfinite(data)].min()
        #cmax = data[np.isfinite(data)].max()
        #print cmin, cmax
        #ax.cax.set_xticks([int(cmin),0,int(cmax)])
        
        #Annotate
        #title = titles[i]
        at = AnchoredText(title,prop=dict(size=10,weight='bold'), frameon=True, loc=9) #upper center
        at.patch.set_boxstyle("round, pad=0.0, rounding_size=0.2")
        ax.add_artist(at)
        
        
    if cmode == 'single':
        ax.cax.colorbar(im)
    
    # Plot arrows
    #ax = grid[1] 
    #ax.quiver(C[::n,::n], R[::n,::n], ux[::n,::n], uy[::n,::n], pivot='mid') #radial
    ax = grid[2]
    ax.quiver(C[::n,::n], R[::n,::n], ux[::n,::n], Z[::n,::n], pivot='mid') #EW
    ax = grid[3]
    ax.quiver(C[::n,::n], R[::n,::n], Z[::n,::n], uy[::n,::n], pivot='mid') #NS
    plt.show()



def clim(fig=1, subplot=0, vmin=-1, vmax=1):
    """change colobar limits of currently active figure"""
    #bewhere colorbars with get_axes method
    ax = plt.figure(fig).axes[subplot]
    plt.sca(ax)
    plt.sci(ax.images[0])
    plt.clim((vmin,vmax))
    
 
def variogram(Interferogram, array=None):
    """ Plot the variogram for an array """
    self = Interferogram
    if array == None:
        phs = rp.tools.load_half(2)
 
 
def correlate(pixX=None, pixY=None, rows=None, fit_order=None, ax=None,
              one2one=False, showeq=False, showr2=False, ms=0.5, markevery=1,
              minX=None):
    """Correlate values of two coregistered 2D arrays.
    #NOTE: minX is to avoid spurious lower elevation points (lakes & errors)
    Example
    -------
    Make plot of dem height vs unwrapped phase, if binning specify pixel box
    size. eg =[700,1000] for the salar de atacama subregion. order=n where n is
    the order of best-fit polynomial
    >>>InSAR.Plot.correlate(phs,hgt,rows=[700,1000],order=1)"""
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if rows: # Only plot subset of rows
        pixX = pixX[rows[0]:rows[1],:]
        pixY = pixY[rows[0]:rows[1],:]

    vecX = pixX.flatten()
    vecY = pixY.flatten()
    masked = np.isnan(vecX) + np.isnan(vecY) # Mask NaN values
    
    if minX:
        masked += (vecX < minX)
    
    vecX = vecX[~masked]
    vecY = vecY[~masked]
    
    line_points = ax.plot(vecX,vecY,'k.',markersize=ms, markevery=markevery)
    #ax.set_title('DEM vs Phase Plot', fontsize=14, fontweight='bold')
    #ax.set_xlabel('Phase')
    #ax.set_ylabel('DEM height (m)')
    #plt.hold(True)
    
    # Plot 1:1 line
    if one2one:
        xmin,xmax = ax.get_xlim()
        points = np.linspace(xmin,xmax,20)
        one2one = ax.plot(points,points,'b-',label='1:1') #NOTE:add legend
    
    # Plot best-fit line
    if fit_order: 
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(vecX,vecY)
        coef = (slope,intercept)
        rsq = r_value**2
        #coef = np.polyfit(vecPHS,vecHGT,order)
        fit = np.poly1d(coef)
        xmin,xmax = ax.get_xlim()
        points = np.linspace(xmin,xmax,20)
        line_fit = ax.plot(points,fit(points),'-')
        if showeq: # Display linear fit equation
            text = 'y = %.4fx + %.4f' % (coef[0],coef[1])
            ax.text(0.2, 0.90, text, 
                horizontalalignment='left',
                verticalalignment='center',
                transform = ax.transAxes) 
        if showr2:
            text = r'$R^2 = {:.2f} $'.format(rsq)
            ax.text(0.2, 0.95, text, 
                horizontalalignment='left',
                verticalalignment='center',
                transform = ax.transAxes) 
            
    plt.draw()

    
def geo_correlate(arrayX, arrayY, order=1):
    """DEPRECATED:  just use correlate()Pass cropped arrays on same grid from gmtsample() output. ideally arrays
    should have already been converted from LOS to Vertical"""
    if type(arrayX) is str:
        ds = gdal.Open(arrayX)
        nLon = ds.RasterXSize
        nLat = ds.RasterYSize
        band = ds.GetRasterBand(1)
        arrayX = band.ReadAsArray(0,0,nLon,nLat)
        arrayX[arrayX==0] = np.NaN # Maybe not needed...

    if type(arrayY) is str:
        ds = gdal.Open(arrayY)
        nLon = ds.RasterXSize
        nLat = ds.RasterYSize
        band = ds.GetRasterBand(1)
        arrayY = band.ReadAsArray(0,0,nLon,nLat)
        arrayY[arrayY==0] = np.NaN

    vecX = arrayX.flatten() #necessary?
    vecY = arrayY.flatten()
    masked = np.isnan(vecX) + np.isnan(vecY) # Mask NaN values
    vecX = vecX[~masked]
    vecY = vecY[~masked]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line_points = ax.plot(vecX,vecY,'m.',markersize=0.4,markevery=2,label='pixel value')
    #ax.set_title('Overlapping Pixel Correlation', fontsize=14, fontweight='bold')
    #ax.set_xlabel('t282 stack velocity (cm/yr)')
    #ax.set_ylabel('t10 stack velocity (cm/yr)')
    plt.hold(True)

    # Plot 1:1 line
    xmin,xmax = ax.get_xlim()
    points = np.linspace(xmin,xmax,20)
    one2one = ax.plot(points,points,'b-',label='one2one')
    
    # Plot best-fit line
    '''
    if order: 
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(vecX,vecY)
        fit_sp = np.poly1d(slope,intercept)
        line_fit = ax.plot(points,fit_sp(points),'b--',label='scipy fit')
        ax.text(0.95, 0.9, str(r_value), 
                horizontalalignment='right',
                verticalalignment='center',
                transform = ax.transAxes)
    '''
    plt.legend()
    plt.draw()


def side_by_side_old(array1, array2, array3, fixCaxis=True):
	""" DEPRECATED Plot two arrays with imshow, keep the same colorbar for each, third
	array has it's own colorbar scale (eg. for residual)"""
	origin = 'upper'
	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	ax1.set_title('Rmp')
	im1 = ax1.imshow(array1,origin=origin,cmap=plt.cm.jet)
	fig.colorbar(im1)
	
	ax2 = fig.add_subplot(132)
	ax2.set_title('Synthetic')
	im2 = ax2.imshow(array2,origin=origin,cmap=plt.cm.jet)
	
	if fixCaxis:
		# Set the same colorscale bounds as im1
		norm = mpl.colors.Normalize(vmin=np.nanmax(array1), vmax=np.nanmax(array1))
	im2.set_norm(norm) 
	fig.colorbar(im2)
	
	ax3 = fig.add_subplot(133)
	ax3.set_title('Residual')
	im3 = ax3.imshow(array3,origin=origin,cmap=plt.cm.jet)
	fig.colorbar(im3)
	
	plt.draw()


def side_by_side(array1, array2, array3):
    """ Use ImageGrid for aligned image plots """
    fig = plt.figure()
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (1, 3),
                    direction="row",
                    axes_pad = 0.05,
                    add_all=True,
                    label_mode = "1",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="all", #"single"
                    cbar_size="10%",
                    cbar_pad=0.05,
                    )
    #grid[0].set_xlabel("X")
    #grid[0].set_ylabel("Y")
    #grid2[0].set_xticks([-2, 0])
    #grid2[0].set_yticks([-2, 0, 2])

    #NOTE: could find global min/max from three arrays here
    norm = mpl.colors.Normalize(vmin=np.nanmin(array1),
                                vmax=np.nanmax(array1))
    for ax,data in zip(grid,[array1,array2,array3]):
        im = ax.imshow(data,origin='lower',norm=norm,cmap=plt.cm.jet)

    cax = grid.cbar_axes[0]
    cax.colorbar(im)
    cax.toggle_label(True)
    # Add internal titles to arrays if desired
    #for ax, im_title in zip(grid2, ["(a)", "(b)", "(c)"]):
    #    t = add_inner_title(ax, im_title, loc=2)
    #    t.patch.set_ec("none")
    #    t.patch.set_alpha(0.5)
    plt.show()


def showmap(path, half=2,dims=None, annotate=True, ax=None, vmin=None, vmax=None, nGrid=1.0):
    """Use Basemap module to plot a georeferenced interferogram, return basemap instance"""
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    if path.endswith('unw'):
        geo = roipy.data.Geogram(path)
        geotrans = geo.Geotrans
        data = tools.load_half(geo,half=half)
    else:
        data, geotrans, proj = tools.load_gdal(path)
    
    nLat, nLon = data.shape
    x0 = geotrans[0] #top left longitude
    dx = geotrans[1] #pixel width
    xR = geotrans[2]
    y0 = geotrans[3] #top left latitude
    yR = geotrans[4]
    dy = geotrans[5] #pixel height
    lowLeft = (x0, y0 + dy*nLat) #NOTE: make sure dy is correct sign
    upRight = (x0 + dx*nLon, y0)
    
    # NOTE: 'cyl' projection implies (lon,lat) = (x,y), see Basemap docs otherwise
    bmap = Basemap(projection='cyl',
                llcrnrlat=lowLeft[1],
                urcrnrlat=upRight[1],
                llcrnrlon=lowLeft[0],
                urcrnrlon=upRight[0],
                resolution='i',
                suppress_ticks=False,
                ) #res=c,l,i,h,f
    #lons = xO + dx * np.arange(nLon)
    #lats = y0 + dy * np.arange(nLat)
    im = bmap.imshow(np.flipud(data),vmin=vmin,vmax=vmax,ax=ax) #NOTE: have to flip it up and down

    if annotate:
        bmap.drawcountries()
        #latVec = np.arange(bmap.llcrnrlat, bmap.urcrnrlat, nGrid)
        #lonVec = np.arange(bmap.llcrnrlon, bmap.urcrnrlon, nGrid)
        #print latVec
        #print lonVec
        #bmap.drawparallels(latVec,color='black',linewidth=1,dashes=[2,2],labels=[1,0,0,1],labelstyle='+/-',fmt='%.1f')
        #bmap.drawmeridians(lonVec,color='black',linewidth=1,dashes=[2,2],labels=[1,0,0,1],labelstyle='+/-',fmt='%.1f')
        
        #Add colorbar
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05) #"top"
        #plt.colorbar(im, cax=cax)
        bmap.colorbar(im)
        
    plt.show()
    return bmap
    
    
# =============================================================================#
# Requires ts instance
# =============================================================================#
def browseSet(Set, convert2cm=False, units=None):
    """ Simple GUI to browse through all interferograms in a roipy.data.Set"""
    self = Set
    pair = str(self.Pairs[0][0]) + ' ' + str(self.Pairs[0][1])
    ig = self.Igrams[pair]
    #ig = self.Igrams.values()[0]
    #dates = ig.Rsc['PAIR']
    name = ig.Name
    data = tools.load_half(ig,2, convert2cm=convert2cm)
    
    # Make plot window with slider to scan through IGs
    fig = plt.figure()
    fig.suptitle(os.path.dirname(ig.Path), fontsize=14, fontweight='bold')
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.set_title(pair)
    im = plt.imshow(data,cmap=plt.cm.jet)
    cb = plt.colorbar()
    
    if units:
        cb.set_label(units)

    # Platform-independent Slider to go through IGs
    axcolor = 'lightgoldenrodyellow'
    axIG = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    sIG = Slider(axIG,
                 'IG#',
                 0,
                 len(self.Pairs)-1,
                 #valinit=1, set below
                 valfmt='%i',
                 closedmin=True,
                 closedmax=True,
                 dragging=True)


    def update(val):
        val = int(val)
        pair = str(self.Pairs[val][0]) + ' ' + str(self.Pairs[val][1])
        ig = self.Igrams[pair]
        #ig = self.Igrams.values()[sIG.val-1]
        #dates = ig.Rsc['PAIR']
        data = tools.load_half(ig,2, convert2cm=convert2cm)
        im.set_data(data)
        im.autoscale() # autoscale colorbar to new data
        ax.set_title(pair)
        plt.draw()

    def onpress(event):
        if event.key not in ('n', 'p'): return
        if event.key=='n':
            newval = sIG.val + 1
        else:
            newval = sIG.val - 1
        if newval < sIG.valmin: newval = sIG.valmax
        if newval > sIG.valmax: newval = sIG.valmin
        sIG.set_val(newval) # update() automatically called

    fig.canvas.mpl_connect('key_press_event', onpress)
    sIG.on_changed(update)
    sIG.set_val(0)


def browse(Timeseries):
    """Browse through 3 products simultaneously """
    self = Timeseries.Set
    pair = str(self.Pairs[0][0]) + ' ' + str(self.Pairs[0][1])
    ig = self.Igrams[pair]
    #ig = self.Igrams.values()[0]
    #dates = ig.Rsc['PAIR']
    #name = ig.Name
    data = tools.load_half(ig,2)
    
    # Make plot window with slider to scan through IGs
    fig = plt.figure()
    fig.suptitle(os.path.dirname(ig.Path), fontsize=14, fontweight='bold')
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.set_title(pair)
    im = plt.imshow(data,cmap=plt.cm.jet)
    cb = plt.colorbar()

    # Platform-independent Slider to go through IGs
    axcolor = 'lightgoldenrodyellow'
    axIG = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    sIG = Slider(axIG,
                 'IG#',
                 0,
                 len(self.Pairs)-1,
                 #valinit=1, set below
                 valfmt='%i',
                 closedmin=True,
                 closedmax=True,
                 dragging=True)


    def update(val):
        val = int(val)
        pair = str(self.Pairs[val][0]) + ' ' + str(self.Pairs[val][1])
        ig = self.Igrams[pair]
        #ig = self.Igrams.values()[sIG.val-1]
        #dates = ig.Rsc['PAIR']
        data = tools.load_half(ig,2)
        im.set_data(data)
        im.autoscale() # autoscale colorbar to new data
        ax.set_title(pair)
        plt.draw()

    def onpress(event):
        if event.key not in ('n', 'p'): return
        if event.key=='n':
            newval = sIG.val + 1
        else:
            newval = sIG.val - 1
        if newval < sIG.valmin: newval = sIG.valmax
        if newval > sIG.valmax: newval = sIG.valmin
        sIG.set_val(newval) # update() automatically called

    fig.canvas.mpl_connect('key_press_event', onpress)
    sIG.on_changed(update)
    sIG.set_val(0)
   


def browse3_new(Timeseries,prefix1='d_',prefix2='ramp_',prefix3='rd_'):
    """ Plot two arrays with imshow, keep the same colorbar for each, third
    array has it's own colorbar scale (eg. for residual)"""
    print('use arrow keys to toggle back and forth')
    # Make plot window with slider to scan through IGs
    #print 'changed'
    self = Timeseries.Set
    dates = self.PairsString[0]
    ig = self.Igrams[dates]
    
    #data = tools.load_binary_old(Timeseries, ig, prefix=prefix1)
    #synthetic = tools.load_binary_old(Timeseries,ig, prefix=prefix2)
    #residual = tools.load_binary_old(Timeseries,ig, prefix=prefix3)
    data = np.load(Timeseries.RunDir +'/'+ prefix1 + ig.Name[:-3] + 'npy')
    synthetic = np.load(Timeseries.RunDir +'/'+ prefix2 + ig.Name[:-3] + 'npy')
    residual = np.load(Timeseries.RunDir + '/'+ prefix3 + ig.Name[:-3] + 'npy')
    
    cmap = 'viridis'
    
    fig = plt.figure(figsize=(11,8.5))
    bigtitle = fig.suptitle(dates, fontsize=14, fontweight='bold')
    ax1 = fig.add_subplot(131)
    #ax1.subplots_adjust(left=0.25, bottom=0.25)
    ax1.set_title(prefix1)
    im1 = plt.imshow(data, cmap=cmap)
    cb1 = fig.colorbar(im1)

    ax2 = fig.add_subplot(132)
    ax2.set_title(prefix2)
    im2 = ax2.imshow(synthetic,cmap=cmap)
    # Set the same colorscale bounds as im1
    #norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
    #                            vmax=data[np.isfinite(data)].max())
    #im2.set_norm(cb1.norm) #alternative to manually accessing norm as above
    fig.colorbar(im2)

    ax3 = fig.add_subplot(133)
    ax3.set_title(prefix3)
    im3 = ax3.imshow(residual,cmap=cmap)
    im3.set_norm(cb1.norm)
    fig.colorbar(im3)

    # Platform-independent Slider to go through IGs
    axcolor = 'lightgoldenrodyellow' # other options?
    axIG = plt.axes([0.25, 0.015, 0.65, 0.03], axisbg=axcolor)
    sIG = Slider(axIG, 'IG#', 0, len(self.Pairs), valinit=0, valfmt='%i',
                 closedmin=True, closedmax=True, dragging=True)
    #print dir(sIG)
    def update(val):
        val = int(val)
        dates = str(self.Pairs[val][0]) + ' ' + str(self.Pairs[val][1]) #just take pairsstring?
        ig = self.Igrams[dates]
        #ig = self.Igrams.values()[sIG.val-1]
        #dates = ig.Rsc['PAIR']
        #data = tools.load_half(ig,2)
        
        #dates = self.PairsString[val]
        #print dates
        #ig = self.Igrams[dates]
        #print ig

        data = np.load(Timeseries.RunDir +'/'+ prefix1 + ig.Name[:-3] + 'npy')
        synthetic = np.load(Timeseries.RunDir +'/'+ prefix2 + ig.Name[:-3] + 'npy')
        residual = np.load(Timeseries.RunDir + '/'+ prefix3 + ig.Name[:-3] + 'npy')
        
        im1.set_data(data)
        im2.set_data(synthetic)
        im3.set_data(residual)
        
        norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
                                    vmax=data[np.isfinite(data)].max())
        im1.autoscale() # autoscale colorbar to new data
        #im2.set_norm(norm) #match synthetic colorbar with Rmp
        im2.autoscale()
        im3.autoscale()
        #im3.set_norm(norm)
        bigtitle.set_text(dates)
        #ax2.set_title(self.Pairs[int(sIG.val)-1])
        plt.draw()
        
    def onpress(event):
        if event.key not in ('right', 'left'): return
        if event.key=='right': newval = sIG.val + 1
        else:  newval = sIG.val - 1
        
        if newval < sIG.valmin: newval = sIG.valmax
        if newval > sIG.valmax: newval = sIG.valmin
        sIG.set_val(newval) # update() automatically called
    
    fig.canvas.mpl_connect('key_press_event', onpress)
    sIG.on_changed(update)
    sIG.set_val(0)
    #plt.show()
    
  

  


def browse3_old(Timeseries,prefix1='Def',prefix2='MskDef',prefix3='RmpMskDef'):
    """ Plot two arrays with imshow, keep the same colorbar for each, third
    array has it's own colorbar scale (eg. for residual)"""
    # Make plot window with slider to scan through IGs
    #print 'changed'
    self = Timeseries.Set
    dates = self.PairsString[0]
    ig = self.Igrams[dates]
    data = tools.load_binary_old(Timeseries, ig, prefix=prefix1)
    synthetic = tools.load_binary_old(Timeseries,ig, prefix=prefix2)
    residual = tools.load_binary_old(Timeseries,ig, prefix=prefix3)
    
    fig = plt.figure(figsize=(11,8))
    bigtitle = fig.suptitle(dates, fontsize=14, fontweight='bold')
    ax1 = fig.add_subplot(131)
    #ax1.subplots_adjust(left=0.25, bottom=0.25)
    ax1.set_title(prefix1)
    im1 = plt.imshow(data, cmap=plt.cm.jet)
    cb1 = fig.colorbar(im1)

    ax2 = fig.add_subplot(132)
    ax2.set_title(prefix2)
    im2 = ax2.imshow(synthetic,cmap=plt.cm.jet)
    # Set the same colorscale bounds as im1
    #norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
    #                            vmax=data[np.isfinite(data)].max())
    #im2.set_norm(cb1.norm) #alternative to manually accessing norm as above
    fig.colorbar(im2)

    ax3 = fig.add_subplot(133)
    ax3.set_title(prefix3)
    im3 = ax3.imshow(residual,cmap=plt.cm.jet)
    im3.set_norm(cb1.norm)
    fig.colorbar(im3)

    # Platform-independent Slider to go through IGs
    axcolor = 'lightgoldenrodyellow' # other options?
    axIG = plt.axes([0.25, 0.015, 0.65, 0.03], axisbg=axcolor)
    sIG = Slider(axIG, 'IG#', 0, len(self.Pairs), valinit=0, valfmt='%i',
                 closedmin=True, closedmax=True, dragging=True)
    #print dir(sIG)
    def update(val):
        val = int(val)
        dates = str(self.Pairs[val][0]) + ' ' + str(self.Pairs[val][1]) #just take pairsstring?
        ig = self.Igrams[dates]
        #ig = self.Igrams.values()[sIG.val-1]
        #dates = ig.Rsc['PAIR']
        #data = tools.load_half(ig,2)
        
        #dates = self.PairsString[val]
        #print dates
        #ig = self.Igrams[dates]
        #print ig
            
        data = tools.load_binary_old(Timeseries, ig, prefix=prefix1)
        synthetic = tools.load_binary_old(Timeseries, ig, prefix=prefix2)
        residual = tools.load_binary_old(Timeseries, ig, prefix=prefix3)
        im1.set_data(data)
        im2.set_data(synthetic)
        im3.set_data(residual)
        
        norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
                                    vmax=data[np.isfinite(data)].max())
        im1.autoscale() # autoscale colorbar to new data
        #im2.set_norm(norm) #match synthetic colorbar with Rmp
        im2.autoscale()
        #im3.autoscale()
        im3.set_norm(norm)
        bigtitle.set_text(dates)
        #ax2.set_title(self.Pairs[int(sIG.val)-1])
        plt.draw()
        
    def onpress(event):
        if event.key not in ('n', 'p'): return
        if event.key=='n': newval = sIG.val + 1
        else:  newval = sIG.val - 1
        
        if newval < sIG.valmin: newval = sIG.valmax
        if newval > sIG.valmax: newval = sIG.valmin
        sIG.set_val(newval) # update() automatically called
    
    fig.canvas.mpl_connect('key_press_event', onpress)
    sIG.on_changed(update)
    sIG.set_val(0)
    #plt.show()


def igrams_sharing_date(ts, date, nrow=1, prefix='RmpMskDef', cbar='each',
                normalize=True, ylabeldate=True, fwidth=17.0, fheight=11.0,
                mask=False):
    """ Plot all igrams containing the specified date """
    row, col = np.where(ts.Set.Pairs == date)
    selections = row
    length = ts.Set.Length
    width = ts.Set.Width
    
    igrams = float(selections.size)
    ncol = np.ceil(igrams/nrow).astype(np.int)
    fig = plt.figure(figsize=(fwidth,fheight))
    fig.suptitle('Interferograms with date={}'.format(date), fontsize=14, fontweight='bold')
    
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (nrow, ncol),
                    direction="row",
                    axes_pad = 0.75,
                    add_all=True,
                    label_mode = 'all', #'all', 'L', '1'
                    share_all = True,
                    cbar_location="top", #top,right
                    cbar_mode='each', #each,single,None
                    cbar_size=0.1,#"7%",
                    cbar_pad=0.0#,"1%",
                    )
    if mask:
        badpix = np.load(mask)
    for i,ignum in enumerate(selections):
        ax = grid[i]
        ig = ts.Set[ignum]
        print(ig.Name)
        data = tools.load_binary_old(ts, ig, prefix=prefix)
        if mask:
            data[badpix] = np.nan
        im = ax.imshow(data, cmap=plt.cm.jet)#, norm=master_norm)
        ax.cax.colorbar(im)
        cmin = int(data[np.isfinite(data)].min())
        cmax = int(data[np.isfinite(data)].max())
        ax.cax.set_xticks([cmin,0,cmax])
        
        if ylabeldate:
            ax.set_ylabel(ig.Rsc['DATE12'])
        
        ax.text(0.9,0.95,str(ignum),
                #weight='bold',
                ha='right',
                va='top',
                bbox=dict(facecolor='white'),
                transform=ax.transAxes)

        ax.tick_params(labelbottom=1,labeltop=0,labelleft=0,labelright=1,
                        bottom=1,top=1,left=1,right=1)
    
    #don't show grid frames without data...
    Nextra = grid.ngrids - ts.Set.Nig
    if Nextra > 0:
        for ax in grid[-Nextra:]:
            #print axtop
            
            ax.set_visible(False)
            ax.cax.set_visible(False)


def map_subsets(Timeseries, file='pix2subsets_svd.mat'):
    """ Map rank of B matrix for each pixel. Number of subsets given by:
    L = (N+1) - Rank(A) where N is the number of dates in inversion """
    self = Timeseries
    path = os.path.join(Timeseries.OutDir, file)
    ss = tools.load_mat(path, self.Set.Length, self.Set.Width)
    ss = ss.reshape((self.Set.Length, self.Set.Width),order='F')
    ss[ss==0] = np.nan
    pcolor(ss, '# Subsets per pixel')
    return ss


def map_date_coverage(Timeseries, prefix='RmpMskDef', platform=None, ax=None):
    """ Map number of dates in timeseries inversion of each pixel """
    #NOTE: also stored in pix2dates_svd.mat
    self = Timeseries
    datatype = np.dtype('<f4')
    width = self.Set.Width
    length = self.Set.Length
    size = width*length
    
    if platform == None:
        platform = '' #for plot titile
        iglist = self.Set.PairsString
    elif platform == 'ERS':
        iglist = self.Set.query('PLATFORM', 'ERS1') + self.Set.query('PLATFORM', 'ERS2')
    elif platform == 'ERS1':
        iglist = self.Set.query('PLATFORM', 'ERS1')
    elif platform == 'ERS2':
        iglist = self.Set.query('PLATFORM', 'ERS2')
    elif platform == 'Envisat':
        iglist = self.Set.query('PLATFORM', 'Envisat')
    ndates = len(iglist)
    print('{} total dates'.format(ndates))
    
    #NOTE: old method, based on matlab output... more efficient, but hard to
    #distinguish ERS from Envisat
    #ndates = np.zeros((length,width),dtype=datatype)
    #path = os.path.join(Timeseries.OutDir,'pix2cumdef_svd.mat')
    #cumdef = tools.load_mat(path, length, width)
    #ndates = np.sum(np.isfinite(cumdef),0) #ndates b npixels
    #ndates = np.reshape(ndates,(length,width),order='F')
    
    # New method: make an array of lists... VERY SLOW, but works
    #NOTE: think of ways to speed this up!
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    dates = np.empty(size, dtype=np.object)
    filler(dates, dates)
    for name in iglist:
        ig = self.Set.Igrams[name]
        data = tools.load_binary_old(self, ig, prefix=prefix).flatten()
        indGood, = np.isfinite(data).nonzero() #NOTE: common needed b/c returns tuple
        for pix in indGood:
            if ig.Rsc['DATE1'] not in dates[pix]:
                dates[pix].append(ig.Rsc['DATE1'])
            if ig.Rsc['DATE1'] not in dates[pix]: 
                dates[pix].append(ig.Rsc['DATE2'])
    len_func = np.frompyfunc(len,1,1)
    dqmap = len_func(dates)
    dqmap = np.reshape(dqmap,(length, width)).astype(int)
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    # Discretized colorbar b/c usually aren't too many dates
    #pcolor(dqmap, '{0} {1} Dates Per Pixel'.format(self.Set.Track, platform))
    discrete_jet = cmap_discretize(plt.cm.jet, ndates)
    im = ax.imshow(dqmap, cmap=discrete_jet)
    cbar = fig.colorbar(im, ticks=list(range(ndates+1)))
    ax.set_title('{0} {1} Dates Per Pixel'.format(self.Set.Track, platform))
    plt.show()
    
    return dqmap


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = np.linspace(0,1.,N)
    # N+1 indices
    indices = np.linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = np.array(cdict[key])
        I = scipy.interpolate.interp1d(D[:,0], D[:,1]) #NOTE: requires scipy
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = np.zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap('colormap',cdict,1024)



def map_interferogram_coverage(Timeseries, prefix='RmpMskDef', platform=None, ax=None):
    """ Map number of interferograms coherent per pixel"""
    #NOTE: could also just go off pix2cumdef matrix
    self = Timeseries
    datatype = np.dtype('<f4')
    width = self.Set.Width
    length = self.Set.Length
    dqmap = np.zeros((length,width),dtype=datatype)
    
    if platform == None:
        platform = '' #for plot titile
        iglist = self.Set.PairsString
    elif platform == 'ERS':
        iglist = self.Set.query('PLATFORM', 'ERS1') + self.Set.query('PLATFORM', 'ERS2')
    elif platform == 'ERS1':
        iglist = self.Set.query('PLATFORM', 'ERS1')
    elif platform == 'ERS2':
        iglist = self.Set.query('PLATFORM', 'ERS2')
    elif platform == 'Envisat':
        iglist = self.Set.query('PLATFORM', 'Envisat') 
    print('{} total interferograms'.format(len(iglist)))
    for ig in iglist:
        data = tools.load_binary_old(self, self.Set.Igrams[ig], prefix=prefix)
        indGood = np.isfinite(data)
        dqmap += indGood
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    # Image with discrete colormap
    im = plt.imshow(dqmap,cmap=plt.cm.jet) #NOTE: what pcolor does
    ax.set_title('{0} {1} Interferograms Per Pixel'.format(self.Set.Track, platform))
    cb = plt.colorbar()
    plt.show()
    
    return dqmap

# NOTE: not working with current ipython --pylab...
'''
def point_lasso(fignum=1,xd=None,yd=None):
    """ Get point coordinates by drawing a lasso around them... very cool """
    print 'Left-click to draw a boundary around points you want to isolate'
    sys.stdout.flush()
    
    data = [Datum(*xy) for xy in np.random.rand(100, 2)]
    fig = plt.figure()
    ax = fig.add_subplot(111, xlim=(0,1), ylim=(0,1), autoscale_on=False)
    lman = LassoManager(ax, data)
    
    #fig = plt.figure(fignum)
    #ax = fig.add_subplot(111, autoscale_on=False)
    #ax = fig.get_axes()[0] #assumes only 1 axes instance
    #lines = ax.get_lines()[0]
    #xd = lines.get_xdata()
    #yd = lines.get_ydata()
    #data = [Datum(x,y) for x,y in zip(xd,yd)]
    #lman = LassoManager(ax, data)
    plt.show()
'''

def vmap_ginput(Timeseries, fignum=1):
    """ Interactively plot timeseries from a figure """
    fig = plt.figure(fignum)
    ax = fig.gca()
    
    # Set up second figure for profiles
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_title('Timeseries for pixel:')
    
    done = False
    while not done:
        print('Select a timeseries point by left-clicking mouse on image')
        points = np.asarray(fig.ginput(),dtype=int)
        ax.plot(points[:,0],points[:,1],'k.',scalex=False,scaley=False)
        index = points[0][::-1] #extract tuple, swap order
        timeseries_separate(Timeseries, index=index, ax=ax2, annotate=True)
        
        # tap any key to stop, mouse click selects a new point
        print('Left-click in image to continue, Press a key in the image to exit')
        done = plt.waitforbuttonpress()
        ax2.clear()
        ax.lines.pop() #remove point marker
        #fig2.show()


def basemap_ginput(fignum, geofile, interp='nn',nsamp=500):
    """Extract profile of arbitrary straight line from a map overlay"""
    print('Select start and end points by left-clicking')
    sys.stdout.flush()
    fig = plt.figure(fignum)
    ax = fig.gca()
    
    # Get user-selected points
    points = np.asarray(fig.ginput(2))
    start = points[0]
    end = points[1]
    
    # Draw profile line on plot
    start=overlay['start']
    end=overlay['end']
    ax.plot([start[0],end[0]], [start[1],end[1]], 'r.-') 
    ax.text(start[0],start[1],overlay['plabel'][0],fontsize=10, fontweight='bold')
    ax.text(end[0],end[1],overlay['plabel'][1],fontsize=10, fontweight='bold')
    
    # Profile in separate plot
    fig2 = plt.figure()
    ax2 = fig.add_subplot(111)
    x0, y0 = latlon2image(geo, start[0], start[1])
    x1, y1 = latlon2image(geo, end[0], end[1])
    length = int(np.hypot(x1-x0, y1-y0)) #sample each pixel line passes through
    x = np.linspace(x0, x1, length)
    y = np.linspace(y0, y1, length)
    print(x0,y0)
    zi = phase[y.astype(np.int), x.astype(np.int)] #check transpose?
    
    # Plot profile w/ kilometers on x-axis
    length_km = latlon2range.latlon2range(start[1],start[0],end[1],end[0])
    x_km = np.linspace(0,length_km,zi.size)
    #ax2 = fig.add_subplot(212)
    ax2.plot(x_km,zi,'r.')
    #profile = ax2.plot(zi) #x-axis in terms of pixels
    atL = AnchoredText(overlay['plabel'][0],prop=dict(size=10,weight='bold'), frameon=True, loc=2)
    atL.patch.set_boxstyle("round, pad=0.0, rounding_size=0.2")
    ax2.add_artist(atL)
    atR = AnchoredText(overlay['plabel'][1],prop=dict(size=10,weight='bold'), frameon=True, loc=1)
    atR.patch.set_boxstyle("round, pad=0.0, rounding_size=0.2")
    ax2.add_artist(atR)
    #ax2.set_ylabel(r'$V_{los}$ (cm/yr)')
    #ax2.set_xlabel(r'$Distance$ (km)')
    if pLabel: ax2.set_ylabel('velocity (cm/yr)')
    ax2.set_xlabel('distance (km)')
    


def profile_ginput(fignum=None, interp='nn', nsamp=500, showline=True):
    """ Use ginput to draw a profile line between two arbitrary points on
    an imshow() figure"""
    #NOTE: if georeferenced used latlon2range for axes in km?
    print('Select start and end points by left-clicking')
    sys.stdout.flush() #NOTE: ensures print statement shows up before function exits
    
    if fignum == None:
        fig = plt.gcf()
    else:
        fig = plt.figure(fignum)
    ax = fig.gca()
    im = ax.get_images()[0]
    z = im.get_array()
    
    # Have user select profile points
    points = np.asarray(fig.ginput(2))
    x = points[:,0]
    y = points[:,1]
    
    # Draw profile line on plot
    if showline:
        start = points[0]
        end = points[1]
        ax.plot([start[0],end[0]], [start[1],end[1]], 'k.-',scalex=0,scaley=0) 
        ax.text(start[0],start[1],"A",fontsize=14, fontweight='bold')
        ax.text(end[0],end[1],"A'",fontsize=14, fontweight='bold')
        plt.draw()
    
    # Extract image values from the profile line
    #http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array
    if interp == 'nn':
        length = int(np.hypot(x[1]-x[0], y[1]-y[0]))
        xsamp = np.linspace(x[0],x[1],length).astype(int)
        ysamp = np.linspace(y[0],y[1],length).astype(int)
        zsamp = z[ysamp, xsamp]
        ax.plot(xsamp,ysamp,'k.-',lw=2)
        
    elif interp == 'cubic':
        xsamp = np.linspace(x[0],x[1],nsamp)
        ysamp = np.linspace(y[0],y[1],nsamp)
        zsamp = scipy.ndimage.map_coordinate(z,np.vstack((xsamp,ysamp)))
    
    # Open a new figure with the profile
    #NOTE: would be nice to show distance
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(zsamp,'k.-')
    
    #label A & A'
    atL = AnchoredText("A",prop=dict(size=10,weight='bold'), frameon=True, loc=2)
    #atL.patch.set_boxstyle("round, pad=0.0, rounding_size=0.2")
    ax.add_artist(atL)
    atR = AnchoredText("A'",prop=dict(size=10,weight='bold'), frameon=True, loc=1)
    #atR.patch.set_boxstyle("round, pad=0.0, rounding_size=0.2")
    ax.add_artist(atR)
    
    plt.title('profile')
    plt.show()
    
    return (xsamp,ysamp,zsamp)


def profile(array, val, ax='x'):
    """ Make a simple profile along a row or column of array in radar coordinates"""
    fig = plt.figure()
    im = plt.imshow(array,cmap=plt.cm.jet)
    cb = plt.colorbar()
    if ax == 'x':
        plt.axvline(val,c='k')
        profile = array[:,val]
    else:
        plt.axhline(val,c='k')
        profile = array[val,:] 
    
    fig = plt.figure()
    plt.plot(profile,'b-',lw=2) #just show line, not points 'b.-' show points
    plt.axhline(0,c='k')
    plt.title('Profile {0}={1}'.format(ax,val))
    #plt.ylabel('deformation')
    plt.xlabel('Pixel #')
    plt.grid()
    #plt.show()
    
    # return datapoints!
    return profile


def profile_asterix(data, center=None, nprofiles=5, clim=None):
    """ Draw series of profiles through center point. Default center for t6089
    x is column number, y is row # center=(311,68), nprofiles=5, clim=(-0.5,1)"""
    fig = plt.figure(figsize=(17,11))
    ax1 = fig.add_subplot(121)
    im = plt.imshow(data,cmap=plt.cm.jet)
    cb = plt.colorbar()
    if clim:
        plt.clim(clim)    
    plt.xlabel('col #')
    plt.ylabel('row #')
    
    ax2 = fig.add_subplot(122)
    plt.axhline(0,c='k')
    plt.title('profiles')
    plt.ylabel('deformation')
    plt.xlabel('pixel')
    
    nrow, ncol = data.shape 
    if center:
        r0, c0 = center
    else:
        r0 = nrow/2
        c0 = ncol/2
        
    #profiles = {}
    #colors = ['b', 'g', 'r', 'c', 'm','y','k']
    slopes = np.linspace(0, np.pi, nprofiles)
    for i,rad in enumerate(slopes[:-1]): #don't repeat 0 & pi    
        # Add profile line to interferogram
        #print i, rad, colors[i]
        #special case division by zeros
        if rad == 0: #could also do m=np.inf
            start = (0, r0)
            end = (ncol, r0)
        elif rad == np.pi/2:
            start = (c0, 0)
            end = (c0, nrow)
        else:
            m = np.tan(rad)
            leftIntercept = r0 + m*-c0  #NOTE: imshow takes care of axes flipping automatically!
            rightIntercept = r0 + m*(ncol-c0)
            start = (0, leftIntercept)
            end = (ncol, rightIntercept)
        ax1.plot([start[0],end[0]], [start[1],end[1]], scalex=0, scaley=0)
        
        # Add profile to adjacent plot
        #NOTE: mean, probably more representative
        length = np.floor(np.hypot(start[0]-end[0], start[1]-end[1]))  #sample each pixel line passes through
        cols = np.linspace(start[0], end[0]-1, length) #NOTE end-2 to make sure indexing works
        rows = np.linspace(start[1], end[1]-1, length)
        
        # Radial-plot
        radii = np.hypot(cols-c0, rows-r0)
        
        # East Positive (to check for E-W symmetry)
        if rad == np.pi/2: #special case for vertical profile
            radii[np.where(rows>r0)] *= -1
        else:
            radii[np.where(cols<c0)] *= -1
        
        # North Positive
        #if rad == 0:
        #    radii[np.where(cols<c0)] *= -1
        #else:
        #    radii[np.where(rows>r0)] *= -1
        
        # not sure why there are indexing errors:
        good = (rows <= data.shape[0]) & (cols <= data.shape[1])
        rows = rows[good]
        indrows = rows.astype(np.int)
        cols = cols[good]
        indcols = cols.astype(np.int)
        pPoints = data[indrows, indcols]
        ax2.plot(radii[good], pPoints, marker='.')
        
    #ax1.plot(c0,r0, marker='s', mec='k', mew=2, mfc='none', scalex=0, scaley=0)
    ax1.plot(c0,r0,'ko', ms=2, scalex=0, scaley=0)
        
    

def profile_swath(data, val, npix=3, ax=0):
    """
    plot an average profile taken from a swath of pixels
    npix=4 means take rows or cols +/-4 from the target line
    ax=0 is x-axis, 1 is y-axis
    val is central line of swath
    """
    nrow, ncol = data.shape
    fig = plt.figure(figsize=(11,8.5))
    #im = plt.imshow(data, cmap=plt.cm.jet)
    
    if ax == 1:
        plt.axvline(color='gray',linestyle='dashed')
        data = data[:,val-npix:val+npix]
    else:
        plt.axhline(color='gray',linestyle='dashed')
        data = data[val-npix:val+npix,:] 
    
    swath_mean = np.mean(data, axis=ax)
    #swath_median = np.median(data, axis=ax)
    swath_std = np.std(data, axis=ax)
    
    ax = fig.add_subplot(111)
    ax.plot(swath_mean, 'k-', lw=2, label='mean')
    ax.plot(swath_mean + swath_std, 'k--', label='+/- 1std')
    ax.plot(swath_mean - swath_std, 'k--')
    
    plt.show()
    
    return swath_mean
    

def profile_crosshair(data):
    """ averaged row & column profiles in radar coordinates"""
    nrow, ncol = data.shape
    fig = plt.figure(figsize=(8,11))
    ax = fig.add_subplot(111)
    im = plt.imshow(data, cmap=plt.cm.jet)
    #cb = plt.colorbar()
    plt.suptitle('mean profiles in radar coordinates', fontsize=14, fontweight='bold')
    plt.ylabel('row #')
    plt.xlabel('col #')
    
    # create new axes on the right and on the top of the current axes.
    divider = make_axes_locatable(ax)
    ax_x = divider.append_axes("top", size=1.2, pad=0.2, sharex=ax)
    ax_y = divider.append_axes("right", size=1.2, pad=0.2, sharey=ax)
    
    # Row profile (Top)
    mean_rows = np.mean(data, axis=0)
    median_rows = np.median(data, axis=0)
    std_rows = np.std(data, axis=0)
    
    ax_x.plot(mean_rows, 'k-', lw=2, scalex=False, scaley=False)
    ax_x.plot(mean_rows + std_rows, 'k--', label='+/- 1std', scalex=False, scaley=False)
    ax_x.plot(mean_rows - std_rows, 'k--', scalex=False, scaley=False)
    ax_x.axhline(c='r')
    for tl in ax_x.get_xticklabels():
        tl.set_visible(False)
    ax_x.set_yticks([-0.5, 0, 0.5])
    ax_x.set_ylabel('cm/yr')
    
 
    # Column Profile (Right)
    mean_cols = np.mean(data, axis=1)
    median_cols = np.median(data, axis=1)
    std_cols = np.std(data, axis=1)

    ax_y.plot(mean_cols, np.arange(nrow), 'k-', lw=2, scalex=False, scaley=False) #errorbar
    ax_y.plot(mean_cols + std_cols, np.arange(nrow), 'k--', label='+/- 1std', scalex=False, scaley=False)
    ax_y.plot(mean_cols - std_cols, np.arange(nrow), 'k--', scalex=False, scaley=False)
    ax_y.axvline(c='r')
    for tl in ax_y.get_yticklabels():
        tl.set_visible(False)
    ax_y.set_xticks([-0.5, 0, 0.5])
    #ax_y.set_xlim(-1,1)
    ax_y.set_xlabel('cm/yr')



def profile_collapsed(data, center, yex=(-0.5,1.0), nbins=100, markevery=10, ax=None):
    """ collapsed profile from arbitrary center point in array in radar coordinates
    NOTE: at Bob's suggestion, should normalize bins by area
    """
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
        
    nrow, ncol = data.shape
    r, c = center
    distance = np.zeros_like(data)
    RR, CC = np.indices(data.shape)
    distance = np.hypot(RR-r, CC-c)
    #pcolor(distance)
    #pcolor(array)
    
    # Plot Mean & Std deviation envelopes based on bins
    print(nbins)
    data = data.flatten()
    bins = np.linspace(distance.min(), distance.max(), nbins)
    indBins = np.digitize(distance.flatten(), bins)
    #pcolor(indBins.reshape(nrow, ncol))
    
    ave = np.zeros_like(bins)
    med = np.zeros_like(bins)
    std = np.zeros_like(bins)
    for i,bin in enumerate(bins):
        defvals = data[indBins == i+1]
        goodvals = defvals[np.isfinite(defvals)]
        ave[i] = np.mean(goodvals)
        std[i] = np.std(goodvals)
        med[i] = np.median(goodvals)
    
    #plt.plot(distance, data, 'k.', alpha=0.1) #actual data points!
    plt.plot(bins, ave, 'k-', lw=2, label='mean')
    plt.plot(bins, med, 'b-', label='median')
    plt.plot(bins, ave+std,'k--', label='+/- 1std')
    plt.plot(bins, ave-std, 'k--')
    plt.ylabel('deformation (cm)')
    plt.xlabel('# pixels from deformation center')
    plt.title('collapsed profile, radar coordinates')
    plt.ylim(yex)
    plt.axhline(c='k')
    plt.legend()
    

def plot_stack(Timeseries, stack, cumDef, cumTime):
    """ Plot output of timeseries stacking """
    self = Timeseries
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = plt.imshow(stack,cmap=plt.cm.jet)
    cb = plt.colorbar()
    cb.set_label('cm / {0:.2f}yr'.format(self.Set.Timespan))
    plt.title('simple stack')
    
    # Subplots --> cumulative deformation, cumulative time, average rate
    fig = plt.figure(figsize=(17,11))
    titlestr = '{0} Stack   {1} Interferograms   {2} : {3} '.format(self.Set.Track,
                                                               self.Set.Nig,
                                                               self.Set.Dates[0],
                                                               self.Set.Dates[-1])
    plt.suptitle(titlestr, fontweight='bold', fontsize=12)
    ax = fig.add_subplot(131)
    im = plt.imshow(cumDef,cmap=plt.cm.jet)
    cb = plt.colorbar()
    cb.set_label('cm')
    plt.title('cumulative deformation')
    
    ax = fig.add_subplot(132)
    im = plt.imshow(cumTime,cmap=plt.cm.jet)
    cb = plt.colorbar()
    cb.set_label('years')
    plt.title('cumulative time')
    
    ax = fig.add_subplot(133)
    im = plt.imshow(stack,cmap=plt.cm.jet)
    cb = plt.colorbar()
    cb.set_label('cm / {0:.2f}yr'.format(self.Set.Timespan))
    plt.title('average velocity')
    
    plt.show()


def plot_bil(Interferogram, cmap=plt.cm.bwr):
    """ Plot output of timeseries stacking """
    ar1, ar2 = roipy.tools.load_bil(Interferogram)
    
    # Subplots --> cumulative deformation, cumulative time, average rate
    fig = plt.figure(figsize=(11,8.5))
    
    plt.suptitle(Interferogram.Name, fontweight='bold', fontsize=12)
    ax = fig.add_subplot(121)
    im = plt.imshow(ar1,cmap=plt.cm.bwr)
    cb = plt.colorbar()
    #cb.set_label('cm')
    plt.title('Array 1')
    
    ax = fig.add_subplot(122)
    im = plt.imshow(ar2,cmap=plt.cm.bwr)
    cb = plt.colorbar()
    #cb.set_label('years')
    plt.title('Array 2')
    
    plt.show()


def stack(Timeseries, prefix='RmpMskDef', clim='auto', timeWeight=True, **kwargs): #, units='pixels' ax=None):
    """ Make a stack (cumulative deformation)/(cumulative time) """
    # NOTE: moved to timeseries method
    '''
    self = Timeseries
    datatype = np.dtype('<f4')
    width = self.Set.Width
    length = self.Set.Length
    weight = np.zeros((length,width),dtype=datatype)
    stack = np.zeros((length,width),dtype=datatype)
    for ig in self.Set:
        data = tools.load_binary_old(self, ig, prefix=prefix)
        indGood = np.isfinite(data)
        weight[indGood] += -float(ig.Rsc['TIME_SPAN_YEAR']) #uplift positive
        stack[indGood] += data[indGood]
    
    if not timeWeight: weight = 1
    stack = stack / weight
    pcolor(stack, '{0} Stack'.format(self.Set.Track), ig=ig, **kwargs)#, units='pixels' ax=None)
    if not clim == 'auto':
        plt.clim(clim)
    #NOTE: do data extraction under rp.tools....!
    return stack
    '''

def timespans_dq():
    """ Color points in timespan plot based on data quality"""
    #date point gets colored value based on average standard deviation of scenes
    #using that date. NOTE: Mask out regions of know signal or else large timespan
    #scenes show high variance.
    print('Not Implemented')
    
def date_variance(Timeseries, signalmask=None, prefix='RmpMskDef', ax=None):
    """ Calculate the average variance and standard deviation of all interferograms
    with a common date """
    self = Timeseries.Set
    vars = np.zeros(self.Ndate)
    stds = np.zeros_like(vars)
    nigs = np.zeros_like(vars, dtype='i4')
    
    for i,date in enumerate(self.Dates):
        row, col = np.where(self.Pairs == date)
        datevar = np.zeros(row.size)
        
        for j,ind in enumerate(row):
            unw = tools.load_binary_old(Timeseries, self[ind], prefix=prefix)
            unw_ma = np.ma.masked_array(unw, np.isnan(unw))
            if signalmask:
                signal = np.load(signalmask)
                unw_ma[signal] = ma.masked
            datevar[j] = np.mean(unw_ma)
            
        nigs[i] = row.size
        vars[i] = np.mean(datevar)
        stds[i] = np.std(datevar)
        
    dates = np.array(pltdate.num2date(self.DatesSerial))
    
    if ax==None:
        fig = plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    #ax.plot(dates, vars, 'bo', markersize=18)
    #ax.errorbar(dates, vars, stds, fmt='bo', markersize=18, ecolor='0.1')
    
    # Dots Colored by number of interferograms
    sc = ax.scatter(dates, vars, s=75, c=nigs, cmap=plt.cm.jet_r, zorder=20) #'s' can also be array for variable sizing
    ax.errorbar(dates, vars, stds, fmt=None, ecolor='0.1') #None--> just error bars
    
    #plot average variance & std range
    avevar = np.mean(vars)
    stdvar = np.std(vars)
    plt.axhline(avevar, color='k', linestyle='--', linewidth=2)
    plt.axhspan(avevar-stdvar, avevar+stdvar, linestyle='dashed', facecolor='y', alpha=0.3)
    plt.axhspan(avevar-2*stdvar, avevar+2*stdvar, linestyle='dotted', facecolor='y', alpha=0.3)
    
    #label number of interferograms containing date at top of plot
    #NOTE: alternative is to color by #IGs
    '''
    trans = ax.get_xaxis_transform() #text specified in (Data, Axes) coords
    for d,n in zip(dates,nigs.astype('S4')):
        ax.text(d,0.95,n,
                    ha='center',
                    va='top',
                    bbox=dict(facecolor='white'),
                    transform=trans)
    
    ax.set_title('Date Quality Assessment')
    ax.set_ylabel('Average Variance')
    '''
    
    # format the ticks
    ax.set_ylabel('variance')
    months = pltdate.MonthLocator()
    ax.xaxis.set_minor_locator(months)
    #more precise date string for the x axis locations in the
    ax.fmt_xdata = pltdate.DateFormatter('%Y-%m-%d')
    datemin = datetime.date(dates[0].year-1, 1, 1)
    datemax = datetime.date(dates[-1].year+1, 1, 1)
    ax.set_xlim(datemin, datemax)
    fig.autofmt_xdate() #rotation=90
    
    #NOTE: for some reason must come after autofmt_xdate
    cbar = plt.colorbar(sc) #could do horizontal position
    cbar.set_label('# of interferograms')
    
    #return dates with variance greater than 1 std to filter from time series
    #Set up another function for this (like match_date... get_iglist())
    baddates = np.where(vars > avevar+stdvar)[0]
    inds = []
    for date in baddates:
        inds.append(np.where(self.Pairs == self.Dates[date])[0])
    alligs = np.unique(np.hstack(inds))
    
    igList = []
    for ig in alligs:
        igList.append(self.PairsString[ig])
    
    return igList
    plt.show()


def timespans(Timeseries, ax=None):
    """Make a plot showing the timespans of all interferograms in a set"""
    self = Timeseries
    
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    
    ers = None
    envi = None
    alos = None
    for i, pair in enumerate(self.Set.PairsString): #loop in order
        ig = self.Set.Igrams[pair]
        date1 = ig.Date1
        date2 = ig.Date2
        IGspan = [date2, date1] #plot command handles datetime objects
        counter = [i+1,i+1]
        # Separate ERS & Envisat by color
        if ig.Rsc['PLATFORM'] == 'Envisat':
            envi, = ax.plot(IGspan, counter, 'r.-',label='Envisat')
        elif ig.Rsc['PLATFORM'].startswith('ERS'):
            #b = plt.plot(IGspan, counter, 'b.-')
            ers, = ax.plot(IGspan, counter, 'b.-',label='ERS')
        elif ig.Rsc['PLATFORM'] == 'ALOS':
            alos, = ax.plot(IGspan, counter, 'b.-',label='ALOS')
    
    # Avoid 100's of legend labels
    handles = [p for p in (ers,envi,alos) if not p==None]
    labels = [p.get_label() for p in handles]
    
    ax.set_ylim(0,i+2)
    dates = [pltdate.num2date(int(date)) for date in self.Set.DatesSerial]
    months = pltdate.MonthLocator()
    ax.xaxis.set_minor_locator(months)
    #more precise date string for the x axis locations in the lower right:
    ax.fmt_xdata = pltdate.DateFormatter('%Y-%m-%d')
    datemin = datetime.date(dates[0].year-1, 1, 1)
    datemax = datetime.date(dates[-1].year+1, 1, 1)
    ax.set_xlim(datemin, datemax)
    fig.autofmt_xdate()
    
    plt.grid(True)
    plt.title('{0} Interferogram Timespans'.format(self.Set.Track))
    plt.xlabel('Date')
    plt.ylabel('Interferogram #')
    #plt.legend(colors,names,loc='lower right') #does this work?
    #plt.legend(proxy_legend, ('Envisat','ERS'), loc='upper left')
    #plt.legend(, ('Envisat','ERS'), loc='upper left')
    plt.legend(handles,labels, loc='upper left')
    plt.show()




def summary2(Timeseries, signalmask=None, prefix='RmpMskDef'):
    """ Visual summary of a set of interferograms """
    self = Timeseries.Set
    timespans = np.array([abs(ig.Timespan) for ig in self]) #prob don't need abs()
    baselines = np.array([ig.Bperp for ig in self])
    if signalmask:
        customMask = np.load(signalmask)
        indSignal = (customMask==1)
    mean = np.zeros(len(self))
    std = np.zeros(len(self))
    #NOTE: call calc_statistics() routine
    for i,ig in enumerate(self):
        #unw = tools.load_half(ig)
        unw = tools.load_binary_old(Timeseries, ig, prefix=prefix)
        unw_ma = np.ma.masked_array(unw, np.isnan(unw))
        if signalmask:
            unw_ma[indSignal] = ma.masked
        mean[i] = np.mean(unw_ma)
        std[i] = np.std(unw_ma)
    
    fig = plt.figure(figsize=(11,8))
    ax = fig.add_subplot(141)
    marker= 'bo'
    ax.plot(mean, marker)
    #ax.errorbar(np.arange(self.Nig), mean, std,fmt=marker)
    ax.set_xlabel('Interferogram #')
    ax.set_ylabel('Mean Phase')
    
    ax = fig.add_subplot(142)
    ax.plot(timespans, std, marker)
    ax.set_xlabel('Timespan')
    #ax.set_ylabel('Mean Phase')
    
    ax = fig.add_subplot(143)
    ax.plot(baselines,std, marker)
    ax.set_xlabel('Perp Basline')
    #ax.set_ylabel('Mean Phase')
    
    ax = fig.add_subplot(144)
    ax.hist(std)
    ax.set_xlabel('DEM phase fit')
    
    plt.suptitle('{} Summary'.format(prefix))
    plt.savefig('summary.pdf')


def summary(Timeseries, signalmask=None, prefix='RmpMskDef'):
    """mean phase w/ standard deviation for rmp removed & filtered
    signal box is for t6089"""
    self = Timeseries.Set
    timespans = np.array([abs(ig.Timespan) for ig in self]) #prob don't need abs()
    baselines = np.array([ig.Bperp for ig in self])
    if signalmask:
        customMask = np.load(signalmask)
        indSignal = (customMask==1)
    mean = np.zeros(len(self))
    std = np.zeros(len(self))
    #NOTE: call calc_statistics() routine
    for i,ig in enumerate(self):
        #unw = tools.load_half(ig)
        unw = tools.load_binary_old(Timeseries, ig, prefix=prefix)
        unw_ma = np.ma.masked_array(unw, np.isnan(unw))
        if signalmask:
            unw_ma[indSignal] = ma.masked
        mean[i] = np.mean(unw_ma)
        std[i] = np.std(unw_ma)
    
    fig = plt.figure(figsize=(11,8))
    ax = fig.add_subplot(131)
    marker= 'o'
    ax.errorbar(np.arange(self.Nig), mean, std,fmt=marker,ecolor='0.75')
    ax.set_xlabel('Interferogram #')
    ax.set_ylabel('Mean Phase')
    ax.axhline(0.0, c='k')
    ax.axhline(5.62/2, c='k')
    
    ax = fig.add_subplot(132)
    #timespans in ascending order, NOTE, sorting automatically happens b/c scatter data
    #t_sort = timespans.argsort()
    #ax.errorbar(timespans.sort(), mean[t_sort], std[t_sort], fmt='b.')
    ax.errorbar(timespans, mean, std,fmt=marker,ecolor='0.75') #gray color
    ax.set_xlabel('Timespan')
    #ax.set_ylabel('Mean Phase')
    
    ax = fig.add_subplot(133)
    #baselines in ascending order
    #b_sort = baselines.argsort()
    #ax.errorbar(baselines.sort(), mean[b_sort], std[b_sort], fmt='b.')
    ax.errorbar(baselines, mean, std,fmt=marker,ecolor='0.75')
    ax.set_xlabel('Perp Basline')
    #ax.set_ylabel('Mean Phase')
    
    #ax = fig.add_subplot(144)
    #baselines in ascending order
    #b_sort = baselines.argsort()
    #ax.errorbar(baselines.sort(), mean[b_sort], std[b_sort], fmt='b.')
    #ax.errorbar(baselines, mean, std,fmt=marker, ecolor='0.75')
    #ax.set_xlabel('Perp Basline')
    
    plt.suptitle('{} Summary'.format(prefix))
    plt.savefig('summary.pdf')
   
'''    
def vmap_timeseries_new(Timeseries, cumdef=None, clim='auto'):
    """ Use pix2cumdef.mat array from matlab routines """
    self = Timeseries
    width = self.Set.Width
    length = self.Set.Length
    print width, length
    if cumdef == None:
        path = os.path.join(Timeseries.OutDir,'pix2cumdef_svd.mat')
        cumdef = tools.load_mat(path, length, width)
    
    tsVel = np.zeros(length*width)
    slope = np.zeros_like(tsVel)
    rsq = np.zeros_like(tsVel)
    stderr = np.zeros_like(tsVel)
    for i,pix in enumerate(tsVel):
        slope[i], rsq[i], stderr[i] = tools.calc_fit(Timeseries, cumdef, pix)
    
    tsVel = np.reshape(tsVel, (self.Set.Width, self.Set.Length))
    pcolor(tsVel, '{0} Best-fit Linear Velocity'.format(self.Set.Track))
    if not clim == 'auto':
        plt.clim(clim)
    
    return tsVel
'''

def vmap_timeseries(Timeseries, result='pix2avevel_svd.mat', clim='auto'):
    """ Read pix2avevel.mat from matlab routines """
    self = Timeseries
    path = os.path.join(self.OutDir, result)
    tsVel = -1 * tools.load_mat(path, self.Set.Length, self.Set.Width)
    #tsVel = np.fliplr(tsVel.reshape((self.Set.Length, self.Set.Width),order='F'))
    
    if self.Set.Orbit == 'descending':
        tsVel = np.rot90(tsVel,2)
    
    pcolor(tsVel, '{0} Best-fit Linear Velocity'.format(self.Set.Track))
    if not clim == 'auto':
        plt.clim(clim)
    
    return tsVel


def colorbar_only(vmin,vmax,outname='colorbar.png',figsize=(4,1),
                  cbsize=[0.05,0.5,0.9,0.2],cmap=None, label='cm/yr',
                  orient='horizontal',extend='both',transparent=0,show=False):
    """ Save colorbar as an image w/ white or transparent background.
    Default colormap = jet
    figsize in inches
    cbsize in [left, bottom, %length, %height]
    """
    print(vmin, vmax)
    if orient == 'vertical':
        figsize = (1,4)
        cbsize = [0.05,0.05,0.1,0.9]
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(cbsize)
    
    if cmap == None:
        cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    
    cb = mpl.colorbar.ColorbarBase(ax,
                                cmap=cmap,
                                norm=norm,
                                extend=extend,
                                orientation=orient,
                                )
    cb.set_label(label)
    
    #Show & save
    if show:
        plt.show()
    
    plt.savefig(outname,
                transparent=transparent,
                 bbox_inches='tight', #doesn't work for ps output...
                 )
    print('output {}'.format(outname))





'''
def timeseries(Timeseries, cumdef=None, sub=(1515,135),order=1):
    """ Plot the timeseries for a specified pixel """
    
    self = Timeseries.Set
    length = self.Length
    width = self.Width
    if cumdef == None:
        path = os.path.join(Timeseries.OutDir,'pix2cumdef_svd.mat')
        cumdef = tools.load_mat(path, length, width)
    
    dates = np.array([pltdate.num2date(int(date)) for date in self.DatesSerial])
    datatype = np.dtype('<f4')
    pixel = np.ravel_multi_index(sub,(length, width))
    #deformation = cumdef[pixel,:]
    deformation = cumdef[:,pixel]
    indDates = np.isfinite(deformation)
    X = dates[indDates]
    Y = -1 * deformation[indDates].transpose()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X,Y,'k.',markersize=0.5, markevery=1)
    titlestr = 'Cumulative Deformation for Pixel={0}'.format(sub)
    ax.set_title(titlestr, fontsize=14, fontweight='bold')
    ax.set_xlabel('Date')
    ax.set_ylabel('Displacement (cm)')
    plt.hold(True)
    
    if order: # Plot best-fit line
        coef = np.polyfit(X,Y,order)
        fit = np.poly1d(coef)
        xmin,xmax = ax.get_xlim()
        points = np.linspace(xmin,xmax,20)
        line_fit = ax.plot(points,fit(points),'--')
        equation = 'y = %.3fx + %.3f' % (coef[0],coef[1])
        ax.text(0.95, 0.95, equation, 
                horizontalalignment='right',
                verticalalignment='center',
                transform = ax.transAxes) 
    #plt.draw()
'''




def calc_timeseries(Timeseries, cumdef=None, index=None, showloc=False):
    """ Extract cumulative deformation values for plotting """
    #NOTE: should move to rp.tools...
    self = Timeseries.Set
    length = self.Length
    width = self.Width
    if cumdef == None:
        path = os.path.join(Timeseries.OutDir,'pix2cumdef_svd.mat')
        cumdef = tools.load_mat(path)
    
    #NOTE: descending tracks are read in flipped up-down,
    #ascending tracks are read in flipped left-right
    if self.Orbit == 'descending':
        indmat = (length - index[0], index[1])
    else:
        indmat = (index[0], width - index[1])
    
    dates = np.array(pltdate.num2date(self.DatesSerial))
    #datatype = np.dtype('<f4')
    pixel = np.ravel_multi_index(indmat,(length,width),order='F') #numpy >1.6
    #pixel = (sub[0] * width) + sub[1] #python linear index
    #pixel = (index[0]-1) * length + index[1] #matlab linear index (columnwise)
    #print width, length, sub, pixel
    deformation = cumdef[:,pixel]
    indDates = np.isfinite(deformation)
    X = dates[indDates]
    Y = -deformation[indDates].transpose()
    
    #show pixel location on velocity map
    if showloc:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        path = os.path.join(Timeseries.OutDir, 'pix2avevel_svd.mat')
        tsVel = -1 * tools.load_mat(path, self.Length, self.Width)
        tsVel = np.fliplr(tsVel.reshape((self.Length, self.Width),order='F'))
        # NOTE: commented out 10/2/13 b/c taken care of in load_mat? -SH
        #if self.Orbit == 'descending':
        #    tsVel = np.rot90(tsVel,2)
        #pcolor(tsVel, '{0} Best-fit Linear Velocity'.format(self.Track))
        #NOTE: must use plt.imshow to interact w/ figure in ipython...
        im = plt.imshow(tsVel, cmap=plt.cm.jet)
        cb = plt.colorbar(im)
        ax.plot(index[1],index[0], marker='s', mec='k', mew=2, mfc='none', scalex=0, scaley=0)
    
    return X,Y,indDates


def timeseries(Timeseries, index=('row','col'), cumdef=None, order=1, annotate=None,
               save=False, style='default', biannual=False, showloc=False, fit='linear',
               showfit=False, errorbar=False, yerr=1, ax=None):
    """ Plot the timeseries for a specified pixel default is uturuncu t10
    data is the pix2cumdef array
    """
    #NOTE: vertical line or rectangle for date range on timeseries plot:
    #eq = datetime.date(2009,12,6)
    #plt.axvline(eq)
    #eq1 = datetime.date(2009,12,6)
    #eq2 = datetime.date(2009,12,19)
    #plt.axvspan(eq1,eq2,facecolor='g',alpha=0.5)
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    
    self = Timeseries.Set
    
    if style != 'default': #e.g. poster vs. paper font sizes
        set_style(style)
    
    X,Y,indDates = calc_timeseries(Timeseries, cumdef, index, showloc)
    X_ = self.DatesSerial[indDates]


    
    if errorbar:    
        ax.errorbar(X,Y,yerr=yerr, fmt='o') #plus or mis 0.5 yunit errorbars
    else:
        ax.plot(X,Y,'bo') #makersize=5 overwrites whatever st by mpl.rc() command
    
    if annotate:
        title = '{0} pixel={1}'.format(self.Track, index)
        ax.set_title(title)
        #ax.set_xlabel('Date',fontsize=14) #set_style takes care of this
        ax.set_ylabel('Displacement (cm)') #alternatively set in rcParams
    #minor tick dates by month
    months = pltdate.MonthLocator()
    ax.xaxis.set_minor_locator(months)
    #more precise date string for the x axis locations in the
    ax.fmt_xdata = pltdate.DateFormatter('%Y-%m-%d')
    #plt.hold(True)
    ax.grid(True)
    
    #Linear Fit 1
    if fit == 'linear':    
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X_,Y)
        poly = np.poly1d((slope,intercept))
        r_sq = r_value**2
        if showfit: 
            lineFit = ax.plot(X_, poly(X_), 'k-',label='linear')
            #equation = r'$y = %.3fx + %.3f$' % (slope,intercept)
            slope_cmyr = slope*365 #cm/day to cm/yr
            print(slope_cmyr,r_sq)
            slopeText = r'$m = {:.1f} cm/yr$'.format(slope_cmyr) #1 decimal point sig fig
            fitQual = r'$R^2 = {:.2f} $'.format(r_sq)
            ax.text(0.95, 0.10, slopeText,  #.95 upper right
                    horizontalalignment='right',
                    verticalalignment='center',
                    transform = ax.transAxes) 
            ax.text(0.95, 0.05, fitQual,  #.9 upper right
                    horizontalalignment='right',
                    verticalalignment='center',
                    transform = ax.transAxes)
    '''
    #Linear Fit 2
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X_[:8],Y[:8])
    poly = np.poly1d((slope,intercept))
    r_sq = r_value**2
    lineFit = ax.plot(X_, poly(X_), 'k-',label='linear')
    #equation = r'$y = %.3fx + %.3f$' % (slope,intercept)
    slope_cmyr = slope*365
    print slope_cmyr,r_sq
    #Linear Fit 3
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X_[8:],Y[8:])
    poly = np.poly1d((slope,intercept))
    r_sq = r_value**2
    lineFit = ax.plot(X_, poly(X_), 'k-',label='linear')
    #equation = r'$y = %.3fx + %.3f$' % (slope,intercept)
    slope_cmyr = slope*365
    print slope_cmyr,r_sq
    '''
    
    if fit == 'quadratic':
        print('fitting 2nd-order polynomial')
        coef = np.polyfit(X_,Y,2)
        poly = np.poly1d(coef)
        xmin,xmax = ax.get_xlim()
        points = np.linspace(xmin,xmax,20)
        line_fit = ax.plot(points,poly(points),'k-',label='quadratic')
        if showfit:
            #plot equation
            slopeText = r'$y = %3.2ex^2 %3.2ex + %3.2e$' % (coef[0],coef[1],coef[2])
            ax.text(0.95, 0.15, slopeText,  #.95 upper right
                    horizontalalignment='right',
                    verticalalignment='center',
                    transform = ax.transAxes)
    
    # Plot date labels every other year for long timseries
    if biannual:
        annual = pltdate.YearLocator(1)
        biannual = pltdate.YearLocator(2)
        ax.xaxis.set_major_locator(biannual)
        ax.xaxis.set_minor_locator(annual)
    
    fig.autofmt_xdate() #Angle dates #NOT working if ax passed to another function...?
    plt.ylabel(r'$\Delta LOS$ (cm)')
    #plt.legend(loc=2)
    #plt.draw()
    if save:
        save_figure('timeseries{0}'.format(sub))#,outdir='/home/scott/figures/AGU_2011')


''' 
def save_figure(name,nfig=None,fmt='eps',alpha=False):
    """Save figure to file"""
    outdir = '/home/scott/figures/AGU_2011/timeseries/'#NOTE: change to initialization attribute
    path = os.path.join(outdir,name)
    savename ='{0}.{1}'.format(path,fmt)
    plt.savefig(savename,
                format=fmt,
                transparent=False,
                bbox_inches='tight',
                pad_inches=0.05)
    print 'saved: {0}'.format(savename)
    
def set_style(style='presentation'):
    """change matplotlib plotting defaults for font sizes for posters or papers"""
    if style == 'presentation': 
        #print mpl.rcParams
        font = {'size': 16}#,
                #'weight':'bold'} bold is too bold
        mpl.rc('font',**font)
        mpl.rcParams['lines.markersize'] = 8 #not working??? must reload script in pylab
        print "changing plot settings to 'presentation'"

def reset_style():
    mpl.rcdefaults() #whatever is in the matplotlib rc file
'''



def timeseries_separate(Timeseries, index=(1515,135), data=None, order=1,
                        annotate=False, errorbar=False, yerr=1, ax=None):
    """ ERS dates plotted as green triangles, ENVISAT dates plotted as blue circles"""
    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    
    self = Timeseries.Set
    X,Y,indDates = calc_timeseries(Timeseries, data, index)
    
    dates_string = self.Dates.astype('S8') 
    dates_serial = self.DatesSerial
    X_serial = dates_serial[indDates]
    X_string = dates_string[indDates]
    
    ersIGs = self.query('PLATFORM','ERS1') + self.query('PLATFORM','ERS2')
    envisatIGs = self.query('PLATFORM','Envisat')
    ersList = set(' '.join(ersIGs).split())
    envisatList = set(' '.join(envisatIGs).split())
    commonList = ersList.intersection(envisatList) #dates with both ERS & ENVISAT acquisitions
    #print commonList
    
    x_ers = []; y_ers = [] #NOTE: alternatively could change point color after plotting
    x_envisat = []; y_envisat = []
    x_common= []; y_common = []
    for date,plotdate,value in zip(X_string,X,Y):
        if date in commonList:
            x_common.append(plotdate)
            y_common.append(value)
        elif date in ersList:
            x_ers.append(plotdate)
            y_ers.append(value)
        else:
            x_envisat.append(plotdate)
            y_envisat.append(value)
        

    
    if errorbar:
        ax.errorbar(x_envisat,y_envisat, yerr=yerr,fmt='bo',label='Envisat')
        ax.errorbar(x_ers,y_ers,yerr=yerr,fmt='g^',label='ERS')
        ax.errorbar(x_common,y_common,yerr=yerr,fmt='r*', label='Both')
    else:
        ax.plot(x_envisat,y_envisat,'bo',markersize=5,label='Envisat')
        ax.plot(x_ers,y_ers,'g^',markersize=5,label='ERS')
        ax.plot(x_common,y_common,'r*',markersize=5,label='Both')
    
    
    ax.legend(loc=2,numpoints=1) #upper left legend
    
    #minor tick dates by month
    months = pltdate.MonthLocator()
    ax.xaxis.set_minor_locator(months)
    #more precise date string for the x axis locations in the
    ax.fmt_xdata = pltdate.DateFormatter('%Y-%m-%d')
    
    if annotate:
        title = '{} Time Series, Pixel={}'.format(self.Track,index)
        ax.set_title(title)
        ax.set_xlabel('Date')
        ax.set_ylabel('Displacement (cm)')
    #ax.hold(True)
    ax.grid(1)
    
    if order: # Plot best-fit line
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X_serial,Y)
        fit = np.poly1d((slope,intercept))
        r_sq = r_value**2
        lineFit = ax.plot(X_serial, fit(X_serial), 'k-')
        #equation = r'$y = %.3fx + %.3f$' % (slope,intercept)
        slope_cmyr = slope*365
        slopeText = r'$m = %3.1f cm/yr$' % slope_cmyr #One significant figure
        fitQual = r'$R^2 = %.3f$' % r_sq
        ax.text(0.95, 0.10, slopeText,  #.95 upper right
                horizontalalignment='right',
                verticalalignment='center',
                transform = ax.transAxes) 
        ax.text(0.95, 0.05, fitQual,  #.9 upper right
                horizontalalignment='right',
                verticalalignment='center',
                transform = ax.transAxes)

    fig.autofmt_xdate()
    plt.show()



def hist(data, signalmask=None, nbins=50):
    """ Plot histogram of all values in an interferogram """
    #unw = tools.load_binary_old(Timeseries, ig, prefix=prefix)
    unw_ma = np.ma.masked_array(data, np.isnan(data))
    if signalmask:
        signal = np.load(signalmask)
        unw_ma[signal] = ma.masked
    stats={}
    min = np.min(unw_ma)
    max = np.max(unw_ma)
    mean = np.mean(unw_ma)
    std = np.std(unw_ma)
    var = std**2
    
    clim = (mean-2*std, mean+2*std)

    histdata = data[np.isfinite(data)].flatten()
    fig = plt.figure()
    plt.hist(histdata,bins=nbins)
    plt.axvline(x=clim[0], color='r')
    plt.axvline(x=clim[1],color='r')
    #plt.title('Histogram of {0} {1}'.format(ig.Name, prefix))
    plt.xlabel('value')
    plt.ylabel('# pixels')    
    
    # for rescaling images to demphasize outliers & bring out middle values/background
    return clim


def baseline_perp(Set):
    """ bargraph of baselines for each interferogram in set"""
    fig = plt.figure()
    self = Set
    x = np.arange(self.Nig)
    y = [ig.Bperp for ig in self]
    y.sort()
    black_grey = np.tile(('0','0.5'), self.Nig)
    ax = plt.bar(x,y,color=black_grey, edgecolor='k')
    plt.title('{0} Perpendicular Baselines'.format(self.Track))
    plt.ylabel('B_perp (m)')
    plt.xlabel('Interferogram #')
    
    
def baseline_time(Set):
    """ bargraph of time separation for each ig in set"""
    #NOTE: may want to convert to days...
    fig = plt.figure()
    self = Set
    x = np.arange(self.Nig)
    y = [abs(ig.Timespan) for ig in self]
    y.sort()
    black_grey = np.tile(('0','0.5'), self.Nig)
    ax = plt.bar(x,y,color=black_grey, edgecolor='k')
    plt.title('{0} Temporal Baselines'.format(self.Track))
    plt.ylabel('B_temporal (years)')
    plt.xlabel('Interferogram #')
    
    
def baseline_both(Set):
    """ Sort by increasing temoral baseline, after Casu et al 2008 fig2"""
    #Make use of mpl.axis_grid1 and numpy recarray & sort by timespan
    print('Low on the priority list...')
    
    
def baseline_plot(Set):
    """ python version of the classic baseline plot. assumes set of
    interferograms already exists & just connects the dots"""
    fig = plt.figure()
    self = Set
    #Make use of mpl.axis_grid1 and numpy recarray & sort by timespan
    #dates = np.array([pltdate.num2date(int(date)) for date in self.DatesSerial])
    #pltdates = dt.datetime.strptime(self.Baselines.keys(), '%Y%m%d') 
    dates = pltdate.num2date(self.DatesSerial)
    baselines = [self.Baselines[date.strftime('%Y%m%d')] for date in dates]
    ax = plt.plot(dates, baselines, 'ko', markersize=8)
    plt.hold(True)
    
    # NOTE: find date pairs & connect the dots, brute force see LineCollection for alternative
    # NOTE: should add text labels to points
    d = [0,0]
    b = [0,0]
    for i,serdate in enumerate(self.DatesSerial):
        d[0] = pltdate.num2date(serdate)
        b[0] = baselines[i]
        for match in tools.match_date(self, serdate): #best in tools or data method?
            d[1] = pltdate.num2date(match)
            b[1] = self.Baselines[d[1].strftime('%Y%m%d')]
            plt.plot(d, b,'b-')
    plt.title('{0} Baseline Plot'.format(self.Track))
    plt.ylabel('B_temporal (years)')
    plt.xlabel('SAR Acquisition Date')


def data_frames(ts, fwidth=17.0, fheight=11.0, nrow=4, normalize=True,
                    show=True, outname='proofsheet.pdf', cbar='each',
                    prefix='RmpMskDef', ylabeldate=False):
    """ Make a 'proof sheet' of all frames in time series default 17x11 paper
    format. cbar can be 'single' (normalized to min/max of set?)"""
    #NOTE: add ability to layout just a subregion...
    #NOTE: implement colorbar=single
    #NOTE: add ability to do a query on timeseries instance, and only show matches
    
    #path = os.path.join(ts.OutDir,'pix2cumdef_svd.mat')
    length = ts.Set.Length
    width = ts.Set.Width
    #cumdef = -1.0 * tools.load_mat(path, length, width)
    #base = np.zeros((1,length*width)) 
    #cumdef = np.vstack((base,cumdef))
    #want first date to have zeros instead of nans
    #cumdef[0,:] = 0.0
    
    #scale all images to final deformation limits
    #data = cumdef[-1]
    #master_norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
    #                                   vmax=data[np.isfinite(data)].max())
    
    # Landscape-style layout, just two columns b/c length is long!
    igrams = float(ts.Set.Nig)
    ncol = np.ceil(igrams/nrow).astype(np.int)
    fig = plt.figure(figsize=(fwidth,fheight))
    fig.suptitle('{} Deformation (cm)'.format(ts.Set.Track), fontsize=14, fontweight='bold')
    
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (nrow, ncol),
                    direction="row",
                    axes_pad = 0.25,
                    add_all=True,
                    label_mode = 'all', #'all', 'L', '1'
                    share_all = True,
                    cbar_location="top", #top,right
                    cbar_mode='each', #each,single,None
                    cbar_size=0.1,#"7%",
                    cbar_pad=0.0#,"1%",
                    )

    for i,ig in enumerate(ts.Set):
        ax = grid[i]
        print(ig.Name)
        #data = np.fliplr(cumdef[i].reshape((length,width),order='F')) #NOTE: 'F' needed b/c matlab
        #data = tools.load_half(ig) #NOTE: load Rmp/Msk etc
        data = tools.load_binary_old(ts, ig, prefix=prefix)
        im = ax.imshow(data, cmap=plt.cm.jet)#, norm=master_norm)
        ax.cax.colorbar(im)
        cmin = np.nanmin(data)
        cmax = np.nanmax(data)
        ax.cax.set_xticks([cmin,0,cmax])
        
        #ax.set_title(str(i))
        #ax.set_title(ig.Rsc['DATE12'], fontsize=10, weight='bold')
        #datepairs = ts.Set.match_date(date).size
        if ylabeldate:
            ax.set_ylabel(ig.Rsc['DATE12'])
        
        #at = AnchoredText(str(i), prop=dict(size=8, weight='bold'), loc=1, #pad=0.2,
        #                  borderpad=0.0, frameon=True)
        #at.patch.set_boxstyle("round, pad=0.1, rounding_size=0.2") #rounded box
        #ax.add_artist(at)
        ax.text(0.9,0.95,str(i),
                #weight='bold',
                ha='right',
                va='top',
                bbox=dict(facecolor='white'),
                transform=ax.transAxes)
        
        #no ticks
        ax.tick_params(labelbottom=0,labeltop=0,labelleft=0,labelright=0,
                        bottom=0,top=0,left=0,right=0)
    
    #single colorbar based on last frame
    #ax.cax.colorbar(im)
    #cbar = grid.cbar_axes[0].colorbar(im)
    #cbar.set_label_text('Unwrapped Phase')

    #don't show grid frames without data...
    Nextra = grid.ngrids - ts.Set.Nig
    if Nextra > 0:
        for ax in grid[-Nextra:]:
            #print ax
            ax.set_visible(False)
            ax.cax.set_visible(False)
    
    plt.savefig(outname)
    outpath = os.path.abspath(outname)          
    print('saved {}'.format(outpath))        
         
         
            
def synthetic_frames_asc(ts, fwidth=11.0, fheight=8.0, nrow=2, normalize=True,
                               show=True):
    """ Make a 'proof sheet' of all frames in time series """
    #NOTE: add ability to layout just a subregion...
    
    path = os.path.join(ts.OutDir,'pix2cumdef_svd.mat')
    length = ts.Set.Length
    width = ts.Set.Width
    cumdef = -1.0 * tools.load_mat(path, length, width)
    #base = np.zeros((1,length*width)) 
    #cumdef = np.vstack((base,cumdef))
    #want first date to have zeros instead of nans
    cumdef[0,:] = 0.0
    
    #scale all images to final deformation limits
    data = cumdef[-1]
    master_norm = mpl.colors.Normalize(vmin=data[np.isfinite(data)].min(),
                                       vmax=data[np.isfinite(data)].max())
    
    # Landscape-style layout, just two columns b/c length is long!
    dates = float(ts.Set.Ndate)
    ncol = np.ceil(dates/nrow).astype(np.int)
    fig = plt.figure(figsize=(fwidth,fheight))
    fig.suptitle('{} Synthetic Reconstruction'.format(ts.Set.Track), fontsize=14, fontweight='bold')
    
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (nrow, ncol),
                    direction="row",
                    axes_pad = 0.25,
                    add_all=True,
                    label_mode = '1', #'all'
                    share_all = True,
                    cbar_location="right",
                    cbar_mode='single', #each
                    cbar_size=0.1,#"7%",
                    cbar_pad=0.1#,"1%",
                    )

    for i,date in enumerate(ts.Set.Dates):
        ax = grid[i]

        data = np.fliplr(cumdef[i].reshape((length,width),order='F')) #NOTE: 'F' needed b/c matlab
        im = ax.imshow(data,cmap=plt.cm.jet, norm=master_norm)
        ax.set_title(date, fontsize=10, weight='bold')
        datepairs = ts.Set.match_date(date).size
        at = AnchoredText(datepairs, prop=dict(size=8, weight='bold'), loc=1, pad=0.2,
                          borderpad=0.5, frameon=True)
        #at.patch.set_boxstyle("round, pad=0.1, rounding_size=0.2") #rounded box
        ax.add_artist(at)
        #no ticks
        ax.tick_params(labelbottom=0,labeltop=0,labelleft=0,labelright=0,
                        bottom=0,top=0,left=0,right=0)
    #colorbar based on last frame
    cbar = grid.cbar_axes[0].colorbar(im)
    cbar.set_label_text('Cumulative Displacement (cm)')

    #don't show grid frames without data...
    Nextra = grid.ngrids - ts.Set.Ndate
    if Nextra > 0:
        for ax in grid[-Nextra:]:
            #print ax
            ax.set_visible(False)



 # BELOW ROUTINES ARE FOR A PROFILE THROUGH A GEOREFERENCED MAP #
 # ------------------------------------------------------------ #
def project2utm(Point, epsg=32719, printout=True):
    """ Transform point object in lat/lon to UTM 192 WGS84 (epsg=32719)"""
    oldCoords = Point.GetPoint()
    
    srs = Point.GetSpatialReference()
    trs = osr.SpatialReference()
    trs.ImportFromEPSG(epsg)
    latlon2utm = osr.CoordinateTransformation(srs, trs)
    Point.Transform(latlon2utm)

    if printout:
        newCoords = Point.GetPoint()
        print(oldCoords, newCoords)
    #easting = newGeometry.GetX()
    #northing = newGeometry.GetY()
    #return easting, northing
    #return xyzTuple
    #return Point


def get_distance(point1, point2, transObject):
    """Distance between two points in km"""
    point1.Transform(transObject)
    point1.GetPoint()
    point2.Transform(transObject)
    point2.GetPoint()
    distance = point1.Distance(point2) / 1000.0 #in kilometers
    print(distance)
    return distance


def make_profile(path, rasters=['stack282','dem_t282'],unit='degrees'):
    """Need Path to ESRI format shape file output by ProfileFromLine in QGIS
    plot is a list of shapefile attribute raster values to put into profile
    assume profile distance measured in degrees & therefore approximate distance
    in meters"""
    ds = ogr.Open(path)
    layer = ds.GetLayer()
    if unit == 'degrees': #otherwise utm
        srs = layer.GetSpatialRef()
        trs = osr.SpatialReference()
        trs.ImportFromEPSG(32719) #WGS84, UTM 19South
        latlon2utm = osr.CoordinateTransformation(srs, trs)

    n = layer.GetFeatureCount()
    X = np.zeros(n); Y = np.zeros(n) 
    Z = np.zeros(n); Z1 = np.zeros(n)


    # Most precise (other than saving file in UTM to begin with!
    # Get first point & then loop through rest
    feature_p = layer.GetNextFeature()
    point_p = feature_p.GetGeometryRef()
    point_p.Transform(latlon2utm)
    #NOTE: it is essential that when using geometry methods eg Distance, the
    #parent feature object must not be destroyed or reassigned before calling
    #the method
    feature = layer.GetNextFeature()
    while feature:
        i = feature.GetFID()    
        point_p = feature_p.GetGeometryRef()
        point = feature.GetGeometryRef()
        point.Transform(latlon2utm)

        X[i] = point.Distance(point_p) / 1000.0 #distance in km
        Z[i] = feature.GetField(rasters[0])
        Z1[i] = feature.GetField(rasters[1])

        feature_p = feature
        feature = layer.GetNextFeature()
    
    # LATLON2KM APPROXIMATION
    """
    deg2km = 105.0 #~105km/deg at 22S
    X_approx = np.zeros(n)
    layer.ResetReading()
    feature = layer.GetNextFeature()
    while feature:
        i = feature.GetFID()    
        X_approx[i] = feature.GetField('accumDst') * deg2km  #distance in km
        #Z[i] = feature.GetField(rasters[0])         #LOS velocity cm/yr
        #Z1[i] = feature.GetField(rasters[1])        #elevation in m
        feature = layer.GetNextFeature()
    """

    # MAKE PLOT
    X = np.cumsum(X) #cumulative sum distances
    xlabel = 'Along profile distance (km)'
    #X = np.arange(n) #pixels
    #xlabel = 'Along profile distance (km)'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.set_title("t2282 Profile B-B'", fontsize=14, fontweight='bold')

    veloz = ax.plot(X,Z,'bo-',markersize=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('LOS velocity (cm/yr)')
    plt.hold(True)

    ax2 = plt.twinx()
    elev = ax2.plot(X,Z1,'ko-',markersize=2)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Elevation (m)')
    plt.legend([veloz,elev],['velocity','elevation'],loc='lower left')

    # Best to just set interactively after plot creation in pylab
    #ax.set_xlim(0,110)
    ax.set_ylim(-1.5,1.5)
    ax2.set_ylim(3000,6000)
    
    #header = 'longitude latitude aveVelocity(cm/yr) elevation(m)' #only in version 2.0
    np.savetxt('profile.out',np.vstack((X,Y,Z,Z1)).transpose(),fmt='%.4f')#,header=header)
    # NOTE: easy to reproduce plot loading data with np.loadtxt()    
    
    
