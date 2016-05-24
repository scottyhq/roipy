"""
I/O tools for interacting with ROI_PAC data arrays
"""

from __future__ import print_function
import matplotlib.pyplot as plt #for saving np array as an image
#import scipy.stats as ss
import numpy as np
import numpy.ma as ma
from scipy.io import loadmat
from scipy.io import netcdf_file as netcdf
from time import sleep

import shutil
import os
import subprocess
import zipfile
import sys
import re
try:
    import h5py
except:
    print("install 'h5py' python module for reading matlab .mat files")

import scipy.stats

import roipy.data #NOTE: should keep this separate? cross referencing...
#import roipy.plot #can't cross reference? since roipy.plot imports roipy.tools
import roipy.models #check that models doesn't import tools?

#from osgeo import osr, ogr, gdal
import osr
import ogr
import gdal
#NOTE: must manually add location of projection files for gdal? or set environment variable...
#/home/scott/Enthought/Canopy_64bit/User/share/gdal/gcs.csv
#sys.path.append('/home/scott/epd-7.2/share/gdal')


# ============================================================================ #
# ARRAY INPUT
# ============================================================================ #
def load_grd(path):
    """ scipy to load GMT file (really netCDF type... or build gdal with netcdf.."""
    grdFile = netcdf(path,'r')
    try: #method if output by GDAL
        print('trying method 1...')
        shape = grdFile.variables['dimension'][::-1]
        data = grdFile.variables['z'][::-1]
        data = np.reshape(data,shape) #NOTE:order='F' b/c from matlab
    except:
        print('using method 2...')
        #NOTE: byte reading is off here...
        shape = (grdFile.dimensions['y'], grdFile.dimensions['x'])
        data = grdFile.variables['z'][::-1]
        data = np.reshape(data,shape)

    return data

def calc_fit(tsInstance, cumdef=None, pixel=None, fit='linear'):
    """Make average  """
    self = tsInstance.Set
    length = self.Length
    width = self.Width
    #cumdef array is flipped lr when read in so
    # what a pain in the ass...
    #indmat = (index[0], width - index[1])
    #dates = np.array(pltdate.num2date(self.DatesSerial))
    #pixel = np.ravel_multi_index(indmat,(length,width),order='F') #numpy >1.6
    deformation = cumdef[:,pixel]
    indDates = np.isfinite(deformation)
    #X = dates[indDates]
    Y = -deformation[indDates].transpose()
    X_ = self.DatesSerial[indDates]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X_,Y)
    rsq = r_value**2

    return slope, rsq, std_err



def make_synthetic(Geogram, center=(-67.22,-22.27), source='mogi', source_params=dict(xoff=1e2,yoff=1e3,d=25e3,dV=30e6,nu=0.25,output='cart'), incpath='.'):
    """ create synthetic deformation source on same georeferenced grid as interferogram
    NOTE: should have a simple way to do this for both radar coords & georeferenced coords
    Example:
        geo = rp.data.Geogram('geo_stack89_8rlks.unw')
        ux,uy,uz,ulos = rp.tools.make_synthetic(geo)
        rp.plot.synthetic_grid(ux,uy,uz,ulos)
    """
    # Extract grid from interferogram
    Lon, Lat = get_grid(Geogram) #lat/lon grid
    #Lat = np.flipud(Lat) # not sure why this does anything
    X = latlon2range(Lat, center[0]*np.ones_like(Lon), Lat, Lon,output='lon') #in meters
    Y = latlon2range(center[1]*np.ones_like(Lat), Lon, Lat, Lon,output='lat')
    #R = latlon2range(center[1]*np.ones_like(Lat), center[0]*np.ones_like(Lon), Lat, Lon,output='hypot')

    # Calculate synthetic
    if source == 'mogi':
        ux,uy,uz = roipy.models.calc_mogi(X,Y,**source_params)

    # Project to los - separate function?
    if os.path.exists(incpath):
        #inc = rp.data.Geogram(incpath)
        #look, head = load_bil(inc)
        look, head = load_bil_file(incpath)
        look = np.radians(look)
        head = np.radians(head)
        ulos = -1 * (np.sin(look)*np.sin(head)*ux + np.sin(look)*np.cos(head)*uy -np.cos(look)*uz)
        return ux,uy,uz,ulos

    # NOTE: incorporate plotting command to show map grid of 4 components
    return ux,uy,uz



def dem_phase_set(Set, hgtPath = '/home/scott/data/insar/t6089/aux_files/radar_32rlks.hgt'):
    """ """
    rsq = np.zeros(Set.Nig)
    hgt = roipy.data.Interferogram(hgtPath)
    elev = load_half(hgt, half=2)
    for i,ig in enumerate(Set):
        #unw = roi_py.data.Interferogram(unwPath)
        phase = load_half(ig, half=2)
        vecPHS = phase.flatten()
        vecHGT = elev.flatten()
        masked = np.isnan(vecPHS) + np.isnan(vecHGT) # Mask NaN values
        vecPHS = vecPHS[~masked]
        vecHGT = vecHGT[~masked]
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(vecHGT,vecPHS)
        rsq[i] = r_value**2

    fig = plt.figure()
    plt.plot(rsq, 'bo')
    plt.axhline(0.2)
    plt.xlabel('Interferogram #')
    plt.ylabel('r^2')
    plt.title('DEM vs Phase fits')


def load_roifile(Interferogram):
    """ Call appropriate load() function based on file suffix """
    self = Interferogram
    ext = os.path.splitext(self.Name)[1]
    if ext in ['rmg','unw','cor','hgt','msk']:
        data = load_bil(self)
    elif ext in ['slc', 'int', 'amp']:
        data = load_cpx(self)
    elif ext in ['dem']:
        data = load_binary(self.Path,dtype='int16')
    else: # assume single array r4 format
        data = load_r4(self.Path) # assume single array r4
    return data


def load_r4(path, length=None, width=None):
    """ Load a ROI_PAC r4 format file """
    data = np.fromfile(path, dtype='f4')
    data = np.flipud(data.reshape((length,width))) # N-up, W-left orientation
    return data


def load_bil_file(path):
    ig = roipy.data.Interferogram(path)
    mag,phs = load_bil(ig)
    return mag,phs


def load_bil(Interferogram, dtype='f4', nanval=0):
    """Load amplitude and phase array from BIL files output by ROI_PAC
    (.cor, .unw, .hgt, .msk)

    Inputs
    Interferogram = roi_py.data.Interferogram instance OR path to file
    TODO: dims = (width,length) tuple to override .rsc dimensions
    geo = True if .unw is georeferenced to maintain correct array orientation

    Outputs
    (leftArray, rightArray) = 2D numpy arrays

    Example
    look,head = incidence file
    amplitude,phase = .unw file
    slope,elevation = .hgt file
    """
    self = Interferogram
    width = int(self.Rsc['WIDTH'])
    length = int(self.Rsc['FILE_LENGTH'])
    fullwidth = width*2

    data = np.fromfile(self.Path, dtype=dtype, count=length*fullwidth)
    data = np.reshape(data,(length,fullwidth))
    data[data==nanval] = np.NaN # Assumes NaNs are 0
    #print data.dtype

    #leftArray = np.ascontiguousarray(data[:,0:width])
    mag = data[:,0:width]
    mag = orient_array(self, mag)

    #rightArray = np.ascontiguousarray(data[:,width:fullwidth])
    phs = data[:,width:fullwidth]
    phs = orient_array(self, phs)

    return mag, phs


def load_cpx(Interferogram, dtype='complex64', nanval=0):
    """ Based on Roiview's ReadData() function
    """
    self = Interferogram
    data = np.fromfile(self.Path, dtype=self.DataType)
    data = np.reshape(data,(self.Length,self.Width))

    # NOTE: could customize here to read any ROI_PAC format!
    # Also add nanval when loading data
    #data = data[np.arange(0, self.Length, 1)]

    #mag = np.hypot(data.real, data.imag) #same as abs()?
    #phs = np.arctan2(data.imag, data.real) #same as angle()?
    mag = np.abs(data)
    phs = np.angle(data)

    mag[mag==nanval] = np.NaN
    phs[phs==nanval] = np.NaN

    # rotate descending (see roi_browser.py)

    return mag, phs


def load_ascii(path):
    """Read ascii format array into a numpy array """
    data = np.loadtxt(path)
    return data


def load_mat(path, length=None, width=None):
    """Read .mat output from matlab routines into an array
    Inputs
    path = path to any .mat file to be opened

    Outputs
    data = 2D numpy array

    Examples:

    #. Load cumulative deformation array

    >>> cumdef = load_mat('pix2cumdef.mat')

    #. Load velocity map & reshape to correct dimensions
    >>> velmap = load_mat('pix2avevel.mat',set.Length,set.Width)
    """
    #print path, length, width

    #Latest version of matlab
    try:
        print(path)
        f = h5py.File(path)
        data = list(f.values())[0]
        print('loading with h5py')

        # Sparse Matrix
        if isinstance(data, h5py.highlevel.Group):
            data = f[prefix]['data'].value
            row_ind = f[prefix]['ir'].value
            row_ptr = f[prefix]['jc'].value

            sparseData = scipy.sparse.csc_matrix((data,row_ind, row_ptr))
            data = np.array(sparseData.todense())
            data = np.ascontiguousarray(data) #Not sure this does anything
            if length and width:
                print(length, width)
                data = data.reshape((length,width),order='F')
        # Full Matrix:
        else:
            data = np.array(data)
            data = np.ascontiguousarray(data)
            #NOTE: careful if reorienting done in some other commands...
            if length and width:
                print(length,width)
                data = np.fliplr(data.reshape((length, width),order='F'))
        f.close()

    # Older versions of matlab
    except Exception as e:
        print(e)
        print('loading with scipy.io.loadmat')
        print('WARNING: function specific to matlab timeseries ouput')
        dict = loadmat(path)
        #data = dict.popitem()[1]
        try:
            data = dict['pix2cumdef']
        except:
            data = dict['pix2avevel']
        data = data.T #transpose to match orientation from h5py
        #NOTE: reorienting done in other commands... maybe should do it here.

    return data

def load_binary_old(tsInstance, igInstance, path=None, prefix=None,dtype='<f4',orient=None):
    """Read .bin output from matlab routines into an array.

    Matlab saves to Column-Major (F) binary format but Numpy saves to
    Row-Major (C) binary by default. This routine automatically converts
    the binary array to C-order in memory for efficient subsequent operations

    Inputs
    prefix = matlab output prefix as a string. ['Def', 'MskDef', 'RmpMskDef']
    path = full path to a binary file

    Outputs
    data = 2D numpy array

    Examples

    #. Load a particular .bin output file

    >>> dem = ig.load_binary('t282.dem') #eg from get_SRTM.pl
    >>> RmpMskDef = ig.load_binary('RmpMskDef')
    """
    self = igInstance
    if prefix:
        path = os.path.join(tsInstance.ProcDir, prefix + '_' + self.Name[:-3] + 'bin')
    #print path
    width = int(self.Rsc['WIDTH'])
    length = int(self.Rsc['FILE_LENGTH'])
    data = np.fromfile(path, dtype=dtype)

    if self.Georeferenced:
        data = np.flipud(np.transpose(data.reshape((width, length)))) # FOR georeferenced files saved in .bin matlab format
    else:
        data = np.fliplr(np.transpose(data.reshape((width, length)))) #added fliplr.. why?
    #data = np.transpose(data.reshape((width, length)))
    #NOTE: geo_ files need 'ORBIT_DIRECTION added
    if 'ORBIT_DIRECTION' not in self.Rsc:
        self.Rsc['ORBIT_DIRECTION'] = None
    if (self.Rsc['ORBIT_DIRECTION'] == 'descending'):
        #print 'reorienting'
        data = np.rot90(data,2)

    data[data==0] = np.NaN #NOTE: use NaNs or no?
    return np.ascontiguousarray(data)


def load_binary(path=None, dtype='<f4'):
    """Read .bin output from matlab routines into an array
    """
    data = np.fromfile(path, dtype=dtype)
    data = np.transpose(data)
    #data[data==0] = np.NaN #NOTE: use NaNs or no?
    #data = np.ascontiguousarray(data)
    return data



def load_half(Interferogram, half=2, path=None, dtype='f4', convert2cm=False):
    """Read a single array from a .rmg format file into memory.

    .rmg format files contain two side-by-side equal dimension arrays that
    are stored in Row-Major (C) format and therefore have interleaved data
    stored as alternating matrix rows.

    Inputs
    path = path to any .rmg format file (.rmg, .unw, .hgt)
    half = (1,2) load either first (1) or second (2) array

    Outputs
    data = 2D numpy array

    Examples

    #. Load the amplitude array in radar coods

    >>> amp = load_half('rect_090822-073109.unw,1)

    #2. Load the phase array from a georeferenced interferogram

    >>> geophs = load_half('geo_stack_32rlks.unw,2,True)
    """

    self = Interferogram
    width = int(self.Rsc['WIDTH'])
    length = int(self.Rsc['FILE_LENGTH'])

    if path is None:
        path = self.Path
    data = np.zeros(self.Shape, dtype=dtype)
    with open(path,'rb') as f:
        rows = list(range(self.Length + 1 - half))
        if half == 2:
            #NOTE: if fromfile has a 'skip' argument, could go easier...
            junk = np.fromfile(f, dtype=dtype, count=self.Width)
        for i in rows:
            data[i,:] = np.fromfile(f, dtype=dtype, count=self.Width)
            junk = np.fromfile(f, dtype=dtype, count=self.Width)

    data[data==0] = np.NaN
    data = orient_array(self,data)

    if convert2cm:
        data = data * self.Phs2cm

    return data



def load_ma(dataFile, maskFile=None):
    """ open phase and mask arrays in .npy format and return a numpy masked array"""
    data = np.load(dataFile)
    # Boolean mask array, 1=bad data, 0=good data
    if maskFile is None:
        mask = ma.nomask
    else:
        mask = np.load(maskFile)
    masked_array = ma.array(data, mask=mask)
    return masked_array



def load_gdal(inFile, data=True, geotrans=True, proj=True):
    """ Read any raster data format supported by GDAL into a numpy array
    returns data,geo,dim where data=numpy array or RGB array, geo=geotransform values,
    dim=array dimensions

    Inputs
    inFile = (str) path to georeferenced file

    Outputs
    (data,geotrans,proj)
    data = np.array containing data
    geotrans = tuple of geotransform info (upper left x, w-e pixel res, rotation, top left y, rotation, n-s pixel resolution)
    proj = EPSG integer identifier

    Example

    #. Read in raster array from a GDAL accepted format

    >>> data,geo,dim = load_gdal(<path-to-raster>)
    """
    try:
        ds = gdal.Open(inFile)
        if geotrans:
            geotrans = ds.GetGeoTransform()
        else:
            geotrans = None
        if proj:
            proj = ds.GetProjection()
        else:
            proj = None

        if data:
            nLon = ds.RasterXSize
            nLat = ds.RasterYSize

            #single band image
            if ds.RasterCount == 1:
                #data = ds.GetRasterBand(1).ReadAsArray()
                band = ds.GetRasterBand(1)
                data = band.ReadAsArray()
                bad = band.GetNoDataValue()
                data[data == bad] = np.NaN
                #data = band.ReadAsArray(0,0,nLon,nLat) #nlon, nlat only necessary if reading partial array
            elif ds.RasterCount == 3:
                r = ds.GetRasterBand(1).ReadAsArray()
                g = ds.GetRasterBand(2).ReadAsArray()
                b = ds.GetRasterBand(3).ReadAsArray()
                data = np.dstack((r,g,b))
        else:
            data = None

        print('Sucessfully read in file:\n{0}'.format(inFile))
        return data, geotrans, proj

    except:
        print('Unable to import {0} with GDAL!'.format(inFile))
        raise


def load_rsc(path=None):
    """Read from unw.rsc file & store metadata as an attribute dictionary
    named rsc"""
    metadata = {}
    # remove blank lines, retaining only lines with key, value pairs
    with open(path, 'r') as rsc:
        allLines = [line for line in rsc.readlines() if line.strip()]
        for line in allLines:
            var, value = line.split()
            metadata[var] = value

    return metadata


# ============================================================================ #
# ARRAY OUTPUT
# ============================================================================ #
def save_rsc(RscDict, filename):
    """ write Interferogram.Rsc dictionary to a file"""
    #NOTE: load_rsc can only handle two columns of data, remove triples:
    try:
        del RscDict['PAIR']
        del RscDict['PAIR_SERIAL']
    except:
        pass #keys not in dictionay

    #NOTE: figure out how to get nice column alignment
    with open(filename, 'w') as rsc:
        for item in list(RscDict.items()):
            rsc.write("{0}           {1}\n".format(*item))


def save_ma(dataFile, maskFile, mask_array):
    """ save a masked array as a data array and logical array in .npy format. Prepend an optional prefix if desired """
    np.save(dataFile, mask_array.data)
    np.save(maskFile, mask_array.mask)


def save_r4(name, array, nanvalue=0):
    """save numpy array to ROI_PAC .r4 format """
    array[np.isnan(array)] = 0
    array = np.flipud(array)
    array = array.flatten()
    array = array.astype('f4')
    array.tofile(name)
    print(name)


def save_bil(IG, outpath, amp, phs, dtype='f4', nanvalue=np.NaN):
    """ save a manipulated array back to god'ole .unw format (band-interleaved) """
    # Assign nan values
    self = IG
    (length, width) = phs.shape
    fullwidth = width*2
    phs[np.isnan(phs)] = nanvalue
    amp[np.isnan(amp)] = nanvalue
    amp = orient_array(self,amp)
    phs = orient_array(self,phs)

    # join amplitude and phas
    join = np.zeros((length, fullwidth))
    join[:,0:width] = amp
    join[:,width:fullwidth] = phs

    # Save as a
    data = join.astype(dtype)
    data.tofile(outpath)#, dtype=self.DataType)#, count=length*fullwidth) only fromfile
    #print data.dtype
    #save associated rsc file
    save_rsc(self.Rsc, outpath + '.rsc')
    #print 'save_bil Output: {0}'.format(os.path.abspath(outpath))

    return os.path.abspath(outpath)


def save_envi(Interferogram, data, outname=None):
    """write to flat binary format a la ENVI"""
    self = Interferogram
    if outname==None:
        outname = self.Path[:-4] + '.bin'
    #data = np.flipud(data) #.unw data must flip for proper GDAL reading, .bin shouldn't flip
    data.tofile(outname)  #keep dtype or change to specific int16? etc

    # Write associated header
    save_envi_header(Interferogram, data, outname)

    return outname


def save_envi_header(self, nparray, outname=None):
    """Save a numpy array as a binary matric with associated ENVI .hdr header
    origin is UPPER LEFT pixel"""
    # Save array to file with float32 binary numberd
    if outname==None:
        headerFile = self.Path[:-4] + '.bin.hdr'
    else:
        headerFile = outname + '.hdr'
    if nparray.dtype == 'int16':
        code = 2 #16-bit signed integer
    else:
        code = 4 #32-bit floating point

    self.enviParams = {'N_LON':self.Rsc['WIDTH'],
                       'N_LAT':self.Rsc['FILE_LENGTH'],
                       'UL_LON':self.Rsc['X_FIRST'],
                       'UL_LAT':self.Rsc['Y_FIRST'],
                       'D_LON':self.Rsc['X_STEP'],
                       'D_LAT':self.Rsc['Y_STEP'].strip('-'),#strip negative
                       'DTYPE':code} #http://www.brockmann-consult.de/beam/doc/help/general/BeamDimapFormat.html

    headerENVI ="""ENVI
description =   {{Generated with igplot.py}}
samples =       {N_LON}
lines =         {N_LAT}
bands =         1
file type =     ENVI Standard
data type =     {DTYPE}
interleave =    bsq
byte order =    0
sensor type =   Unknown
wavelength =    Unknown
map info =      {{Geographic Lat/Lon, 1.0000, 1.0000, {UL_LON}, {UL_LAT}, {D_LON}, {D_LAT}, WGS-84, units=Degrees}}
    """.format(**self.enviParams)
        # Write header to file header same as data file with .hdr ending
    with open(headerFile, "w") as f: #NOTE: automatically closes file
        f.write(headerENVI)


def save_image(Interferogram, outname=None, data=None, vmin=None, vmax=None,
               cmap=plt.cm.jet, format='png', worldfile=False, nodata=None, phs2cm=True):
    """Save a 2D numpy array as a png with optional associated world file
    png has transparent background by default"""
    self = Interferogram

    if type(data) != np.ndarray:
        #amp, data = load_bil(self)
        data = load_half(self,half=2)
    if not outname:
        outname = self.Name + '.' + format

    #NOTE: this works, but changes array dimensions through dpi &
    #fig = plt.figure()
    #plt.axis('off')
    #im = plt.imshow(data,cmap=plt.cm.jet)
    #plt.savefig(savename, transparent=True)
    if phs2cm:
        data = data * self.Phs2cm

    if nodata:
        data[data==nodata] = np.nan

    # Array elements=image pixels, transparent pixels are np.NaNs
    plt.imsave(outname, data, vmin=vmin, vmax=vmax, cmap=cmap, format=format)

    if worldfile:
        # [dx,dy,rotx?,roty?,xUL,yUL]
        coords = (self.Rsc['X_STEP'], self.Rsc['Y_STEP'],'0','0',self.Rsc['X_FIRST'],self.Rsc['Y_FIRST'])
        with open(outname + 'w', 'w') as worldfile:
            worldfile.write('\n'.join(coords))

    return os.path.abspath(outname)


def save_gdal(Interferogram=None, nparray=None, outfile=None, geotrans=None, proj=4326,
              nanval=None, format='GTiff', half=2):
    """ Write numpy array to any raster format supported by GDAL.

    Inputs
    Interferogram = roipy.data.Interferogram instance
    nparray = np.array associated with Interferogram
    outfile = str output name
    geotrans = tuple of (upper left x, w-e pixel res, rotation, top left y, rotation, n-s pixel resolution)
    proj = EPSG integer identifier
    format = str format string recognized by GDAL (eg. ENVI, GTiff, GMT)
    NOTE: GMT support not necessarily built into EPD python

    Examples:

    #. Save .unw as a GeoTiff in WGS84 lat/lon

    >>> ig = rp.data.Interferogram(<path-to-ig>)
    >>> ref = (-67.8716,0.00083,0.0,-21.1199,0.0,-0.000833)
    >>> rp.tools.save_gdal(ig,'example.tif',ref,4326,'Gtiff')
    """
    self = Interferogram

    suffi = {'GTiff':'tif',
             'ENVI':'bin'}
    try:
        suffix = suffi[format]
    except:
        suffix = ''

    if nparray == None:
        nparray = load_half(self, half=half)
    if outfile == None:
        outfile = self.Name + '.' + suffix
    try:
        driver = gdal.GetDriverByName(format)
        nd = driver.Register()
    except:
        print("unrecognized format, check 'gdalinfo --formats'")

    # DATA
    rows,cols = nparray.shape
    #note order: cols, rows, nbands, dtype (default GDT_Byte=uint8)
    outDataset = driver.Create(outfile, cols, rows, 1, gdal.GDT_Float32)
    outBand = outDataset.GetRasterBand(1)
    if nanval:
        print('setting nans to {}'.format(nanval))
        outBand.SetNoDataValue(nanval) #or np.NaN?
        nparray[np.isnan(nparray)] = nanval
    #order: array, xoffset, yoffset (integers)... if successful result=0
    result = outBand.WriteArray(nparray, 0, 0)

    # GEOTRANS
    if geotrans == None:
        geotrans = get_geotrans(self.Path)
    result = outDataset.SetGeoTransform(geotrans)

    # PROJECTION
    if proj:
        outSR = osr.SpatialReference()
        result = outSR.ImportFromEPSG(proj)
        outWkt = outSR.ExportToWkt()
        result = outDataset.SetProjection(outWkt)

    #Write to disk & Close
    outDataset.FlushCache()
    del outDataset

    #print 'Done saving georeferenced file:\n{0}'.format(outfile)
    return outfile


def get_geotrans(geofile):
    '''(upper left x, w-e pixel res, rotation, top left y, rotation, n-s pixel resolution)'''
    try:
        self = roipy.data.Geogram(geofile)
        ulx = float(self.Rsc['X_FIRST'])
        dx = float(self.Rsc['X_STEP'])
        xrot = 0.0
        uly = float(self.Rsc['Y_FIRST'])
        yrot = 0.0
        #NOTE: carefl with positive & negative conventions...
        dy = float(self.Rsc['Y_STEP'])
    except:
        print('ERROR Unable to load geotransform information!')
        raise

    geotrans = (ulx, dx, xrot, uly, yrot, dy)
    #print geotrans

    return geotrans


def gdal2grd(gdalFile, outName=None):
    """ Since EPD gdal doesn't have netCDF support, must call external command"""
    if outName == None:
        outName = gdalFile + '.grd'
    cmd = 'gdal_translate -of GMT {0} {1}'.format(gdalFile, outName)
    print(cmd); os.system(cmd)
    return outName


def grd2gdal(grdFile, outName=None):
    """ Since EPD gdal doesn't have netCDF support, must call external command
    intermediary command,,, or use scipy load_netcdf (see load_gmt()"""
    if outName == None:
        outName = grdFile + '.bin'
    cmd = 'gdal_translate -of ENVI {0} {1}'.format(grdFile, outName)
    print(cmd); os.system(cmd)
    return outName


def save_kmz(Geogram, data=None, outname=None, cmap=None, vmin=None, vmax=None,
             colorbar=False, nodata=None, phs2cm=True):
    """Write geocoded array as a png and save google earth kmz image overlay.
    cmap is 'jet by default. vmin & vmax autoscaled unless specified.
    colorbar=path to colorbar.png file to use as a google earth static overlay"""
    print('saving kmz file...')
    self = Geogram
    kmzname = outname
    kmlname = kmzname[:-3] + 'kml'

    pngpath = save_image(self, data=data, outname=outname[:-3] + 'png', cmap=cmap, vmin=vmin, vmax=vmax, nodata=nodata, phs2cm=True)
    pngname = os.path.basename(pngpath)


    self.Rsc['IMG'] = pngname #don't use abspath here!
    self.Rsc['CBAR'] = colorbar
    self.Rsc['N_LAT'] = float(self.Rsc['Y_FIRST'])
    self.Rsc['S_LAT'] = self.Rsc['N_LAT'] + (float(self.Length) * float(self.Rsc['Y_STEP']))
    self.Rsc['W_LON'] = float(self.Rsc['X_FIRST'])
    self.Rsc['E_LON'] = self.Rsc['W_LON'] + (float(self.Width) * float(self.Rsc['X_STEP']))
    self.Rsc['KML'] = kmlname

    #cludge fix
    if 'ORBIT_NUMBER' not in self.Rsc:
        self.Rsc['ORBIT_NUMBER'] = 'N/A'

    #NOTE: avoid use of Timespan to just show image automatically...
    headertxt = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
    <Document>
    """

    overlaytxt = """
          <GroundOverlay>
            <name> {PAIR} </name>
            <description>orbit {ORBIT_NUMBER} </description>
            <Icon>
                <href> {IMG} </href>
                <viewBoundScale> 0.75 </viewBoundScale>
            </Icon>
            <LatLonBox>
                <north> {N_LAT} </north>
                <south> {S_LAT} </south>
                <east> {E_LON} </east>
                <west> {W_LON} </west>
                <rotation>0</rotation>
            </LatLonBox>
        </GroundOverlay>
    """.format(**self.Rsc)

    # to put colorbar in lower left change overlay & screen to y=0
    # to resize so that it is always 1/2 widow size set size y=0
    # should save with bold font, high resolution, looks iffy...
    cbartxt = """
    <ScreenOverlay>
    <name>Colorbar</name>
    <Icon>
    <href> {CBAR} </href>
    </Icon>
    <overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>
    <screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>
    <rotationXY x="0" y="0" xunits="pixels" yunits="pixels"/>
    <size x="0" y="0" xunits="pixels" yunits="fraction"/>
    </ScreenOverlay>
    """.format(**self.Rsc)


    footertxt = """
    </Document>
</kml>
    """

    kml = open(kmlname, 'w')
    kml.write(headertxt)
    kml.write(overlaytxt)
    if colorbar:
        kml.write(cbartxt)
    kml.write(footertxt)
    kml.close()

    #zip together image and kml
    kmz = zipfile.ZipFile(kmzname, 'w')
    kmz.write(kmlname,compress_type=zipfile.ZIP_DEFLATED)
    kmz.write(pngname,compress_type=zipfile.ZIP_DEFLATED)
    if colorbar:
        kmz.write(colorbar,compress_type=zipfile.ZIP_DEFLATED)
    kmz.close()

    #NOTE: commented out b/c automatic time slider in GE is annoying, could be
    #cool to show development of deformation though...
    #<TimeSpan>
    #    <begin> {DATE1} </begin>
    #    <end> {DATE2} </end>
    #</TimeSpan>

    return kmzname



# ============================================================================ #
# COORDINATE TRANSFORMATIONS
# ============================================================================ #
def basemap2overlay(bmap, geotrans, x, y):
    """ Convert Basemap coordinates to row/col of image overlay array"""
    lon,lat = bmap(x,y, inverse=True)
    x0 = geotrans[0] #top left longitude
    y0 = geotrans[3] #top left latitude
    dx = geotrans[1] #pixel width
    dy = geotrans[5] #pixel height

    #row = int((lon - x0) / dx)
    #col = int((lat - y0) / dy)
    col = int((lon - x0) / dx)
    row = int((lat - y0) / dy)
    return row, col

def overlay2basemap(bmap, geotrans, row, col):
    """ Convert row/col of image overlay array to lat/lon Basemap coordinates"""
    x0 = geotrans[0] #top left longitude
    y0 = geotrans[3] #top left latitude
    dx = geotrans[1] #pixel width
    dy = geotrans[5] #pixel height

    x = x0 + (row * dx)
    y = y0 + (col * dy)
    lon,lat = bmap(x,y)
    return lon, lat


def load_overlay(path, convert2cm=False):
    """ Load a georeferenced unw or GDAL-supported file for basemap
    if loading a unw file convert2cm reads RSC to convert radians to cm
    """
    #if path.endswith('unw'):
    try:
        geo = roipy.data.Geogram(path)
        if geo.DataType == 'float32':
            image = load_half(geo, half=2, convert2cm=convert2cm)
        elif geo.DataType == 'complex64':
            mag, image = load_cpx(geo)
        geotrans = get_geotrans(path)
    #else:
    except:
        #raise
        image, geotrans, proj = load_gdal(path)

    nLat, nLon = image.shape[:2] #3rd= number of bands
    x0 = geotrans[0] #top left longitude
    dx = geotrans[1] #pixel width
    #xR = geotrans[2] #rotation
    y0 = geotrans[3] #top left latitude
    #yR = geotrans[4] #rotation
    dy = geotrans[5] #pixel height
    LL = (x0, y0 + dy*nLat) #NOTE: make sure dy is correct sign
    UR = (x0 + dx*nLon, y0)
    extent = (LL[0], UR[0], LL[1], UR[1])

    return image, extent


def calc_ramp(array, ramp='quadratic', custom_mask=None):
    """
    Remove a quadratic surface from the interferogram Subtracting
    the best-fit quadratic surface forces the background mean surface
    displacement to be zero

    Note: exclude known signal, unwrapping errors, etc. with custom_mask

    ramp = 'dc','linear','quadratic'

    returns ramp
    """
    X,Y = np.indices(array.shape)
    x = X.reshape((-1,1))
    y = Y.reshape((-1,1))

    # Work with numpy mask array
    phs = ma.masked_invalid(array)

    if custom_mask != None:
        phs[custom_mask] = ma.masked

    d = phs.reshape((-1,1))
    g = ~d.mask
    dgood = d[g].reshape((-1,1))

    if ramp == 'quadratic':
        print('fit quadtratic surface')
        G = np.concatenate([x, y, x*y, x**2, y**2, np.ones_like(x)], axis=1) #all pixels
        Ggood = np.vstack([x[g], y[g], x[g]*y[g], x[g]**2, y[g]**2, np.ones_like(x[g])]).T
        try:
            m,resid,rank,s = np.linalg.lstsq(Ggood,dgood)
        except ValueError as ex:
            print('{}: Unable to fit ramp with np.linalg.lstsq'.format(ex))


    elif ramp == 'linear':
        print('fit linear surface')
        G = np.concatenate([x, y, x*y, np.ones_like(x)], axis=1)
        Ggood = np.vstack([x[g], y[g], x[g]*y[g], np.ones_like(x[g])]).T
        try:
            m,resid,rank,s = np.linalg.lstsq(Ggood,dgood)
        except ValueError as ex:
            print('{}: Unable to fit ramp with np.linalg.lstsq'.format(ex))

    elif ramp == 'dc':
        G = np.ones_like(phs)
        m = np.mean(phs)
        print('fit dc offset')

    ramp = np.dot(G,m)
    ramp = ramp.reshape(phs.shape)

    return ramp



def radar2ground(Interferogram, inc_mean=42.0):
    """ Estimate dimensions of data in km rather than radar coordinates """
    self = Interferogram
    #NOTE: these keys aren't in all .rsc files...
    #corners = np.array([self.Rsc['LOOK_REF1'],
    #                   self.Rsc['LOOK_REF2'],
    #                   self.Rsc['LOOK_REF3'],
    #                   self.Rsc['LOOK_REF4']], dtype='float')
    #look_mean = (corners.min() + corners.max()) / 2
    #inc_mean = look_mean + 5
    slant2groundx = 1e3 * (1 / np.sin(np.deg2rad(inc_mean)))
    slant2groundy = 1e3 * float(self.Rsc['EARTH_RADIUS']) / (float(self.Rsc['EARTH_RADIUS']) + float(self.Rsc['HEIGHT']))
    col2km = float(self.Rsc['RANGE_PIXEL_SIZE']) * slant2groundx
    row2km = float(self.Rsc['AZIMUTH_PIXEL_SIZE']) * slant2groundy

    return col2km, row2km


def radar2latlon(Interferogram, row=None, col=None, KmPerDeg=110.0):
    """ APPROXIMATE lat,lon coordinate of radar pixel"""
    #see latlon2radar
    # call IntSim command line utility
    print('TODO')


def latlon2radar(Interferogram, lat=None, lon=None, KmPerDeg=110.0):
    """ APPROXIMATE lat,lon coordinate of radar pixel"""
    #NOTE: must use affine transformation image is mostly rotated and translated
    #and geo-referencing just accounts for earth's curvature...
    #see scipy.ndimage.interpolation.affine_transformation or matplotlib.transform
    print('TODO')


def latlon2range_cp(pointLat,pointLon,gridLat,gridLon):
    """ Vectorized for numpy to be called from collapsed_profile()
    Example:
        Lon,Lat = rp.tools.get_grid(Geo)
        radial_distances = latlon2range_cp(-22,-67,Lat,Lon)
    """
    req = 6378.136
    rp = 6356.751

    twopi = 2*np.pi
    piover2 = np.pi/2
    dlat1 = np.radians(pointLat)
    dlon1 = np.radians(pointLon)
    dlat2 = np.radians(gridLat)
    dlon2 = np.radians(gridLon)

    epsln = (req/rp)**2
    #NOTE: how to get something that works for both scalar and array inputs?
    #dlat1[np.abs(dlat1) < piover2] = np.arctan(epsln * np.tan(dlat1))
    if np.abs(dlat1) < piover2:
        dlat1 = np.arctan(epsln * np.tan(dlat1))

    index = np.abs(dlat2) < piover2
    dlat2[index] = np.arctan(epsln*np.tan(dlat2))[index]
    # broadcasting
    #dlat2 = np.where(np.abs(dlat2) < piover2, np.arctan(epsln * np.tan(dlat2), dlat2))

    alat = (dlat1+dlat2)/2
    dlat = dlat2-dlat1
    dlon = dlon1-dlon2
    #dlon[dlon < -pi] = dlon+twopi
    #dlon[dlon >  pi] = dlon+twopi
    dlon = np.where(((dlon < -np.pi) & (dlon > np.pi)), dlon+twopi, dlon)
    #index = (dlon < -pi) & (dlon > pi) #NOTE: can't use python 'and' here, need numpy '&'

    # scott's version
    e2 = 1 - (rp/req)**2
    q = 1 - e2*np.sin(alat)**2
    ulat = (req*(1-e2)/q**1.5)*dlat
    vlon = (req*np.cos(alat)/q**0.5)*dlon
    d = np.hypot(ulat, vlon) #euclidean distance in kilometers
    return d


def latlon2range(lat1,lon1,lat2,lon2,output='hypot'):
    """ Translation of matlab script latlon2range.m
        output can be 'lat', 'lon', 'hypot'
%Convert pairs of Earth centric latitudes and longitudes to distance
%in meridional and zonal directions and total range

%Inputs are in degrees and output is in meters
%
%       INPUTS
%       lat1  =   centric latitude of first point (point or vector)
%       lon1  =   longitude of first point
%       lat2  =   centric latitude of second point
%       lon2  =   longitude of second point
%       req   =   equatorial radius, km
%       rp    =   polar radius, km

%       OUTPUTS
%       ulat  =   range in zonal direction
%       vlon  =   range in meridional direction
%       d  =   total range
"""
    req = 6378136.0
    rp = 6356751.0
    #pi = 3.1415926535898

    twopi = 2*np.pi
    piover2 = np.pi/2
    dlat1 = np.radians(lat1)
    dlon1 = np.radians(lon1)
    dlat2 = np.radians(lat2)
    dlon2 = np.radians(lon2)

    epsln = (req/rp)**2


    ind = (np.abs(dlat1) < piover2)
    dlat1[ind] = np.arctan(epsln * np.tan(dlat1[ind]))

    ind = (np.abs(dlat2) < piover2)
    dlat2[ind] = np.arctan(epsln * np.tan(dlat2[ind]))

    alat = (dlat1+dlat2)/2
    dlat = dlat2-dlat1
    dlon = dlon1-dlon2


    ind = np.logical_or(dlon < -np.pi, dlon >  np.pi)
    dlon[ind] += twopi

    # scott's Euclidean distance version
    e2 = 1 - (rp/req)**2
    q = 1 - e2*np.sin(alat)**2
    ulat = (req*(1-e2)/q**1.5)*dlat
    vlon = (req*np.cos(alat)/q**0.5)*dlon

    if output == 'lat':
        d = ulat
    elif output == 'lon':
        d = -vlon
    else:
        d = np.hypot(ulat, vlon) #euclidean distance in km


    return d


def distance_pyproj(lat1,lon1,lat2,lon2):
    """
    distance using pyproj, which draws on PROJ.4 library
    NOTE: only works between single points, not numpy arrays...
    """
    geod = pyproj.Geod(ellps='WGS84')
    angle1,angle2,distance = geod.inv(lon1, lat1, lon2, lat2)

    #angle1 = bearing
    #angle2 = back azimuth

    return distance



def distance_haversine(point1, point2, radius=6371e3):
    """ calculate the haversine distance between points (lon,lat)
    accurate over short distances (~3m accuracy over 1km)

    usage: d = distance_haversine((lon1,lat1), (lon2,lat2))

    distance formulas from:
    http://www.movable-type.co.uk/scripts/latlong.html
    """
    #NOTE: can probably use gdal or google earth api's as well, GE apparently uses Re=6,378,137m
    lon1, lat1 = np.radians(point1)
    lon2, lat2 = np.radians(point2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2)
    distance = 2 * radius * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    return distance


def distance_cosines(point1, point2, radius=6371e3):
    """ spherical law of cosines
    accurate over large distances (~3m accuracy over 1km using/ float64) """
    lon1, lat1 = np.radians(point1)
    lon2, lat2 = np.radians(point2)

    dlon = lon2 - lon1
    distance = radius * (np.arccos(np.sin(lat1)*np.sin(lat2) +
                                   np.cos(lat1)*np.cos(lat2)*np.cos(dlon)))

    return distance


def distance_cyl(point1, point2, radius=6371e3):
    """ equirectangular projection (pythogorean theorem w/ effective cartesian
    coordinates. Fine for 'small' distances (depends also on bearing & latitude)
    Basemap 'cyl' projection
    """
    #NOTE: km2mi = 1/1.609344
    #km2Nmi = 1/1.852 #nautical mile
    lon1, lat1 = np.radians(point1)
    lon2, lat2 = np.radians(point2)

    x = (lon2 - lon1) * np.cos((lat1+lat2)/2)
    y = lat2 - lat1

    distance = np.hypot(x,y) * radius

    return distance


def bearing(point1, point2):
    """ Calculate initial bearing (forward azimuth) of great circle arc
    bearing from true North (theta=0)"""
    lon1, lat1 = np.radians(point1)
    lon2, lat2 = np.radians(point2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    y = np.sin(dlon)*np.cos(lat2)
    x = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(dlon)
    bearing = np.arctan2(y, x) #in radians (-pi, pi)

    bearingDeg = (np.degrees(bearing) + 360) % 360
    #NOTE back azimuth conversion backDeg = (bearing + 180) % 360

    return bearingDeg


def cross_track(start, end, point, radius=6371e3, dist='haversine'):
    """ shortest distance from point to great circle path (aka 'cross-track error')
    start = p1, end = p2, point = p3 (lon,lat) degrees
    NOTE: negative implies to the west of the line
    """
    if dist=='haversine':
        d13 = distance_haversine(start, point,radius=radius)
    elif dist=='cyl':
        d13 = distance_cyl(start,point,radius=radius)

    theta13 = np.radians(bearing(start,point))
    theta12 = np.radians(bearing(start,end))

    distance = radius * np.arcsin(np.sin(d13/radius)*sin(theta13-theta12))

    return distance


def distance_vicenty(point1, point2, radius=6371e3):
    """ generalized haversine for ellipsoid, expect to be most accurate
    most accurate (~1mm accuracy over 1km)
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    print('Use latlon2range function')


# ============================================================================ #
# COMMON OPERATIONS
# ============================================================================ #
def orient_array(Interferogram, array):
    '''Called from load() routines to put arrays in North=up, West=left, upper left
    pixel = [0,0] orientation'''
    #NOTE: add special case for georeferenced arrays (or is this okay?)
    self = Interferogram
    if 'ORBIT_DIRECTION' in self.Rsc:
        if self.Rsc['ORBIT_DIRECTION'] == 'ascending':
            #print 'ascending'
            array = np.flipud(array)
        elif self.Rsc['ORBIT_DIRECTION'] == 'descending':
            #print 'descending'
            array = np.fliplr(array)
    return array


def fit_surface(array, type='quad', show=False):
        """ fit a quadratic surca"""
        x,y = np.indices(array.shape)
        xf = x.flatten()
        yf = y.flatten()
        df = array.flatten()
        g, = np.isfinite(df).nonzero() #NOTE: not sure if comma is needed..

        G = np.array([np.ones(np.size(g)), xf[g], yf[g], xf[g]**2, xf[g]*yf, yf[g]**2])
        m, res, rank, singval = np.linalg.lstsq(G.T, df)
        print('fit coefficients:', m)
        surface = np.dot(G.T,m)
        if show:
            plt.figure()
            plt.imshow(surface.reshape(array.shape))
            plt.colorbar()
        return surface


def calc_statistics(Set, signalmask=None):
    """ Assign statistics to unaltered interferograms. NOTE: could do this for
    each stage in analysis in timeseries.py routine. signalmask is a saved .npy
    file where values of 1 makr known signal"""
    #add routine to calculate semivariogram...
    self = Set
    if signalmask:
        custommask = np.load(signalmask)
        signal = (custommask==1)
    for ig in self:
        ig.Stats = {}
        unw = load_half(ig)
        unw_ma = np.ma.masked_array(unw, np.isnan(unw))
        if signalmask:
            unw_ma[signal] = ma.masked
        ig.Stats['min'] = np.min(unw_ma)
        ig.Stats['max'] = np.max(unw_ma)
        ig.Stats['mean'] = np.mean(unw_ma)
        ig.Stats['std'] = np.std(unw_ma)
        ig.Stats['var'] = ig.Stats['std']**2


def geotrans2grid(geotrans,data):
    """ Return X, Y meshgrids for 2D modeling """
    nx = data.shape[1]
    ny = data.shape[0]

    ullon = geotrans[0]
    dlon = geotrans[1]
    ullat = geotrans[3]
    dlat = geotrans[5]
    lons = ullon + dlon*np.arange(nx)
    lats = ullat + dlat*np.arange(ny)
    Xi,Yi = np.meshgrid(lons,lats)

    return Xi,Yi



def get_grid(geo, center=True):
    """Return meshgrid of lat lon corresponding to unw
    Coordinates are center points by default, otherwise lower left corner of pixels
    Input = Geogram object
    NOTE: may not work for northern hemisphere igs
    NOTE: modelling better suited for UTM grid
    """
    xmin = float(geo.Rsc['X_FIRST'])
    nx = int(geo.Rsc['WIDTH'])
    xmax = xmin + (nx * float(geo.Rsc['X_STEP']))


    ymax = float(geo.Rsc['Y_FIRST'])
    ny = int(geo.Rsc['FILE_LENGTH']) #negative...
    ymin = ymax + (ny * float(geo.Rsc['Y_STEP']))

    if center:
        shiftx = float(geo.Rsc['X_STEP']) / 2
        shifty = float(geo.Rsc['Y_STEP']) / 2
    else:
        shiftx = 0
        shifty = 0

    lons = np.linspace(xmin+shiftx, xmax+shiftx, nx)
    lats = np.linspace(ymax+shifty, ymin+shifty, ny)

    Xi,Yi = np.meshgrid(lons,lats)

    return Xi,Yi



def image2latlon(geo,x,y):
    """ for a georeferenced array convert row,col location to lat,lon """
    x0 = float(geo.Rsc['X_FIRST'])
    y0 = float(geo.Rsc['Y_FIRST'])
    dx = float(geo.Rsc['X_STEP'])
    dy = float(geo.Rsc['Y_STEP'])
    lon = x0 + (x * dx)
    lat = y0 + (y * dy)
    return lon,lat


def latlon2image(geo, lon, lat):
    """ for a georeferenced array convert lat,lon location to row,col """
    x0 = float(geo.Rsc['X_FIRST'])
    y0 = float(geo.Rsc['Y_FIRST'])
    dx = float(geo.Rsc['X_STEP'])
    dy = float(geo.Rsc['Y_STEP'])
    x = int((lon - x0) / dx)
    y = int((lat - y0) / dy)
    return x,y


def export_latex_table_agu(Set, baseline_list='list.out', outdir='.', runTeX=False):
    """Multipage table in AGU format. baseline_list from dolist3.m output.
    Line spacing is determined by document class. runTeX=1 to automatically
    convert .tex output into pdf"""
    self = Set

    # Load baselines if not yet loaded
    if not hasattr(self,'Baselines'):
        self.load_baselines(baseline_list)
        #self.assign_baselines() #within load_baselines

    texFile = os.path.join(self.Track + '_agu.tex')
    tex = open(texFile, 'w')

    # Write Header
    tex.write(r"""%% Automatically generated LaTeX table
\documentclass[draft,grl]{AGUTeX}
\usepackage{longtable}
\begin{document}
\begin{article}
\end{article}

\begin{center}
\begin{longtable}{l cccc}
\caption{%s Interferograms}\\
\hline \hline
\bf{Dates} & \bf{Platform} & \bf{Timespan} (yr) & \bf{Baseline$_\perp$} (m) \\
\hline
\endfirsthead

\bfseries \tablename\ \thetable{} -- continued from previous page \\
\hline
\bf{Dates} & \bf{Platform} & \bf{Timespan} (yr) & \bf{Baseline$_\perp$} (m) \\
\hline
\endhead

Continued on next page...\\
\endfoot

\hline \hline
\endlastfoot

""" % self.Track)

    #Write details for each file
    #NOTE: could change date format if desired...
    for ig in self:
        tex.write('{PAIR} & {PLATFORM} & {0:.3f} & {BASELINE_PERP} \\\\ \n'.format(ig.Timespan,**ig.Rsc))

    # Write Footer
    tex.write("""
\label{tab:alligrams}
\end{longtable}
\end{center}
\end{document}
""")

    print('Wrote {} interferogram stats'.format(len(self)))
    print('Output: {0}'.format(texFile))

    if runTeX:
        os.system('latex {}'.format(texFile))
        os.system('dvipdf {}'.format(texFile.replace('tex','dvi')))


def get_stats(Interferogram, histogram=False):
	self = Interferogram
	self.amp_stats = {}
	self.phs_stats = {}
	amp,phs = load_bil(Interferogram)

	for d,data in zip([self.amp_stats, self.phs_stats],[amp,phs]):
		d['min'] = np.nanmin(data)
		d['max'] = np.nanmax(data)
		d['mean'] = np.nanmean(data)
		d['std'] = np.nanstd(data)
		d['var'] = np.nanvar(data)
		print(d)


def export_latex_table(Set, baseline_list='list.out', outdir='.', runTeX=False):
    """Interferogram table, max length 1 page"""
    self = Set

    # Load baselines if not yet loaded
    if not hasattr(self,'Baselines'):
        self.load_baselines(baseline_list)
        #self.assign_baselines() #within load_baselines

    texFile = os.path.join(self.Track + '.tex')
    tex = open(texFile, 'w')

    #Write Header
    tex.write(r"""%% Automatically generated LaTeX table
\documentclass[11pt]{amsart}
\usepackage{geometry}
\geometry{letterpaper}
\usepackage{amssymb}
\begin{document}

\begin{table}
\caption{%s Interferograms}
\centering
\begin{tabular}{l cccc}
\hline \hline
\bf{Dates} & \bf{Platform} & \bf{Timespan} (yr) & \bf{Baseline$_\perp$} (m) \\
\hline
""" % self.Track)

    #Write details for each file
    #NOTE: could change date format if desired...
    for ig in self:
        tex.write('{PAIR} & {PLATFORM} & {0:.3f} & {BASELINE_PERP} \\\\ \n'.format(ig.Timespan,**ig.Rsc))

    # Write Footer
    tex.write("""\hline
\end{tabular}
\label{tab:alligrams}
\end{table}
\end{document}
""")
    print('Wrote {} interferogram stats'.format(len(self)))
    print('Output: {0}'.format(texFile))

    if runTeX:
        os.system('latex {}'.format(texFile))
        os.system('dvipdf {}'.format(texFile.replace('tex','dvi')))


def export_latex_table_long(Set, baseline_list='list.out', outdir='.', runTeX=False):
    """Interferogram table, max length 1 page"""
    self = Set

    # Load baselines if not yet loaded
    if not hasattr(self,'Baselines'):
        self.load_baselines(baseline_list)
        #self.assign_baselines() #within load_baselines

    texFile = os.path.join(self.Track + '_long.tex')
    tex = open(texFile, 'w')

    #Write Header
    tex.write(r"""%% Automatically generated LaTeX table
\documentclass[11pt]{amsart}
\usepackage{geometry}
\geometry{letterpaper}
%%\pagestyle{empty} %%remove line numbers

\usepackage{amssymb}
\usepackage{setspace}
\usepackage{longtable}

\begin{document}
%%\thispagestyle{empty} %%remove line numbers
\onehalfspacing

\begin{center}
\begin{longtable}{l cccc}
\caption{%s Interferograms}\\
\hline \hline
\bf{Dates} & \bf{Platform} & \bf{Timespan} (yr) & \bf{Baseline$_\perp$} (m) \\
\hline
\endfirsthead

\bfseries \tablename\ \thetable{} -- continued from previous page \\
\hline
\bf{Dates} & \bf{Platform} & \bf{Timespan} (yr) & \bf{Baseline$_\perp$} (m) \\
\hline
\endhead

Continued on next page...\\
\endfoot

\hline \hline
\endlastfoot
""" % self.Track)

    #Write details for each file
    #NOTE: could change date format if desired...
    for ig in self:
        tex.write('{PAIR} & {PLATFORM} & {0:.3f} & {BASELINE_PERP} \\\\ \n'.format(ig.Timespan,**ig.Rsc))

    # Write Footer
    tex.write("""
\label{tab:alligrams}
\end{longtable}
\end{center}
\end{document}
""")
    print('Wrote {} interferogram stats'.format(len(self)))
    print('Output: {0}'.format(texFile))

    # Run system LaTeX commands to create pdf
    if runTeX:
        os.system('latex {}'.format(texFile))
        os.system('dvipdf {}'.format(texFile.replace('tex','dvi')))


def array2kmz(ig, data, transfile='aux_files/geomap_32rlks.trans',
              outname='geo_array.unw',cmap=None,vmin=-1,vmax=1,
              colorbar=False):
    """ Steps needed to create a google earth overlay from a numpy array
    specify png colorbar to attach via colorbar="""
    geopath = geocode(ig, transFile=transfile, data=data, outname=outname)
    Geo = roipy.data.Geogram(geopath)
    phs = load_half(Geo)
    save_kmz(Geo,
             data=phs,
             outname=outname[:-3] + 'kmz',
             cmap=cmap,
             vmin=vmin,
             vmax=vmax,
             colorbar=colorbar)


def extents2kml(ig):
    ''' save google earth kml polygon outlining the extents of the interferogram'''
    name = '{0}_extents.kml'.format(ig.Rsc['TRACK'])
    with open(name,'w') as kml:
        kml.write('''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>{TRACK}_extents.kml</name>
	<Style id="sn_ylw-pushpin">
		<IconStyle>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
		<LineStyle>
			<color>ffff0000</color>
			<width>3</width>
		</LineStyle>
		<PolyStyle>
			<fill>0</fill>
		</PolyStyle>
	</Style>
	<StyleMap id="msn_ylw-pushpin">
		<Pair>
			<key>normal</key>
			<styleUrl>#sn_ylw-pushpin</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#sh_ylw-pushpin</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="sh_ylw-pushpin">
		<IconStyle>
			<scale>1.3</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
		<LineStyle>
			<color>ffff0000</color>
			<width>3</width>
		</LineStyle>
		<PolyStyle>
			<fill>0</fill>
		</PolyStyle>
	</Style>
	<Placemark>
		<name>{TRACK}</name>
		<styleUrl>#msn_ylw-pushpin</styleUrl>
		<Polygon>
			<tessellate>1</tessellate>
			<outerBoundaryIs>
				<LinearRing>
					<coordinates>
						{LON_REF1},{LAT_REF1},0 {LON_REF2},{LAT_REF2},0 {LON_REF4},{LAT_REF4},0 {LON_REF3},{LAT_REF3},0 {LON_REF1},{LAT_REF1},0
					</coordinates>
				</LinearRing>
			</outerBoundaryIs>
		</Polygon>
	</Placemark>
</Document>
</kml>
'''.format(**ig.Rsc)
        )
    print('wrote:', name)


def geocode(Interferogram, transFile, data=None, outname=None, kml=False):
    """call geocode.pl, which makes rect_lookup & rsc
    NOTE: outname must end in .unw"""
    print('running geocode.pl')
    self = Interferogram
    inname = self.Name

    #NOTE: output to cwd, get rid of 32rlks, append geo prefix
    if outname == None:
        outname = 'geo_' + self.Name.replace('_32rlks','')


    if type(data) is np.ndarray:
        amp = load_half(self,1)
        inname = save_bil(self, 'tmp.unw', amp, data)

    cmd = 'geocode.pl {0} {1} {2}'.format(transFile, inname, outname)
    print(cmd)
    os.system(cmd)

    #NOTE: not working b/c perl script has multiple returns throughout script
    #NOTE: could try subprocess instead of os.system()...
    # eg whether or not keywords exist
    #sleep(10)
    #geogram = roipy.data.Geogram(outname)
    #if kml:
    #    unw2png(geogram)
    #return geogram

    return outname


def unw2png(Geogram, outname=None):
    """call unw2png.pl, which makes .png and .kml files for overlay"""
    self=Geogram
    if outname==None:
        outname = self.Name + '.png'
    cmd = 'unw2png.pl {0} {1}'.format(self.Name, outname)
    print(cmd)
    os.system(cmd)


def georef_timeseries(Timeseries):
    """Save reconstructed deformation for each date in the timeseries
    """
    print('Work in progress')
    os.mkdir('tmp_geotimeseries')


def georef_timeseries_smooth(Timeseries):
    """ Instread"""
    print('Work in progress')

def animation_kml_timeseries(Timeseries):
    """save a kml with output from georef_timeseries """
    print('work in progress')


def georeference(Interferogram, trans, rsc, data=None, outname='rect_result.unw',
                 dtype='f4', nanval=np.nan, nlooks=1):
    """
    NOTE: use geocode instead!
    call ROI_PAC rect_lookup.pl command to georeference a file or array.
    NOTE: outname must start with 'rect_' !!!

    Keyword arguments:

        .. tabularcolumns:: |l|L|

        =============   ====================================================
        Keyword         Description
        =============   ====================================================
        Interferogram   :class:`~roipy.data.Interferogram` instance
        amp             (np.array) amplitude values
        phs             (np.array) phase or deformation values
        transfile       (str) path to .trans file for georeferencing
        =============   ====================================================


    Outputs:
    geo_<>.unw,  geo_<>.unw.rsc

    Examples:

    #. view an interferogram in Google Earth:

    >>> ig = rp.data.Interferogram(path)
    >>> outpath = rp.tools.georeference(ig, './geomap.trans')
    >>> ig_g = rp.data.Interferogram(outpath)
    >>> rp.tools.save_kmz(ig_g)

    #. georeference a stack

    >>> ts = rp.timeseries.Timeseries(set89, 'all')
    >>> stack = rp.plot.stack(ts)
    >>> ig = set89[0]
    >>> rp.tools.georeference(ig, 'geomap.trans', 'geo.rsc', stack')
    """
    #NOTE: can just call geocodem.pl and avoid copying the geo.rsc file at the end

    self = Interferogram
    #save numpy array as unw with rsc
    #if (type(amp) is np.ndarray) and (type(phs) is np.ndarray):
    #    path = save_bil(self, outname, amp, phs, dtype, nanval)
    #    name = outname
    if type(data) is np.ndarray:
        amp = load_half(self,1)
        path = save_bil(self, outname, amp, data, dtype, nanval)
        #NOTE: save_bil must finish before calling rect_lookup
        #sleep(10) Doesn't do it...
        name = outname
    else:
        path = self.Path #if georeferencing any old interferogram
        name = self.Name

    # NOTE: Only needed for wrapped interferograms
    # NOTE: should be able to
    #Get metadata for transfile & write rect_lookup.in file
    transRSC = load_rsc(trans + '.rsc')
    transRSC['TRANS'] = trans
    transRSC['NUMCOL'] = int(transRSC['XMAX']) - int(transRSC['XMIN']) + 1
    transRSC['NUMROW'] = int(transRSC['YMAX']) - int(transRSC['YMIN']) + 1
    transRSC['INPUT'] = path
    transRSC['INPUTWIDTH'] = self.Width
    transRSC['INPUTLENGTH'] = self.Length
    transRSC['OUTPUT'] = name.replace('rect_','geo_')
    rectlkup = open('rect_lookup.in','w')
    rectlkup.write(
    """
    Longitude-Latitude Lookup File Name    (-)   = {TRANS}			! Lookup file
    Lookup File Dimensions                 (-)   = {WIDTH} {FILE_LENGTH}	! across ,  down
    Input Image File Name                  (-)   = {INPUT}			! file to be resampled
    Input File Dimensions                  (-)   = {INPUTWIDTH} {INPUTLENGTH}     	! across , down
    Output Image File Name                 (-)   = {OUTPUT}			! output filename
    Lookup File Start Sample for Output    (-)   = {XMIN}      		! pixel across
    Number of Samples for Output           (-)   = {NUMCOL}     		! number of pixels
    Lookup File Start Line for Output      (-)   = {YMIN}    		! start line
    Number of Lines for Output             (-)   = {NUMROW}      		! number of lines
    Skip Every N Lookup points Across-Down (-)   = 1 1        			! across , down
    File Type                              (-)   = RMG         			! [RMG,CPX]
    Interpolation Method                   (-)   = Bilinear    			! [Bilinear,Sinc,NN]
    """.format(**transRSC))
    rectlkup.close()

    # georeference & create .rsc file for georefed array
    #cmd = ['rect_lookup', 'rect_lookup.in']
    #subprocess.call(cmd)
    cmd = 'rect_lookup rect_lookup.in'
    success = os.system(cmd) #automatically waits for exit status
    #if success:
    #print 'success' NOTE: Doesn't seem to work for the perl scripts...
    #self.Georeferenced = True
    # NOTE: see geocode.pl for creation of geo_*.rsc
    #cmd = 'cp {self.Rsc} {OUTPUT}.rsc'.format(transRSC)
    #os.system(cmd4)
    shutil.copyfile(rsc, transRSC['OUTPUT'] + '.rsc')

    # downsample result if requested
    if nlooks != 1:
       transRSC['OUTPUT'] = lookdown(transRSC['OUTPUT'], nlooks)

    # Return georeferenced interferogram object
    geo = roipy.data.Geogram(transRSC['OUTPUT'])
    print(os.path.abspath(transRSC['OUTPUT']))
    return geo


def lookdown(path, nlooks=2):
    """ downsample to hald original resolution """
    cmd = 'look.pl {0} {1}'.format(path, nlooks)
    success = os.system(cmd) #return code doesn't seem to work...
    #if success:
    found = re.findall('_\d+rlks.unw',path)
    if found:
        suf = found[0]
        num = suf.strip('_rlks.unw')
        newnum = str(int(num)*nlooks)
        newsuf = suf.replace(num, newnum)
        outname = path.replace(suf, newsuf)
    else:
        outname = path[:-4] + '_{}rlks.unw'.format(nlooks)

    return outname


def get_cart2los(geoincidenceFile):
    """ Convert cartesian arrays into line-of-sight """
    # NOTE: careful of sign changes for ascending & descending tracks!
    # NOTE: assumes standarard WGS84lat/lon format of file
    geo = roipy.data.Geogram(geoincidenceFile)
    look,head = load_bil(geo)

    look = np.deg2rad(look)
    head = np.deg2rad(head)

    # Convention is negative --> reduction in line of sight, therefore uplift
    EW2los = np.sin(head) * np.sin(look)
    NS2los = np.cos(head) * np.sin(look)
    Z2los = -np.cos(look)

    # Change convention to uplift --> positive
    cart2los = -np.dstack([EW2los, NS2los, Z2los])

    return cart2los


def match_date(Set, date):
    """Return a list of SERIAL dates that pair with a given date """
    self = Set
    row, col = np.where(self.PairsSerial == date)
    matches = self.PairsSerial[row, ~col]
    return matches
    #print '\n'.join(iglist)


def match_igrams(Set, date):
    """Return an array of SERIAL PAIRS interferograms containing a certain date (integer)"""
    self = Set
    row, col = np.where(self.PairsSerial == date)
    matches = self.PairsSerial[row]
    return matches


def crop_grds(grd1, grd2, bounds=None):
    """Crop grd2 to grd1. take subsection defined by boiunds if given """
    out1 = grd1[:-4] + '_crop.grd'
    if bounds:
        cmd = 'grdsample {0} -G{1} -R{2}/{3}/{4}/{5}r'.format(grd1,out1,*bounds)
        print(cmd); os.system(cmd)
    else:
        cmd = 'grdsample {0} -G{1}'.format(grd1,out1)
        print(cmd); os.system(cmd)

    out2 = grd2[:-4] + '_crop.grd'
    cmd = 'grdsample {0} -G{1} -R{2}'.format(grd2,out2,out1)
    print(cmd); os.system(cmd)

    return out1, out2
