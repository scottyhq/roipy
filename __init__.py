"""
Object-oriented package for post-processing ROI_PAC data in pylab.

**Example Usage:**

>>> import roipy as rp
>>> path89 = 'examples/t6089/stack'
>>> set89 = rp.data.Set(path89)
>>> rp.plot.browseSet(set89)
 
Tab-complete to see available methods & attributes!
 
Contains modules for organization of data structures and common tasks:

  :mod:`roipy.data` 
    defines the :class:`~roipy.data.Interferogram` class and the :class:`~roipy.data.Set` class. 
    These are objects with useful attributes such as an attached RSC dictionary.
  
  :mod:`roipy.tools` 
    load & save numpy arrays to various formats, do geocoding stuff with gdal, etc.
  
  :mod:`roipy.plot` 
    convenient matplotlib plots of numpy arrays & basemap maps
  
  :mod:`roipy.timeseries`
    methods for calculating SBAS timeseries on a :class:`~roipy.data.Set` of interferograms

  :mod:`roipy.noise`
    methods for adding & annalyzing different types of noise

  :mod:`roipy.models`
    methods for analytic forward model synthetics

test making a table
    
    .. tabularcolumns:: |l|L|
    
    ================        ====================================================
    Keyword                 Description      
    ================        ==================================================== 
    Interferogram           rp.data.Interferogram
    amp                     (np.array) amplitude values
    phs                     (np.array) phase or deformation values
    transfile               (str) path to .trans file for georeferencing
    ================        ==================================================== 

"""

__all__ = ['data','tools','plot','timeseries','noise','models']
__version__ = 1.0
__author__ = 'sth54@cornell.edu'
__lastupdate = '1/12/12'

#Python 2
#import data, tools, plot, timeseries, noise, models

# Python 3
from roipy import data, tools, plot, timeseries, noise, models

