*************
configuration
*************

:mod:`roipy`
====================

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

new testing link :func:`~roipy.tools.load_bil`
