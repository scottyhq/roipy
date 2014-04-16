.. _installing:

**********
Installing
**********

Dependencies
============

**Requirements**


matplotlib 1.0.0 (or later, `download <http://sf.net/projects/matplotlib/>`__)

python 2.4 (or later but not python3)
    matplotlib requires python 2.4 or later (`download <http://www.python.org/download/>`__)

numpy 1.2.1 (or later)
    array support for python (`download <http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103>`__)

**Optional libraries**

h5py 2.0.0
    Python interface to HDF5 library (`download <http://code.google.com/p/h5py/>`__)
    only needed for matlab file I/O e.g. :func:`~roipy.tools.load_mat`

gdal 1.9
    geospatial data abstraction library. (`download <http://trac.osgeo.org/gdal/wiki/DownloadSource/>`__)
    NOTE: you must build from source to have python bindings.


Installation
============
1) Download source code `here <http://www.geo.cornell.edu/eas/gstudent/sth54/pubs/roipy/>`__. 

2) Make sure your python version and dependencies are up to date

3) Add the roipy directory to your PYTHONPATH environment varible::
	
	export PYTHONPATH=/home/software/roipy

4) from the ipython prompt type ``import roipy as rp`` and start using roipy




