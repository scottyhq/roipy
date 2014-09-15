#roipy

python tools for interactive post-processing of roi_pac files. 


### Warning
These scripts are not the cleanest. It's a collection of utilities I've been creating and using over the years to work with ROI_PAC data, so you're on your own to figure out cryptic errors etc. 

Nevertheless, there is some useful stuff here, but you'll probably only find it useful if you're in our lab group and have access to other matlab scripts etc.

### InSAR timeseries
To get the necessary scripts installed on the computer run:

```
git clone https://github.com/scottyhq/roipy.git
```

Make sure the parent folder you're installing to is in the `$PYTHONPATH` search path. For example, assuming you now have `/home/me/github/roipy`:

```
export PYTHONPATH=$PYTHONPATH:/home/me/github
```

Open an iPython terminal (just type `ipython`) and run the following commands to set-up the matlab codes for running an InSAR timeseries. Note that this assumes you've run `stackum.pl` and have a directory of rectified interferograms that follow a standard naming scheme (e.g. a bunch of `rect_070827-061016_32rlks.unw` files in a directory called `stack`):

```python
import roipy as rp
set = rp.data.Set('stack')
ts = rp.timeseries.Timeseries(set,'timeseries')
ts.prep_matlab()
```

You should now have a folder called 'timeseries' which contains a file called `defaults.m` that you need to run timeseries scripts in Matlab. You might need to edit paths in defaults.m.


There are a bunch of plotting routings in `rp.plot` and other utilities in `rp.tools`. Feel free to explore, but you might run into deadends!
