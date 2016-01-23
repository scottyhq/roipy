"""
Interferogram and Set classes for working with ROI_PAC output:
attributes in UpperCase
methods in lower_case
"""

import sys
import numpy as np
import glob
import datetime as dt
import os.path
import re



class Interferogram():
    """General interferogram object (from any platform: ERS, Envisat, ALOS)
    Designed for BIL/RMG files output by ROI_PAC that have associated rsc's
    
    Example
    -------
    >>> import roi_py as rp
    >>> ig = rp.data.Interferogram(path/to/interferogram)
    >>> print(ig)
    """

    def __init__(self, filepath):
        """Assign attributes & call setup method()"""
        self.Path = filepath
        self.DataDir, self.Name = os.path.split(filepath)
        self.Rsc = self.load_rsc()
        self.Rsc['PLATFORM'] = 'unspecified'
        self.Georeferenced = False
        self.Width = int(self.Rsc['WIDTH'])
        self.Length = int(self.Rsc['FILE_LENGTH'])
        self.Shape = (self.Length, self.Width)
        self.Size = self.Width * self.Length
        self.Wavelength = float(self.Rsc['WAVELENGTH'])
        self.Phs2cm = self.Wavelength / (4*np.pi) * 100

        self.ProcDir = os.path.join(os.path.dirname(self.DataDir), 'clean_stack')
        if not os.path.isdir(self.ProcDir):
            os.mkdir(self.ProcDir)
        self.setup()
        
        

    def __str__(self):
        message ='''
Path: {PATH}
Track: {TRACK}
Platform: {PLATFORM}
Date1: {DATE1}
Date2: {DATE2}
Timespan: {TIME_SPAN_YEAR}
Length: {FILE_LENGTH}
Width: {WIDTH}
Dtype: {DATA_TYPE}'''.format(**self.Rsc)

        return message

    def assign_dtype(self):
        roi_dtypes = dict(unw='float32',
                          cor='float32',
                          hgt='float32',
                          msk='float32',
                          rmg='float32',
                          r4='float32',
                          dem='int16',
                          amp='complex64',
                          int='complex64',
                          slc='complex64')
        suffix = self.Path[-3:]
        if suffix in roi_dtypes:
            self.DataType = roi_dtypes[suffix]
        else:
            self.DataType = 'Unknown'
        self.Rsc['DATA_TYPE'] = self.DataType


    def setup(self):
        #self.associate_files() # Not really useful?
        self.add2rsc()
        self.assign_dtype()
        try:
            self.pix2km()
        except:
            print('pix2km info not in rsc file...')
        
    def pix2km(self, inc_mean=None):
        """ ESTIMATE dimensions of data in km rather than radar coordinates """
        #incidence = {'t6089':41.1,
        #             't2282':22.8,
        #             't2010':20.2,
        #             't2003':20.0,
        #             't2318':23.1,
        #             'unspecified':22.8}
        #print(self.Name)
        #inc_mean = incidence[self.Rsc['TRACK']]
        
        # NOTE: could make more accurate version of this that uses incidence file...
        # NOTE: mean value from incidence file much more accurate, this works
        # for envisat, but 'LOOK_REF' keys not in ERS2 metadata...
        if inc_mean==None:
            try:
                corners = np.array([self.Rsc['LOOK_REF1'],
                           self.Rsc['LOOK_REF2'],
                           self.Rsc['LOOK_REF3'],
                           self.Rsc['LOOK_REF4']], dtype='float')
                #look_mean = (corners.min() + corners.max()) / 2.0
                look_mean = corners.mean()
                inc_mean = look_mean + 4.0 #estimate... 3-5 based on experience
            except:
                print('Unable to estimate mean incidence angle from LOOK_REF entries in .rsc')

        slant2groundx = 1e-3 * (1 / np.sin(np.deg2rad(inc_mean)))
        slant2groundy = 1e-3 * float(self.Rsc['EARTH_RADIUS']) / (float(self.Rsc['EARTH_RADIUS']) + float(self.Rsc['HEIGHT']))
        scalex = float(self.Rsc['RANGE_PIXEL_SIZE']) * slant2groundx
        scaley = float(self.Rsc['AZIMUTH_PIXEL_SIZE']) * slant2groundy
        
        
        self.Rsc['RANGE_GROUND_PIXEL_SIZE'] = scalex #km per pixel in x
        self.Rsc['AZI_GROUND_PIXEL_SIZE'] = scaley #km per pixel in y
        self.Rsc['MEAN_INCIDENCE'] = inc_mean
        self.col2km = scalex
        self.row2km = scaley
        self.radius2km = np.hypot(scalex,scaley)
        

    def add2rsc(self):
        """Add custom header lines to rsc metadata stored in memory"""
        try:
            #regular expression search for 't' followed by 4 numbers
            self.Rsc['TRACK'] = re.search("t\d{4}",self.Path).group()
        except:
            self.Rsc['TRACK'] = 'unspecified'
        
        try:
            self.get_dates() #conversion from string to serial dates
        except:
            print('Warning: Unable to extract date information from rsc')
        
        if 'TIME_SPAN_YEAR' not in self.Rsc:
            self.fix_timespan()
        self.Rsc['PATH'] = self.Path


    def associate_files(self):
        """Add dictionary of associate ROI_PAC files to the Interferogram
        object"""
        #NOTE: maybe not that useful...
        self.Files = {}
        self.Files['rsc'] = self.Path + '.rsc'
        self.Files['msk'] = self.Path[:-4] + '.msk'
        self.Files['inc'] = self.Path.replace(self.Rsc['DATE12'],'incidence')
        self.Files['hgt'] = os.path.join(self.DataDir, 'radar_%srlks.hgt' % self.Rsc['RLOOKS'])     

        
    def get_dates(self):
        """Extract dates from the metadata and convert them to serial dates"""
        #self.Date1 = dt.datetime.strptime(self.Rsc['FIRST_FRAME_SCENE_CENTER_TIME'],'%Y%m%d%H%M%S%f')                    
        #NOTE: should store DATE1 & DATE2 in .rsc file
        date1 = self.Rsc['DATE12'].split('-')[0] #most recent date listed first
        if date1.startswith('9'):
            self.Date1 = dt.datetime.strptime('19' + date1, '%Y%m%d')
        else:
            self.Date1 = dt.datetime.strptime('20' + date1, '%Y%m%d') 
        self.Rsc['DATE1'] = self.Date1.strftime('%Y%m%d')
        self.Rsc['DATE1_SERIAL'] = str(self.Date1.toordinal())

        date2 = self.Rsc['DATE12'].split('-')[1]
        #Convert from 2-digit years to 4 digit years
        if date2.startswith('9'):
            self.Date2 = dt.datetime.strptime('19' + date2, '%Y%m%d')
        else:
            self.Date2 = dt.datetime.strptime('20' + date2, '%Y%m%d') 
        self.Rsc['DATE2'] = self.Date2.strftime('%Y%m%d')
        self.Rsc['DATE2_SERIAL'] = str(self.Date2.toordinal())

        self.Rsc['PAIR'] = self.Rsc['DATE1'] + ' ' + self.Rsc['DATE2']
        self.Rsc['PAIR_SERIAL'] = self.Rsc['DATE1_SERIAL'] + ' ' + self.Rsc['DATE2_SERIAL']
        
        #note: must come after add2rsc() b/c some rsc files missing timespan:
        self.Timespan = abs(float(self.Rsc['TIME_SPAN_YEAR']))
        self.ID = self.Rsc['DATE1'] + ' ' + self.Rsc['DATE2']
        

    def load_rsc(self):
        """Read from unw.rsc file & store metadata as an attribute dictionary
        named rsc"""
        metadata = {}
        rscFile = open(self.Path + '.rsc', 'r')
        # remove blank lines, retaining only lines with key, value pairs
        #if line.strip() gets rid of blank line cases
        allLines = [line for line in rscFile.readlines() if line.strip()]
        for line in allLines:
            try:
                var, value = line.split()
            except:
                print('Unable to import from .rsc:',line)
            metadata[var] = value
        rscFile.close()

        return metadata


    def load_bil(self, dims=None):
        """Load amplitude and phase array from BIL files (eg .unw) output by ROI_PAC

        Input
        -----
        self = roi_py.data.Interferogram instance
        dims = (width,length) tuple to override .rsc dimensions
        geo = True if .unw is georeferenced to maintain correct array orientation

        Output
        ------
        (leftArray, rightArray) = 2D numpy arrays

        Examples
        ------
        look,head = incidence file
        amplitude,phase = .unw file
        slope,elevation = .hgt file
        """
        #Convenience method copied from roi_py.tools
        if dims:
            width = dims[0]
            length = dims[1]
        else:
            width = int(self.Rsc['WIDTH'])
            length = int(self.Rsc['FILE_LENGTH'])

        fullwidth = width*2
        data = np.fromfile(self.Path, dtype=self.DataType, count=length*fullwidth)
        data = np.reshape(data, (length,fullwidth))  #NOTE: order=F(fortran-like) fill-by column not row! as matlab does
        data[data==0] = np.NaN # Assume nan=0
        leftArray = np.ascontiguousarray(data[:,0:width])
        rightArray = np.ascontiguousarray(data[:,width:fullwidth])

        # Orient to N up, W left
        if self.Georeferenced:
            leftArray = np.rot90(leftArray.transpose(),1)
            rightArray = np.rot90(rightArray.transpose(),1)
        else:
            leftArray = np.rot90(leftArray,2)#N up, W left 
            rightArray = np.rot90(rightArray,2)

        
        return leftArray, rightArray


    def fix_timespan(self):        
        """Write the missing 'TIME_SPAN_YEAR' keyword to the .rsc file, get the
        timespan from extracting the second date from the DATE12 keyword"""
        deltaT = self.Date2 - self.Date1
        self.Rsc['TIME_SPAN_YEAR'] = str(deltaT.days / 365.242)
        print('adding TIME_SPAN_YEAR to %s.rsc' % self.Name)
        
        rsc = open(self.Path + '.rsc', 'a')
        rsc.write('TIME_SPAN_YEAR       %s\n' % self.Rsc['TIME_SPAN_YEAR'])
        rsc.close()



class Geogram(Interferogram):
    """ Class for georectified files, for now only difference Georeferenced=True
    NOTE: just keep this the same as Interferogram & add flag for georeferenced...
    """
    def __init__(self, filepath):
        self.Georeferenced = True
        self.Path = filepath
        self.DataDir, self.Name = os.path.split(filepath)
        #self.DataType = np.dtype('<f4')
        self.Rsc = self.load_rsc()
        self.Rsc['PLATFORM'] = 'unspecified'
        self.Width = int(self.Rsc['WIDTH'])
        self.Length = int(self.Rsc['FILE_LENGTH'])
        self.Shape = (self.Length, self.Width)
        self.Wavelength = float(self.Rsc['WAVELENGTH'])
        self.Phs2cm = self.Wavelength / (4*np.pi) * 100

        self.ProcDir = os.path.join(os.path.dirname(self.DataDir), 'clean_stack')
        if not os.path.isdir(self.ProcDir):
            os.mkdir(self.ProcDir)
        self.setup()
        
    
    def setup(self):
        #self.associate_files()
        self.assign_dtype()
        self.add2rsc()
        #self.Rsc['PLATFORM'] = 'unspecified'
        self.get_geotrans()
        

    def get_geotrans(self):
        '''(upper left x, w-e pixel res, rotation, top left y, rotation, n-s pixel resolution)'''
        try:
            #self = roipy.data.Geogram(geofile)
            ulx = float(self.Rsc['X_FIRST'])
            dx = float(self.Rsc['X_STEP'])
            xrot = 0.0
            uly = float(self.Rsc['Y_FIRST'])
            yrot = 0.0
            dy = float(self.Rsc['Y_STEP'])
            self.Geotrans = (ulx, dx, xrot, uly, yrot, dy)
        except:
            print('WARING: Unable to load geotransform information!')
            raise

# ---------------------------------------------------------------------------- #        
class Set():
    """Class for representing a set of rectified interferograms.

    Example                                                                                                                                                                                                                                                                                                   
    --------
    >>> track227 = Set(<path/to/rectified/227/interferograms>) 
    >>> print(track227)
    """    
    
    
    def __init__(self, input, pattern='rect*unw'):
        """ input can either be a path or list of file strings"""
        if type(input) is str:
            self.DataDir    = input.rstrip(os.path.sep) #no trailing slash
            self.Igrams     = self.load_directory(pattern=pattern)
        else:
            self.DataDir = os.path.dirname(os.path.abspath(input[0]))
            #print(self.DataDir)
            self.Igrams = self.load_paths(input)
        self.Nig        = len(self.Igrams)
        self.Width      = self.Igrams.values()[0].Width
        self.Length     = self.Igrams.values()[0].Length
        self.Shape       = (self.Length, self.Width)
        self.Omissions  = {}
        self.Tandems    = []
        try:
            self.Track = re.search("t\d{4}",self.DataDir).group()
            #print('got track')
        except:
            self.Track = 'unspecified'
        self.setup()
        
        
    def __str__(self):
        message ='''
DataDir: {DataDir}
Track: {Track}
Orbit: {Orbit}
Interferograms: {Nig}
Dates: {Ndate}
Length: {Length}
Width: {Width}
'''.format(**self.__dict__)
        return message
    
    
    def __getitem__(self,i):
        """ PairsString are Igrams dictionary keys in chronological order"""
        return self.Igrams[self.PairsString[i]]
    
    
    def __len__(self):
        return len(self.Igrams)


    def __iter__(self):
        self.index = 0
        return self


    def next(self):
        """Iterate through interferograms in chronological order """
        if self.index < len(self.Igrams):
            item = self.Igrams[self.PairsString[self.index]]
            self.index += 1
            return item
        else:
            raise StopIteration


    def setup(self):
        """ Associate important timeseries arrays with the interferogram
            ORDER MATTERS... eg self.Nig reset before get_interferograms()"""
        #self.associate_files() #redo timeseries w/ python
        self.Nig    = len(self.Igrams)
        self.Pairs, self.PairsSerial    = self.get_interferograms()
        self.Dates, self.DatesSerial    = self.get_dates()
        self.TimeIntervals              = self.get_time_intervals()
        self.TimeIndex, self.TimeIndexString    = self.get_time_index()
        self.ChronList                  = self.get_chronlist()
        self.Timespan                   = (self.DatesSerial[-1] - self.DatesSerial[0]) / 365.24

        # Reset parameters in case some interferograms are omitted...
        self.Ndate  = len(self.Dates)
        self.Names  = [ig.Name for ig in self]
        try:
            self.Orbit = self[0].Rsc['ORBIT_DIRECTION']
        except:
            self.Orbit = 'undefined'
        

    def load_directory(self, pattern):
        """Create a dictionary of Interferogram.General() instances from all rect*unw
        files in the specificed path
        
        Usage
        -----
        igrams = Interferogram.load_set(os.getcwd())"""
        igrams = {} 
        paths = glob.glob(os.path.join(self.DataDir, pattern))
        for path in paths:
            ig = Interferogram(path)
            igrams[ig.Rsc['PAIR']] = ig

        return igrams


    def load_paths(self, paths):
        """Instead of globbing directory, load only specified files"""
        igrams = {}
        for path in paths:
            ig = Interferogram(path)
            igrams[ig.Rsc['PAIR']] = ig

        return igrams

    
    def get_dates(self):
        """Assign lists of metadata values as instance attributes"""
        recentDates = [int(ig.Rsc['DATE1']) for ig in self]
        laterDates = [int(ig.Rsc['DATE2']) for ig in self]
        Dates = np.unique(np.array(recentDates + laterDates,dtype='int32'))
        # NOTE: could use datetime module to convert self.DATES list or:
        recentDates = [int(ig.Rsc['DATE1_SERIAL']) for ig in self]
        laterDates = [int(ig.Rsc['DATE2_SERIAL']) for ig in self]
        DatesSerial = np.unique(np.array(recentDates + laterDates,dtype='int32'))

        return Dates, DatesSerial


    def get_time_intervals(self):
        """Subtract serial dates to get day intervals"""
        TimeIntervals = np.diff(self.DatesSerial)
        #self.TimeIntervals = [str(dt) for dt in dtList]
        
        return TimeIntervals
    
    ''' 
    def Array2ListoStrings(self):
        """ Just like it sounds """
        #TODO
        strlist = []
    '''
    
    def get_interferograms(self):
        """Assign list of date pairs as Set() attributes in both YYMMDD format
        and serial format"""
        #NOTE Doesn't Work
        #pairs = np.zeros((self.Nig,2),dtype='int32')
        #pairs[:,0] = [int(ig.Rsc['DATE1']) for ig in self]
        #pairs[:,1] = [int(ig.Rsc['DATE2']) for ig in self]
        #tricky way to sort entire matrix by first column
        #Pairs = pairs[pairs[:,0].argsort(),:] #Problem... does not do secondary sort of second column
        
        #NOTE complicated
        #dates1 = [int(ig.Rsc['DATE1']) for ig in self]
        #dates2 = [int(ig.Rsc['DATE2']) for ig in self]
        #sorter = np.array([str(d1) + str(d2) for d1, d2 in zip(dates1,dates2)])
        #pairs = np.array([dates1,dates2]).T #transpose
        #inds = np.argsort(sorter) #sort list by another list
        #pairs_sorted = np.take(pairs, inds) #only1D array
        #Pairs = pairs[inds] 

        #NOTE Best, use structured array
        pairs = np.zeros(self.Nig, dtype=[('date1','int'),('date2','int')])
        pairs['date1'] = [int(ig.Rsc['DATE1']) for ig in self.Igrams.values()]
        pairs['date2'] = [int(ig.Rsc['DATE2']) for ig in self.Igrams.values()]
        sorted = np.sort(pairs, order=['date1','date2']) #first date1, then date2
        #NOTE: could rewrite other scripts to access as Pairs['date1']
        Pairs = np.array([sorted['date1'],sorted['date2']]).T
        self.PairsString = [str(c1) + ' ' + str(c2) for c1,c2 in Pairs]
    
        pairsSer = np.zeros(self.Nig, dtype=[('date1','int'),('date2','int')])
        pairsSer['date1'] = [int(ig.Rsc['DATE1_SERIAL']) for ig in self.Igrams.values()]
        pairsSer['date2'] = [int(ig.Rsc['DATE2_SERIAL']) for ig in self.Igrams.values()]
        sorted = np.sort(pairsSer, order=['date1','date2']) #first date1, then date2
        PairsSerial = np.array([sorted['date1'],sorted['date2']]).T
        self.PairsSerialString = [str(c1) + ' ' + str(c2) for c1,c2 in PairsSerial]

        return Pairs, PairsSerial 
        
        

    def get_time_index(self):
        """Generate the time_index assigning an integer to each SAR acquistion
        in chronological order. Used to construct B matrix for SBAS inversion"""
        TimeIndex = np.zeros((self.Nig,2),dtype='int32')
        TimeIndexString = ''
        chronDates = list(self.DatesSerial)
        dates1 = self.PairsSerial[:,0]
        dates2 = self.PairsSerial[:,1]
        for i,d1,d2 in zip(np.arange(self.Nig),dates1,dates2):
            ind1 = chronDates.index(d1)
            ind2 = chronDates.index(d2)
            TimeIndex[i,:] = [ind1, ind2]
            TimeIndexString += '{0} {1}\n'.format(ind1+1,ind2+1) #for prep_matlab

        return TimeIndex, TimeIndexString


    #NOTE: probably can remove this method
    def get_chronlist(self):
        """Get a chronologically sorted list of IG files"""
        ChronList = [] #for matlab...
        for pair in self.PairsString:
            ChronList.append("'" + self.Igrams[pair].Name + "'")
        return ChronList


    def load_baselines(self, path='list.out'):
        """ Read perpendicular baselines from list.out from dolist.m output,
        store as a dictionary {'date': Bperp}"""
        #vectorized numpy approach:
        record = np.loadtxt(path, dtype={'names':('date','bperp'),'formats':('S6','i4')})
        pre2000 = np.char.startswith(record['date'], '9')
        #NOTE: numpy doesn't like to change #of characters for each element in array
        longDates = np.char.add('20', record['date'])
        longDates[pre2000] = np.char.add('19', record['date'][pre2000])
        self.Baselines = dict(zip(longDates,record['bperp']))
        
        #self.assign_baselines()
    
    
    def assign_baselines_new():
        """Get baseline values from rscs """
        print('TODO')
    
    def assign_baselines(self):
        """ Use self.Baselines dictionary to assign perp baseline as ig attribute"""
        #NOTE: could also add to .rsc file...
        for ig in self:
            date1 = ig.Rsc['DATE1']
            date2 = ig.Rsc['DATE2']
            ig.Bperp = abs(self.Baselines[date1] - self.Baselines[date2])
            ig.Rsc['BASELINE_PERP'] = ig.Bperp
             


    def merge_tandems(self):
        """
        NOTE: not sure this is necessary
        Combine SAR acquisitions that differ by only 1 day into a single date
        ASSUMPTION: one day's deformation is NOT resolvable"""
        dtList = [int(dt) for dt in self.TimeIntervals]
        indOnes = [i for i,dt in enumerate(dtList) if dt==1]
        dates = [self.Dates[i] for i in indOnes]
        dates_ser = [self.DatesSerial[i] for i in indOnes]
        
        for d1,d2 in zip(dates,dates_ser):
            # Works for lists not numy arrays
            #self.Dates.remove(d1)
            #self.DatesSerial.remove(d2)
            
            # Works for numpy arrays
            self.Dates = np.delete(self.Dates, self.Dates==d1)
            self.DatesSerial = np.delete(self.DatesSerial, self.DatesSerial==d2)
        
        self.get_time_intervals()
        self.get_time_index()
        self.Tandems = dates
        print('Removed %s from timeseries' % self.Tandems)


        
    def query(self, attr=None, val=None, outattr=None, IG=None):
        """Search metadata for all interferograms and return a list of results

        Examples
        --------
        >>>query(IG=#)                          return metadata for the specified IG#
        >>>query(<metadata field>)              return value of specified metadata field for all IGs
        >>>query(<metadata field>, <value>)     return list of IGs sharing common metadata attribute
        >>>query('PLATFORM','ERS1')             return list of ERS1 interferograms
        >>>query('PLATFORM','ERS2','DATE12')    return list of ordered DATE12 values for ERS2 IGs
        >>>query(IG=23)                         return metadata for interferogram #23 in the timeseries
        """
        if attr and val:
            LIST = [ig.Rsc['PAIR'] for ig in self if ig.Rsc[attr] == val]
            LIST.sort()
            if outattr:
                outLIST = [self.Igrams[x].Rsc[outattr] for x in LIST]
                inds = np.argsort(LIST) #sort list by another list
                LIST = list(np.take(outLIST,inds))
            #print('\n'.join(LIST))
            return LIST

        if attr and not val:
            LIST = [ig.Rsc['PAIR'] + '  ' + ig.Rsc[attr] for ig in self]
            LIST.sort()
            return LIST
        
        if IG:
            if type(IG) is int:
                IG = self.Pairs[IG]
            metadata = vars(self.Igrams[IG])
            if outattr:
                metadata = metadata[outattr]
            return metadata     



    def match_date(self, date):
        """Return a list of dates that pair with a given date """
        row, col = np.where(self.Pairs == date)
        matches = self.Pairs[row,~col]
        return matches


    def match_igrams(Set,date):
        """Return an array of all interferograms containing a certain date (integer)"""
        self = Set
        row, col = np.where(self.Pairs == date)
        matches = self.Pairs[row]
        return matches

    
    def omit(self, date=False, IG=False, range=False):
        """Omit interferograms from set. Can do one of the following:
        1) omit all interferograms containing a certain date
        2) omit specific listed interferograms
        3) omit all interferograms OUTSIDE a specific date range

        IMPORTANT! Run this method before running prep_matlab()

        Input
        -----
        date = list of dates in quotes, eg '20010812'
        IG = list of interferogram data identifiers, eg ['20010812-20100723']
        range = tuple of start & end dates, eg ('20010812', '20100723)'"""

        iglist = []
        if IG:
            iglist = IG #list of interferograms
        
        elif date:
            iglist1 = self.query('DATE1',date)
            iglist2 = self.query('DATE2',date)
            iglist = iglist1 + iglist2
        
        elif range:
            keep = []
            for ig in self: 
                if (ig.Rsc['DATE2'] <= range[0]) or (ig.Rsc['DATE1'] >= range[1]):
                    iglist.append(ig.Rsc['PAIR'])
        
        for IG in iglist:
            self.Omissions[IG] = self.Igrams[IG]
            del self.Igrams[IG]
        
        print('Omitted the following interferograms:')        
        print('\n'.join(iglist))

        # Set up timeseries for new subset of interferograms
        self.setup()
        



    def remit(self, date=False, IG=False):
        """Put an omitted date back into the timeseries """
        iglist = []
        if IG:
            if type(IG) is str:
                iglist.append(IG)
            else:
                iglist = IG 
        if date:
            for IG in self.Omissions.keys():
                if date in IG.split():
                    iglist.append(IG)
        for IG in iglist:
            self.Igrams[IG] = self.Omissions[IG]
            del self.Omissions[IG]
        self.setup() #re-run setup() routine without ommitted files
        print('Returning the following interferograms:')        
        print('\n'.join(iglist))




        
