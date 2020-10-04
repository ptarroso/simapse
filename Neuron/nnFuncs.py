#!/usr/bin/env python
'''
SIMAPSE - simulation maps for ecological niche modelling
Version 1.01 beta
Copyright (C) 2010  Pedro Tarroso

Please cite: 
"Tarroso, P., Carvalho, S. & Brito, J.C. (2012) Simapse - Simulation
Maps for Ecological Niche Modelling. Methods in Ecology and Evolution
doi: 10.1111/j.2041-210X.2012.00210.x"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from os import curdir, path

#from matplotlib.backends.backend_agg import FigureCanvasAgg
#from matplotlib import figure, cm, ticker
#import numpy as np

# SOME MESSAGES TO DISPLAY
WARNING_ROWS = 'WARNING: Found %s rows in the file %s! Data was trimmed for '+\
               '%s rows.'
WARNING_COLS = 'WARNING: Found %s columns in the file %s! Data was trimmed ' +\
               'for %s columns.'
ERROR_HEADER = 'ERROR: Please check the raster %s! The number of cols or ' +\
               'rows found in the raster does not correspond to values on ' +\
               'the header.'

conn = None
def sendout(func, *args, **kwargs):
    '''Sends information to console or gui'''
    if conn == None:
        if func == 'put':
            print args[0]
        elif func <> 'put':
            pass
    elif conn <> None:
        if func == 'put':
            text = args[0]
            conn.display_msg(text)
        elif func == 'progress':
            conn.progress_bar(*args, **kwargs)

def read_ascii(filename, output = 0):
    '''Reads data from ascii raster file and the corner's coordinates
       The output option may be:
            0 to export just the array (list of lists)
            1 to export the array, ncols, nrows, xllcorner, yllcorner, cellsize, nodata
            2 to export all parameters without the array

        NOTE: if the columns or the lines found in the raster body will be
              trimmed when higher than the ncols or nrows values of the header.
              In the case of lower values, an Exception is raised.'''
    def readbody():
        try:
            MyData = [[float(value) for value in row.split()] for row in contents[6:]]
            rRows, rCols = len(MyData), len(MyData[0])
            if rRows > nrows:
                MyData = MyData[:nrows]
                print WARNING_ROWS % (rRows, filename, nrows)
            if rCols > ncols:
                MyData = [row[:ncols] for row in MyData]
                print WARNING_COLS & (rCols, filename, ncols)
            if rRows < nrows or rCols < ncols:
                raise Exception
            return MyData
        except Exception:
            print ERROR_HEADER % filename 

    #tem q fazer um print do cabecalho com um %14string para ficar igual
    myfile = open(filename, "r")
    contents = myfile.readlines()
    myfile.close()

    #reads header and convert data to variables
    ncols = int(contents[0].split()[-1])
    nrows = int(contents[1].split()[-1])
    xllcorner = float(contents[2].split()[-1])
    yllcorner = float(contents[3].split()[-1])
    cellsize = float(contents[4].split()[-1])
    nodata = float(contents[5].split()[-1])

    if output == 0:
        #Create a list with data
        MyData = readbody()    
        return MyData
    elif output == 1:
        #Create a list with data
        MyData = readbody()  
        return MyData, ncols, nrows, xllcorner, yllcorner, cellsize, nodata
    elif output == 2:
        #reads header and convert data to variables
        return ncols, nrows, xllcorner, yllcorner, cellsize, nodata

def calc2dlists(list1, list2, func):
    '''Applies a function 'func' to each value of two 2d lists.'''
    newlist = [[func([v1,v2]) for v1,v2 in zip(l1,l2)] for l1,l2 in zip(list1,list2)]
    return newlist

def transpose(list1):
    '''Transpose a 2D list (list of lists)'''
    cols = len(list1[0])
    newlist = [[x[y] for x in list1] for y in xrange(cols)]
    return newlist

def sum2d(list2d):
    '''Sums all rows of a 2d list (list of lists).
       All lists must be of same size.'''
    list2d = list2d[:]
    newlist = list2d.pop(0)
    for row in list2d:
        newlist = [x + y for x,y in zip(newlist, row)]
    return newlist

def stvar(values, avg, std, nodata=None):
    '''Standardize a variable (list of lists) based on the average and
       standard deviation. The formula applied is:

            Z = (x - Average) / Standard Deviation

        values - Raster values (list of lists)
        avg    - Average value for raster
        std    - Standard deviation for raster'''
    avg, std = float(avg), float(std)
    newraster, results = [], []
    for row in values:
        temp = []
        for value in row:
            if value == nodata and nodata <> None:
                temp.append(nodata)
            else:
                result = (value - avg) / std
                temp.append(result)
                results.append(result)
        newraster.append(temp)
    stmax, stmin = max(results), min(results)
    return newraster, stmax, stmin

def unstd(values, stats):
    '''Un-standardize a value or a list of values based on stats'''
    avg, std, rmax, rmin, stmax, stmin = stats
    if type(values) == list:
        results = [(x * std) + avg for x in values]
    else:
        results = (values * std) + avg
    return results

def readstats(stsfile):
    '''Reads sts file and return list with values.'''
    stsfile = open(stsfile, 'r')
    contents = stsfile.readlines()
    stats = contents[0].split(';')
    stats = [float(x) for x in stats]
    return stats

def reducelist(lst):
    '''Returns an item or the first list of items inside a list of lists'''
    if type(lst) is list:
        if type(lst[0]) is list and len(lst[0]) == 1:
            result = reducelist(lst[0])
        elif type(lst[0]) is not list and len(lst) > 1:
            result = lst
        else: 
            result = lst[0]
    else:
        result = lst
    return result

def check(test):
    '''Checks if all elements in a list are equal'''
    N = len(test)
    result = min([test[x] == test[x+1] for x in xrange(N-1)])
    return result


class spatial_functions:
    '''Several spatial functions that works in a given extent
       Output dir is current directory and may be overwritten in self.out_dir'''
    from random import choice
    def __init__(self, ncols, nrows, xllcorner, yllcorner, cellsize, nodata):
        #nodata_list is a list with raster size where nodata = 1 and else = 0

        self.ncols = ncols
        self.nrows = nrows
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.nodata = nodata

        self.nodata_list = None
        self.out_dir = curdir

        self.Row, self.Col = 0,0 #When passing functions to rastercalc to access current row, col

    def write_ascii(self, OutData, filename):
        '''Writes an ascii raster with header
           Data to save must be a list of lists'''

        myfile = open(filename, "w")
        header = ["ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "nodata_value"]
        myfile.write('%-14s %s\n' % (header[0], self.ncols))
        myfile.write('%-14s %s\n' % (header[1], self.nrows))
        myfile.write('%-14s %s\n' % (header[2], self.xllcorner))
        myfile.write('%-14s %s\n' % (header[3], self.yllcorner))
        myfile.write('%-14s %s\n' % (header[4], self.cellsize))
        myfile.write('%-14s %s\n' % (header[5], self.nodata))

        separator = ' '
        format = '%5.3f'

        for row in OutData: 
            myfile.write(separator.join([format % value for value in row]))
            myfile.write('\n')
        myfile.close()
        return

    def create_nodata_list(self, raster):
        #Creates a nodata list of the same size of raster where
        # nodata values are 1 (True) and else 0 (False)
        if type(raster) == str:
            temp = read_ascii(raster)
        elif type(raster) == list:
            temp = raster
        self.nodata_list = [map(lambda x: x == self.nodata, line) for line in temp]

    def ExtractValues(self, indata, raster_values, order = None, APRatio = 1,
                      outdir = None):
        '''Extract values from the intersection of data points
           and a list of rasters filenames
           indata  -  Input x,y data for presences (a list or a text file)
           rasters -  A dictionary of raster values 
           order   -  A list of rasters names to extract by order
           args    -  Absence/Presence ratio'''

        if outdir == None:
            outdir = self.out_dir

        rasters = raster_values.keys()
        if order:
            rasters = order

        if self.nodata_list == None:
            self.create_nodata_list(raster_values[rasters[0]])
        
        if type(indata) is list: 
            #Assumes that it has absences or pseudo-absences already
            coordinates = indata
        else:           
            # reads from file if input data is not as list                   
            coordinates, MyData = self.read_presencedata(indata, APRatio, outdir) 
        
        variables = []
        counter = 1
        for raster in rasters:
            RasterValue = raster_values[raster]
            
            if type(indata) <> list: 
                text = "Reading values for raster \'%s\'" % (path.basename(raster))
                sendout('put', text)
            
            temp_array = [RasterValue[y][x] for x,y in coordinates]
            variables.append(temp_array)
            sendout('progress', counter, len(rasters), msg=raster)
            counter += 1
        MyVariables = [reducelist([variables[y][x] for y in range(len(rasters))]) for x in range(len(coordinates))]

        if type(indata) == list:
            return MyVariables
        else:
            self.MyData, self.MyVariables, self.coordinates = MyData, MyVariables, coordinates
            return MyData, MyVariables, coordinates

    def read_presencedata(self, indata, APRatio = 1, outdir = None):
            '''Reads the presence data file and checks if it has absence or abundance data.
               Converts the input coordinates to raster coordinates (pixel location).

               Indata - Text file with presence data
               APRatio - Absence / Presence ratio to create pseudo absences'''

            if outdir == None:
                outdir = self.out_dir

            myfile = open(indata, "r")
            contents = myfile.readlines()
            myfile.close()
            data = map(lambda x: (map(float, x.split(';'))), contents[1:])
            
            ndata = len(data) #TODO Add 'no data to read' error

            abundance_data, absence = False, 0
            MyData, coordinates = [], []

            for line in data:
                MyData.append([line[0]])
                x = int((line[1] - self.xllcorner) / self.cellsize)
                y = int(self.nrows - ((line[2] - self.yllcorner) / self.cellsize))
                coordinates.append((x, y))
                if line[0] == 0: 
                    absence += 1
                elif line[0] <> 1:
                    abundance_data = True

            if absence == 0 and abundance_data == False:
                sendout('put', "Randomly generating pseudo-absences...")
                pa = int(ndata * APRatio)
                # NO limit set!!! TODO: add a warning if more than 1000 pseudoabsences!
                #if pa > 1000:
                #    pa = 1000
                a_data, a_coordinates, a_MyData = self.pseudo_absences(coordinates, pa)
                data = data + a_data
                coordinates = coordinates + a_coordinates
                MyData = MyData + a_MyData                
                self.save_presences(data, outdir)
                text = "%s pseudo-absences for input data were created" % (pa)
                sendout('put', text)
            elif abundance_data == True:
                sendout('put', "Abundance data found! AUC will be unavailable")
            else:
                text = "%s absences found. No pseudo-absences were generated." % (absence)
                sendout('put', text)

            self.abundance = abundance_data

            return coordinates, MyData

    def savecache(self, raster_values, rasterstats, outdir= None):
        '''Saves all rasters values into new ascii raster and also the raster
           statistics into a 'raster'.sts file'''

        if outdir == None:
            outdir = self.out_dir

        for raster in raster_values:
            outraster = '%s/%s.txt' % (outdir, raster)
            outrstats = '%s/%s.sts' % (outdir, raster)
            self.write_ascii(raster_values[raster], outraster)
            rstats = open(outrstats, 'w')
            text = map(str, rasterstats[raster])
            rstats.write(';'.join(text))
            rstats.close()

    def pseudo_absences(self, coordinates, pa_number):
        '''Creates pseudo absences in the raster area, excluding nodata.
           The number of pseudo-absences is based on the number of presences
           (multiplied by the A/P ratio, if given), until a maximum of 1000

            pseudo_absences(coordinates, data, APRatio, *kwargs) -->  a_data, a_coordinate, a_MyData

            coordinates -  list of coordinates with presence
            pa_number   -  number of pseudo-absences to create'''

        def absences(coordinates):
            # Extracts a list of all raster pixels without presence
            y = [n for n in xrange(self.nrows) for r in xrange(self.ncols)]
            x = [x for x in xrange(self.ncols)] * self.nrows
            allCoord = map(None, x, y)
            lastPresence = 0
            for presence in coordinates:
                if presence <> lastPresence: #TODO May add a duplicate presence error!!!
                    allCoord.remove(presence)
                lastPresence = presence
            return allCoord

        testRaster = self.nodata_list
        absences = absences(coordinates)

        if pa_number > len(absences):
            text = "Number of pseudo absences higher than available locations for absences!\n\
                    Limiting the number to the available locations..."
            sendout('put', text)            
            pa_number = len(absences)

        a_data, a_coordinates, a_MyData = [], [], []
        
        p = 0
        while p < pa_number: 
            xy = self.choice(absences)
            if testRaster[xy[1]][xy[0]] <> 1:
                a_coordinates.append(xy)
                a_data.append([0, (((xy[0] * self.cellsize) + self.xllcorner) + (self.cellsize / 2)), 
                                (((self.nrows - xy[1])* self.cellsize) + self.yllcorner) - (self.cellsize / 2)])
                a_MyData.append([0]) 
                p += 1
                del xy
            else:
                continue
            sendout('progress', p+1, pa_number)

        return a_data, a_coordinates, a_MyData

    def save_presences(self, data, outdir):
        #save new presence/absence data file
        myfile = open(outdir + '/new_data.txt', "w")
        myfile.write('Pres;X;Y\n')
        for line in data:
            export = '%s;%s;%s\n' % tuple(line)
            myfile.write(export)
        myfile.close()

    def rastercalc(self, func, rasterdic, outname = None, order = None):
        '''Computes a final map by apllying a function 'func' to the same
           pixel in all rasters. The function input is a list of all rasters'
           values for each pixel and, if outname is None, returns the new 
           raster dic. Otherwise, it will write the raster with outname.
            
            func      - Function to apply
            outname   - Name for the output ascii raster (adds the self.out_dir to the name)
            rasterdic - Dictionary of all rasters'''

        FinalModel = [[0.0] * self.ncols for x in xrange(self.nrows)]

        rasters = rasterdic.keys()
        if order:
            rasters = order
        #Computes the final model map

        for row in xrange(self.nrows):
            for col in xrange(self.ncols):
                self.Row, self.Col = row, col
                dataArray = [rasterdic[raster][row][col] for raster in rasters]
                if self.nodata_list[row][col] == 1:
                    FinalModel[row][col] = self.nodata
                else:
                    FinalModel[row][col] = reducelist(func(dataArray))

        if outname == None:
            return FinalModel
        else:
            self.write_ascii(FinalModel, outname)
            return

    def modelstats(self, rasters, outdir = None, output = True, sufix=''):
        '''Computes the final average and standard deviation models and saves as ascii raster.
           When output is True, returns average and standard deviation rasters.'''

        if outdir == None:
            outdir = self.out_dir

        #Create a dictionary of rasters
        rasterdic = {}
        for raster in rasters:
            rasterdic[raster] = read_ascii(raster, 0)

        N = float(len(rasters))
        
        #Averages all maps
        average = lambda x: float(sum(x)) / N
        avg = self.rastercalc(average, rasterdic)

        #Calculate standard deviation of all maps:
        stdev = lambda x: ((sum([(value - avg[self.Row][self.Col])**2 for value in x])) / N)**0.5
        std = self.rastercalc(stdev, rasterdic)

        if sufix != "":
            sufix = '_' + sufix

        average_file = outdir + "/average" + sufix + ".txt"
        std_file = outdir + "/std" + sufix + ".txt"
        self.write_ascii(avg, average_file)
        self.write_ascii(std, std_file)

        if output == True:
            return avg, std
        else:
            return

    def emptyraster(self, nodata=True):
        '''Creates a raster with 0.0 values
           If nodata is True, the new raster will have nodata values
           equivalent to nodata_list'''
        try:
            if nodata==True and self.nodata_list == None:
                raise Exception('There is no information for nodata. \n\
                                 Create nodata_list.')
            ncols, nrows = self.ncols, self.nrows
            raster = [[0.0] * ncols for x in nrows]
            if nodata ==True:
                ndlist = self.nodata_list
                f = lambda value, check: check == True and self.nodata or value
                raster = [[f(v, c) for v,c in zip(x,y)] for x,y in zip(raster,ndlist)]
            return raster
        except Exception, e:
            print e

class subsets():
    from random import randint, choice, shuffle
    def __init__(self, MyData, MyVariables, coordinates):

        self.MyData = MyData
        self.MyVariables = MyVariables
        self.coordinates = coordinates

    def random_data(self, Percentage):
        '''Creates two random datasets for test and train. Train dataset is a 
           'Percentage' of the orignal dataset'''
        MyData, MyVariables, coordinates = self.MyData, self.MyVariables, self.coordinates

        totalNumb = len(MyData)
        testNumb = int(totalNumb * (float(Percentage) / 100.00))

        value, valueTest = 0, 0
        #Test if there are only zeros or only ones and repeates 
        #until both have presences and absences
        while value == 0 or value == totalNumb or valueTest == 0 or valueTest == testNumb: 
            RndData, RndVariables = MyData[:], MyVariables[:]
            RndTestData, RndTestVariables = [], []
            trainCoord, testCoord = coordinates[:], []
            for item in range(testNumb):
                index = self.randint(0, len(RndData) - 1)
                RndTestData.append(RndData.pop(index))
                RndTestVariables.append(RndVariables.pop(index))
                testCoord.append(trainCoord.pop(index))
            value = sum([x[0] for x in RndData])
            valueTest = sum([x[0] for x in RndTestData])

        self.RndData = RndData
        self.RndVariables = RndVariables 
        self.RndTestData = RndTestData
        self.RndTestVariables = RndTestVariables
        self.RndCoordinates = [trainCoord, testCoord]

        return [self.RndData, self.RndVariables], \
               [self.RndTestData, self.RndTestVariables],\
                self.RndCoordinates

    def repeatData(self, Percentage, Repetitions):
        '''Creates random repetitions of input Data and Variables'''
        for i in xrange(Repetitions):
            self.random_data(Percentage)
            yield self.RndData, self.RndVariables, self.RndTestData, self.RndTestVariables

    def bootstrapData(self, Bsize, Btstrps, Percentage, **kwargs):
        '''Creates a bootstrap sample of Data with Bsize length

            Bsize - percentage of the original size
            Btstrps - Number of bootstraps to generate
            Percentage - Percantage to calculate test and train data

           Output is a generator yielding:

                BdataTrain, BvariablesTrain, BdataTest, BvariablesTest'''
        #from random import randint

        if kwargs == {}:
            self.random_data(Percentage)
            NTrain, NTest = len(self.RndData), len(self.RndTestData)
            BsizeTrain = int(NTrain * (float(Bsize) / 100.00))
            BsizeTest  = int(NTest * (float(Bsize) / 100.00))
        else:
            NTrain, NTest = kwargs['NTrain'], kwargs['NTest']
            BsiseTrain, BsizeTest = kwargs['BsiseTrain'], kwargs['BsizeTest']
        
        try:
            if Bsize < 1 :
                raise Exception("Small_bootstrap")
            value, valueTest = 0, 0

            for bootstrap in xrange(Btstrps):
                BdataTrain, BvariablesTrain, BdataTest, BvariablesTest = [], [],[], []
                for item in xrange(BsizeTrain):
                    index = self.randint(0, NTrain-1)
                    BdataTrain.append(self.RndData[index])
                    BvariablesTrain.append(self.RndVariables[index])
                for item in xrange(BsizeTest):
                    index = self.randint(0, NTest-1)
                    BdataTest.append(self.RndTestData[index])
                    BvariablesTest.append(self.RndTestVariables[index])
                if check(BdataTrain) == True or check(BdataTest) == True:
                    data = self.bootstrapData(Bsize, 1, Percentage,
                                              NTrain = NTrain, NTest = NTest,
                                              BsiseTrain = BsizeTrain, BsizeTest = BsizeTest)
                    BdataTrain, BvariablesTrain, BdataTest, BvariablesTest = data.next()
                    #TODO Add a maximum amount of trials to get bootstraps
                    #TODO or to add a clause where there is a minimum size
                    #TODO to do bootstraps (percentage data and bsize can
                    #TODO not result in less of 4 points to test!)
                yield BdataTrain, BvariablesTrain, BdataTest, BvariablesTest

        except Exception, msg:
            if str(msg) == "Small_bootstrap":
                sendout('put', "Bootstrap sample to small! Increse bootstrap sample size.")

    def kfoldData(self, k):
        '''Creates k Folds subsets of Data'''
        try:
            if k > len(self.MyData) :
                raise Exception("nfolds")
            elif k == len(self.MyData):
                sendout('put', "ATENTION: Leave-one-out Cross Validation (k-folds = number of data points")

            ndata = len(self.MyData)
            Items = ndata/k
            Extra = k - (ndata % k)
            NormalItems = Items * Extra

            sortList = range(ndata)
            self.shuffle(sortList)
            self.RndData = [self.MyData[x] for x in sortList]
            self.RndVariables = [self.MyVariables[x] for x in sortList]
            self.RndCoordinates = [self.coordinates[x] for x in sortList]

            for i in xrange(k):
                trainData = self.RndData[:]
                trainVariables = self.RndVariables[:]
                if i < Extra:
                    testData = [self.RndData[x+(Items*i)] for x in xrange(Items)]
                    testVariables = [self.RndVariables[x+(Items*i)] for x in xrange(Items)]
                    trainData[Items*i:Items*(i+1)] = []
                    trainVariables[Items*i:Items*(i+1)] = []
                elif i >= Extra:
                    testData = [self.RndData[x+(NormalItems+(Items+1)*(i-Extra))] for x in xrange(Items+1)]
                    testVariables = [self.RndVariables[x+(NormalItems+(Items+1)*(i-Extra))] for x in xrange(Items+1)]
                    trainData[NormalItems+(Items+1)*(i-Extra):NormalItems+(Items+1)*(1+i-Extra)] = []
                    trainVariables[NormalItems+(Items+1)*(i-Extra):NormalItems+(Items+1)*(1+i-Extra)] = []
                yield trainData, trainVariables, testData, testVariables

        except Exception, msg:
            if str(msg) == "nfolds":
                sendout('put', "Decrease the number of folds.")

class roc():
    def __init__(self, real, pred):
        '''ROC, Precision/Recall and AUC calculation. Algorithm is based on:
            Fawcett, T., 2004. ROC Graphs: Notes and Practical Considerations
            for Researchers. http://home.comcast.net/~tom.fawcett/public_html/papers/ROC101.pdf'''

        try:
            #Count all real Positives and Negatives
            P, N = real.count(1), real.count(0)
            if P==0 or N==0:
                raise Exception("zero_values")

            self.counter= (P,N)

            #Initiate False and True Positives counter
            FP, TP = 0.0, 0.0
            points = []
            
            #Sort the real values by the predicted values
            combined = map(lambda x,y: [x,y], real,pred)
            combined.sort(lambda x,y:cmp(x[1],y[1]), reverse=True)

            temp = -1E400
            #For all available values
            for i in xrange(len(combined)):
                if combined[i][1] != temp and combined[i][1] <= 1:
                    points.append((FP,TP))
                    temp = combined[i][1]
                if combined[i][0] == 1:
                    TP += 1.0
                else:
                    FP += 1.0
                if len(combined) == 1:
                    points.append((FP,TP))
            points.append((FP,TP))
            self.points = points
            self.RocPoints = None

        except Exception, e:
            if str(e) == "zero_values":
                #If there are only zeros or ones is not possible to compute the Roc and AUC (divide by zero) 
                #One solution may be this exception handles this by giving a auc value of 1.
                print "Error: Real values are only presences or absences. ROC cannot be computed."

    def roc(self):
        P,N = self.counter
        RocPoints = self.points[:]
        #(1 - Specificity (FP/N), Sensivity (TP/P))
        RocPoints = [(x/N, y/P) for x,y in RocPoints]
        self.RocPoints = RocPoints

    def precision_recall(self):
        '''Calculates precision/recall values'''
        P,N = self.counter
        #(Precision (TP / (TP+FP)), Recall(TP/P))
        f = lambda x,y: x + y > 0 and (y/(x+y)) or 1
        pr = [(f(x,y), y/P) for x,y in self.points]
        self.PRPoints = pr

    def auc(self, Points = None):
        '''Calculates AUC value.'''
        if Points == None and self.RocPoints <> None:
            Points = self.RocPoints
        elif Points == None and self.RocPoints == None:
            self.roc()
            Points = self.RocPoints
        
        #TODO Find a better way to define how to calculate AUC for ROC and PR
        x, y = 0, 1
        if Points[0][0] > Points[0][1]:
            x,y = 1, 0

        area = 0
        for position in range(len(Points)-1):
            deltaX = abs(Points[position][x] - Points[position + 1][x])
            deltaY = (Points[position][y] + Points[position + 1][y])/2
            area += deltaX * deltaY
        return area

    def process_all(self):
        self.roc()
        self.precision_recall()
        aucROC = self.auc(self.RocPoints)
        aucPR = self.auc(self.PRPoints)
        rocplot = zip(*self.RocPoints)
        prplot = zip(*self.PRPoints)
        return rocplot, prplot, aucROC, aucPR

class profiler():
    '''Produces the range of each variable and stores the real and the
       standardized values.'''
    def __init__(self, varstats, order = None, r=100):
        self.__percent(varstats, int(r))
        self.r = r
        self.order = order

    def __percent(self, varstats, r):
        '''Stores a dictionary of variables that contains standardize and real
           values based on the max and min values of each raster.
           NOTE: Assumes standardization of variables.'''
        percent_vars, real_vars = {}, {}
        for var in varstats:
            avg, std, rmax, rmin, stmax, stmin = varstats[var]
            stdif = stmax - stmin
            rdif = rmax - rmin
            percent_vars[var] = [((x * stdif) / r) + stmin for x in range(r+1)]
            real_vars[var] = [((x * rdif) / r) + rmin for x in range(r+1)]

        self.pvars = percent_vars
        self.rvars = real_vars

    def onewayprofile(self, raster, func):
        '''Calculates the profile for a raster variable.'''
        pvars = self.pvars
        r = self.r
        pos = self.order.index(raster)
        N = len(self.order)

        # For each raster creates a temporary list for varsur and profile to append results
        profile = []
        for value in xrange(r+1):
            #Creates the profile data with all values 0.0 except for raster
            temp = [0.0] * N
            temp[pos] = pvars[raster][value] 
            profile.append(reducelist(func(temp)))
   
        return profile

    def twowayprofile(self, raster1, raster2, func):
        '''Creates a two-way profile: one raster variable against others'''
        pvars = self.pvars
        r = self.r
        pos1 = self.order.index(raster1)
        pos2 = self.order.index(raster2)
        N = len(pvars.keys())

        # For each raster creates a temporary list for varsur and profile to append results
        twowayprof = []
        for row in xrange(r, -1, -1):
            temp = [0.0] * N
            temp[pos2] = pvars[raster2][row] 
            twowayprof.append([])
            for col in xrange(r+1):
                temp[pos1] = pvars[raster1][col] 
                twowayprof[r-row].append(reducelist(func(temp)))

        return twowayprof

    def varsurface(self, raster, func):
        '''Processes the variation surface for one raster, i.e., the
          variation of prediction of one raster variable throughout all the
          others' raster variables range.'''

        pvars = self.pvars
        r = self.r
        pos = lambda x: self.order.index(x)
        N = len(pvars.keys())

        # Create a list of 'others' raster variables
        others = self.order[:]
        others.remove(raster)

        # For each raster creates a temporary list for varsur and profile to append results
        varsur = []
        for row in xrange(r, -1, -1):
            temp = [0.0] * N
            varsur.append([])
            for other in others:
                temp[pos(other)] = pvars[other][row]

            for col in xrange(r+1):
                temp[pos(raster)] = pvars[raster][col]
                varsur[r-row].append(reducelist(func(temp)))

        return varsur

