#!/usr/bin/env python3
'''
SIMAPSE - simulation maps for ecological niche modelling
Version 2.00 beta
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


from datetime import date
import logging
import logging.handlers as handlers

from .nnFuncs import calc2dlists

class recorder():
    '''Just a simple I/O manager to save and recall results files.
       Note:
          Levels represent diferent files (different variables).
          Names represent diferent data in the same file (different models).'''
    headeronly = False
    def __init__(self, outdir, fileend, levels = None, sep=';', 
                 append=False, prename='Model_'):
        '''Opens the recfile to start writing results
              outdir -> Output directory
              fileend-> A string to place at the end of filename
              levels -> Levels of data (if None, filename = prefix.txt)
              sep    -> separator of values
              append -> If True, opens the file in append mode.'''
        self.sep       = sep
        self.levels    = levels
        self.names     = [] #Keeps record of available data
        self.mem       = [] #Keeps track of memorized data
        self.prename   = str(prename)
        self.memfuncs  = ['Average', 'StDev']
        if outdir[-1] not in ['/', '\\']: outdir = outdir + '/'
        mode = 'w'
        if append: mode = 'a'
        #Just creates the empty file or deletes the old one
        if levels is not None:
            recfile = outdir + '%s' + '_' + fileend + '.txt'
            for level in levels:
                newfile = open(recfile % level, mode)
                newfile.close()
            self.file = recfile
        else:
            recfile = outdir + fileend + '.txt'
            newfile = open(recfile, mode)
            newfile.close()
            self.file = recfile

    def addheader(self, values, level = None, headertxt = 'Variable values'):
        '''Adds one line header to the file.'''
        txt = self.sep.join([str(x) for x in values])
        txt = headertxt + self.sep + txt
        if level is not None:
            filename = self.file % level
        else:
            filename = self.file
        recfile = open(filename, 'a')
        recfile.write(txt + '\n')
        recfile.close()

    def getheader(self, level):
        '''Returns level header.'''
        #TODO: Memorize headers?
        self.headeronly = True
        header, other = self.recall(level, None)
        self.headeronly = False
        return header        
       
    def write(self, values, level = None, name = False, mem=True):
        '''Writes the list of values in the file. When values are a list of 
           lists, each of the lists will be written in a diferent line.
            level  -> level where the value should be written (if 'None'
                      assumes no levels)
            name   -> is any string to be written in the first column
            mem    -> Memorizes the sum and the sum of squares of the values'''
        if type(values[0]) == list:
            if mem == True: self.memorize(values, level)
            if name:
                self.write([self.prename + name], level, False, mem=False)
                if name not in self.names  + self.memfuncs:
                    self.names.append(name)
            for value in values:
                self.write(value, level, False, mem=False)
        else:
            txt = self.sep.join([str(x) for x in values])
            if name:
                if name not in self.names + self.memfuncs:
                    self.names.append(name)
                txt = self.prename + str(name) + self.sep + txt
            if level is not None:
                filename = self.file % level
            else:
                filename = self.file
            recfile = open(filename, 'a')
            recfile.write(txt + '\n')
            recfile.close()
            if mem: self.memorize(values, level)

    def write_levels(self, diclevel, name = False, mem = True):
        '''Writes all values from a dictionary to file. Dictionary keys must
           be the levels of class'''
        for level in diclevel:
            self.write(diclevel[level], level, name, mem=mem)

    def recall(self, level , name):
        '''Reads results files and returns a list'''
        if level == None:
            datafile = open(self.file, 'r')
        else:
            datafile = open(self.file % level, 'r')
        foundata, data2d = False, False
        header, datalist = None, None
        line = 0
        while 1:
            data = datafile.readline()
            if not data: break
            data = data.strip().split(self.sep)
            if line == 0:
                header = [float(x) for x in data[1:]]
                line = 1
                if self.headeronly: break
                continue
            if data[0] == self.prename + name:
                if not foundata: foundata = True
                if len(data) > 1:
                    datalist = [float(x) for x in data[1:]]
                    break
                else:
                    data2d = True
                    break
        if data2d:
            datalist = []
            while 1:
                data = datafile.readline()
                if not data: break
                data = data.strip().split(';')
                if len(data) == 1: break
                datalist.append([float(x) for x in data])
        return header, datalist   

    def memorize(self, values, level=None):
        '''Memorizes a list of level's values that can be retrieved.
           If level already exists, it is summed to existing values.'''
        v = self.__dict__
        if level in self.mem:
            mem_s, mem_ss = v[level][0], v[level][1]
            if type(values[0]) == list:
                #from nnFuncs import calc2dlists
                s = calc2dlists(mem_s, values, sum)
                temp = [[y**2 for y in x] for x in values]
                ss = calc2dlists(mem_ss, temp, sum)
            else:
                s = [x + y for x,y in zip(mem_s, values)]
                ss = [x + y**2 for x,y in zip(mem_ss, values)]
            v[level][0], v[level][1] = s, ss
        else:
            v[level] = []
            v[level].append(values)
            if type(values[0]) == list: 
                v[level].append([[y**2 for y in x] for x in values])
            else:
                v[level].append([x**2 for x in values])
            self.mem.append(level)

    def finalize(self, level=None):
        '''Calculates the average and standard deviation of the level.'''
        v = self.__dict__
        N = len(self.names)
        stdfunc = lambda x: ((x[1] - (x[0]**2 / N )) / (N - 1))**0.5
        if N == 1: stdfunc = lambda x: 0 # to avoid division by 0
        if level in self.mem:
            mem_s, mem_ss = v[level][0], v[level][1]
            if type(mem_s[0]) == list:
                #from nnFuncs import calc2dlists
                avg = [[y/N for y in x] for x in mem_s]
                std = calc2dlists(mem_s, mem_ss, stdfunc)
            else:
                avg = [x/N for x in mem_s]
                std = [stdfunc([x,y]) for x,y in zip(mem_s, mem_ss)]
            v[level][0], v[level][1] = avg, std
            self.write(avg, level = level, name = 'Average', mem=False)
            self.write(std, level = level, name = 'StDev', mem=False)

    def getnames(self, level):
        '''Reads all names in the file if the list is empty.
           Note:
            Always ignores the first line (header)'''
        if self.names == []:
            sep = self.sep
            size = len(self.prename)
            if level == None:
                datafile = open(self.file, 'r')
            else:
                datafile = open(self.file % level, 'r')
            header = 1
            while 1:
                data = datafile.readline()                
                if header == 1: 
                    header = 0
                    continue
                if not data: break
                if data[:size] == self.prename: 
                    data = data.strip().split(sep)
                    self.names.append(data[0][size:])
        return self.names

    def getmem(self, level):
        '''Returns the list in memory for the level.
           If finalized retunrs [avg, std], otherwise
           it returns [sum, sum**2].'''
        return self.__dict__[level]

    def getavg(self):
        '''Returns the averages. If there are more than one level, returns a 
           dictionary, else it returns a list.'''
        if self.levels:
            avg = {}
            for lvl in self.levels:
                avg[lvl] = self.__dict__[lvl][0]
        else: avg = self.__dict__[None][0]
        return avg

    def getstd(self):
        '''Returns the standard deviations. If there are more than one level,
           returns a dictionary, else it returns a list.'''
        if self.levels:
            std = {}
            for lvl in self.levels:
                std[lvl] = self.__dict__[lvl][1]
        else: std = self.__dict__[None][1]
        return std

class html():
    '''Creates a simple html with text, links and tables.'''
    def __init__(self, name):
        '''Constructor with the name for the project.'''
        html = '<html>\n'
        html += '<head><title>%s</title></head>\n' % name
        html += '<body>\n'
        self.html = html
        self.htmldoc = True
        self.htmlbody = True
        self.htmlparagraph = False
        self.htmltable = False

    def title(self, title):
        '''A text in title style.'''
        if self.htmlbody == True:
            html = '<h1>%s<br></h1>\n' % title
            self.html += html
            self.hline()

    def chapter(self, chapter):
        '''Text in chapter style.'''
        if self.htmlbody == True:
            html = '<h2>%s<br></h2>\n' % chapter
            self.html += html

    def hline(self):
        '''Adds a horizontal line.'''
        html = '<hr>\n'
        self.html += html

    def paragraph(self):
        '''Starts or ends a paragraph.'''
        if self.htmlparagraph == False:
            html = '<p>\n'
            self.html += html
            self.htmlparagraph = True
        else:
            html = '</p>\n'
            self.html += html
            self.htmlparagraph = False

    def body(self):
        '''Starts or ends body tag.'''
        if self.htmlbody == True:
            html = '</body>\n'
            self.html += html
            self.htmlbody = False
        else:
            html = '<body>\n'
            self.html += html
            self.htmlbody = True         

    def textline(self, text):
        '''Adds a line of text to the project.'''
        if self.htmlparagraph == False:
            self.paragraph()
        if self.htmlbody == True:
            html = self.txtref(text)
            self.html += html

    def txtref(self, text):
        '''Returns a string with html code for a text line.'''
        html = '%s<br>\n' % text
        return html

    def bline(self):
        '''Adds a blanck line to the project'''
        if self.htmlbody == True:
            html = self.bref()
            self.html += html

    def bref(self):
        '''Returns a string with html code for a blanck line.'''
        html = '<br>\n'
        return html

    def endhtml(self):
        '''Ends the html tag.'''
        if self.htmlbody == True:
            self.body()
        if self.htmldoc == True:
            html = '</html>'
            self.html += html
            self.htmldoc = False

    def write(self, filename):
        '''Exports the project to a file.'''
        if self.htmldoc == True:
            self.endhtml()
        if filename[-5:] is not '.html':
            filename += '.html'
        htmlfile = open(filename, 'w')
        for line in self.html:
            htmlfile.write(line)
        htmlfile.close()

    def newtable(self, rows=2, cols=2, border=0, padding=2, spacing=2):
        '''Starts a new table of rows*cols size.'''
        if self.htmlbody == True:
            if self.htmltable == False:
                if self.htmlparagraph == True:
                    self.paragraph()
                html = '<table border=\"%s\" cellpadding=\"%s\" cellspacing=\"%s\">\n<tbody>\n' % (border, padding, spacing)
                self.html += html
                self.virtualtable = [[None for y in range(cols)]for x in range(rows)]
                self.htmltable = True

    def closetable(self):
        '''Closes the table and saves to project.'''
        if self.htmlbody == True:
            if self.htmltable == True:
                html = ''
                for row in self.virtualtable:
                    html += '<tr>\n'
                    for col in row:
                        if col == None:
                            col = '<br>'
                        html += '<td>%s</td>\n' % col
                    html += '</tr>\n'
                html += '</tbody>\n</table>\n'
                self.html += html
                self.htmltable = False
                del self.virtualtable
                #self.paragraph()

    def cell(self, row, col, text):
        '''Adds information to an individual cell in a open table.'''
        if self.htmltable == True:
            vt = self.virtualtable
            row = row - 1
            col = col - 1
            if row <= len(vt) and col <= len(vt[0]):
                vt[row][col] = text

    def image(self, imagefile, link = True, width = 0, height = 0, alt = 'Image file'):
        '''Adds an image to the project.'''
        if self.htmlbody == True:
            html = self.imgref(imagefile, link, width, height, alt)
            self.html += html

    def imgref(self, imagefile, link = True, width = 0, height = 0, alt = 'Image file'):
        '''Returns a string with html code for a image.'''
        html = ''
        if link == True:
            html += '<a href=\"%s\">' % imagefile
        dim = ''
        if width > 0:
            dim += 'width: %spx;' % width
        elif height > 0:
            dim += 'height: %spx;' %height
        html += '<img style=\"%s\" alt=%s\" src=\"%s\">' % (dim, alt, imagefile)
        if link == True:
            html += '</a>'
        html += '<br>\n'
        return html

    def href(self, link, text = None):
        'Returns a string with html code for a link.'
        if text == None:
            text = link
        link = '<a href=\"%s\">%s</a>' % (link, text)
        return link

class htmlreport(html):
    '''Creates an html report for SIMAPSE results.'''
    def __init__(self):
        '''Constructor with project name'''
        html.__init__(self, 'SIMAPSE Results')
        self.title('SIMAPSE - Simulation Maps in Ecology')

    def summary(self, nnoptions):
        '''Adds a summary of the models and neural networks details.'''
        #from datetime import date

        if nnoptions['method'] == "Random repetition":
            details = '%s repetitions' % nnoptions['repetitions']
        elif nnoptions['method'] == "Cross validation":
            details = '%s k-folds' % nnoptions['repetitions']
        elif nnoptions['method'] == "Bootstrapping":
            details = '%s bootstraps of %s%% of the original size' % (nnoptions['repetitions'], nnoptions['bsize'])

        totaliter = nnoptions['iterreport'] * nnoptions['iterinter']

        text1 = ['Model created on %s' % date.today(),
                 'Data file: %s' % nnoptions['file_data'],
                 'Raster folder: %s' % nnoptions['dir_rasters'],
                 'Output folder: %s' % nnoptions['out_dir']]
        text2 = ['Neural Network Scheme: %s' % nnoptions['scheme'],
                 'Learning Rate: %s' % nnoptions['lrate'],
                 'Momentum: %s' % nnoptions['momentum'],
                 'Iterations: %s (%s auc reports x %s internal iterations)' % (totaliter, nnoptions['iterreport'], nnoptions['iterinter'])]
        text3 = ['Subsampling Method: %s with %s' % (nnoptions['method'], details),
                 'Test Percentage: %s' % nnoptions['percentage'],
                 'Pseudoabsences / presence ratio: %s' % nnoptions['apratio'],
                 'AUC Thresholds: %s for train and %s for test' % (nnoptions['auctrain'], nnoptions['auctest'])]

        self.chapter('Summary')
        self.paragraph()
        for line in text1:
            self.textline(line)
        self.bline()
        for line in text2:
            self.textline(line)
        self.bline()
        for line in text3:
            self.textline(line)
        self.hline()
        self.paragraph()

    def model(self):
        '''Adds the model results that include the average and standard
           deviation models, the variable importance and ROC and Precision
           Recall curves.'''
        self.chapter('Model results')
        self.paragraph()
        text = 'The image below shows the average model and the standard deviation of all the sucssesfuly built models.'
        self.textline(text)
        self.image('Results_map.png', True, 750)
        text = 'The following plots show the average variable importance for the model with standard deviation and the ROC and Precision/Recall curve with the respective AUC values.'
        self.textline(text)
        self.image('Results_variables_roc.png', True, 750)
        ref = self.href('log.txt', 'log file')
        text = 'The sum of squared partial derivatives for each variable indicates the influence of each variable in the model. This plot represents the average score from all models built with the standard deviation. You can find more details of variables scores in the %s.' % ref 
        self.textline(text)
        text = 'The ROC and Precision-Recall curves indicate the fitness of the model and are represented in blue and green, respectively. The Area Under the Curve (AUC) is calculated for each curve, and ranges from 0.0 to 1.0.' 
        self.textline(text)
        self.hline()
        self.paragraph()

    def variables(self, variables):
        '''Automaticaly adds tables with the variables and links to data.'''
        nvars = len(variables)
        ncols = 3 #only 3 plots per line
        nrows = int(nvars / ncols)
        
        if nvars % ncols != 0:
            nrows += 1
        self.chapter('Variables results')
        self.paragraph()
        text = 'The SIMAPSE produce variable plots to infer the behaviour of each variable in the final model. All values are averaged from the models that were successfully built.'
        self.textline(text)
        self.bline()
        text = 'The profile plots represent the variable predictivity when all other variables are set to 0. The range of values of each variable (x-axis) is plotted against the output of the model (y-axis), showing the usage of the variable by the model.'
        self.textline(text)
        self.addtable(variables, 'profile', nrows, ncols)
        self.bline()
        text = 'The graphs of the partial derivatives represent the network sensitivity to the variation of each input on the presence / absence list. The values of each variable (x-axis) are represented according to the partial derivatives of each input (y-axis).'
        self.textline(text)
        self.addtable(variables, 'PaD', nrows, ncols)
        self.bline()
        text = 'The variation surface of each variable repesents the behaviour of each variable within the range of the others variables. The variables range (x-axis) is ploted against the the range of all others variables (between 0 and 100% - y-axis) and the predictiviy is given by the color (0.0 to 1.0, from red to green, respectively).'
        self.textline(text)
        self.addtable(variables, 'varsurface', nrows, ncols)
        self.hline()
        self.paragraph()

    def addtable(self, variables, vartype, nrows, ncols):
        '''Creates a table with variables.'''
        text  = self.txtref   
        href  = self.href
        image = self.imgref
        bline = self.bref
        newtb = self.newtable
        clstb = self.closetable
        cell  = self.cell
        nvars = len(variables)
        newtb(nrows,ncols)
        for i in range(nvars):
            row = int((i / 3 ) + 1)
            col = int((i - (3 * row)) + 1)
            varname = variables[i]
            txt = '%s%s%s' % (text(varname + ' - ' + href('%s_%s.txt' % (varname, vartype), 'Data file')), 
                              image('%s_%s.png' % (varname, vartype), True, 300),
                              bline())
            cell(row, col, txt)
        clstb()


class logger(object):
    '''Setups and manages a logger'''

    format = logging.Formatter("%(message)s")

    def __init__(self, namestr, logFile=None, level=logging.DEBUG):
        ''' Initiates a memory handler if no file is defined'''
        self.log = logging.getLogger(namestr)
        self.log.setLevel(level)

        self.logFile = logFile

        if logFile is None:
            self._memHandler()
        else:
            self._fileHandler()

    def _memHandler(self):
        '''Defines a memory handler'''

        self.memhandler = handlers.MemoryHandler(capacity=1024)
        self.memhandler.setFormatter(self.format)
        self.log.addHandler(self.memhandler)

    def _fileHandler(self):
        '''Sets a file handler'''

        self.filehandler = logging.FileHandler(self.logFile)
        self.filehandler.setFormatter(self.format)
        self.log.addHandler(self.filehandler)

    def _removeMemHandler(self):
        '''Removes the memory handler after flushing to file'''
        if self.memhandler:
            self.memhandler.setTarget(self.filehandler)
            self.memhandler.flush()
            self.log.removeHandler(self.memhandler)
            self.memhandler = None

    def setLogFile(self, logFile):
        '''Creates a file handler. 
           Flushes and remove the memory handler'''

        self.logFile = logFile
        self._fileHandler()
        self._removeMemHandler()
