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


import random
import threading
import os.path as path

from .nnRecorder import logger
from . import nnManager

# CODE VARS
CODE_TEXT          = 'message'
CODE_ERROR         = 'error'
CODE_PROGRESS      = 'progress'
CODE_MODIFY_BUTTON = 'button'
CODE_WRITELOG      = 'writelog'
CODE_SHOWRESULTS   = 'showresults'
CODE_GRAPHS        = 'graphs'

class simargs(object):
    # Set some default values
    simargs = {}
    simargs['file_data']     = '.'
    simargs['dir_rasters'] = '.'
    simargs['out_dir']       = '.'
    simargs['project_dir']   = '.'
    simargs['method']        = 'Random repetition'
    simargs['repetitions']   = 5
    simargs['bsize']         = 100
    simargs['iterreport']    = 250
    simargs['iterinter']     = 10
    simargs['lrate']         = 0.1
    simargs['momentum']      = 0.1
    simargs['hiddenlyrs']    = '3'
    simargs['apratio']       = 1
    simargs['percentage']    = 25
    simargs['aucfilter']     = False
    simargs['auctrain']      = 0.800
    simargs['auctest']       = 0.800
    simargs['abundvar']      = False
    simargs['burnin']        = 25

    def __init__(self):
        pass

    def getVar(self, var):
        '''Returns the value of the variable.'''
        try:
            if var in self.simargs.keys():
                return self.simargs[var]
            else:
                raise NameError('Variable %s does not exist.')
        except:
            pass
    
    def setVar(self, var, value):
        '''Sets the value of the variable or creates a new variable'''
        self.simargs[var] = value

class connection(simargs):
    def __init__(self):
        #threading.Thread.__init__(self)
        self.manager = nnManager.Manager(self)
        self.gui = False
        self.logger = logger('Simapse log')
        self.LOG = self.logger.log
        self.queue = False

    def _initiateLogFile(self):
        '''initiates a file handler for log after defining output folder'''
        self.logger.setLogFile(path.join(self.simargs['out_dir'], 'log.txt'))

    def display_msg(self, *args, **kwargs):
        ''' Displays a message.
            If GUI: Wraper for queue.put_nowait() but with code for message
            If CL: print message on terminal'''
        if self.gui:
            self.queue.put_nowait([CODE_TEXT, args, kwargs])
        else:
            print(args[0]) #TODO add a streamer to the logger
        self.LOG.info(args[0])

    def modify_button(self, *args, **kwargs):
        ''' Modifies buttons state in GUI
            Wraper for queue.put_nowait() but with code for modify button '''
        if self.gui:
            self.queue.put_nowait([CODE_MODIFY_BUTTON, args, kwargs])

    def progress_bar(self, *args, **kwargs):
        ''' Modifies progress bar in GUI
            Wraper for queue.put_nowait() but with code for update progress bar'''
        if self.gui:
            self.queue.put_nowait([CODE_PROGRESS, args, kwargs])

    def showResults(self, *args, **kwargs):
        ''' Wraper for queue.put_nowait() to show results'''
        if self.gui:
            self.queue.put_nowait([CODE_SHOWRESULTS, args, kwargs])
        else:
            print('Check results in output folder')

    def processGraph(self, graph_object):
        "Process graphs. If in GUI mode, GUI thread processes."
        if self.gui:
            self.queue.put_nowait([CODE_GRAPHS, graph_object, {}])
        else:
            graph_object()

    def put_message(self, msg):
        """ A wrapper to queue.put_nowait() """
        self.queue.put_nowait(msg)

    def get_message(self):
        """ A wrapper to queue.get_nowait() """
        return self.queue.get_nowait()

    def processor(self, target, name=None):
        varsdic = self.simargs
        if self.logger.logFile is None:
            self._initiateLogFile()
        if self.gui:
            # Start in a different thread to avoid a non-responding gui
            processor_thread = threading.Thread(None, self.wrap_processor, name, (target, varsdic))
            processor_thread.start()
        else:
            print(name)
            self.wrap_processor(target, varsdic)

    def wrap_processor(self, target, varsdic):
        ''' A wrapper for the manager.model '''
        try:
            target(**varsdic)
        except (Exception) as e:
            self.LOG.error(e)
            print('error:', e)

