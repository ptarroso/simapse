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

from tkinter import Tk, Toplevel, Label, filedialog
import queue as qu
import threading
from math import ceil
from . import Connector

class EMNN():
    def __init__(self, conn):
        conn.gui = True
        queue = self.queue = qu.Queue()
        conn.queue = queue
        self.conn = conn
   
    def startgui(self):
        self.root = root = Tk()
        root.title('SIMAPSE - Simulation Maps for Ecological Niche Modelling')
        root.geometry('695x445+375+115')
        root.tk.call('encoding', 'system', 'utf-8')
        root.configure(bg='#d9d9d9')
        root.resizable(width='false', height='false')
        self.gui()
        self.root.after(100, self.periodicUpdate)
        root.mainloop()

    def gui(self):
        from tkinter import Button, Entry, Frame, Label, Checkbutton, \
                            Scrollbar, Text, IntVar, StringVar, OptionMenu

        self.dir_project = ''
        self.folds = 0
        self.func = ""
        self.aucfilter = IntVar()

        self.lightgray = '#d9d9d9'
        self.darkgray  = '#d3d3d3'

        self.normaltext = ("Helvetica", -10)
        self.boldtext   = ("Helvetica", -10, "bold")
        self.bigtext    = ("Helvetica", -12, "bold")
        self.smalltext  = ("Helvetica", -9)

        self.butData = Button(self.root)
        self.butData.place(in_=self.root,x=5,y=5)
        self.butData.configure(height=1, width=10,
                               text="Data File",
                               font = self.normaltext,
                               highlightbackground = self.lightgray,
                               command=self.dataOpen)

        self.butRasterFolder = Button(self.root)
        self.butRasterFolder.place(in_=self.root,x=5,y=35)
        self.butRasterFolder.configure(height=1, width = 10,
                                       text="Raster Folder",
                                       font = self.normaltext,
                                       highlightbackground = self.lightgray,
                                       command=self.rasterOpen)

        self.butOutFolder = Button(self.root)
        self.butOutFolder.place(in_=self.root,x=5,y=65)
        self.butOutFolder.configure(height=1, width = 10,
                                    text="Out Folder",
                                    font = self.normaltext,
                                    highlightbackground = self.lightgray,
                                    command=self.outOpen)

        self.entData = Entry(self.root)
        self.entData.place(in_=self.root,x=100,y=9)
        self.entData.configure(textvariable="file_data",
                               font = self.normaltext,      
                               width=97,
                               background=self.darkgray,
                               relief="flat",
                               highlightbackground = self.lightgray)
        self.entData.insert(0, '')

        self.entOutFolder = Entry(self.root)
        self.entOutFolder.place(in_=self.root,x=100,y=69)
        self.entOutFolder.configure(textvariable="out_dir",
                                    font = self.normaltext,
                                    width=97,
                                    background=self.darkgray,
                                    relief="flat",
                                    highlightbackground = self.lightgray)
        self.entOutFolder.insert(0, '')

        self.entRasterFolder = Entry(self.root)
        self.entRasterFolder.place(in_=self.root,x=100,y=39)
        self.entRasterFolder.configure(textvariable="dir_rasters",
                                       font = self.normaltext,
                                       width=97,
                                       background=self.darkgray,
                                       relief="flat",
                                       highlightbackground = self.lightgray)
        self.entRasterFolder.insert(0, '')

        self.activeMethod = StringVar(self.root)
        self.activeMethod.set("Random repetition")
        self.butMETHOD = OptionMenu(self.root, self.activeMethod, 
                                    "Random repetition",
                                    "Cross validation",
                                    "Bootstrapping",
                                    command = lambda x:self.displayMethodFrame(x))
        self.butMETHOD.place(in_=self.root, x=4, y=97,
                                 height= "25", width="120")
        self.butMETHOD.configure(font = self.smalltext,
                                 background = self.lightgray) 
        self.displayMethodFrame(self.activeMethod.get())

        self.activeOption = StringVar(self.root)
        self.activeOption.set("Network structure")
        self.butOPTION = OptionMenu(self.root, self.activeOption, 
                                    "Network structure",
                                    "Options",
                                    command = lambda x:self.displayOptionFrame(x))
        self.butOPTION.place(in_=self.root, x=4, y=182,
                                 height= "25", width="120")
        self.butOPTION.configure(font=self.smalltext,
                                 background = self.lightgray) 
        self.displayOptionFrame(self.activeOption.get())

        self.Progress_frame = Frame(self.root)
        self.Progress_frame.place(in_=self.root, x=5, y=423)
        self.Progress_frame.configure(borderwidth="2",
                                      relief='sunken',
                                      height="20",
                                      width="105", 
                                      bg='white')

        self.Progress_bar = Label(self.Progress_frame)
        self.Progress_bar.place(x=0, y=0)
        self.Progress_bar.configure(font = self.smalltext)

        self.Progress_info = Label(self.root)
        self.Progress_info.place(x=110, y=425)
        self.Progress_info.configure(font=self.smalltext,
                                     bg = self.lightgray)

        self.frameButtons = Frame(self.root)
        self.frameButtons.place(in_=self.root,x=5,y=336)
        self.frameButtons.configure(borderwidth="2",
                                    bg = self.lightgray,
                                    relief="raise",
                                    height="84",
                                    width="260")

        self.butREAD = Button(self.root)
        self.butREAD.place(in_=self.frameButtons,x=5,y=5, width = 70)
        self.butREAD.configure(font= self.bigtext,
                               highlightbackground = self.lightgray,
                               text="READ",
                               command=self.read)

        self.butRUN = Button(self.root)
        self.butRUN.place(in_=self.frameButtons,x=80,y=5, width = 70)
        self.butRUN.configure(font= self.bigtext,
                              highlightbackground = self.lightgray,
                              text="RUN",
                              command=self.returnVar,
                              state='disabled')

        self.butRESULTS = Button(self.root)
        self.butRESULTS.place(in_=self.frameButtons,x=155,y=5, width = 80)
        self.butRESULTS.configure(font= self.bigtext,
                                  highlightbackground = self.lightgray,
                                  text="RESULTS",
                                  command=self.returnResults,
                                  state='disabled')

        self.butPROJECT = Button(self.root)
        self.butPROJECT.place(in_=self.frameButtons,x=5,y=45, width = 70)
        self.butPROJECT.configure(font= self.boldtext,
                                  highlightbackground = self.lightgray,
                                  text="PROJECT",
                                  command=self.project,
                                  state='disabled')

        self.frameText = Frame(self.root)
        self.frameText.place(in_=self.root,x=270,y=100)
        self.frameText.configure(height=320, width=400)

        self.scrollbar = Scrollbar(self.root)

        self.textFrame = Text(self.frameText, wrap='word')
        self.textFrame.configure(font= self.normaltext,
                                 height="24",
                                 width="65")
        self.textFrame.place(x=0,y=0)
        self.textFrame.configure(yscrollcommand = self.scrollbar.set)

        self.scrollbar.configure(command=self.textFrame.yview,
                                 highlightbackground = self.lightgray)
        self.scrollbar.place(x=675,y=100, height=320)

    def displayMethodFrame(self, method):
        ''' Displays individual frames for the subsetting method'''
        from tkinter import Button, Entry, Frame, Label

        if 'frameMethodSelection' in self.__dict__: self.frameMethodSelection.destroy()

        self.varupdate()
        c = self.conn.simargs

        self.frameMethodSelection = Frame(self.root)
        self.frameMethodSelection.place(in_=self.root,x=5,y=122)
        self.frameMethodSelection.configure(borderwidth="2", 
                                            relief="raise",
                                            height="60",
                                            width="260",
                                            bg = self.lightgray)

        if method == "Random repetition":
            self.labRepetitions = Label(self.root)
            self.labRepetitions.place(in_=self.frameMethodSelection,x=2,y=10)
            self.labRepetitions.configure(font = self.normaltext,
                                          borderwidth="1", 
                                          justify='left',
                                          anchor= 'e',
                                          bg = self.lightgray,
                                          text="Number of repetitions:") 

            self.entRepetitions = Entry(self.root)
            self.entRepetitions.place(in_=self.frameMethodSelection,x=125,y=10)
            self.entRepetitions.configure(textvariable="repetitions", width="7",
                                          font = self.normaltext,
                                          highlightbackground = self.lightgray)
            self.entRepetitions.delete(0, 'end')
            self.entRepetitions.insert('end', c['repetitions'])
            
        elif method == "Cross validation":
            self.labRepetitions = Label(self.root)
            self.labRepetitions.place(in_=self.frameMethodSelection,x=2,y=10)
            self.labRepetitions.configure(font = self.normaltext, 
                                          bg = self.lightgray,
                                          text="Number of folds:") 

            self.entRepetitions = Entry(self.root)
            self.entRepetitions.place(in_=self.frameMethodSelection,x=100,y=10)
            self.entRepetitions.configure(textvariable="repetitions",
                                          width="7",
                                          highlightbackground = self.lightgray,
                                          font = self.normaltext)
            self.entRepetitions.delete(0, 'end')
            self.entRepetitions.insert('end', c['repetitions'])

        elif method == "Bootstrapping":
            self.labRepetition = Label(self.root)
            self.labRepetition.place(in_=self.frameMethodSelection,x=2,y=5)
            self.labRepetition.configure(borderwidth="1",
                                         text="Number of Bootstraps:",
                                         bg = self.lightgray,
                                         font = self.normaltext)

            self.entRepetitions = Entry(self.root)
            self.entRepetitions.place(in_=self.frameMethodSelection,x=125,y=5)
            self.entRepetitions.configure(textvariable="repetitions", 
                                          width="7",
                                          highlightbackground = self.lightgray,
                                          font = self.normaltext)
            self.entRepetitions.delete(0, 'end')
            self.entRepetitions.insert('end', c['repetitions'])

            self.labBsize = Label(self.root)
            self.labBsize.place(in_=self.frameMethodSelection,x=2,y=30)
            self.labBsize.configure(borderwidth="1", 
                                    text="Bootstraps Sample Size:",
                                    bg = self.lightgray,
                                    font = self.normaltext)

            self.entBsize = Entry(self.root)
            self.entBsize.place(in_=self.frameMethodSelection,x=125,y=30)
            self.entBsize.configure(textvariable="Bsize",
                                    width="7",
                                    highlightbackground = self.lightgray,
                                    font = self.normaltext)
            self.entBsize.delete(0, 'end')
            self.entBsize.insert('end', c['bsize'])

    def displayOptionFrame(self, option):
        ''' Displays individual frames for the subsetting method'''
        from tkinter import Button, Entry, Frame, Label, Checkbutton

        if 'frameOptionSelection' in self.__dict__: self.frameOptionSelection.destroy()

        self.varupdate()
        c = self.conn.simargs

        self.frameOptionSelection = Frame(self.root)
        self.frameOptionSelection.place(in_=self.root,x=5,y=207)
        self.frameOptionSelection.configure(borderwidth="2", 
                                            relief="raise", 
                                            height="125",
                                            width="260",
                                            bg = self.lightgray)

        if option == "Network structure":
            self.labMaxiter = Label(self.root)
            self.labMaxiter.place(in_=self.frameOptionSelection,x=190,y=5)
            self.labMaxiter.configure(borderwidth="1",
                                      font = self.normaltext,
                                      text="Internal",
                                      bg = self.lightgray)
            
            self.labNNN = Label(self.root)
            self.labNNN.place(in_=self.frameOptionSelection,x=95,y=5)
            self.labNNN.configure(borderwidth="1",
                                  font = self.normaltext,
                                  text="Reported",
                                  bg = self.lightgray)

            self.labTI = Label(self.root)
            self.labTI.place(in_=self.frameOptionSelection,x=5,y=25)
            self.labTI.configure(borderwidth="1",
                                 font = self.normaltext,
                                 text="Total iterations =",
                                 bg = self.lightgray)

            self.entITERReport = Entry(self.root)
            self.entITERReport.place(in_=self.frameOptionSelection,x=88,y=25)
            self.entITERReport.configure(textvariable="AUCReport",
                                         width="10",
                                         font = self.normaltext,
                                         highlightbackground = self.lightgray)
            self.entITERReport.delete(0, 'end')
            self.entITERReport.insert('end', c['iterreport'])

            self.times = Label(self.root)
            self.times.place(in_=self.frameOptionSelection, x=160, y=25)
            self.times.configure(text="x",
                                 font = self.normaltext,
                                 bg = self.lightgray)

            self.entITERInter = Entry(self.root)
            self.entITERInter.place(in_=self.frameOptionSelection,x=180, y=25)
            self.entITERInter.configure(textvariable="maxiter",
                                        width="10",
                                        font = self.normaltext,
                                        highlightbackground = self.lightgray)
            self.entITERInter.delete(0, 'end')
            self.entITERInter.insert('end', c['iterinter'])

            self.labEta = Label(self.root)
            self.labEta.place(in_=self.frameOptionSelection,x=5,y=55)
            self.labEta.configure(borderwidth="1",
                                  font = self.normaltext,
                                  text="Learning Rate",
                                  bg = self.lightgray)

            self.butHINT = Button(self.root)
            self.butHINT.place(in_=self.frameOptionSelection,x=65, y=75, 
                               height = 23, width=20)
            self.butHINT.configure(font= self.smalltext,
                                   text="H",
                                   command=self.hint,
                                   state='disabled',
                                   highlightbackground = self.lightgray)

            self.labMomentum = Label(self.root)
            self.labMomentum.place(in_=self.frameOptionSelection,x=88,y=55)
            self.labMomentum.configure(borderwidth="1",
                                       font = self.normaltext,
                                       text="Momentum",
                                       bg = self.lightgray)

            self.entLRATE = Entry(self.root)
            self.entLRATE.place(in_=self.frameOptionSelection,x=5,y=75)
            self.entLRATE.configure(textvariable="eta",
                                    width="8",
                                    font = self.normaltext,
                                    highlightbackground = self.lightgray)
            self.entLRATE.delete(0, 'end')
            self.entLRATE.insert('end', c['lrate'])

            self.entMomentum = Entry(self.root)
            self.entMomentum.place(in_=self.frameOptionSelection,x=90,y=75)
            self.entMomentum.configure(textvariable="momentum",
                                       width="8",
                                       font = self.normaltext,
                                       highlightbackground = self.lightgray)
            self.entMomentum.delete(0, 'end')
            self.entMomentum.insert('end', c['momentum'])

            self.labNNS = Label(self.root)
            self.labNNS.place(in_=self.frameOptionSelection,x=165,y=55)
            self.labNNS.configure(borderwidth="1",
                                  font = self.normaltext,
                                  text="Hidden Layers",
                                  bg = self.lightgray)

            self.entNNShape = Entry(self.root)
            self.entNNShape.place(in_=self.frameOptionSelection,x=160,y=75)
            self.entNNShape.configure(textvariable="HiddenLyr",
                                      width="14",
                                      font = self.normaltext,
                                      highlightbackground = self.lightgray)
            self.entNNShape.delete(0, 'end')
            self.entNNShape.insert('end', c['hiddenlyrs'])
            
        elif option == "Options":

            self.labPercentage = Label(self.root)
            self.labPercentage.place(in_=self.frameOptionSelection,x=2, y=5)
            self.labPercentage.configure(borderwidth="1", 
                                         text="Test %:",
                                         font = self.normaltext,
                                         bg = self.lightgray)

            self.entPercentage = Entry(self.root)
            self.entPercentage.place(in_=self.frameOptionSelection,
                                     x=45, y=5, 
                                     width = 30)
            self.entPercentage.configure(textvariable="Percentage",
                                         font = self.normaltext,
                                         highlightbackground = self.lightgray)
            self.entPercentage.delete(0, 'end')
            self.entPercentage.insert('end', c['percentage'])

            self.labAPRatio = Label(self.root)
            self.labAPRatio.place(in_=self.frameOptionSelection,x=80, y=5)
            self.labAPRatio.configure(borderwidth="1",
                                      text="Pseudoabsences/Presences:",
                                      font = self.normaltext,
                                      bg = self.lightgray)

            self.entAPRatio = Entry(self.root)
            self.entAPRatio.place(in_=self.frameOptionSelection,x=220,y=5)
            self.entAPRatio.configure(textvariable="apratio", width="4",
                                      font = self.normaltext,
                                      highlightbackground = self.lightgray)
            self.entAPRatio.delete(0, 'end')
            self.entAPRatio.insert('end', c['apratio'])

            self.labBurnin = Label(self.root)
            self.labBurnin.place(in_=self.frameOptionSelection,x=2, y=30)
            self.labBurnin.configure(borderwidth="1",
                                     text="Burn-in iterations:",
                                     font = self.normaltext,
                                     bg = self.lightgray)

            self.entBurnin = Entry(self.root)
            self.entBurnin.place(in_=self.frameOptionSelection,x=90,y=30)
            self.entBurnin.configure(textvariable="burnin", width="8",
                                     font = self.normaltext,
                                     highlightbackground = self.lightgray)
            self.entBurnin.delete(0, 'end')
            self.entBurnin.insert('end', c['burnin'])

            self.chkAucFilter = Checkbutton(self.root)
            self.chkAucFilter.place(in_=self.frameOptionSelection,x=2,y=60)
            self.chkAucFilter.configure(font = self.normaltext,
                                        text="Filter with AUC threshold",
                                        variable=self.aucfilter,
                                        bg = self.lightgray,
                                        command = self.aucstate)

            self.labAUCTrain = Label(self.root)
            self.labAUCTrain.place(in_=self.frameOptionSelection,x=5,y=85)
            self.labAUCTrain.configure(borderwidth="1",
                                       font = self.normaltext,
                                       text="training data",
                                       bg = self.lightgray)

            self.entAUCTrain = Entry(self.root)
            self.entAUCTrain.place(in_=self.frameOptionSelection,x=70,y=85)
            self.entAUCTrain.configure(textvariable="valueAuc",
                                       width="8",
                                       font = self.normaltext,
                                       highlightbackground = self.lightgray)
            self.entAUCTrain.delete(0, 'end')
            self.entAUCTrain.insert('end', c['auctrain'])

            self.labAUCTest = Label(self.root)
            self.labAUCTest.place(in_=self.frameOptionSelection,x=130,y=85)
            self.labAUCTest.configure(borderwidth="1",
                                      font = self.normaltext,
                                      text="testing data",
                                      bg = self.lightgray)

            self.entAUCTest = Entry(self.root)
            self.entAUCTest.place(in_=self.frameOptionSelection,x=195, y=85)
            self.entAUCTest.configure(textvariable="valueAucTest",
                                      width="8",
                                      font = self.normaltext,
                                      highlightbackground = self.lightgray)
            self.entAUCTest.delete(0, 'end')
            self.entAUCTest.insert('end', c['auctest'])
        
            self.aucstate()

    def varupdate(self):
        extract = self.extract
        c = self.conn.simargs

        c['file_data']     = extract('entData', c['file_data'])
        c['dir_rasters']   = extract('entRasterFolder', c['dir_rasters'])
        c['out_dir']       = extract('entOutFolder', c['out_dir'] )
        c['method']        = extract('activeMethod', c['method'])
        c['iterreport']    = int(extract('entITERReport', c['iterreport']))
        c['iterinter']     = int(extract('entITERInter', c['iterinter']))
        c['lrate']         = float(extract('entLRATE', c['lrate']))
        c['momentum']      = float(extract('entMomentum', c['momentum']))
        c['hiddenlyrs']    = extract('entNNShape', c['hiddenlyrs'])
        c['apratio']       = float(extract('entAPRatio', c['apratio']))
        c['percentage']    = int(extract('entPercentage', c['percentage']))
        c['burnin']        = int(extract('entBurnin', c['burnin']))
        c['auctrain']      = float(extract('entAUCTrain', c['auctrain']))
        c['auctest']       = float(extract('entAUCTest', c['auctest']))
        c['repetitions']   = int(extract('entRepetitions', c['repetitions']))
        c['bsize']         = int(extract('entBsize', c['bsize']))
        c['aucfilter']     = bool(extract('aucfilter', int(c['aucfilter'])))

    def extract(self, test, default):
        if test in self.__dict__.keys():
            value = self.__dict__[test].get()
        else:
            value = default
        return value

    def read(self):
        import time
        self.varupdate()
        self.conn.processor(self.conn.manager.read_all, 'read')

    def returnVar(self):
        # Updates variables
        self.varupdate()
        self.conn.processor(self.conn.manager.model, 'run')

    def returnResults(self):
        self.varupdate()
        self.conn.processor(self.conn.manager.results, 'results')

    def project(self):
        self.varupdate()
        project_dir = filedialog.askdirectory()
        self.conn.simargs['project_dir'] = project_dir
        self.conn.processor(self.conn.manager.project, 'project')

    def hint(self):
        self.varupdate()
        self.conn.processor(self.conn.manager.hint, 'hint')

    def dataOpen(self):
        self.entData.delete(0, 'end')
        file_data = filedialog.askopenfilename(filetypes=[("text files","*.txt"), ("allfiles","*")])
        self.entData.insert('end', file_data)
        
    def rasterOpen(self):
        self.entRasterFolder.delete(0, 'end')
        dir_rasters = filedialog.askdirectory()
        self.entRasterFolder.insert('end', dir_rasters)
        
    def outOpen(self):
        self.entOutFolder.delete(0,'end')
        out_dir = filedialog.askdirectory()
        self.entOutFolder.insert('end', out_dir)
        
    def update_text(self, string_txt):
        txt = string_txt + " \n"
        self.textFrame.insert('end', txt)
        self.textFrame.yview('end')

    def processGraph(self, graph_object):
        '''Just a wraper to call the graphics creation object'''
        graph_object()

    def periodicUpdate(self):
        """Executes periodic checks to GUI:
            - if there are new messages and displays when true"""
        try:
            while 1:
                code, args, kwargs = self.queue.get_nowait()
                if code == Connector.CODE_TEXT:
                    self.update_text(*args)
                elif code == Connector.CODE_PROGRESS:
                    self.progress(*args, **kwargs)
                elif code == Connector.CODE_MODIFY_BUTTON:
                    self.modify_but(*args)
                elif code == Connector.CODE_SHOWRESULTS:
                    self.showResults(*args)
                elif code == Connector.CODE_GRAPHS:
                    self.processGraph(args)
                else:
                    self.update_text('Unknown message...')
                self.queue.task_done()
                self.root.update()
        except qu.Empty:
            pass
        self.root.after(100, self.periodicUpdate)

    def modify_but(self, state, buttonlist):
        if buttonlist == 'all':
            buttonlist = ['READ', 'RUN', 'RESULTS', 'PROJECT', 'HINT', 'METHOD', 'OPTION']
        for button in buttonlist:
            but = "self.but%s.configure(state=\'%s\')" % (button, state)
            exec(but)

    def abundance(self):
        self.entAUCTrain.configure(state='disabled')
        self.entAUCTest.configure(state='disabled')

    def aucstate(self):
        if self.aucfilter.get(): state = 'normal'
        else: state = 'disabled'
        self.entAUCTrain.configure(state=state)
        self.entAUCTest.configure(state=state)

    def progress(self, value, maxvalue, **kwargs):
        '''Shows the progress bar in GUI
           args are:
            Value    - Value of the current progress
            MaxValue - Value where progress terminates
           kwargs are
            color    - Color for the progess bar
            msg      - message to display'''
        color = 'blue'
        if 'color' in kwargs: color = kwargs['color']
        percent = ceil((value * 100) / maxvalue)
        self.Progress_bar.configure(font = self.smalltext,
                                    foreground="white",
                                    background=color) #"#0000ff"
        if percent != 100:
            width = int((percent / 100.000) * 20)
            text = '%s%%' % percent
            self.Progress_bar.configure(text = text, width = width)
            if 'msg' in kwargs:
                self.Progress_info.configure(text = kwargs['msg'])
        elif percent == 100:
            self.Progress_bar.configure(text = "",
                                        width = 0,
                                        relief="flat",
                                        background="#ece9d8")
            self.Progress_info.configure(text="")

    def showResults(self, figures):
        from tkinter import Button
        figures = [self.entOutFolder.get() + "/" + x for x in figures]
        ResultsWindow = mytoplevel(self.root, figures, self)
        ResultsWindow.title('Results')
        butRW = Button(ResultsWindow, text = 'CLOSE', command = ResultsWindow.destroy)
        butRW.pack()

    def update_queue(self):
        try:
            while self.queue.qsize():
                msg = '%s\n' % self.queue.get(0)
                self.textFrame.insert('end', msg)
                self.textFrame.yview('end')
        except Queue.Empty:
            pass

class mytoplevel(Toplevel):
    def __init__(self, master, im, gui):
        try:
            Toplevel.__init__(self)
            from PIL import ImageTk
            Toplevel.__init__(self, master)
            self.images = [ImageTk.PhotoImage(file = x) for x in im]
            for image in self.images:
                self.image_label = Label(self, image=image)
                self.image_label.pack()
        except (ImportError) as e:
            msg = "Probably you don't have the PIL module but all the data is saved on the output folder. Try again after installing the Python module."
            gui.textFrame.insert('end', msg)
            gui.textFrame.yview('end')



