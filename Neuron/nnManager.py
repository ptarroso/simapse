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

from os import path, mkdir, remove, listdir
import threading

import nnFuncs
from nnEngine import NN, savenet, loadnet, sigm, dsigm
from nnRecorder import recorder, htmlreport
import nnGraphs


### SOME VARIABLES ###
CROSS_METHOD = "Cross validation"
BTSTRP_METHOD = "Bootstrapping"
RANDOM_METHOD = "Random repetition"
report = htmlreport()
EXTRA_MODULES = nnGraphs.EXTRA_MODULES

### SOME MESSAGES ###
SUMMARY_MSG   = '\n Summary:' +\
                '\n   Number of models: %s' +\
                '\n   Re-sampling technique: %s' +\
                '\n   Neural Network: %s (%s momentum and %s learning rate)' +\
                '\n   Iterations = %s (%s internal and %s reports)' 
SUMMARYAUC_MSG= '\n   AUC threshold: %s for train and %s for test'
CHOSEN_ITER   = "Iteration chosen: %s"
SOME_FAIL_AUC = "\nNot all the models could meet the AUC theshold. Those " +\
                "models will be removed from the final results. To try to " +\
                "avoid this in the future, change the hidden structure, " +\
                "increase the number of neurons, or add more iterations if " +\
                "the network doesn\'t seem fully trained. The AUC threshold" +\
                " may also be set to very high values."
ALL_FAIL_AUC  = "\nNone of the models could meet the AUC " +\
                "theshold. To try to avoid this in the future, change the " +\
                "hidden structure, increase the number of neurons, or add " +\
                "more iterations if the network doesn't seem fully trained." +\
                " The AUC threshold may also be set to very high values."
SINGLE_MODEL  = "Single model produced."
MODEL_FAIL    = '\nModel %s failed to achieve the AUC threshold.'
NFOLDS_ERR    = "ATENTION! The number of folds cannot exceed the number of " +\
                "available points!"
FOLDS_MSG     = "Each fold contains %s points."
NPOINTS_MSG   = "%s points for train and %s points for test."
READ_MSG      = "Reading data..."
KFOLDSMT_MSG  = "K-fold Cross Validation"
BTSRTP_MSG    = "Bootstraps (Subsets are %s%% of the original size)"
RAND_MSG      = "Random subsampling"
SENSIT_MSG    = "Sensitivity analysis for model %s"
FINALMAPS_MSG = "\nPress Results to produce the final maps."
NETAUCINF_MSG = "Net %3s -> Train: error - %2.3f AUC - %2.3f | Test: error - %2.3f AUC - %2.3f"
NETINF_MSG    = "Net %3s -> Train: error - %2.3f | Test: error - %2.3f"
PREPRSLT_MSG  = '\nPreparing final results...'
READPROJ_MSG  = "Reading project rasters..."
HINTS_MSG     = "\nHints for learning rate value:"
BACKUP_MSG    = 'Backup of old data in output folder done!'
MODELNO_MSG   = "\nModel no. %s"
READDONE_MSG  = "\nReading files done!"
BURNIN_MSG    = "Burn in phase..."
CALCVARS_MSG  = "Creating plots for %s"
CALCAUC_MSG   = "Calculating ROC, PR and AUCs for final averaged model"
SHOWMAPS_MSG  = "Opening maps window"
COMPMODEL_MSG = "Computing model %s of %s (%s fail)"
PROJECT_MSG   = "Projecting model with network %s"
NO_MOD_MSG    = "\nNo extra modules found. Please install the modules to " +\
                "create and display the images or check the text files of " +\
                "the results in the results folder."

def filematch(pathname, ext, fullpath=True):
    files = listdir(pathname)
    if type(ext) == str:
        ext = [ext] 
    match = [x for x in files if x.split('.')[-1] in ext]
    if match and fullpath:
        match = [pathname + path.sep + x for x in match]
    return match          

def varname(rasterfile):
    '''Extracts the name of a file without extension
       from a full path string'''
    filename = path.basename(rasterfile)
    var = filename.split('.')[0]
    return var

class Manager():
    def __init__(self, conn):
        self.conn = conn
        nnFuncs.conn = conn

    def read_all(self, dir_rasters, file_data, percentage, out_dir = None,
                 repetitions = None, method = None, apratio = 1, **kwargs):
        '''Reads all rasters from raster directory and extracts data from text
           file containing Presences/Absences. Also generates Pseudo Absences
           when needed.

            dir_rasters - Directory where the rasters can be found
            file_data   - Text file with presence data (and absence data
                          if available)
            percentage  - percentage value to divide available points in
                          train and test datasets
            repetitions - Run repetitions
            method      - Data subset method
            apratio     - Absence/Presence ratio to create pseudoabsences
                          when needed

           NOTE: It will allways standardize the variables if the directory
                 'standardvars' is not found inside rasters directory.'''
        #Disable buttons while reading
        self.conn.modify_button('disable', 'all')

        #Standardize check
        #TODO: check if the number and name of rasters in std folder are correct
        std_dir = dir_rasters + '/standardvars'
        if path.isdir(std_dir):
            standard = False
            dir_rasters = std_dir
        else:
            standard = True
            mkdir(std_dir)
            
        rasters, raster_values, rasterstats, header = self.read_rasters(dir_rasters, standard)
    
        self.conn.display_msg(READ_MSG)

        # Initializes nnFuncs.spatial_functions with data from the first raster
        spfuncs = nnFuncs.spatial_functions(*header)

        allData, allVariables, DataCoords = spfuncs.ExtractValues(file_data, raster_values, rasters, apratio, out_dir)

        self.conn.abundvar = spfuncs.abundance
        self.totaldata = (allData, allVariables, DataCoords)
        self.testNumb = testNumb = int(len(allData) * (percentage / 100.00))
        self.trainNumb = trainNumb = len(allData) - testNumb
        self.rasters = rasters
        self.raster_values = raster_values
        self.rasterstats = rasterstats

        #Checks some possible errors and add messages to gui
        # ESTE ERRO PODE PASSAR PARA UM ERROR HANDLER!!!
        try:
            if method == CROSS_METHOD:
                folds = int(len(allData) / repetitions)
                if folds == 0:
                    raise Exception(NFOLDS_ERR)
                else:
                    self.conn.display_msg(FOLDS_MSG % folds)
            else:
                self.conn.display_msg(NPOINTS_MSG % (trainNumb, testNumb))
            self.conn.modify_button('normal','all')

        except Exception, msg:
            self.conn.display_msg(str(msg))
            self.conn.modify_button('normal', ['READ', 'HINT', 'METHOD', 'OPTION'])
        if standard:
            spfuncs.savecache(raster_values, rasterstats, std_dir)

        self.conn.display_msg(READDONE_MSG)
        self.spfuncs = spfuncs

    def read_rasters(self, dir_rasters, standard=True):
        '''Reads all rasters from directory 'dir_rasters'.
           The output is a list of rasters, a dictionary of raster values and a
           dictionary of raster statistics (max, min, mean, mean and standard
           deviation) and a raster header

                dir_rasters - Directory where rasters can be found
                standard    - Boolean for standardize the raster values (the
                              raster_values will be standardized. If it is 
                              False then it will look for rasters and stats
                              files.'''
        rasters_list = filematch(dir_rasters, 'txt')
        rasters_list.sort()
        #TODO Add an error of 'no rasters in directory'
        header = nnFuncs.read_ascii(rasters_list[0], 2)
        na = header[-1]
        raster_values, rstats = {}, {}
        rasters = []
        if type(standard) == dict:
            rstats = standard

        for rasterfile in rasters_list:
            raster = path.basename(rasterfile)[:-4]
            temp_raster = nnFuncs.read_ascii(rasterfile)
            if standard == False:
                stsfile = '%s/%s.sts' % (dir_rasters, raster)
                rstats[raster] = nnFuncs.readstats(stsfile)
            elif standard == True:
                temp_list = [value for x in temp_raster for value in x if value <> na]
                N =  float(len(temp_list))
                avg = sum(temp_list) / N
                std = (sum([(x - avg)**2 for x in temp_list]) / N)**0.5
                rmax, rmin = max(temp_list), min(temp_list)
                temp_raster, stmax, stmin = nnFuncs.stvar(temp_raster, avg, std, na)
                rstats[raster] = (avg, std, rmax, rmin, stmax, stmin)
            elif type(standard) == dict:
                avg, std, rmax, rmin, stmax, stmin = rstats[raster]
                temp_raster, nstmax, nstmin = nnFuncs.stvar(temp_raster, avg, std, na)
            raster_values[raster] = temp_raster   
            rasters.append(raster)

        return rasters, raster_values, rstats, header

    def model(self,  hiddenlyrs, lrate, momentum, iterinter, out_dir, 
              method = RANDOM_METHOD, repetitions = 10, aucfilter = False,
              percentage = 50, iterreport = 1, **kwargs):
        '''Neural Networks with repetition

            hiddenlyrs   - Number of Hidden Layers in the network
            lrate        - Learning Rate parameter for Neural Networs
            momentum     - Momentum paramenter for Neural Networks
            iterinter    - Iteration parameter number for Neural Networks
            out_dir      - Output directory for models
            method       - method for subseting data 
            repetitons   - Repetions of Neural Networks (No of models)
            aucfilter    - True for Neural Networks with AUC calculations
                           False for Neural Networks without AUC calculation
            percentage   - percentage of training and test data
            iterreport   - No of AUC reports

            kwargs:
            auctrain     - Value of AUC threshold for train data
            auctest      - Value of AUC threshold for test data'''

        self.conn.modify_button('disable', 'all')

        #Some needed variables
        rasters = self.rasters
        raster_values = self.raster_values
        ninputs = len(rasters)
        modelFiles = self.modelFiles = []
        showmsg = self.conn.display_msg

        #Initiate graphs and results logs
        varprof     = nnFuncs.profiler(self.rasterstats, rasters)
        profile_log = recorder(out_dir, 'profile', rasters)
        varsur_log  = recorder(out_dir, 'varsurface', rasters)
        pderiv_log  = recorder(out_dir, 'PaD', rasters)
        result_log  = recorder(out_dir, 'results_data')

        #Add headers to logs
        result_log.addheader(rasters + ['Chosen_net', 'TrainError', 'TrainAuc', 'TestError', 'TestAUC'])
        headers_pd = nnFuncs.transpose(self.totaldata[1])
        for r_index in xrange(ninputs):
            raster = rasters[r_index]
            header = varprof.rvars[raster]
            profile_log.addheader(header, raster)
            varsur_log.addheader(header, raster)
            header = nnFuncs.unstd(headers_pd[r_index], self.rasterstats[raster])
            pderiv_log.addheader(header, raster, headertxt='Variables')

        #Cleanup
        del header, headers_pd, r_index         

        ### Creates the network ###
        NeuralShape = []
        NeuralShape.append(ninputs)
        if len(hiddenlyrs) <> 0:
            for item in  hiddenlyrs.split(','): NeuralShape.append(int(item))
        NeuralShape.append(1) # One output only

        net = NN(NeuralShape, iterations=iterinter, LR=lrate, momentum=momentum, 
                 verbosity=0)
       
        ### Creates repeated networks to produce n models ###
        allData, allVariables, DataCoords = self.totaldata
        subsets = nnFuncs.subsets(allData, allVariables, DataCoords)
    
        if method == CROSS_METHOD:
            repmethod = subsets.kfoldData(repetitions)
            samplingName = KFOLDSMT_MSG
        elif method == BTSTRP_METHOD:
            bsize = kwargs['bsize']
            repmethod = subsets.bootstrapData(bsize, repetitions, percentage)
            samplingName = BTSRTP_MSG % (bsize)
        elif method == RANDOM_METHOD:
            repmethod = subsets.repeatData(percentage, repetitions)
            samplingName = RAND_MSG

        #Displays summary
        msg = SUMMARY_MSG % (repetitions, samplingName, NeuralShape, momentum,
                            lrate, iterinter*iterreport, iterinter, iterreport)
        if aucfilter:
            auctrain, auctest = kwargs['auctrain'], kwargs['auctest']
            msg += SUMMARYAUC_MSG % (auctrain, auctest)
        showmsg(msg)

        #Keeps track of the networks repetitions that did not achieve proposed AUC values
        self.failed = 0
        
        for self.rep in xrange(1, repetitions + 1):
            self.conn.display_msg(MODELNO_MSG % (self.rep))

            #show progress bar
            msg = COMPMODEL_MSG % (self.rep, repetitions, self.failed)
            self.conn.progress_bar(self.rep-1, repetitions, msg=msg)

            targets, inputs, targetsTest, inputsTest = repmethod.next()
            #Prepares the net with random weights and burnin
            net.rndWeights()
            if 'burnin' in kwargs:
                net = self.burnin(net, inputs, targets, kwargs['burnin'])

            if aucfilter:
                self.repnet(net, inputs, targets, inputsTest, targetsTest,
                            iterreport, auctrain, auctest)
                if len(self.values) == 0:
                    self.failed += 1
                    continue
                    #TODO: Should update the progress bar!

            else:
                self.repnet(net, inputs, targets, inputsTest, 
                            targetsTest, iterreport)

            self.bestnet()
            self.writeChosenNet(out_dir)
          
            net, details = self.chosennet

            #Sensitivity analysis of the network
            self.conn.display_msg(SENSIT_MSG % (self.rep))
            VarSurfaces, Profiles = {}, {}
            pcounter = 0
            ptotal = ninputs * 2 + len(self.totaldata[1])
            for rst in rasters:
                self.conn.progress_bar(pcounter, ptotal, color='darkgreen')
                VarSurfaces[rst] = varprof.varsurface(rst, net.testnet)
                Profiles[rst] = varprof.onewayprofile(rst, net.testnet)
                pcounter += 1
            deriv = []
            for line in self.totaldata[1]:
                self.conn.progress_bar(pcounter, ptotal, color='darkgreen')
                deriv.append(nnFuncs.reducelist(net.pderiv(line)))
                pcounter += 1

            #Write data to log and calculate variable importance per repetition
            tderiv = nnFuncs.transpose(deriv)
            varimp = []
            for r_index in xrange(ninputs):
                self.conn.progress_bar(pcounter, ptotal, color='darkgreen')
                pderiv_log.write(tderiv[r_index], rasters[r_index], str(self.rep))
                varimp.append(sum([x**2 for x in tderiv[r_index]]))
                pcounter += 1
            result_log.write(varimp + details, name=str(self.rep))
            profile_log.write_levels(Profiles, str(self.rep))
            varsur_log.write_levels(VarSurfaces, str(self.rep))

            #Save model
            file_model = '%s/Model_%s.txt' % (out_dir, self.rep)
            modelFiles.append(file_model)
            self.spfuncs.rastercalc(net.testnet, raster_values, file_model, rasters)

        #Check if there are enough models
        if self._checkModelErrors(repetitions):
            self.conn.modify_button('normal', ['READ', 'RUN', 'HINT', 'METHOD', 'OPTION'])

        else:
            # Finalize the variables logs
            for rst in rasters:
                profile_log.finalize(rst)
                varsur_log.finalize(rst)
                pderiv_log.finalize(rst)
            result_log.finalize()

            self.plog = profile_log
            self.vlog = varsur_log
            self.dlog = pderiv_log
            self.rlog = result_log

            self.conn.display_msg(FINALMAPS_MSG)
            self.conn.modify_button('normal', ['READ', 'RUN', 'RESULTS', 'HINT', 'METHOD', 'OPTION'])

    def _checkModelErrors(self, N):
        '''check if model produced have errors'''       
        fail = self.failed
        if fail == 0:
            return(False)
        if fail > 0 and fail < N-1:
            # Warning: some models failed
            # (assumes aucfilter to fail)
            self.conn.display_msg(SOME_FAIL_AUC)
            return(False)
        elif N-fail == 1:
            # Warning: single model
            self.conn.display_msg(SINGLE_MODEL)
            return(False)
        elif fail == N:
            # Error: all models fail
            # (assumes aucfilter to fail)
            self.conn.display_msg(ALL_FAIL_AUC)
            return(True)


    def burnin(self, net, inputs, targets, burnin):
        '''trains a network for a Burnin Period to adujst the weights'''
        #TODO Adjust to internal iterations
        net.loaddata(inputs, targets)
        iterations = net.iterations
        net.iterations = 1
        self.conn.display_msg(BURNIN_MSG)
        for i in xrange(burnin):
            net.trainnet(0)
        net.iterations = iterations
        return net

    def repnet(self, net, inputs, targets, inputsTest, targetsTest,
               iterreport = 1, auctrain = None, auctest = None):
        '''Trains and tests Neural Networks based on the number of internal
           iterations and AUC/error reports (when needed)'''
        values, nets = [], []
        calculateAUC = False
        if auctrain <> None and auctest <> None: calculateAUC = True

        if calculateAUC:
            # Get a plain list of real values for roc
            real = [x for line in targets for x in line]
            realTest =[x for line in targetsTest for x in line]

        for i in xrange(iterreport):
            net.loaddata(inputs, targets)
            net.trainnet(0)

            output = net.testnet(inputs, 0)
            outputTest = net.testnet(inputsTest, 0)

            error = net.neterror(inputs, targets)[0]
            errorTest = net.neterror(inputsTest, targetsTest)[0]
            if calculateAUC:
                #Get plain list of predicted values for AUC
                pred = [x for line in output for x in line]                    
                predTest = [x for line in outputTest for x in line]
               
                roc = nnFuncs.roc(real, pred)
                auc = roc.auc()

                rocTest = nnFuncs.roc(realTest, predTest)
                aucTest = rocTest.auc()

                if auc >= auctrain and aucTest >= auctest:
                    values.append([i, error, auc, errorTest, aucTest])
                txt = NETAUCINF_MSG % (i, error, auc, errorTest, aucTest)
                self.conn.display_msg(txt)

            else:
                values.append([i, error, 0, errorTest, 0])
                txt = NETINF_MSG % (i, error, errorTest)
                self.conn.display_msg(txt)
            
            nets.append(net)

        self.values = values
        self.nets   = nets

    def bestnet(self):
        ''' Finds best net and deletes the others. The best net is 
            defined by the sorting value (test error).'''
        best = None
        values = self.values
        nets   = self.nets
        if len(values) > 0:
            #sort by the index 3 (4th element) of the list
            #TODO Add different sorting options
            order = sorted(range(len(values)), key =lambda x: values[x][3])
            values = [values[x] for x in order]
            nets = [nets[x] for x in order]
            best = [nets[0], values[0]]
            self.conn.display_msg("\n" + CHOSEN_ITER % best[1][0])

        self.chosennet = best

    def writeChosenNet(self, outdir):
        '''Writes the chosen net in the output folder.'''
        net, detail = self.chosennet
        filenet = '%s/net%s_rep%s.net' % (outdir, detail[0], self.rep)
        savenet(net, filenet)

    def results(self, repetitions, out_dir, **kwargs):
        '''Processes all results to produce final maps of average and standard
           deviation. Also produces the final graphs for profiles, partial 
           derivatives and variable surfaces'''
        self.conn.display_msg(PREPRSLT_MSG)
        self.conn.modify_button('disable', 'all')

        if EXTRA_MODULES == False:
            msg = NO_MOD_MSG
            self.conn.display_msg(msg)
            self.conn.modify_button('normal', 'all')
            return
        
        graphs       = None #nnFuncs.graphs()
        ninputs      = len(self.rasters)
        raster_names = self.rasters
        plog         = self.plog
        vlog         = self.vlog
        dlog         = self.dlog
        rlog         = self.rlog
        spfuncs      = self.spfuncs
        ptotal       = ninputs + 2
        pcounter     = 0

        #Creates variable plots
        for rst in raster_names:
            msg = CALCVARS_MSG % rst
            self.conn.progress_bar(pcounter, ptotal, msg=msg)
            varsur = vlog.getmem(rst)[0]
            header = vlog.getheader(rst)

            varsur = nnGraphs.throwgraph(nnGraphs.varsur2d, varsur,
                                         header, rst, out_dir)
            self.conn.processGraph(varsur)

            sctplot = nnGraphs.throwgraph(nnGraphs.scatterplot,
                                          dlog.getheader(rst),
                                          dlog.getmem(rst)[0],
                                          rst, out_dir,
                                          std=dlog.getmem(rst)[1])
            self.conn.processGraph(sctplot)

            lnplot = nnGraphs.throwgraph(nnGraphs.linesplot, 
                                         plog.getheader(rst), 
                                         plog.getmem(rst)[0], 
                                         rst, out_dir, 
                                         std=plog.getmem(rst)[1])
            self.conn.processGraph(lnplot)
            pcounter += 1

        # Calculate ROC/PR & AUC for the average model with total data
        self.conn.progress_bar(pcounter, ptotal, msg=CALCAUC_MSG)
        average, std = spfuncs.modelstats(self.modelFiles, out_dir)
        real = [x for line in self.totaldata[0] for x in line]
        pred = spfuncs.ExtractValues(self.totaldata[2], {'average':average})

        if self.conn.abundvar:
            plotvalues = [real, pred]
        else:
            roc = nnFuncs.roc(real, pred)
            plotvalues = roc.process_all()
        pcounter += 1
        # Show maps
        self.conn.progress_bar(pcounter, ptotal, msg=SHOWMAPS_MSG)
        Maps = nnGraphs.throwgraph(nnGraphs.MapsGraph, average, std, 
                                   "Results_map", out_dir, spfuncs.nodata)
        self.conn.processGraph(Maps)

        Analysis = nnGraphs.throwgraph(nnGraphs.AnalysisGraph, 
                                       rlog.getavg()[:ninputs],
                                       plotvalues, "Results_variables_roc",
                                       out_dir, dnames = raster_names, 
                                       dstd=rlog.getstd()[:ninputs], 
                                       corr=self.conn.abundvar)
        self.conn.processGraph(Analysis)


        self.conn.showResults(["Results_map.png","Results_variables_roc.png"])
        
        pcounter +=1

        #Make html report
        kwargs['scheme'] = '%s,%s,%s' % (ninputs, kwargs['hiddenlyrs'], 1)
        kwargs['repetitions'] = repetitions
        kwargs['out_dir'] = out_dir
        report.summary(kwargs)
        report.model()
        report.variables(raster_names)
        report.write(out_dir + '/report.html')

        self.conn.progress_bar(pcounter, ptotal)
        self.conn.modify_button('normal', 'all')

    def project(self, out_dir, project_dir, **kwargs):
        '''Projects all saved models by loading the trained neural network
           to the new raster set found in the \'project_dir\'
           Standardization of projection rasters is processed with the values
           average and standard deviation of the original rasters (for train).'''
        #from filelist import filelist
        self.conn.modify_button('disable', 'all')
        self.conn.display_msg(READPROJ_MSG)
        
        networks = filematch(out_dir, 'net')
        N = len(networks)
        prj_rasters = filematch(project_dir, 'txt')
        prj_raster_values = {}
        prjFiles = []

        rasters_prj, raster_values_prj, rasterstats_prj, header_prj = self.read_rasters(project_dir, self.rasterstats)
        spfuncs_prj = nnFuncs.spatial_functions(*header_prj)
        spfuncs_prj.create_nodata_list(raster_values_prj[rasters_prj[0]])

        ncounter = 1
        for network in networks:
            rep = network.split('_')[-1][:-4]
            msg = PROJECT_MSG % rep.split('rep')[-1]
            self.conn.progress_bar(ncounter, N, msg=msg)
            net = loadnet(network)
            file_prj = '%s/Project%s.txt' % (out_dir, rep)
            prjFiles.append(file_prj)
            spfuncs_prj.rastercalc(net.testnet, raster_values_prj, file_prj, self.rasters)
            ncounter += 1

        average_prj, std_prj = spfuncs_prj.modelstats(prjFiles, out_dir, sufix="_prj")

        Maps = nnGraphs.throwgraph(nnGraphs.MapsGraph, average_prj, std_prj, 
                                   "Projection_map", out_dir, spfuncs_prj.nodata)
        self.conn.processGraph(Maps)

        self.conn.showResults(["Results_map.png","Projection_map.png"])
        
        self.conn.modify_button('normal', 'all')

    def hint(self, hiddenlyrs = 3, repetitions = 5, percentage = 50,
                 method = RANDOM_METHOD, iterinter = 25, **kwargs):
        '''Gives a Learning Rate hint based on the network scheme, momentum
           and internal iterations. The results is a percentage of the maximum
           value of error change'''

        self.conn.modify_button('disable', 'all')

        LR = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 
              0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

        ninputs = len(self.rasters)
        NeuralShape = []
        NeuralShape.append(ninputs)
        if len(hiddenlyrs) <> 0:
            for item in  hiddenlyrs.split(','): NeuralShape.append(int(item))
        NeuralShape.append(1) # One output only

        # Create network
        net = NN(NeuralShape, iterations=iterinter)
        
        if repetitions > 5: repetitions == 5 # Maximum allowed of repetitions for hint
        allData, allVariables, DataCoords = self.totaldata
        subsets = nnFuncs.subsets(allData, allVariables, DataCoords)

        hint_list = []
        total = len(LR)*repetitions
        progress = 1
        for learning_rate in LR:
            grad_error = 0
            net.LearningRate = learning_rate

            if method == CROSS_METHOD:
                repmethod = subsets.kfoldData(repetitions)
                samplingName = KFOLDSMT_MSG
            elif method == BTSTRP_METHOD:
                bsize = kwargs['bsize']
                repmethod = subsets.bootstrapData(bsize, repetitions, percentage)
                samplingName = BTSRTP_MSG % (bsize)
            elif method == RANDOM_METHOD:
                repmethod = subsets.repeatData(percentage, repetitions)
                samplingName = RAND_MSG

            for item in xrange(repetitions):
                targets, inputs, ttest, itest = repmethod.next()
                net.rndWeights()
                initial_error = net.neterror(inputs, targets)[0]
                net.trainnet(0)
                end_error = net.neterror(inputs, targets)[0]
                grad_error = grad_error + (initial_error - end_error)
                self.conn.progress_bar(progress, total)
                progress += 1
            repmethod.close()
            if grad_error >= 0:
                hint_list.append(grad_error / repetitions)
            else:
                hint_list.append(0)

        hint = LR[hint_list.index(max(hint_list))]
        hint_percent = map(lambda x: int((x / max(hint_list)) * 75), hint_list)

        i = 0
        msg = ""
        for item in LR:
            value = map(lambda x: x*hint_percent[i], '*')
            msg += "\n%.5f : %s" % (item, value[0])
            i = i + 1
        self.conn.display_msg(HINTS_MSG)
        self.conn.display_msg(msg)
        self.conn.modify_button('normal', ['READ', 'RUN', 'PROJECT', 'HINT', 'METHOD', 'OPTION'])




