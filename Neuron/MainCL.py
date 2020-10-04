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

import argparse

# Text Variables
HELP_FILE_DATA   = "File containing the species presence/absence " +\
                   "(1 or 0) and coordinates (file format: presence;X;Y)."
HELP_DIR_RASTERS = "Directory containing t.he independent variables " +\
                   "as ASCII Rasters files."
HELP_OUT_DIR     = "Directory to write all results."
HELP_METHOD      = "Method for resampling the data. 1 - Random " +\
                   "repetition (default), 2 - K-fold Cross " +\
                   "validation, 3 - Bootstrapping."
HELP_REPETITIONS = "Number of repetitions (default is 5). If the " +\
                   "subsampling is K-fold cross validation, it is " +\
                   "the number of folds to divide the data."
HELP_BSIZE       = "Percentage of the resample size for each " +\
                   "repetition. Default is 100%%, meaning that each "+\
                   "resample has the same number of elements as the "+\
                   "original dataset."
HELP_ITERREPORT  = "Number of iterations when the network error is " +\
                   "reported (default 250)."
HELP_ITERINTER   = "Number of internal iterations between reporting "+\
                   "iterations (default 10; Total iterations = ir*ii)."
HELP_BURNIN      = "Number of begining iterations to discard."
HELP_LRATE       = "Learning rate value (default 0.1)."
HELP_MOMENTUM    = "Momentum (default 0.1)."
HELP_HIDDENLYRS  = "Structure of the hidden layers (default 3). " +\
                   "Number of neurons per layer separated by commas."
HELP_APRATIO     = "Pseudo absences ratio. The default value of 1 " +\
                   "creates as many pseudo absences as presences."
HELP_PERCENTAGE  = "Percentage of the data to test the network."
HELP_AUCFILTER   = "Logical value for filtering iterations with AUC " +\
                   "value. Default is 0, which means that AUC is not used."
HELP_AUCTRAIN    = "AUC threshold for training results."
HELP_AUCTEST     = "AUC thresholf for test results."
HELP_PROJECT_DIR = "Directory with rasters to create a projection."
HELP_ONLY_PROJ   = "Logical value for building models or only to " +\
                   "project already built models using a different " +\
                   "set of rasters. Default is 0: models are built " +\
                   "before projection."

TITLE           = "SIMAPSE- Simulation Maps for Ecological Niche " +\
                  "Modelling. Version 1.1"
CITATION        = "Please cite: Tarroso, P, Carvalho, S & Brito, " +\
                  "JC (2012) Simapse - Simulation Maps for " +\
                  "Ecological Niche Modelling. Methods in Ecology " +\
                  "and Evolution, doi:10.1111/j.2041-210X.2012.00210.x" 


startGUI = False

parser = argparse.ArgumentParser(description=TITLE, epilog=CITATION)
parser.add_argument("file_data", nargs='?',
                    help=HELP_FILE_DATA)
parser.add_argument("dir_rasters", nargs='?',
                    help=HELP_DIR_RASTERS)
parser.add_argument("out_dir", nargs='?',
                    help=HELP_OUT_DIR)
parser.add_argument("-s", "--method", type=int, choices=[1, 2, 3], 
                    default=1, help=HELP_METHOD)
parser.add_argument("-r", "--repetitions", type=int, default = 5,
                    help=HELP_REPETITIONS)
parser.add_argument("-bs", "--bsize", type=int, default = 100,
                    help=HELP_BSIZE)
parser.add_argument("-ir", "--iterreport", type=int, default = 250, 
                    help=HELP_ITERREPORT)
parser.add_argument("-ii", "--iterinter", type=int, default = 10, 
                    help=HELP_ITERINTER)
parser.add_argument("-b", "--burnin", type=int, default=25, 
                    help=HELP_BURNIN)
parser.add_argument("-lr", "--lrate", type=float, default = 0.1, 
                    help=HELP_LRATE)
parser.add_argument("-m", "--momentum", type=float, default = 0.1,
                    help=HELP_MOMENTUM)
parser.add_argument("-hl", "--hiddenlyrs", type=str, default = "3",
                    help=HELP_HIDDENLYRS)
parser.add_argument("-pa", "--apratio", type=float, default = 1,
                    help=HELP_APRATIO)
parser.add_argument("-tp", "--percentage", type=int, default = 25,
                    help=HELP_PERCENTAGE)
parser.add_argument("-af", "--aucfilter", type=bool, default=False,
                    help=HELP_AUCFILTER)
parser.add_argument("-atr", "--auctrain", type=float, default=0.800, 
                    help=HELP_AUCTRAIN)
parser.add_argument("-ate", "--auctest", type=float, default=0.800, 
                    help=HELP_AUCTEST)
parser.add_argument("-prj", "--project_dir", type=str, 
                    help=HELP_PROJECT_DIR)
parser.add_argument("-p", "--only_project", type=bool, default=False,
                    help=HELP_ONLY_PROJ)
args = parser.parse_args()

if [args.file_data, args.dir_rasters, args.out_dir] == [None, None, None]:
    startGUI = True

def startCL(conn):
    if args.method == 1:
        args.method = "Random repetition"
    elif args.method == 2:
        args.method = "Cross validation"
    elif args.method == 3:
        args.method = "Bootstrapping"

    print(args.method)
    conn.simargs.update(args.__dict__)
    conn.processor(conn.manager.read_all, 'read rasters')

    if not args.only_project:
        conn.processor(conn.manager.model, 'run model')
        conn.processor(conn.manager.results, 'creating results')

    if args.project_dir:
        conn.processor(conn.manager.project, 'project')
    

