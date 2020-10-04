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
import Connector

import MainCL

__author__ = "Pedro Tarroso"
__copyright__ = "2012, Pedro Tarroso"
__credits__   = "Contributors: Silvia Carvalho\
                 Jose Carlos Brito" 
__license__ = "GPL"
__version__ = "1.1"


#Create Connector
conn = Connector.connection()

#Check if there are arguments (Command line or GUI)
if MainCL.startGUI:
    import MainGUI
    GUI = MainGUI.EMNN(conn)
    GUI.startgui()
else:
    MainCL.startCL(conn)
