#!/usr/bin/env python3
'''
SIMAPSE - simulation maps for ecological niche modelling
Version 2.00
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

from os import curdir

from .nnFuncs import transpose, read_ascii, sendout

try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib import figure, cm, ticker
    import numpy as np
    EXTRA_MODULES = True
except ImportError:
    EXTRA_MODULES = False

class throwgraph():
    def __init__(self, func, *args, **kwargs):  
        self.func = func
        self.args = args
        self.kwargs = kwargs
    def __call__(self):
        self.func(*self.args, **self.kwargs)
        

def varsur2d(varsur, labels, varname = None, outdir = None, dpi = 300):
    '''Outputs one variation surface plot per variable'''
    if outdir == None: outdir = curdir

    txt = 'Variation surface: %s' % varname
    N = len(labels)

    #Create Figure and Canvas
    fig  = figure.Figure()
    canvas = FigureCanvasAgg(fig)
    plt = fig.add_subplot(111)
    X = np.array([labels] * N)
    Y = np.array(transpose([range(N)] * N))
    Z = varsur[:]
    Z.reverse()
    im = plt.pcolor(X, Y, Z, cmap=cm.RdYlGn, shading='auto')
    fig.colorbar(im)
    outfile = '%s/%s_varsurface.png' % (outdir, varname)
    fig.savefig(outfile, dpi= dpi)
    fig.clf()


def varsur3d(self, varsur, labels, varname = None, outdir = None, dpi = 300):
    '''Outputs one variation surface plot per variable'''   
    if outdir == None: outdir = curdir

    #Create Figure and Canvas
    fig  = figure.Figure()
    canvas = FigureCanvasAgg(fig)
    ax = Axes3D(fig)
    
    #Create numpy variables  
    X = np.array(labels)
    Y = np.array(range(len(labels)))
    X, Y = np.meshgrid(X, Y)
    Z = np.array(varsur)

    img = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.jet)

    ax.set_xlabel(varname + ' values')
    ax.set_ylabel('Other variables (% of range)')
    ax.set_zlabel('Network output')

    outfile = '%s/%s_varsurface.png' % (outdir, varname)
    fig.savefig(outfile, dpi= dpi)
    fig.clf()

def scatterplot(Xdata, Ydata, varname, outdir = None, dpi = 300, **kwargs):
    '''creates a 2D scatterplot per variable'''
    if outdir == None: outdir = curdir

    #Create figure and canvas
    fig  = figure.Figure()
    canvas = FigureCanvasAgg(fig)
    plt = fig.add_subplot(111)

    std = None
    if 'std' in kwargs.keys():
        std = kwargs['std']

    plt.errorbar(Xdata, Ydata, yerr=std, xerr=None, 
                 fmt='.', c='k', ecolor='lightgrey')

    plt.set_title(varname + " PaD")
    plt.set_ylabel("Predictive effect", fontsize = 10)
    plt.set_xlabel(varname + " values", fontsize = 10)

    #outfile = '%s/%s_PaD.svg' %(outdir, varname)
    outfile = '%s/%s_PaD.png' %(outdir, varname)
    fig.savefig(outfile, dpi= dpi)
    fig.clf()

def linesplot(Xdata, Ydata, varname, outdir = None, dpi = 300, **kwargs):
    '''creates a line plot per variable'''
    if outdir == None: outdir = curdir

    #Create figure and canvas
    fig  = figure.Figure()
    canvas = FigureCanvasAgg(fig)
    plt = fig.add_subplot(111)

    if 'std' in kwargs.keys():
        std = kwargs['std']
        minimum = [x-y for x,y in zip(Ydata,std)]
        maximum = [x+y for x,y in zip(Ydata,std)]
        plt.fill_between(Xdata, minimum, maximum, 
                         facecolor='lightblue', edgecolor='lightblue')
    if 'minmax' in kwargs.keys():
        minimum = kwargs['minmax'][0]
        maximum = kwargs['minmax'][1]
        plt.fill_between(Xdata, minimum, maximum, 
                         facecolor='lightblue', edgecolor='lightblue')           

    plt.plot(Xdata, Ydata)
    plt.axis((min(Xdata), max(Xdata), 0.0, 1.0))

    plt.set_title(varname + " Profile")
    plt.set_ylabel("Predictive effect", fontsize = 10)
    plt.set_xlabel(varname + " values", fontsize = 10)

    #outfile = '%s/%s_profile.svg' %(outdir, varname)
    outfile = '%s/%s_profile.png' %(outdir, varname)
    fig.savefig(outfile, dpi= dpi)
    fig.clf()

def MapsGraph(plotgraph1, plotgraph2, name, 
              outdir = None, nodata = -9999, dpi = 96):
    '''Opens a new plot window with a raster image
       plotgraph1 -> raster file with model to plot
       plotgraph2 -> raster file with standar deviation to plot
       name       -> name for outputfile
       outdir     -> Output folder
       dpi        -> Resolution of output
       nodata     -> No Data value'''
    if outdir == None: outdir = curdir

    if type(plotgraph1) == str: 
        plotgraph1 = read_ascii(plotgraph1, 0)
    if type(plotgraph2) == str:
        plotgraph2 = read_ascii(plotgraph2, 0)

    fig  = figure.Figure()
    fig.set_figheight(2)

    #plots average model
    canvas = FigureCanvasAgg(fig)
    averageplot = fig.add_subplot(121)

    temp1 = [x for y in plotgraph1 for x in y if x != nodata]
    max_value = max(temp1)
    min_value = min(temp1)

    im = averageplot.imshow(plotgraph1, interpolation='bilinear', 
                           cmap=cm.RdYlGn, vmin=min_value)

    fig.colorbar(im)
    averageplot.set_title('Average Model', position = (0.5,0.962))
    averageplot.axis('off')

    #plots standard deviation model

    temp2 = [x for y in plotgraph2 for x in y if x != nodata]
    max_value = max(temp2)
    min_value = min(temp2)

    stdplot = fig.add_subplot(122)

    im2 = stdplot.imshow(plotgraph2, interpolation='bilinear', 
                         cmap=cm.RdYlGn, vmin=min_value)
    stdplot.axis('off')
    fig.colorbar(im2)

    stdplot.set_title('Standard Deviation', position = (0.5,0.968))

    filename = '%s/%s.png' % (outdir, name)
    fig.savefig(filename, dpi= dpi)
    fig.clf()

def AnalysisGraph(derivatives, plotvalues, name,
                  outdir = None, dpi=96, corr=False, **kwargs): 
    '''Plots derivatives and final roc plot'''
    if outdir == None: outdir = curdir
    fig  = figure.Figure()
    fig.set_figheight(4)
    fig.subplots_adjust(bottom= 0.25, wspace=0.3)
    
    canvas = FigureCanvasAgg(fig)
    varBars = fig.add_subplot(121)

    #Some options
    dstd = None
    if 'dstd' in kwargs:
        dstd = kwargs['dstd']
    dnames = None
    if 'dnames' in kwargs:
        dnames = kwargs['dnames']
    
    position = range(len(derivatives))
    width = 0.9
    barplot = varBars.bar(position, derivatives, width,
                          color='b', yerr=dstd, ecolor='k')
    varBars.set_ylabel('Sum of squared partial derivatives', fontsize = 10)
    varBars.set_title('Variable contribution', fontsize = 12)
    varBars.set_xticks([x + width/2 for x in position]) 
    varBars.set_xticklabels(dnames, rotation = 90, fontsize = 8)
    varBars.set_xlim(-width,len(position))
    
    # Plots final Roc and precision/recall curve with
    # respective Aucs with the final averaged prediction map
    secondPlot = fig.add_subplot(122)
    if corr: 
        corrplot(secondPlot, plotvalues)
    else: 
        rocplot(secondPlot, plotvalues)

    fig.savefig(outdir + "/" + name + ".png", dpi= dpi)
    #fig.savefig(outdir + "/" + name + ".svg", dpi= 300)
    fig.clf()
    
    #Displays the exact value for each variable importance
    bar_labels = map(lambda x,y,z: "%s = %.4f (%.4f)" % (x,y,z), dnames, derivatives, dstd)
    sendout('put', "\nSum of squared partial derivatives and standard deviation by variable")
    for item in bar_labels:
        sendout('put', item)

def rocplot(fig, rocvalues):
    roc, pr, aucROC, aucPR = rocvalues
    Final_ROC = fig.plot(roc[0], roc[1], 'b')
    Rinal_PR = fig.plot(pr[0], pr[1], 'g')
    halfline = fig.plot([0,1], [0,1], 'k--')
    fig.axis([0, 1, 0, 1])
    fig.set_xlabel('1-specificity (recall)', fontsize = 10)
    fig.set_ylabel('sensivity (precision)', fontsize = 10)
    fig.set_title('ROC curve (Precison/Recall)', fontsize = 12)
    fig.text(0.6, 0.30, 'AUC Roc:%.3f' % aucROC, fontsize = 10)
    fig.text(0.6, 0.20, 'AUC PR:%.3f' % aucPR, fontsize = 10)

def corrplot(fig, corrvalues):
    real, predicted = corrvalues
    corr_plot = fig.plot(real, predicted, 'o', ms=2.5)
    halfline = fig.plot([0,1], [0,1], 'k--')
    fig.axis([0, 1, 0, 1])
    fig.set_xlabel('Real values', fontsize = 10)
    fig.set_ylabel('Predicted values)', fontsize = 10)
    fig.set_title('Real vs. Predicted', fontsize = 12)
    #fig.text(0.6, 0.30, 'Correlation:%.3f' % corr, fontsize = 10)


