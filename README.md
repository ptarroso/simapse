# SIMAPSE

Simapse is a software for building niche models using artificial neural networks (ANN) written in Python. It provides a graphical interface with all possible configurations, including data resampling methods, ANN setup and other options. It provides a transparent process as all data generated can be found in the chosen ouptut directory in plain text format. Simapse can deal with continuous data (p.e. abundance), binary presence absence and presence-only data sets. In the latter case it randomly generates pseudo-absences from the non-presence pixels of the area defined by the rasters.

Besides the graphical interface, Simapse also offers a command line option fully configurable. It allows to automate the process of niche modeling which is particularly useful when dealing with multiple species.


## Installation

Simapse is very easy to install: download the code as ip, unpack it and run either the 'simapse.py' (Linux, mac os) or the 'start.bat' file (Windows).

### Dependencies

Simapse has very few dependencies. If you don't have already, you will need [Python](https://www.python.org/) installed in your system. Simapse was originally programmed for Python 2 and the current 'Master' brach here is for this Python version. The branch 'port2python3' adds support with a fully working Python 3 version with all functionality. It was a direct port to Python 3, correcting incompatibilities, and does not add any additional functionality.

Simapse has a two other optional dependencies: [matplotlib](https://matplotlib.org/) for producing plots and [Python Imaging Library - Pillow](https://pypi.org/project/Pillow/) for showing the plots. If you don't have these dependencies installed, a warning is shown but you can still use the software. Outputs will be only data and pots can be produced in any other plotting software.

## Usage

Simapse tries to provide a very simple interface for producing the models, however it does not do exhaustive checks on the input data. The inputs are in plain text and include observation and raster data.

### Observation/presence data

The observation data is a text file with 3 fields separated by a semicolon. It follows the format

| Presence |  X  |  Y  |
|:--------:|-----|-----|
|1         |-8.12|40.24|
|1         |-7.34|38.23|
|...       |...  |...  |

The first column has the observational data and might be:
- Presence-only: a series of "1"s
- Presence/absence: ones and zeros
- Continuous: a number with "." as decimal separator.

### Raster data

The raster data follows the plain text ASCII grid format. It consists of a 6 line header and a matrix-like disposition of pixel values.


ncols       |4
nrows       |5
xllcorner   |-8.0
yllcorner   |34.0
cellsize    |1000
NODATA_value|-9999
1 |2 |3 |4
5 |6 |7 |8
9 |10|11|12
13|14|15|-9999
16|17|18|-9999


## Citation

[Tarroso, P., Carvalho, S. B., & Brito, J. C. (2012). Simapse–simulation maps for ecological niche modelling. *Methods in Ecology and Evolution*, 3(5), 787-791.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00210.x)

## Note

The software is fully functional and was exhaustively tested with many data sets. However the code grew fast and lacks a good structure. A rewriting of the whole program is desirable and might happen in the future.
