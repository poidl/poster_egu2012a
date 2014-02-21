With the files in this directory, you can reproduce the figures of the following poster:

Riha, S., Peliz, A., 2012: An Idealized Isopycnal 2-Layer Model of the Alboran Sea.
Presented at EGU General Assembly 2012

Reference of the abstract:
Geophysical Research Abstracts Vol. 14, EGU2012-6292-1, 2012
EGU General Assembly 2012 © Author(s) 2012


The run** directories contain configuration files for Robert Hallberg's
Hallberg Isopycnal Model (HIM). This model is the most important tool for this study, 
and you will have to install it to reproduce the figures.
The source code of HIM (GPL License) is available from GFDL:
http://data1.gfdl.noaa.gov/nomads/forms/him_beta_nomads.html

Runs 53-66 take about an hour (total) to complete on a commodity laptop (Intel i5 CPU M430  
@ 2.27GHz). I found that running the model with OPENMPI on a multi-core single-cpu system 
(i.e. running code designed for distributed-memory systems on a shared-memory system) works well. The other runs (for Fig. 8) will take on the order of a day to complete with the processor mentioned above. The code for Fig. 9) is not included here.

The directory figures_python contains python scripts to generate the figures. The code works
with matplotlib 1.1.0, scipy 0.10.0b2, netCDF4 0.9.8. Other versions may work. 

http://matplotlib.sourceforge.net/users/installing.html
http://www.scipy.org/Installing_SciPy/Linux
http://code.google.com/p/netcdf4-python/wiki/UbuntuInstall


Generate the input files for the model (run**/data/*):
	python prep_run**_**.py


In the run**/Makefiles, change the line

	SRCDIR = /home/stefan/arbeit/him/source/

to point to the HIM source


in run**/init.h change the line

	#define INPUTDIR "/home/stefan/arbeit/him/......"
	
to point to the run directory. Use the script python_utils/search_replace.py if you find it
useful.


If you don't want to use parallelization, # undef PARALLEL_X in run**/init.h.
In this case, change the lines 

	mpirun -n 2 ./run53/HIM < ./run53/temp_in
	
to

	./run53/HIM < ./run53/temp_in
	
in run**_**.sh 


Run runs 53 to 56 for Fig. 3), 4) and 5)

	./run53_56.sh 


Generate the figures with the scripts in figures_python.


To read the content of the netcdf files you can use 'slidenc' (same repository), but note that
the program is work in progress. 
In the current version (on May 5th, 2012) of slidenc, you have to remove the lines

	import pyomReader
	import romsReader
	import himReader
	
in the file '__init__.py' before using it.

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/23dfc47a09f888141e3ac3753bd99439 "githalytics.com")](http://githalytics.com/poidl/poster_egu2012a)
