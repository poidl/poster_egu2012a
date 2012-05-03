import os as os
import fileinput
import sys

original_stdout = sys.stdout
cd=os.getcwd()

output=True

vec=[64,58,66,61,60,62,63];
#vec=[53,54,55,56];
#vec=[81,82,83]

searchstr='#define INPUTDIR'
oldstr='/home/stefan/arbeit/gib/poster2012/'
newstr='/home/stefan/arbeit/gib/trash/'

for ii in vec:
    fstring=cd+'/../run'+str(ii)+'/init.h'

    for line in fileinput.FileInput(fstring,inplace=1):
        if searchstr in line:
            if output:
                print  >> original_stdout,'  OLD: '+line,
            line=line.replace(oldstr,newstr)
            if output:
                print >> original_stdout, '  NEW: '+line,
        print line,
         
    
