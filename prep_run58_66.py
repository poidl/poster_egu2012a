import os as os

vec=[64,58,66,61,60,62,63];
for ii in vec:
    os.system('python ./run'+str(ii)+'/preprocess/triang_ana.py')
    os.system('python ./run'+str(ii)+'/preprocess/initial.py')
    os.system('python ./run'+str(ii)+'/preprocess/sponge.py')

