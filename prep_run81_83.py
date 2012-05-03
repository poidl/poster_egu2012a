import os as os

vec=[81,82,83];
for ii in vec:
    os.system('python ./run'+str(ii)+'/preprocess/triang_ana.py')
    os.system('python ./run'+str(ii)+'/preprocess/initial.py')
    os.system('python ./run'+str(ii)+'/preprocess/sponge.py')
