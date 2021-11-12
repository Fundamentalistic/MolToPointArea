import os
import sys
import json

pdbFilePath = sys.argv[1]
configurationFilePath = sys.argv[2]
stepOfSpace = int(sys.argv[3])

hyperparameters = json.loads(open(configurationFilePath, 'r').read())

min_x = hyperparameters['min_x']
min_y = hyperparameters['min_y']
min_z = hyperparameters['min_z']
max_x = hyperparameters['max_x']
max_y = hyperparameters['max_y']
max_z = hyperparameters['max_z']
x_dif = hyperparameters['x_dif']
y_dif = hyperparameters['y_dif']
z_dif = hyperparameters['z_dif'] 

