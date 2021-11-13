import os
import sys
import json
import numpy

pdbFilePath = sys.argv[1]
configurationFilePath = sys.argv[2]
stepOfSpace = float(sys.argv[3])

# Устанавливаем параметры гиперпространств
hyperparameters = json.loads(open(configurationFilePath, 'r').read())

min_x = hyperparameters['min_x']
min_y = hyperparameters['min_y']
min_z = hyperparameters['min_z']
max_x = hyperparameters['max_x']
max_y = hyperparameters['max_y']
max_z = hyperparameters['max_z']
x_dif = max_x - min_x
y_dif = max_y - min_y
z_dif = max_z - min_z

print(f"min_x: {min_x} min_y: {min_y} min_z: {min_y}")
print(f"max_x: {max_x} max_y: {max_y} max_z: {max_z}")
print(f"x_dif: {x_dif} y_dif: {y_dif} z_dif: {z_dif}")

# Устанавливаем параметры перехода из одной системы координат в другую
x_shift = min_x if min_x < 0 else -min_x
y_shift = min_y if min_y < 0 else -min_y
z_shift = min_z if min_z < 0 else -min_z

# Устанавливаем размеры координатных осей дискретного пространства
x_size = int(x_dif // stepOfSpace)
y_size = int(y_dif // stepOfSpace)
z_size = int(z_dif // stepOfSpace)

print(f"x_size: {x_size}, y_size: {y_size}, z_size: {z_size}")

# Создаем дискретное пространство
space = numpy.full((x_size + 100, y_size + 100, z_size + 100), 1, dtype=numpy.byte)

# Преобразование координат из обычного пространства в дискретное
x_to_descrete_x = lambda x: int((x - x_shift) // stepOfSpace)
y_to_descrete_y = lambda y: int((y - y_shift) // stepOfSpace)
z_to_descrete_z = lambda z: int((z - z_shift) // stepOfSpace)

# Обратное преобразование с потерей
x_discrete_to_x = lambda x: float((x * stepOfSpace) + x_shift)
y_discrete_to_y = lambda y: float((y * stepOfSpace) + y_shift)
z_discrete_to_z = lambda z: float((z * stepOfSpace) + z_shift)

print(f"PDB file {pdbFilePath} reading...")
pdbFile = open(pdbFilePath, 'r')
for line in pdbFile:
    record = line[0:4]
    if record == "ATOM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        x_d = x_to_descrete_x(x)
        y_d = y_to_descrete_y(y)
        z_d = z_to_descrete_z(z)

        space[x_d, y_d, z_d] = 0
print("PDB file reading is complete. Hyperspace calculation...")
hyperspaces = hyperparameters['data']
coordinatesArray = []
for hyperspace in hyperspaces:
    h_start_x = x_to_descrete_x(hyperspace['min_x'])
    h_end_x = x_to_descrete_x(hyperspace['max_x'])
    h_start_y = y_to_descrete_y(hyperspace['min_y'])
    h_end_y = y_to_descrete_y(hyperspace['max_y'])
    h_start_z = z_to_descrete_z(hyperspace['min_z'])
    h_end_z = z_to_descrete_z(hyperspace['max_z'])
    for x_d in range(h_start_x, h_end_x):
        for y_d in range(h_start_y, h_end_y):
            for z_d in range(h_start_z, h_end_z):
                if space[x_d, y_d, z_d] == 0:
                    # print("x_d, y_d, z_d")
                    # print(x_d, y_d, z_d)
                    coordinatesArray.append((x_discrete_to_x(x_d), y_discrete_to_y(y_d), z_discrete_to_z(z_d)))

print("TCL generating...")
tclScript = open('waterPoints.tcl', 'w+')
for point in coordinatesArray:
    tclScript.write("graphics top point  {" + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) +"}\n")
tclScript.close()
print("Complete!")
