import numpy
import json
import sys

minPoint = [-211.50100708007813, -224.50900268554688, 11.904999732971191]
maxPoint = [31.152999877929688, 44.63100051879883, 110.29499816894531]
stepOfSpace = 0.1
pdbFilePath = "/home/main/3rjp/frames/frames.pdb"
configurationFilePath = "/home/main/Water/hyperspaces.json"
regime = sys.argv[1]
resultTclFilePath = sys.argv[2]
NEAR_CONSTANT = 500

max_x = maxPoint[0]
max_y = maxPoint[1]
max_z = maxPoint[2]

min_x = minPoint[0]
min_y = minPoint[1]
min_z = minPoint[2]

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

# Создаем дискретное пространство
space = None
if regime == 'frame_read':
    space = numpy.full((x_size + NEAR_CONSTANT, y_size + NEAR_CONSTANT, z_size + NEAR_CONSTANT), 1, dtype=numpy.byte)
    print(x_size, y_size, z_size)
if regime == 'load_from_disk':
    space = numpy.load(f'{resultTclFilePath}.space.np')

# Преобразование координат из обычного пространства в дискретное
x_to_descrete_x = lambda x: int((x - x_shift) // stepOfSpace)
y_to_descrete_y = lambda y: int((y - y_shift) // stepOfSpace)
z_to_descrete_z = lambda z: int((z - z_shift) // stepOfSpace)

# Обратное преобразование с потерей
x_discrete_to_x = lambda x: float((x * stepOfSpace) + x_shift)
y_discrete_to_y = lambda y: float((y * stepOfSpace) + y_shift)
z_discrete_to_z = lambda z: float((z * stepOfSpace) + z_shift)

if regime == 'frame_read':
    print(f"PDB file {pdbFilePath} reading...")
    pdbFile = open(pdbFilePath, 'r')
    for ind, line in enumerate(pdbFile):
        record = line[0:4]
        if ind % 1000000 == 0:
            print(f"{ind} lines is read", end='\r')
        if record == "ATOM":
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            x_d = x_to_descrete_x(x)
            y_d = y_to_descrete_y(y)
            z_d = z_to_descrete_z(z)

            space[x_d, y_d, z_d] = 0

    numpy.save(f'{resultTclFilePath}.space.np', space)

    print("PDB file reading is complete. Hyperspace calculation...")

tclScript = open(resultTclFilePath, 'w+')
hyperparameters = json.loads(open(configurationFilePath, 'r').read())
hyperspaces = hyperparameters['data']
hlen = len(hyperspaces)
for index, hyperspace in enumerate(hyperspaces):
    h_start_x = x_to_descrete_x(hyperspace['min_x'])
    h_end_x = x_to_descrete_x(hyperspace['max_x'])
    h_start_y = y_to_descrete_y(hyperspace['min_y'])
    h_end_y = y_to_descrete_y(hyperspace['max_y'])
    h_start_z = z_to_descrete_z(hyperspace['min_z'])
    h_end_z = z_to_descrete_z(hyperspace['max_z'])
    for x_d in range(h_start_x, h_end_x):
        print(f'x_d: {x_d} from {h_end_x}. Hyperspace index: {index} of {hlen}', end='\r')
        for y_d in range(h_start_y, h_end_y):
            for z_d in range(h_start_z, h_end_z):
                if space[x_d, y_d, z_d] == 1:
                    tclScript.write(
                        "graphics top point  {" + str(x_discrete_to_x(x_d)) + " " + str(
                            y_discrete_to_y(y_d)) + " " + str(z_discrete_to_z(z_d)) + "}\n")

tclScript.close()
"""
tclScript = open('waterPoints.tcl', 'w+')
for x_d in range(x_size + NEAR_CONSTANT):
    for y_d in range(y_size + NEAR_CONSTANT):
        print(f"{x_d} x is read of {x_size} and {y_d} read of {y_size}", end='\r')
        for z_d in range(z_size + NEAR_CONSTANT):
            if space[x_d, y_d, z_d] == 0:
                tclScript.write(
                    "graphics top point  {" + str(x_discrete_to_x(x_d)) + " " + str(y_discrete_to_y(y_d)) + " " + str(z_discrete_to_z(z_d)) + "}\n")
tclScript.close()
"""
print("Complete!")
