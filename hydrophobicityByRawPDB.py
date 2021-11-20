import numpy
import json
import sys
from PIL import Image

minPoint = [-211.50100708007813, -224.50900268554688, 11.904999732971191]
maxPoint = [31.152999877929688, 44.63100051879883, 110.29499816894531]
stepOfSpace = 0.1
pdbFilePath = "/home/main/3rjp/frames/frames.pdb"
configurationFilePath = "/home/main/Water/hyperspaces.json"
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
space = numpy.full((x_size + NEAR_CONSTANT, y_size + NEAR_CONSTANT, z_size + NEAR_CONSTANT), 1, dtype=numpy.ubyte)
print(x_size, y_size, z_size)

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

numpy.save(f'result.space.np', space)

print("PDB file reading is complete. Hyperspace calculation...")

print("Complete!")
