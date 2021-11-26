import numpy
import os
from PIL import Image

minPoint = [-211.50100708007813, -224.50900268554688, 11.904999732971191]
maxPoint = [31.152999877929688, 44.63100051879883, 110.29499816894531]
stepOfSpace = 0.1
pdbFilePath = "D:\\3rjp\\frames\\frames.pdb"
NEAR_CONSTANT = 500
PROTEIN_ATOMS_LEN = 10000
MAX_ATOM_TO_PLANE = 50

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
print("Initiate protein sequence array")
protein = numpy.zeros((PROTEIN_ATOMS_LEN, MAX_ATOM_TO_PLANE, 2), dtype=numpy.float)
protein_cursors = numpy.zeros(PROTEIN_ATOMS_LEN, dtype=numpy.uint)

# Преобразование координат из обычного пространства в дискретное
x_to_descrete_x = lambda x: int((x - x_shift) // stepOfSpace)
y_to_descrete_y = lambda y: int((y - y_shift) // stepOfSpace)
z_to_descrete_z = lambda z: int((z - z_shift) // stepOfSpace)

# Обратное преобразование с потерей
x_discrete_to_x = lambda x: float((x * stepOfSpace) + x_shift)
y_discrete_to_y = lambda y: float((y * stepOfSpace) + y_shift)
z_discrete_to_z = lambda z: float((z * stepOfSpace) + z_shift)


print(f"PDB file {pdbFilePath} reading...")
protein_is_read = False
pdbFile = open(pdbFilePath, 'r')
for ind, line in enumerate(pdbFile):
    record = line[0:4].replace(" ", "")
    if ind % 1000000 == 0:
        print(f"{ind} lines is read", end='\r')
    if record == "ATOM":

        res_name = line[17:21].replace(" ", "")

        if res_name != "TIP3" and not protein_is_read:

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            x_d = x_to_descrete_x(x)
            y_d = y_to_descrete_y(y)
            z_d = z_to_descrete_z(z)

            protein[x_d][protein_cursors[x_d]][0] = y_d
            protein[x_d][protein_cursors[x_d]][1] = z_d
            protein_cursors[x_d] += 1
            if protein_cursors[x_d] >= MAX_ATOM_TO_PLANE:
                print("COORDINATE VOLUME MORE THEN MAX_ATOM_TO_PLANE")
                quit()

    elif record.find("END") != -1:
        protein_is_read = True
        break

print("PDB file reading is complete. Images writing...")

for i in range(x_size):
    img = None
    try:
        img = Image.open(f"img/{i}.png")
    except:
        continue
    img1 = img.convert('RGB')
    rgbi = Image.new('RGB', img1.size)
    rgbi.paste(img1)
    for iy in range(protein_cursors[i]):
        x = int(protein[i][iy][0])
        y = int(protein[i][iy][1])
        rgbi.putpixel((x, y), (0, 255, 0))
        print(x, y)
    print(f"Saving image: {i}.png.", end=" ")
    rgbi.save(f"img/{i}.png")
    print("Complete")

print("Program complete!")
