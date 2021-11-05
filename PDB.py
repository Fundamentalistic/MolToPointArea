import os

HYDROPHOBIC = 0
CHARGED = 1
POLAR = 2
SPECIAL = 3
PROLINE = 4

def roundBy(val, base):
    return (val // base)*base

class Resolve_Dicts:
    types = {
        'R': CHARGED,
        'H': CHARGED,
        'K': CHARGED,
        'D': CHARGED,
        'E': CHARGED,
        'S': POLAR,
        'T': POLAR,
        'N': POLAR,
        'Q': POLAR,
        'C': SPECIAL,
        'U': SPECIAL,
        'G': SPECIAL,
        'A': HYDROPHOBIC,
        'I': HYDROPHOBIC,
        'L': HYDROPHOBIC,
        'M': HYDROPHOBIC,
        'F': HYDROPHOBIC,
        'W': HYDROPHOBIC,
        'Y': HYDROPHOBIC,
        'V': HYDROPHOBIC,
        'P': PROLINE
    }

    simple_numbers = {
        'A': 3,
        'R': 5,
        'N': 7,
        'D': 11,
        'B': 13,
        'C': 17,
        'E': 23,
        'Q': 29,
        'Z': 31,
        'G': 37,
        'H': 41,
        'I': 43,
        'L': 47,
        'K': 53,
        'M': 59,
        'F': 61,
        'P': 67,
        'S': 71,
        'T': 73,
        'W': 79,
        'Y': 83,
        'V': 89,
        'U': 97,
        'O': 101
    }
    backward_simple_numbers = {
        3: 'A',
        5: 'R',
        7: 'N',
        11: 'D',
        13: 'B',
        17: 'C',
        23: 'E',
        29: 'Q',
        31: 'Z',
        37: 'G',
        41: 'H',
        43: 'I',
        47: 'L',
        53: 'K',
        59: 'M',
        61: 'F',
        67: 'P',
        71: 'S',
        73: 'T',
        79: 'W',
        83: 'Y',
        89: 'V',
        97: 'U',
        101: 'O'
    }
    three_to_simple = {
        'GLY': 'G',
        'ALA': 'A',
        'VAL': 'V',
        'ILE': 'I',
        'LEU': 'L',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'CYS': 'C',
        'MET': 'M',
        'ASP': 'D',
        'ASN': 'N',
        'GLU': 'E',
        'GLN': 'Q',
        'LYS': 'K',
        'ARG': 'R',
        'HIS': 'H',
        'PHE': 'F',
        'TYR': 'Y',
        'TRP': 'W',
        'SEC': 'U',
        'PYR': 'O'
    }
    main_dict = {
        'A': 3,
        'R': 5,
        'N': 7,
        'D': 11,
        'B': 13,
        'C': 17,
        'E': 23,
        'Q': 29,
        'Z': 31,
        'G': 37,
        'H': 41,
        'I': 43,
        'L': 47,
        'K': 53,
        'M': 59,
        'F': 61,
        'P': 67,
        'S': 71,
        'T': 73,
        'W': 79,
        'Y': 83,
        'V': 89,
        'U': 97,
        'O': 101
    }
    coordinates = None
    atype = None
    numerical = None


class PDB:
    content = None
    alphaTrace = []
    protein = []
    water = []
    waterPoints = []
    alphaTraceLen = 0
    waterProbabilityNetwork = []
    test = {}

    minX = 0
    maxX = 0
    minY = 0
    maxY = 0
    minZ = 0
    maxZ = 0

    def __init__(self, path):
        file = open(path, 'r')
        self.content = file.readlines()
        file.close()
        self.__parse_content()

    def __parse_content(self):
        self.AASequence = []
        for line in self.content:
            record = line[0:4]
            if record == "ATOM":
                res_name = line[17:21].replace(" ", "")
                name = line[12:16]
                name = name.replace(' ', '')
                serial = int(line[6:11], 16)
                seq_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if x < self.minX:
                    self.minX = x
                if y < self.minY:
                    self.minY = y
                if z < self.minZ:
                    self.minZ = z

                if x > self.maxX:
                    self.maxX = x
                if y > self.maxY:
                    self.maxY = y
                if z > self.maxZ:
                    self.maxZ = z

                atom = {
                    'res_name': res_name,
                    'serial': serial,
                    'seq_num': seq_num,
                    'x': x,
                    'y': y,
                    'z': z
                }

                if res_name != "TIP3" and name == "CA":
                    self.AASequence.append(res_name)
                    self.alphaTrace.append(atom)
                    self.protein.append(atom)
                elif res_name != "TIP3":
                    self.protein.append(atom)
                else:
                    self.water.append(atom)
                    self.waterPoints.append([x, y, z])

        self.alphaTraceLen = len(self.alphaTrace)

    def getAlphaTrace(self):
        return self.alphaTrace

    def getWater(self):
        return self.water

    def generateWaterProbabilityNetwork(self):
        for index, point in enumerate(self.waterPoints):

            if index % 1000 == 0:
                print(f"Index: {index} Len: {len(self.waterPoints)}")

            x = roundBy(point[0] * 10, 1) / 10
            y = roundBy(point[1] * 10, 1) / 10
            z = roundBy(point[2] * 10, 1) / 10
            probePoint = [x, y, z]

            sx = str(x)
            sy = str(y)
            sz = str(z)

            try:
                self.test[sx][sy][sz] += 1
            except:
                self.test[sx] = {sy: {sz: 1}}
                self.waterProbabilityNetwork.append(probePoint)



    def getWaterProbabilityNetwork(self):
        return self.waterProbabilityNetwork

    def generatePointsTCL(self):

        tclScript = open('waterPoints.tcl', 'w+')
        for point in self.waterProbabilityNetwork:
            tclScript.write("graphics top point  {" + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) +"}\n")
        tclScript.close()
