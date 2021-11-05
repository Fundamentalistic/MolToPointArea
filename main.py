import sys
import os
import torch
from PDB import PDB

path = "D:\\Water\\result3.pdb"
pdb = PDB(path)

water = pdb.getWater()
alphaTrace = pdb.getAlphaTrace()

print(len(pdb.getWater()))
print(len(pdb.getAlphaTrace()))
pdb.generateWaterProbabilityNetwork()
print(pdb.getWaterProbabilityNetwork())
print(len(pdb.getWaterProbabilityNetwork()))
pdb.generatePointsTCL()