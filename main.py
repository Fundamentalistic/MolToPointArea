import sys
import os
import torch
from PDB import PDB

path = "/home/main/Water/result4.pdb"
hyperspacesPath = "/home/main/Water/hyperspaces.json"
pdb = PDB(path)
pdb.read_hyperspaces_configuration_file(hyperspacesPath)
pdb.produce_empty_sets_from_hyperspaces()

water = pdb.get_water()
alphaTrace = pdb.get_alpha_trace()

print(len(pdb.get_water()))
print(len(pdb.get_alpha_trace()))
pdb.generate_water_probability_network()
# print(pdb.get_water_probability_network())
print(len(pdb.get_water_probability_network()))
pdb.generate_points_tcl()
