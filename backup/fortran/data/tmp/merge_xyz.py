#!/bin/python3
from os import listdir
import glob


xyzfilelist = glob.glob('*.xyz')
xyzfilelist.sort()
#print(xyzfilelist)


merged_file = open("amatter.xyz","a")
for filename in xyzfilelist:
    #for line in file.readline()
    file = open(filename,"r")
    print(filename)
    data = file.read()
    file.close()
    merged_file.write(data)


merged_file.close()
