#!/bin/python3
from os import listdir
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='Example script with command-line option')
    parser.add_argument('input_string', type=str, help='Input string from command line')
    args = parser.parse_args()

    input_string = args.input_string
    #print("Input string:", input_string)
    xyzfilelist = glob.glob('*.xyz')
    xyzfilelist.sort()
    #print(xyzfilelist)


    merged_file = open(input_string,"a")
    for filename in xyzfilelist:
        #for line in file.readline()
        file = open(filename,"r")
        print(filename)
        data = file.read()
        file.close()
        merged_file.write(data)

    merged_file.close()

if __name__ == '__main__':
    main()
