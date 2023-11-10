#!/bin/python3
from os import listdir
import glob
import tqdm
import hashlib
import argparse

def main():
    parser = argparse.ArgumentParser(description='Example script with command-line option')
    parser.add_argument('input_string', type=str, help='Input string from command line')
    args = parser.parse_args()

    input_string = args.input_string
    print("Input file:", input_string)
    if input_string.split(".")[1] != "xyzv":
        print("Input file does not contain velocity information")
        exit()

    print("Couting worms...")
    A_count = B_count = S_count = I_count = E_count = frame_count = line_count = 0
    with open(input_string) as infile:
        for line in tqdm.tqdm(infile):
            line_count += 1
            if "A" in line:
                A_count += 1
            if "B" in line:
                B_count += 1
            if "S" in line:
                S_count += 1
            if "I" in line:
                I_count += 1
            if "E" in line:
                E_count += 1
            if "#" in line:
                frame_count += 1
    worm_count = A_count + B_count
    print("Lines: ",line_count,line_count/frame_count)
    print("Worms: ",worm_count, worm_count/frame_count)
    print("Fluid: ",S_count, S_count/frame_count)
    print("Boundary: ",I_count, I_count/frame_count)
    print("Marker: ",E_count, E_count/frame_count)
    print("Frames: ",frame_count)
    
    # particles per frame = line_count/frame_count - 2 (2 for both header lines)

    output_string = input_string.split(".")[0] + "_worms_only.xyzv"
    print("Output file: ",output_string)
    # file1 = open(input_string, 'r')
    # file2 = open(output_string,'a')
    bad_particles = ["S", "I", "E"]
    with open(input_string) as infile, open(output_string, 'w') as outfile:
        for line in tqdm.tqdm(infile,total=line_count):
            if not any(bad_ptc in line for bad_ptc in bad_particles):
                if line.strip() == str(int(line_count/frame_count - 2)):
                    line = str(int(worm_count/frame_count))+"\n"
                outfile.write(line)
    # file1.close()
    # file2.close()
    print(output_string+" written")

if __name__ == '__main__':
    main()