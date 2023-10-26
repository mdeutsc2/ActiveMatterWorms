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
    #print("Input string:", input_string)
    xyzfilelist = glob.glob('*.xyz')
    xyzfilelist.sort()
    #print(xyzfilelist)

    print("reading data...")
    merged_file = open(input_string,"a")
    for filename in tqdm.tqdm(xyzfilelist):
        #for line in file.readline()
        file = open(filename,"r")
        data = file.read()
        file.close()
        merged_file.write(data)

    merged_file.close()
    print("checking data...")
    # Python program to find MD5 hash value of a file
    md5_hash = hashlib.md5()
    with open(input_string,"rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in tqdm.tqdm(iter(lambda: f.read(4096),b"")):
            md5_hash.update(byte_block)
            # print(md5_hash.hexdigest())

    md5_file = input_string.split(".")[0] + ".md5"
    with open(md5_file,"w+") as f:
        f.write(md5_hash.hexdigest())

    print(input_string)
    print(md5_hash.hexdigest()," > ",md5_file)

if __name__ == '__main__':
    main()
