#!/bin/python3
from os import listdir, remove
from os.path import isfile, isdir, join, basename, getsize
import glob
import tqdm
import hashlib
import zipfile
import argparse
import numpy as np

def human_size(nbytes):
    suffixes = ['B','KB','MB','GB','TB','PB']
    i=0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def get_md5(input_string,md5_file_string):
    print("checking data...")
    # Python program to find MD5 hash value of a file
    md5_hash = hashlib.md5()
    with open(input_string,"rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in tqdm.tqdm(iter(lambda: f.read(4096),b"")):
            md5_hash.update(byte_block)
            # print(md5_hash.hexdigest())
    
    #md5_file = input_string.split(".")[0] + ".md5"
    md5_file = md5_file_string + ".md5"
    with open(md5_file,"w+") as f:
        f.write(md5_hash.hexdigest())

    print(input_string)
    print(md5_hash.hexdigest()," > ",md5_file)

def main():
   parser = argparse.ArgumentParser(description='Example script with command-line option')
   parser.add_argument('input_string', type=str, help='Input string from command line')
   parser.add_argument("-i","--intermediate",help="generates an intermediate file for visualization",default=False)
   parser.add_argument("-d","--delete",help="delete raw files after compression",default=True)
   parser.add_argument("-c","--compresslevel",help="compression level for zlib",default=9)

   args = parser.parse_args()

   input_string = args.input_string
   #print("Input string:", input_string)
   xyzfilelist = glob.glob('*.xyz')
   xyzfilelist.sort()
   #print(xyzfilelist)

   # getting data about the sizes of the arrays that should be generated

   with open(xyzfilelist[0]) as f:
      num_ptc_total = int(f.readline().strip())
   print("Num total:\t",num_ptc_total)
   ptc_type,_,_,_,_,_,_ = np.loadtxt(xyzfilelist[0],dtype="str",unpack=True,skiprows=2)
   unique,counts = np.unique(ptc_type,return_counts=True)
   print(unique)
   print(counts)
   filament_count = counts[0] + counts[1]
   solvent_count = counts[-1]
   print(filament_count,solvent_count)

   fpos_data = np.zeros((len(xyzfilelist),filament_count,3))
   fvel_data = np.zeros((len(xyzfilelist),filament_count,3))
   spos_data = np.zeros((len(xyzfilelist),solvent_count,3))
   svel_data = np.zeros((len(xyzfilelist),solvent_count,3))

   for i in tqdm.tqdm(np.arange(len(xyzfilelist))):
      print(xyzfilelist[i])
      ptc_type,px,py,pz,vx,vy,vz = np.loadtxt(xyzfilelist[i],dtype="str",unpack=True,delimiter=" ",skiprows=2)
      px = px.astype(np.float32)
      py = py.astype(np.float32)
      pz = pz.astype(np.float32)
      vx = vx.astype(np.float32)
      vy = vy.astype(np.float32)
      vz = vz.astype(np.float32)
      fpos_data[i,:,0] = px[0:filament_count]
      fpos_data[i,:,1] = py[0:filament_count]
      fpos_data[i,:,2] = pz[0:filament_count]
      fvel_data[i,:,0] = vx[0:filament_count]
      fvel_data[i,:,1] = vy[0:filament_count]
      fvel_data[i,:,2] = vz[0:filament_count]

      sol_locs = np.argwhere(ptc_type == 'S')
      spos_data[i,:,0] = px[sol_locs].flatten()
      spos_data[i,:,1] = py[sol_locs].flatten()
      spos_data[i,:,2] = pz[sol_locs].flatten()
      svel_data[i,:,0] = vx[sol_locs].flatten()
      svel_data[i,:,1] = vy[sol_locs].flatten()
      svel_data[i,:,2] = vz[sol_locs].flatten()

   np.savez_compressed("./filament_data",a=fpos_data, b=fvel_data)
   np.savez_compressed("./solvent_data",a=spos_data, b=svel_data)


   

   #  print("reading data...")
   #  merged_file = open(input_string,"a")
   #  xyzfile_basenames = []
   #  total_size = 0
   #  for filename in tqdm.tqdm(xyzfilelist):
   #      #for line in file.readline()
   #      file = open(filename,"r",encoding='utf-8')
   #      data = file.read()
   #      print(data[0])
   #      file.close()
   #      merged_file.write(data)
   #      xyzfile_basenames.append(basename(filename))
   #      total_size += getsize(filename)

   # merged_file.close()
   #  if args.intermediate == False:
   #      get_md5(input_string,input_string.split(".")[0])
    
   #      print('compressing data...')
   #      # zipping up the raw data files
   #      loczip = input_string.split(".")[0]+".zip"
   #      with zipfile.ZipFile(loczip,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=args.compresslevel) as zf:
   #          for i,t in tqdm.tqdm(enumerate(xyzfilelist),total=len(xyzfilelist)):
   #              zf.write(t,arcname=xyzfile_basenames[i])
   #      print(len(xyzfilelist)," files (",human_size(total_size),") >>> ",loczip,"(",human_size(getsize(loczip)),")")


   #      get_md5(loczip,loczip)
   #      print("Removing duplicate files...?",args.delete)
   #      if args.delete:
   #          for f in xyzfilelist:
   #              remove(f)
    
    
    



if __name__ == '__main__':
    main()
