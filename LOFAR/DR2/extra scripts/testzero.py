#!/usr/bin/python

import os
import glob

files=glob.glob("Cutout_Catalogue-*")

for elem in files:
    size=os.stat(elem).st_size
    print "Size of cat "+elem+" is "+str(size)
    if size is 0:
        print "size is zero!"
