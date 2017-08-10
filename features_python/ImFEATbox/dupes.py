#!/usr/bin/env python

# Syntax: duplicates.py DIRECTORY

import os, sys

top = sys.argv[1]
d = {}

for root, dirs, files in os.walk(top, topdown=False):
    for name in files:
        fn = os.path.join(root, name)
        basename, extension = os.path.splitext(name)

        basename = basename.lower() # ignore case

        if basename in d:
            print(d[basename])
            print(fn)
        else:
            d[basename] = fn
