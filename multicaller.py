#!/usr/bin/env python3

import subprocess
import os
import sys


if os.path.isfile('logger.log'):
    os.remove('logger.log')

logger = open('logger.log', 'a')
dnull = open(os.devnull, 'w')

pid = str(os.getpid())

logger.write('My pid is: '+pid+'\n')

interval = list(range(-80,81))

for k in interval:
    line = './helios.x -f %s -p field:along-x=%4.3f' % (sys.argv[1], k*0.05)
    logger.write(line+'\n')
    logger.flush()
    process = subprocess.Popen(line, stdout=logger, stderr=logger, shell=True)
    process.wait()
