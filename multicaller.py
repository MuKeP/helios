#!/usr/bin/env python3

import subprocess
import os

logger = open('logger.log', 'a')
dnull = open(os.devnull, 'w')

pid = str(os.getpid())

#logger.write('My pid is: '+pid+'\n')
# os.system('disown '+pid)

interval = list(range(-200,201))

#interval = list(range(-240,-120))+list(range(121,241))
#interval = list(range(226,241))
#interval = list(range(-40,41))

for k in interval:
    line = './helios.x -f btransistor.inp -p field:along-x=%4.3f' % (k*0.05)
    logger.write(line+'\n')
    logger.flush()
    process = subprocess.Popen(line, stdout=logger, stderr=logger, shell=True)
    process.wait()
