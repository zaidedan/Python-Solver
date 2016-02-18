"""
This will serve as a place to run unit tests, 
as this is becoming more necessary,

Basic components of the code can be tested here
The idea is that this should be very quick to run, 
and can be included the main run files, and run 
prior to every solve
"""

import poisson2D
import diffusion2D
import time as clocktime

if __name__ == '__main__':
  start = clocktime.time()
  print "running tests ...",
  rate = poisson2D.test()
  assert (rate[0] > 2.3 and rate[1] > 2.3), \
    "testing failed poisson2D: solver error"
  rate = diffusion2D.test()
  assert (rate[0] > 2.3 and rate[1] > 2.3), \
    "testing failed diffusion2D: solver error"
  print "all passed in", '%.2f' % (clocktime.time()-start),"seconds"