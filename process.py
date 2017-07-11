#!/usr/bin/env python

################################################################################
# process.py
#
# by John S. McCaskill and Norman H. Packard
#
################################################################################
# execute with "./activity.py ./genelifeAct.py"  


import sys
import os.path

usage = "usage:  launch.py designfile responsefile"
if __name__ == '__main__':
  if len(sys.argv) != 3:
    print usage
    sys.exit()
  design = sys.argv[1]
  response = sys.argv[2]
  if (not os.path.isfile(design)) or (not os.path.isfile(response)):
    print usage
    sys.exit()

  # get responses into dictionary:
  resp = {}
  with open(response, 'r') as f:
    for x in f:
      xx = x.split()
      res = xx[-1]              # grab response
      xx = xx[:-1]              # clip response
      xxx = ' '.join([str(float(s)) for s in xx])        # make hashable
      if xxx in resp:
        resp[xxx].append(res)
      else:
        resp[xxx] = [res]
  mx = 0
  for x in resp:
    if len(resp[x])>mx:
      mx = len(resp[x])
  for x in resp:
    if len(resp[x])!= mx:
      print "missing responses!"
      sys.exit()

  idx = {}
  for x in resp:
    idx[x] = 0;
  des = []
  with open(design, 'r') as f:
    for x in f:
      x = x.rstrip()
      xx = x.split()
      xxx = ' '.join([str(float(s)) for s in xx]) # make hashable, with canonical spaces
      if xxx in resp:
        y = x + '\t' + resp[xxx][idx[xxx]] # tab for compatibility with HOT
        print y
        idx[xxx] = idx[xxx] + 1
      else:
        print 'warning: ',xxx,'not in response dictionary.'
