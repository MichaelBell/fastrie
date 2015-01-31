#!/usr/bin/env python

n=0
for i in xrange(63,-1,-1):
  n = n << 32
  n = n + i

print "Correct answer (python): %d" % (n%49979693,)
