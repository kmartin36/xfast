#!/usr/bin/python3

from random import random
import sys
from math import exp, log
from subprocess import call

print('N std::map::emplace kfast::insert std::map:lower_bound kfast::lower_bound log(maxrss)/10')

#for i in range(int(sys.argv[1])):
while True:
  count = int(exp(random()*log(80000000/1000) + log(1000)))
  call(['./xfast', str(count)])
