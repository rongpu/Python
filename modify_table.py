from __future__ import print_function, division
from astropy.table import Table, Column
import numpy as np

t = Table()
t['a'] = np.arange(10)
mask = np.zeros(10, dtype=bool)
mask[0:3]=True

t[mask]['a']=99
print(t)

t['a'][mask]=99
print(t)
