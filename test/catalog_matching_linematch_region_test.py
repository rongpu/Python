# Test region selection in catalog_matching_linematch.py

import numpy as np

x = np.random.rand(30)
x = (x-0.5)*3
print(x)
ra2 = x.copy()
ra2 = ra2%360

ra2max = ra2.max()
ra2min = ra2.min()
print((ra2min-ra2max)%360)

if (ra2min-ra2max)%360<1.:  # If the cat2 region crosses with RA=0
    print('Warning: cat2 crosses with RA=0')
    print(x.min()+360)
    print(x.max())
    ra2min = (((ra2+180)%360.).min()-180)%360.
    ra2max = (((ra2+180)%360.).max()-180+360)%360.
    print(ra2min) # should be equal to x.min()+360
    print(ra2max) # should be equal to x.max()