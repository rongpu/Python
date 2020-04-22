# Also downsize image

import os, glob
import numpy as np
import cv2

input_dir = '/Users/rongpu/Downloads/dr9e/coadd/'
output_dir = '/Users/rongpu/Downloads/dr9e/coadd_brighter/'

image_list = glob.glob(os.path.join(input_dir, '*.jpg'))
image_list = sorted(image_list)

for fn in image_list:
    
    output_fn = os.path.join(output_dir, os.path.basename(fn))

    if os.path.isfile(output_fn):
        continue

    img = cv2.imread(fn).astype(np.int32)

    img_small = []
    binsize = 3
    for index in range(3):
        img1 = img[:, :, index].copy()
        # trim edges to enable downsizing
        # trimmed image size need to be multiples of binsize
        trim_size_x = img1.shape[1] % binsize
        trim_size_y = img1.shape[0] % binsize
        img1 = img1[trim_size_y:(img1.shape[0]-trim_size_y), trim_size_x:(img1.shape[1]-trim_size_x)]
        img1 = np.mean(np.mean(img1.reshape((img1.shape[0]//binsize, binsize, img1.shape[1]//binsize,-1)), axis=3), axis=1)
        img_small.append(img1)
    img_small = np.array(img_small)
    img_small = np.moveaxis(img_small, 0, -1)

    img_small = np.clip(img_small * 3, 0, 255).astype(np.uint8)
    cv2.imwrite(output_fn, img_small)
    
