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

    img = cv2.imread(fn).astype(np.uint16)
    img1 = np.clip(img * 2, None, 255).astype(np.uint8)
    cv2.imwrite(output_fn, img1)
    