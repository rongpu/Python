# Resize 6000x4000 to 4K

import cv2
import os
import glob

# video_name = 'video.avi'
video_name = 'video_4k.mp4'

image_fns = sorted(glob.glob('*.JPG'))
image_fns.sort(key=os.path.getmtime)  # Sort by date modified

frame = cv2.imread(image_fns[0])
height, width, layers = frame.shape

width, height = 3840, 2160  # 4K

# video = cv2.VideoWriter(video_name, 0, 1, (width,height))

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'avc1')  # Be sure to use lower case
# avc1 is slightly better than mp4v
video = cv2.VideoWriter(video_name, fourcc, 60.0, (width, height))

for index, image_fn in enumerate(image_fns):

    if index%100==0:
        print('{} / {}'.format(index, len(image_fns)))

    img = cv2.imread(image_fn)

    # Crop size (in downsized pixels)
    # Trim 6000x4000 to 6000x3375 (16:9)
    crop_left, crop_right, crop_top, crop_bottom = 0, 0, 312, 313
    img_resized = img[crop_top:img.shape[0]-crop_bottom, crop_left:img.shape[1]-crop_right]

    # resize image to 4K
    dim = (width, height)
    img_resized = cv2.resize(img_resized, dim, interpolation=cv2.INTER_LANCZOS4)
    # INTER_LANCZOS4 is better than INTER_AREA or the default option, unless downsampling by an integer factor (e.g. 6000x4000 to 3000x2000)

    video.write(img_resized)

cv2.destroyAllWindows()
video.release()
