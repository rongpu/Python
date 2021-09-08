import cv2
import os
import glob

video_name = 'video.mp4'
downsize = True

image_fns = sorted(glob.glob('*.JPG'))
image_fns.sort(key=os.path.getmtime)  # Sort by date modified

frame = cv2.imread(image_fns[0])
height, width, layers = frame.shape

if downsize:
    width, height = width//2, height//2

# video = cv2.VideoWriter(video_name, 0, 1, (width,height))

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'avc1')  # Be sure to use lower case
video = cv2.VideoWriter(video_name, fourcc, 30.0, (width, height))

for index, image_fn in enumerate(image_fns):

    if index%100==0:
        print('{} / {}'.format(index, len(image_fns)))

    img = cv2.imread(image_fn)

    if downsize:
        # downsize image
        dim = (width, height)
        img = cv2.resize(img, dim, interpolation=cv2.INTER_AREA)
        # INTER_AREA is better than INTER_LANCZOS4 when downsampling by an integer factor (e.g. 6000x4000 to 3000x2000)

    video.write(img)

cv2.destroyAllWindows()
video.release()
