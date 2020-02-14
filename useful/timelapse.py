import cv2
import os
import glob

image_folder = 'timelapse'
# video_name = 'video.avi'
video_name = 'video.mp4'

images_fn = sorted(glob.glob(os.path.join(image_folder, '*.JPG')))
frame = cv2.imread(images_fn[0])
height, width, layers = frame.shape

width, height = width//2, height//2

# video = cv2.VideoWriter(video_name, 0, 1, (width,height))

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
video = cv2.VideoWriter(video_name, fourcc, 30.0, (width, height))

for image_fn in images_fn:

    img = cv2.imread(image_fn)

    # resize image
    dim = (width, height)
    img_resized = cv2.resize(img, dim, interpolation = cv2.INTER_AREA)

    video.write(img_resized)

cv2.destroyAllWindows()
video.release()

