import cv2

# Crop size (in downsized pixels)
# Trim 6000x4000 to 6000x3375 (16:9)
crop_left, crop_right, crop_top, crop_bottom = 0, 0, 312, 313

image_fn = "DSC06008.JPG"

img = cv2.imread(image_fn)

# crop
img_resized = img[crop_top:img.shape[0]-crop_bottom, crop_left:img.shape[1]-crop_right]

# resize image to 4k (3840×2160)
dim = (3840, 2160)
img_resized = cv2.resize(img, dim, interpolation=cv2.INTER_LANCZOS4)
# INTER_LANCZOS4 is better than INTER_AREA or the default option, unless downsampling by an integer factor (e.g. 6000x4000 to 3000x2000)

cv2.imwrite('DSC06008_resized.JPG', img_resized, [cv2.IMWRITE_JPEG_QUALITY, 97])
