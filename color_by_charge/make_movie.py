import cv2
import os

######----Variables to Change----#########
image_folder = '.'
video_name = 'charge_test.mp4'
fps = 6
#####---------------------------##########

# Get sorted list of image filenames
images = [img for img in os.listdir(image_folder) if img.endswith(".bmp")]
images = sorted(images, key=lambda img: int(img.split('.')[0]))
print(images)

# Read the first image to get video dimensions
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

# Use 'mp4v' codec for saving as .mp4
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

# Write images to the video
for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
