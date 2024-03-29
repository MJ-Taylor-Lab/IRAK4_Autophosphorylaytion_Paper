import tifffile, ray, settings, time, os
import cv2 as cv # Install as opencv-python
import numpy as np
from scipy import ndimage


@ray.remote
def remove_frame(image_path, frame_path, tiff_compression_level):
    try:
        # Save name
        path = os.path.dirname(image_path)
        file = os.path.basename(image_path)
        file = os.path.splitext(file)[0] + '_darkframe_removed.tif'
        save_image_path = os.path.join(path, file)
        # Import images
        original_image = tifffile.imread(image_path)
        frame_image = tifffile.imread(frame_path)

        # Subtract images
        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)
        # Subtract
        frame_removed_image = []
        for t in frames:
            frame_removed = original_image[t]
            frame_removed = cv.subtract(frame_removed, frame_image)
            frame_removed_image.append(frame_removed)
        # Save image
        tifffile.imsave(save_image_path, frame_removed_image, bigtiff=True, compress=tiff_compression_level, dtype=frame_removed_image[0].dtype)
        time.sleep(5)
    except:
        print("Cannot complete remove")
    return save_image_path


@ray.remote
def median_blur_remove(img, median_filter_size):
    # Remove blur
    med = ndimage.median_filter(img, size=median_filter_size)
    med_rm = img - np.minimum(med, img)

    return med_rm

def remove_median_blur(image_path, old_term, new_term, median_size, tiff_compression_level):
    try:
        # Save name
        path = os.path.dirname(image_path)
        file = os.path.basename(image_path)
        file = os.path.splitext(file)[0]
        file = file.replace(old_term, new_term) + '.tif'
        save_image_path = os.path.join(path, file)

        # Import image
        original_image = tifffile.imread(image_path)

        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)

        # Run in parallel
        ray.init()
        result_ids = []
        for frame in frames:
            result_id = median_blur_remove.remote(original_image[frame], median_size)
            result_ids.append(result_id)
        med_rm_img = settings.parallel.ids_to_vals(result_ids)

        # Save image
        tifffile.imsave(save_image_path, med_rm_img, bigtiff=True, compress=tiff_compression_level, dtype=med_rm_img[0].dtype)
        time.sleep(10)
        ray.shutdown()
    except:
        print("Cannot complete remove")
    return save_image_path


@ray.remote
def tracking_image(median_image_path, old_term, new_term, median_size, tiff_compression_level):
    try:
        # Save name
        path = os.path.dirname(median_image_path)
        file = os.path.basename(median_image_path)
        file = os.path.splitext(file)[0]
        file = file.replace(old_term, new_term) + '.tif'
        save_image_path = os.path.join(path, file)

        # Import image
        original_image = tifffile.imread(median_image_path)

        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)

        # Blur image
        median_imgs = []
        for frame in frames:
            median_img = ndimage.median_filter(original_image[frame], size=median_size)
            median_imgs.append(median_img)

        # Average frames
        frames = original_image.shape[0]
        frames = range(1, frames-1)
        moving_avg_imgs = []
        for frame in frames:
            moving_avg_img = (median_imgs[frame - 1], median_imgs[frame], median_imgs[frame + 1])
            moving_avg_img = np.array(np.mean(moving_avg_img, axis=0))
            moving_avg_imgs.append(moving_avg_img)

        # Save image
        tifffile.imsave(save_image_path, moving_avg_imgs, bigtiff=True, compress=tiff_compression_level, dtype=original_image[0].dtype)
    except:
        print("Cannot complete remove")
    return save_image_path
