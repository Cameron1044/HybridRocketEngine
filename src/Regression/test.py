import cv2
import numpy as np
import math

def generate_grain_geometry(outer_diameter, num_points, point_length, point_base_width):
    # Create a white image of the desired dimensions
    img = np.ones((outer_diameter, outer_diameter), dtype=np.uint8) * 0
    
    # Calculate the center of the image
    center = (img.shape[0] // 2, img.shape[1] // 2)

    # Draw the black circle on the white image
    cv2.circle(img, center, outer_diameter // 2, (0), -1)  # -1 means the circle will be filled
    
    y, x = np.ogrid[-center[0]:img.shape[0]-center[0], -center[1]:img.shape[1]-center[1]]
    mapX = x / outer_diameter
    mapY = y / outer_diameter
    
    for i in range(0, num_points):
        theta = 2 * np.pi / num_points * i
        comp0 = np.cos(theta)
        comp1 = np.sin(theta)

        rect = abs(comp0 * mapX + comp1 * mapY)
        width = point_base_width / 2 * (1 - (((mapX ** 2 + mapY ** 2) ** 0.5) / point_length))
        vect = rect < width
        near = comp1 * mapX - comp0 * mapY > -0.025
        img[np.logical_and(vect, near)] = 255  # Set to white
    
    return img

def apply_gaussian_to_star(img, pixel_radius):
    # Apply Gaussian blur to the extracted star
    blurred_star = cv2.GaussianBlur(img, (2 * pixel_radius + 1, 2 * pixel_radius + 1), 0)

    # Merge the blurred star with the original image to get the combined result
    combined = cv2.bitwise_or(blurred_star, img)

    # # Convert non-black pixels to white
    combined[combined != 0] = 255

    return combined

def find_and_draw_contour(img):
    # Convert the image to grayscale
    gray = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)

    # Find contours
    contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Draw the contour in red
    cv2.drawContours(gray, contours, -1, (0, 0, 255), 1)  # Red color contour

    return gray

# Example usage:
img = generate_grain_geometry(500, 5, 0.4, 0.2)
img = apply_gaussian_to_star(img, 10)
img_with_contour = find_and_draw_contour(img)
cv2.imshow('Contour of Regressed Star', img_with_contour)
cv2.waitKey(0)
cv2.destroyAllWindows()