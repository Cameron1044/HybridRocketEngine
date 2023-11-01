import cv2
import numpy as np
import math

def normalize(value, diameter):
    """Transforms real unit quantities into self.mapX, self.mapY coordinates. For use in indexing into the
    coremap."""
    return value / (0.5 * diameter)

def unNormalize(value, diameter):
        """Transforms self.mapX, self.mapY coordinates to real unit quantities. Used to determine real lengths in
        coremap."""
        return (value / 2) * diameter

def mapToLength(value, diameter, mapDim):
    """Converts pixels to meters. Used to extract real distances from pixel distances such as contour lengths"""
    return diameter * (value / mapDim)

def mapToArea(value, diameter, mapDim):
    """Used to convert sq pixels to sqm. For extracting real areas from the regression map."""
    return (diameter ** 2) * (value / (mapDim ** 2))

def generate_grain_geometry(outer_diameter, num_points, point_length, point_base_width, mapDim):
    point_length = normalize(point_length, outer_diameter)
    point_base_width = normalize(point_base_width, outer_diameter)
    # Creating a white image 1000 pixels long
    img = np.ones((mapDim, mapDim), dtype=np.uint8) * 0
    mapX, mapY = np.meshgrid(np.linspace(-1, 1, mapDim), np.linspace(-1, 1, mapDim))

    # Calculate the center of the image
    center = (img.shape[0] // 2, img.shape[1] // 2)

    # Draw the black circle on the white image
    cv2.circle(img, center, outer_diameter // 2, (0), -1)  # -1 means the circle will be filled
    # cv2.circle(img, center, mapDim // 2, (150, 150, 150), 1)
    
    y, x = np.ogrid[-center[0]:img.shape[0]-center[0], -center[1]:img.shape[1]-center[1]]

    for i in range(0, num_points):
        theta = 2 * np.pi / num_points * i
        comp0 = np.cos(theta)
        comp1 = np.sin(theta)

        rect = abs(comp0 * mapX + comp1 * mapY)
        width = point_base_width / 2 * (1 - (((mapX ** 2 + mapY ** 2) ** 0.5) / point_length))
        vect = rect < width
        near = (comp1 * mapX) - (comp0 * mapY) > -0.025
        img[np.logical_and(vect, near)] = 255  # Set to white
    
    return img

def getCorePerimeter(img, mapDim):
    """Used to calculate the Initial Core Perimeter in pixels """
    # Calculates a countour of the given img
    contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    # Calculates the initial perimeter of the countour
    firstPass = cv2.arcLength(contours[0], True)
    # Douglas-Peucker Algorithm for Shape Approximation to 0.1 % error
    epsilon = 0.001*firstPass
    DPfit = cv2.approxPolyDP(contours[0],epsilon,True)

    # Recalculates the perimeter of the countour
    perimeter = cv2.arcLength(DPfit, True)
    return perimeter

def getCoreArea(img, mapDim):
    """Used to calculate the Initial core area, A_port, in pixels"""
    # Calculates a countour of the given img
    contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    # Calculates the area of the countour
    area = cv2.contourArea(contours[0])
    
    return area

mapDim = 10000
diameter = 6
pointLength = 1
pointWidth = 0.5

# img = generate_grain_geometry(diameter, 5, pointLength, pointWidth, mapDim)

### CIRCLE TEST ###
img = np.ones((mapDim, mapDim), dtype=np.uint8) * 0
center = (img.shape[0] // 2, img.shape[1] // 2)
cv2.circle(img, center, mapDim // 2, (255), -1)
contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
cv2.drawContours(img, contours, -1, (0, 0, 255), 1)  # Red color contour

Perimeter = getCorePerimeter(img, mapDim)
Area = getCoreArea(img, mapDim)
print("\n Perimeter: ", mapToLength(Perimeter, diameter, mapDim))
print("\n Area: ", mapToArea(Area, diameter, mapDim), "\n")

cv2.imshow('Map', img)
cv2.waitKey(0)
cv2.destroyAllWindows()