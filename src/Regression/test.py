from abc import abstractmethod

import cv2
import numpy as np
import math

class Regression():
    def __init__(self, outerDiameter, mapDim=1000):
        self.outerDiameter = outerDiameter
        self.mapDim = mapDim
        self.mapX, self.mapY = np.meshgrid(np.linspace(-1.5, 1.5, mapDim), np.linspace(-1.5, 1.5, mapDim))
        self.r = 1 / 2

    @abstractmethod
    def generate_grain_geometry(self):
        """Generates a grain geometry based on the given parameters"""

    def normalize(self, value):
        """Transforms real unit quantities into self.mapX, self.mapY coordinates. For use in indexing into the
        coremap."""
        return value / (0.5 * self.outerDiameter)

    def unNormalize(self, value):
            """Transforms self.mapX, self.mapY coordinates to real unit quantities. Used to determine real lengths in
            coremap."""
            return (value / 2) * self.outerDiameter

    def mapToLength(self, value):
        """Converts pixels to meters. Used to extract real distances from pixel distances such as contour lengths"""
        return self.outerDiameter * (value / self.mapDim / 1.5)

    def mapToArea(self, value):
        """Used to convert sq pixels to sqm. For extracting real areas from the regression map."""
        return (self.outerDiameter ** 2) * (value / ((self.mapDim/1.5) ** 2))

    def find_contour(self, img):
        # Find contours
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        return contours

    def draw_all_contours(self, base_img, all_contours):
        # Colors for drawing alternating contours
        colors = [(0, 0, 255), (0, 255, 0), (255, 0, 0), (255, 255, 0), (255, 0, 255), (0, 255, 255)]
        
        # Convert the image to a colored one
        colored_img = cv2.cvtColor(base_img, cv2.COLOR_GRAY2BGR)

        # Draw all contours with alternating colors
        for i, contours in enumerate(all_contours):
            cv2.drawContours(colored_img, contours, -1, colors[i % len(colors)], 1)

        return colored_img

    def apply_disk_filter(self, img, radius):
        # Create a disk-shaped kernel
        y, x = np.ogrid[-radius: radius+1, -radius: radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel = np.zeros((2*radius+1, 2*radius+1))
        kernel[mask] = 1
        kernel /= kernel.sum()  # Normalize
        
        # Convolve the image with the kernel
        filtered_img = cv2.filter2D(img, -1, kernel)
        
        return filtered_img

    def apply_threshold(self, prev_img, img, threshold_value):
        if threshold_value < 0 or threshold_value > 1:
            raise ValueError("Threshold value must be between 0 and 1")
        _, thresholded_img = cv2.threshold(img, threshold_value * 255, 255, cv2.THRESH_BINARY)
        # _, thresholded_img0 = cv2.threshold(img, 0 * 255, 255, cv2.THRESH_BINARY)
        # prev_nonzero = np.count_nonzero(prev_img)
        # current_nonzero = np.count_nonzero(thresholded_img)
        # current_nonzero0 = np.count_nonzero(thresholded_img0)
        # area_thresh = current_nonzero - prev_nonzero
        # area_thresh0 = current_nonzero0 - prev_nonzero
        # area_ratio = area_thresh / area_thresh0
        # rdot = self.mapToLength(40*area_ratio)
        # self.r = self.r + rdot
        # print(f"{threshold_value} : {area_thresh} : {area_thresh0} : {area_ratio} : {rdot} : {self.r}")
        return thresholded_img

    # Function to determine if a contour intersects with a circular boundary
    def contour_intersects_boundary(self, contour, center, radius):
        for point in contour:
            distance = np.linalg.norm(point[0] - center)
            if distance > radius:
                return True
        return False

    # Function to determine if a contour intersects with a circular boundary
    def contour_intersects_boundary2(self, contour, center, radius):
        for point in contour:
            distance = np.linalg.norm(point[0] - center)
            if distance < radius:
                return False
        return True

    def getCorePerimeter(self, img):
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

    def getCoreArea(self, img):
        """Used to calculate the Initial core area, A_port, in pixels"""
        # Calculates a countour of the given img
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        # Calculates the area of the countour
        area = cv2.contourArea(contours[0])
        
        return area
    
class StarGeometry(Regression):
    def __init__(self, outer_diameter, num_points, point_length, point_base_width, mapDim=1000):
        super().__init__(outer_diameter, mapDim)
        self.outer_diameter = outer_diameter
        self.num_points = num_points
        self.point_length = point_length
        self.point_base_width = point_base_width

    def generate_grain_geometry(self):
        outer_diameter = self.outer_diameter
        num_points = self.num_points
        point_length = self.point_length
        point_base_width = self.point_base_width

        # Create a white image of the desired dimensions
        buffer_outer_diameter = math.ceil(outer_diameter * 1.2)
        img = np.ones((buffer_outer_diameter, buffer_outer_diameter), dtype=np.uint8) * 0
        
        # Calculate the center of the image
        center = (img.shape[0] // 2, img.shape[1] // 2)
        
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
    
class FinocylGeometry(Regression):
    def __init__(self, outer_diameter, inner_diameter, num_fins, fin_length, fin_width, mapDim=1000):
        super().__init__(outer_diameter, mapDim)
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        self.num_fins = num_fins
        self.fin_length = fin_length
        self.fin_width = fin_width

    def generate_grain_geometry(self):
        innerDiameter = self.normalize(self.inner_diameter)
        finLength = self.normalize(self.fin_length)
        finWidth = self.normalize(self.fin_width)
        numFins = self.num_fins
        mapDim = self.mapDim

        img = np.ones((mapDim, mapDim), dtype=np.uint8) * 0

        # Open up core
        img[self.mapX**2 + self.mapY**2 < (innerDiameter / 2)**2] = 255

        # Add fins
        for i in range(0, numFins):
            theta = 2 * np.pi / numFins * i
            # Initialize a vector pointing along the fin
            vect0 = np.cos(theta)
            vect1 = np.sin(theta)
            # Select all points within half the width of the vector
            vect = abs(vect0*self.mapX + vect1*self.mapY) < finWidth / 2
            # Set up two perpendicular vectors to cap off the ends
            near = (vect1 * self.mapX) - (vect0 * self.mapY) > 0 # Inside of the core
            far = (vect1 * self.mapX) - (vect0 * self.mapY) < finLength # At the casting tube end of the vector
            ends = np.logical_and(far, near)
            # Open up the fin
            img[np.logical_and(vect, ends)] = 255
        return img
    
outer_diameter = 4
finocyl = FinocylGeometry(outer_diameter=outer_diameter, inner_diameter=1.25, num_fins=8, fin_length=1.5, fin_width=0.2)
base_img = finocyl.generate_grain_geometry()
cv2.imshow('All Contours of Regressed Stars', base_img)
cv2.waitKey(0)
cv2.destroyAllWindows()

all_contours = []
contours = finocyl.find_contour(base_img)
all_contours.append(contours)
processed_img = base_img.copy()
outer_diameter = int(finocyl.mapDim/1.5)
outer_radius = outer_diameter // 2
center = (processed_img.shape[0] // 2, processed_img.shape[1] // 2)

edgeFlag = False
i = 0
while True:
    # Apply regression
    disk_filtered_img = finocyl.apply_disk_filter(processed_img, 20)
    thresholded_img = finocyl.apply_threshold(processed_img, disk_filtered_img, 0.36)

    mask = np.zeros_like(thresholded_img)
    cv2.circle(mask, center, outer_diameter // 2, (255), -1)
    masked_thresholded_img = cv2.bitwise_and(thresholded_img, thresholded_img, mask=mask)

    # Check termination condition
    contours = finocyl.find_contour(thresholded_img)

    if not edgeFlag and (not contours or any(finocyl.contour_intersects_boundary(cont, center, outer_radius) for cont in contours)):
        edgeFlag = True

    if edgeFlag and (not contours or all(finocyl.contour_intersects_boundary2(cont, center, outer_radius) for cont in contours)):
        break

    contours_masked = finocyl.find_contour(masked_thresholded_img)

    # Append contours
    all_contours.append(contours_masked)

    # Prepare for next iteration
    processed_img = thresholded_img.copy()

    img = masked_thresholded_img.copy()
    # if len(img.shape) == 2:
    img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    cv2.circle(img, center, outer_diameter // 2, (0, 0, 255), 3)
    cv2.drawContours(img, contours_masked, -1, (0, 255, 0), 3)
    cv2.imshow('All Contours of Regressed Stars', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

colored_img = finocyl.draw_all_contours(base_img, all_contours)

cv2.circle(colored_img, center, outer_diameter // 2, (0, 0, 255), 3)
cv2.imshow('All Contours of Regressed Stars', colored_img)
cv2.waitKey(0)
cv2.destroyAllWindows()
