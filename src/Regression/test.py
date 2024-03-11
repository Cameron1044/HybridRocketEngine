from abc import abstractmethod

import cv2
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt



class Regression():
    def __init__(self, outerDiameter, mapDim=1000, threshold=0.0, diskFilterRadius=20):
        self.outerDiameter = outerDiameter
        self.mapDim = mapDim
        self.gridSize = 1.5
        self.pixelDiameter = int(self.mapDim / self.gridSize)
        self.mapX, self.mapY = np.meshgrid(np.linspace(-self.gridSize, self.gridSize, mapDim), np.linspace(-self.gridSize, self.gridSize, mapDim))
        self.threshold = threshold
        self.diskFilterRadius = diskFilterRadius
        self.rdot = self.calcR()

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
        return self.outerDiameter * (value / self.pixelDiameter)

    def mapToArea(self, value):
        """Used to convert sq pixels to sqm. For extracting real areas from the regression map."""
        return (self.outerDiameter ** 2) * (value / ((self.pixelDiameter) ** 2))

    def find_contour(self, img):
        # Find contours
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        return contours

    def draw_all_contours(self, base_img, all_contours):
        # Check if the image is already in color (3 channels)
        if len(base_img.shape) == 3 and base_img.shape[2] == 3:
            colored_img = base_img.copy()
        else:
            # Convert the image to a colored one if it's in grayscale
            colored_img = cv2.cvtColor(base_img, cv2.COLOR_GRAY2BGR)

        # Colors for drawing alternating contours
        colors = [(0, 0, 255), (0, 255, 0), (255, 0, 0), (255, 255, 0), (255, 0, 255), (0, 255, 255)]
        
        # Draw all contours with alternating colors
        for i, contours in enumerate(all_contours):
            cv2.drawContours(colored_img, contours, -1, colors[i % len(colors)], 1)

        return colored_img

    def apply_disk_filter(self, img):
        radius = self.diskFilterRadius + 1
        # Create a disk-shaped kernel
        y, x = np.ogrid[-radius: radius+1, -radius: radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel = np.zeros((2*radius+1, 2*radius+1))
        kernel[mask] = 1
        kernel /= kernel.sum()  # Normalize
        
        # Convolve the image with the kernel
        filtered_img = cv2.filter2D(img, -1, kernel)
        
        return filtered_img

    def apply_threshold(self, img):
        if self.threshold < 0 or self.threshold > 1:
            raise ValueError("Threshold value must be between 0 and 1")
        _, thresholded_img = cv2.threshold(img, self.threshold * 255, 255, cv2.THRESH_BINARY)
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

    def getCorePerimeter(self, contours):
        """Used to calculate the Initial Core Perimeter in pixels """
        # Calculates a countour of the given img
        # contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        # Calculates the initial perimeter of the countour
        firstPass = cv2.arcLength(contours[0], True)
        # Douglas-Peucker Algorithm for Shape Approximation to 0.1 % error
        epsilon = 0.001*firstPass
        DPfit = cv2.approxPolyDP(contours[0],epsilon,True)

        # Recalculates the perimeter of the countour
        perimeter = cv2.arcLength(DPfit, True)

        return perimeter

    def getCoreArea(self, contours):
        """Used to calculate the Initial core area, A_port, in pixels"""
        # Calculates a countour of the given img
        # contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        # Calculates the area of the countour
        area = cv2.contourArea(contours[0])
        
        return area
    
    def calcR(self):
        def count_white_pixels(img):
            return np.sum(img[1, :] == 255)
        half_black_white = np.hstack((np.zeros((self.mapDim, self.mapDim//2), dtype=np.uint8), 255*np.ones((self.mapDim, self.mapDim//2), dtype=np.uint8)))

        white_pixel_count1 = count_white_pixels(half_black_white)

        blurred_img = self.apply_disk_filter(half_black_white)
        thresholded_img = self.apply_threshold(blurred_img)

        white_pixel_count2 = count_white_pixels(thresholded_img)
        return white_pixel_count2 - white_pixel_count1
    
    def showImage(self, img):
        cv2.imshow('image', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()
    
    def runLoop(self, rStop=-1, fuelImagePath=None):
        base_img = self.generate_grain_geometry()
        # self.showImage(base_img)
        base_img = self.generate_grain_geometry()
        real_fuel_img = None
        if fuelImagePath:
            # Attempt to load the real image of burned fuel
            real_fuel_img = cv2.imread(fuelImagePath)
            if real_fuel_img is None:
                print(f"Warning: Unable to open '{fuelImagePath}'. Check the file path and integrity.")
            else:
                # Resize the real image to match the dimensions of the base_img
                real_fuel_img = cv2.resize(real_fuel_img, (self.mapDim, self.mapDim))

        all_contours = []
        contours = self.find_contour(base_img)
        all_contours.append(contours)
        processed_img = base_img.copy()
        outer_diameter = self.pixelDiameter
        outer_radius = int(outer_diameter/2)
        center = (processed_img.shape[0] // 2, processed_img.shape[1] // 2)

        data_columns = [
            "r", "pixels", "area", "perimeter"
        ]
        df = pd.DataFrame(columns=data_columns)

        edgeFlag = False
        r = 0
        while True:
            # Apply regression
            disk_filtered_img = self.apply_disk_filter(processed_img)
            # self.showImage(disk_filtered_img)
            thresholded_img = self.apply_threshold(disk_filtered_img)
            # self.showImage(thresholded_img)

            mask = np.zeros_like(thresholded_img)
            cv2.circle(mask, center, outer_radius, (255), -1)
            masked_thresholded_img = cv2.bitwise_and(thresholded_img, thresholded_img, mask=mask)

            # Check termination condition
            contours = self.find_contour(thresholded_img)

            if not edgeFlag and (not contours or any(self.contour_intersects_boundary(cont, center, outer_radius) for cont in contours)):
                break
                edgeFlag = True

            if edgeFlag and (not contours or all(self.contour_intersects_boundary2(cont, center, outer_radius) for cont in contours)):
                break

            contours_masked = self.find_contour(masked_thresholded_img)

            # Append contours
            all_contours.append(contours_masked)

            perimeter = self.getCorePerimeter(contours_masked)
            area = self.getCoreArea(contours_masked)
            r = r + self.rdot

            new_data = {
                "r": self.mapToLength(r),
                "pixels": r,
                "area": self.mapToArea(area),
                "perimeter": self.mapToLength(perimeter)
            }
            
            # Append the new data to the dataframe
            df.loc[len(df)] = new_data

            # Prepare for next iteration
            processed_img = thresholded_img.copy()

            img = masked_thresholded_img.copy()
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            cv2.circle(img, center, outer_radius, (0, 0, 255), 3)
            cv2.drawContours(img, contours_masked, -1, (0, 255, 0), 3)
            if rStop > 0 and self.mapToLength(r) >= rStop:
                break
            # self.showImage(img)

        df.to_csv("src/Regression/burnback_table.csv", index=False)

        colored_img = self.draw_all_contours(base_img, all_contours)

        if real_fuel_img is not None:
            colored_img = self.draw_all_contours(real_fuel_img, all_contours)  # Pass real_fuel_img instead of base_img
        else:
            colored_img = self.draw_all_contours(base_img, all_contours)
        cv2.circle(colored_img, center, outer_diameter // 2, (0, 0, 255), 3)
        self.showImage(colored_img)

        return df
    
class StarGeometry(Regression):
    def __init__(self, outer_diameter, num_points, point_length, point_base_width, mapDim=1000, threshold=0.36, diskFilterRadius=20):
        super().__init__(outer_diameter, mapDim, threshold, diskFilterRadius)
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
    def __init__(self, outer_diameter, inner_diameter, num_fins, fin_length, fin_width, mapDim=1000, threshold=0.36, diskFilterRadius=20):
        super().__init__(outer_diameter, mapDim, threshold, diskFilterRadius)
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
    
def ToMetric(value, conversion, n=1):
    ## To take in a value and the current unit it is in to change it into metric
    ## INPUTS: value - the value we want to convert
    ##         unit  - the unit that will be changed to metric
    ##         n     - Fuel Regression constant used to calculate a
    ##
    ## OUTPUTS: The Value will have it's units changed into metric / SI units
    ##          If the unit is already in metric, nothing will happen
    conversionDict = {
        ## For Units that are English
        # Length --> m
        "in": 1/39.37,
        "ft": 1/3.281,
        # Surfrace Area --> m2
        "in^2": 1/1550,
        "ft^2": 1/10.764,
        # Volume --> m3
        "in^3": 1/61020,
        "ft^3": 1/35.315,
        # Pressure --> Pa
        "psi": 6895,
        # Mass --> kg
        "lbm": 1/2.205,
        # Force --> N
        "lbf": 4.448,
        # Density --> kg/m^3
        "lbm/in^3": 27680,
        "lbm/ft^3": 16.01846337396,
        # Molecular Weight --> lb/lbmo
        "lb/lbmol": 1/1000,
        # Temperature --> K
        "R": 5/9 ,

        ## For Units that are Metric
        # Length --> m
        "m": 1,
        "cm": 1/100,
        "mm": 1/1000,
        # Surfrace Area --> m2
        "m^2": 1,
        "cm^2": 1/10000,
        "mm^2": (1e-6),
        # Volume --> m3
        "m^3": 1,
        "cm^3": (1e-6),
        "mm^3": (1e-9),
        "L": 1/1000,
        # Pressure --> Pa
        "Pa": 1,
        "kPa": 1000,
        # Mass --> kg
        "kg": 1,
        "g": 1/1000,
        # Force --> N
        "N": 1,
        "kN": 1000,
        # Density --> kg/m^3
        "kg/m^3": 1,
        # Molecular Weight --> kg/mol
        "g/mol": 1/1000,
        "kg/mol": 1,
        # Temperature --> K
        "K": 1,
        # Gas Constant --> J/(mol*K)
        "J/(mol*K)": 1,

        ## Unitless
        "unitless": 1,

        ## For a Burn coefficient
        "a": (0.0254**(1 + 2*(n)))*(0.453592**(-n))
    }
    return value * conversionDict[conversion]

    
outer_diameter = ToMetric(3.375, 'in')
inner_diameter = ToMetric(1.0, 'in')
fin_length = ToMetric(0.9+0.5, 'in')
fin_width = ToMetric(0.2, 'in')
num_fins = 6

# outer_diameter = ToMetric(2, 'in')
# inner_diameter = ToMetric(0.5, 'in')
# fin_length = ToMetric(1.4, 'in')
# fin_width = ToMetric(0.35, 'in')
# num_fins = 0

mapDim = 2500
finocyl = FinocylGeometry(outer_diameter=outer_diameter, inner_diameter=inner_diameter, num_fins=num_fins, fin_length=fin_length, fin_width=fin_width, mapDim=mapDim, threshold=0.36, diskFilterRadius=40)
df = finocyl.runLoop(rStop=0.003636, fuelImagePath="src/Regression/fuelGrainTop2.png")

# plt.figure()
# plt.plot(df['r'], df['area']*1000)
# plt.xlabel('Radius (m)')
# plt.ylabel('Area (m^2)')
# plt.title('Area vs Radius')

# plt.figure()
# plt.plot(df['r'], df['perimeter']*1000)
# plt.xlabel('Radius (m)')
# plt.ylabel('Perimeter (m)')
# plt.title('Perimeter vs Radius')

# plt.show()