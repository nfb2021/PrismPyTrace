import numpy as np
import math
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR) # deals with a bug in matplotlib
import pandas as pd
import os
from tqdm import trange
from shapely.geometry import LineString, Polygon, MultiLineString, MultiPolygon
import random as rn 
from scipy.optimize import curve_fit
import time

###################################################################################################################################
###################################################################################################################################

'''Huge parts of this code rely on the Shapely package, wich enables straightforward 2D geometry in Python.
This then directly points to a fundamental limit of this code: It cannot be used for three dimensional 
geometric analysis. For further information, consult https://shapely.readthedocs.io/en/stable/index.html

'''

class Prism:

    '''     ^
           / \\
          / ^ \\
         /  |  \\
        /'alpha'\\     <-- lenght of this side is calculated based on 'width' and 'alpha'
       /         \\
       ----------- 
             ^
             |
       This side is defined via 'width'.
       It is parallel to the z-axis of SigRay (or beam axis).
       The angle opposite to this side is 'alpha', so alpha is the full (!) angle.
       'height' is defined as (0.5 * width) / tan(0.5 * alpha).
       Thus, a prism is defined by these two parameters only.
    '''

    def __init__(self, width, alpha):
        self.alpha = alpha
        self.width = width

    def left_side(self, coords):
        '''Returns only the left side of the prism in the form of a 
        shapely.geometry.linestring.LineString object (in principal not required)'''
        x, y = coords
        start = (x, y)
        end = (x + (self.width / 2), y + (self.width / (2*np.tan(math.radians(self.alpha / 2)))))
        return LineString([start, end])

    def right_side(self, coords):
        '''Returns only the right side of the prism in the form of a 
        shapely.geometry.linestring.LineString object (in principal not required)'''
        x, y = coords
        start = (x + (self.width / 2), y + (self.width / (2*np.tan(math.radians(self.alpha / 2)))))
        end = (x + self.width, y)
        return LineString([start, end])

    def bottom(self, coords):
        '''Returns only the bottom of the prism in the form of a 
        shapely.geometry.linestring.LineString object (in principal not required)'''
        x, y = coords
        start = (x , y)
        end = (x + self.width, y)
        return LineString([start, end])
        
    def prism_whole(self, coords):
        '''Returns a whole prism in the form of a 
        shapely.geometry.polygon.Polygon object'''
        x, y = coords
        return Polygon([(x, y), (x + self.width, y), (x + (self.width / 2), y + (self.width / (2*np.tan(math.radians(self.alpha / 2))))), (x, y)])
            
    def get_centroid(self, prism_polygon):
        '''Return the centroid of a shapely.geometry.polygon.Polygon object'''
        return prism_polygon.centroid.x, prism_polygon.centroid.y

    def get_height(self):
        '''Returns the height of the prisms'''
        return self.width / (2 * (np.tan(math.radians(self.alpha / 2))))


class Snell:

    '''Takes one refractive index and a refractive index decrement (e.g. n_air, delta) and calculates the second refractive index from this.
    Using Snell's law and the specified input angle (in degrees), the second angle is calculated accordingly.
    Furthermore, the corresponding critical angle can be calculated.
    All input and output is un units of degrees.
    For further infomration, consult https://en.wikipedia.org/wiki/Prism#Deviation_angle_and_dispersion
    '''

    def __init__(self, angle_in, n_air, delta):
        self.angle_in = math.radians(angle_in)
        self.n_air = n_air
        n_pol = n_air - delta
        self.n_pol = n_pol

    def air_to_polymer(self, angle_in):
        '''Returns the outbound angle (inside the prism, made of polymer) relative to the normal of the interface.'''
        return math.degrees(np.arcsin(np.sin(math.radians(angle_in)) * (self.n_air / self.n_pol) ))

    def polymer_to_air(self, angle_in):
        '''Returns the outbound angle (outside the prism) relative to the normal of the interface.'''
        return math.degrees(np.arcsin(np.sin(math.radians(angle_in)) * (self.n_pol / self.n_air)))

    def critical_angle(self):
        '''Returns the critical angle. If light coming from the material of larger refractive index hits the interface to a material of smaller
        refrative index under the critical angle, then total internal refelction occurs.'''
        return math.degrees(np.arcsin(self.n_air / self.n_pol))



class PrismDistribution(Prism):
    '''This class can be used to obtain different kinds of distributions of prisms, namely
    one-dimensional arrays and two-dimensional forests.
    Note: All returned or input coordinates of prisms are NOT the centroid coordinates, but the lower left corner
    (see schematic above).
    '''

    def __init__(self, n_prisms, separation, width, alpha, forest_length, forest_width):
        Prism.__init__(self, width, alpha)

        self.separation = separation
        self.n_prisms = n_prisms
        self.forest_length = forest_length
        self.forest_width = forest_width

    def export_to_blender(self, coordinates_list, array_type, prism_height = 100):
        root_dir = os.getcwd()
        timestamp = time.time()
        with open(f"{timestamp}_Forest.txt", "w") as f:
            f.write("x;y;z;width;alpha\n")
            for coords in coordinates_list:  # Note: These are the coordinates of the left corner of the prism, NOT the centroid
                x, y = coords
                f.write(f"{x};{y};{0};{self.width};{self.alpha}\n")

        pyfile = os.path.join(root_dir, str(timestamp) + f"_{array_type}.txt")
        with open("Blender_Forest.py", "w") as f:
            f.write(r"# >>> This file contains the code for the prism forest to be implemented in Blender <<<")
            f.write("\n")
            f.write(r"# >>> In order for this to function, you need the library 'pyblend_prism.py'. You can get this file from https://github.com/nfb2021/PrismPyTrace/tree/main/PyBlend <<<")
            f.write("\n")
            f.write(r"# >>> Then,  copy the file to 'C:\Program Files\Blender Foundation\Blender x.x\x.x\scripts\modules\pyblend_prism.py' <<<")
            f.write("\n")
            f.write(r"# >>> Once this is done, open Blender and click on 'General' to open a new work space <<<")
            f.write("\n")
            f.write(r"# >>> Click on one of the two smaller windows on the right size of the blender main window (Should contain 'scene collection' and 'cube') <<<")
            f.write("\n")
            f.write(r"# >>> Using the cursor you can adapt the size of each window accoridng to your liking <<<")
            f.write("\n")
            f.write("# >>> Then, after selecting one of the two windows, press 'Shift' +  'F11' to open the editor <<<")
            f.write("\n")
            f.write(r"# >>> Click on 'open' and select the 'Blender_Forest.py' file <<<")
            f.write("\n")
            f.write(r"# >>> All code required is written below. Just run the script in Blender using 'Alt' + 'P' <<<")
            f.write("\n")
            f.write("\n")
            f.write("\n")
            f.write("from pyblend_prism import Prism\n")
            f.write("import os\n")
            f.write(f"prism = Prism(20, {prism_height}, 60) # initialize class with some standard values\n") 
            f.write("prism.clear_scene()\n")
            f.write(f"with open(r'{pyfile}', 'r') as f:\n")
            f.write("\tlines = f.readlines()\n")
            f.write("\tfor l, line in enumerate(lines):\n")
            f.write("\t\tif l == 0:\n")              
            f.write("\t\t\tcontinue\n\t\t")
            f.write(r"x = float(line.split(';')[0])")
            f.write("\n\t\t")
            f.write(r"y = float(line.split(';')[1])")
            f.write("\n\t\t")
            f.write(r"z = float(line.split(';')[2])")
            f.write("\n\t\t")
            f.write(r"width = float(line.split(';')[3])")
            f.write("\n\t\t")
            f.write(r"alpha = float(line.split(';')[4])")
            f.write("\n\t\t")
            f.write("new_loc = (x, y, 0)")
            f.write("\n\t\t")
            f.write("prism_obj = prism.define_prism(loc = new_loc, angle = alpha, base_width = width)\n")
            f.write("\n\t\t")
            f.write("prism.link_prism(prism_obj)\n")
            f.write("\n")

    def get_h_array(self, start_coords, y_offset = 0):
        '''Returns a horizontal array of prisms, based on the starting coordinate and the earlier specified number of prisms and the separation between two adjacent prisms.
        By default, the y-offset is set to zero. This parameter allows to shift the whole array vertically in the xy-plane.
        The returned array is of type 'list' and contains prisms of type shapely.geometry.polygon.Polygon object.
        '''
        prism_arr = []
        for p in range(self.n_prisms):
            if p == 0:
                xs, ys = start_coords
                if ys != y_offset and y_offset != 0:
                    raise Warning(f"The chosen y-offset of {y_offset} does not match the chosen start coordinates {xs, ys}.")
                coords = start_coords
            else:
                coords = (p * (self.width + self.separation), y_offset)
           
            prism = self.prism_whole(coords)
            prism_arr.append(prism)

        return prism_arr

    def get_v_array(self, start_coords, x_offset = 0):
        '''Returns a vertical array of prisms, based on the starting coordinate and the earlier specified number of prisms and the separation between two adjacent prisms.
        By default, the x-offset is set to zero. This parameter allows to shift the whole array horizontally in the xy-plane.
        The returned array is of type 'list' and contains prisms of type shapely.geometry.polygon.Polygon object.
        '''
        prism_arr = []
        for p in range(self.n_prisms):
            if p == 0:
                xs, ys = start_coords
                if xs != x_offset and x_offset != 0:
                    raise Warning(f"The chosen x-offset of {x_offset} does not match the chosen start coordinates {xs, ys}.")
                coords = start_coords
            else:
                coords = (x_offset, p * (self.width + self.separation))
           
            prism = self.prism_whole(coords)
            prism_arr.append(prism)

        return prism_arr

    def get_random_forest(self,  seed, threshold, blender_export, shiftx = 0, shifty = 0, patches_list = None):
        ''' Main function for obtaining a 2D random distribution of prisms. Prisms are not allowed to overlap and have all the same 
        size, angle and orientation. Is called by function "get_random_patchy_forest" below, therefore all input parameters defined there
        will be used here.
        Works also without calling "get_random_patchy_forest" before. First input parameter is a seed to initialize the 
        random number generator, which is used to distribute the prisms. The threshold paramater defines the number of attempts to place 
        a prism without overlapping with another prism. Thirdly, the blender-export paramter is a boolean and if True generates a secondary
        python script which can be called from within blender to directly render the just generated prism forest. To do so, a text file containing
        the lower left corner coordinates of each prism is exported, too. Furhter information on the implementation with blender are found 
        in the corresponding generated python script's header. Both shift parameters allow to shift the whole forest to different coordinates.
        Lastly, the patches_list paramater, by default "None", allows to split the forest into pre-defined patches with boundaries in between.
        This then enables individual 3D printing of each of the patches, without a pushing the limits of the printer.
        '''

        min_x, max_x = 0, self.forest_length + shiftx - self.width
        min_y, max_y = 0, self.forest_width + shifty - self.get_height()

        coords_list = []
        rn.seed = seed      

        if patches_list == None:    # no patches for printing
            prism_forest = []   # will contain all prisms objects
            # an intial prism for the whole forest is now generated
            # without it no MultiPolygon could be constructed furter below
            coords = rn.uniform(min_x, max_x), rn.uniform(min_y, max_y) 
            init_prism = self.prism_whole(coords)
            prism_forest.append(init_prism)
            coords_list.append(coords)  # will contain all coordinates to be output as .txt file for blender visualization
            current_n_prisms = 1

        else: 
            prism_forest = list(np.copy(patches_list))
            # for patch in patches_list:            # uncommenting those two lines will plot the patches
            #     plt.plot(*patch.exterior.xy, color = "gold")
            current_n_prisms = 0


        # to ensure there is no overlap between prisms a./o. patches' boundaries, a simple approach is used.
        # all prisms as well as the patched themselves, all are shapely.geometry.polygon.Polygon objects, are combined into a
        # single shapely.geometry.polygon.MultiPolygon object. Then, for each iterattion of the while loop a proposal for the next prism's
        # location is made. By checking if the proposed prism would overlap with the MultiPolygon object, this proposal is accepted or rejected.

        super_prism = MultiPolygon(prism_forest).buffer(0) # buffer(0) prevents funny behavior

        iterations = 0
        while current_n_prisms <= self.n_prisms - 1:
            iterations += 1
            print(f"iteration\t{iterations}, prism no.\t{current_n_prisms}", end = "\r")
            coords = rn.uniform(min_x, max_x), rn.uniform(min_y, max_y) # proposal for next prisms's coordinates
            try_prism = self.prism_whole(coords)
            overlap = try_prism.intersection(super_prism)   # contains if prisms and MultiPolygon do or do not overlap
            # in case they overlap, an intersection is found.
            # weirdly, the class used by shapely to store this information is either of type "Polygon" or "MultiPolygon".
            # since both classes need to be dealt with differently, the proposal is simply direclty rejected if it happens to be of
            # the MultiPolygon class

            if overlap.geom_type == 'Polygon' and len(overlap.exterior.coords) == 0:
                prism_forest.append(try_prism)
                coords_list.append(coords)
                super_prism = MultiPolygon(prism_forest).buffer(0)
                current_n_prisms += 1
                iterations = 0

            if iterations == threshold:
                print(f"\n\nThere was not enough space in the forest to place the intended {self.n_prisms} prisms.")
                print(f"Instead, {current_n_prisms} prisms were placed. If this is not what you want, you should increase the forest size accordingly.\n\n")
                break
            
        if blender_export == True:
            self.export_to_blender(coords_list, "Forest", prism_height = 100)

        if patches_list == None:
            return prism_forest

        else:
            _list = [x for x in prism_forest if x.geom_type == "Polygon" and len(x.exterior.coords) != len(patches_list[0].exterior.coords)] # sorting out weird things, just in case
            return _list



    def get_random_patchy_forest(self, seed = 2021, threshold = 1000, nanoscribe_patch_size = (250, 250), blender_export = False):
        ''' This function should be the go-to when dealing with prism forests. All input parameter sdefined here are only passed on
        to the function above. 
        Creates a list wich contains rectangles in the form of shapely.geometry.polygon.Polygon objects, representing the boundary regions of the patches 
        The return of this function then calls the get_random_forest() function above and passes all input parameters as well as said list on.       
        '''
        patchx, patchy = nanoscribe_patch_size

        if patchx >= self.forest_length and patchy >= self.forest_width: # in this case a single patch is larger than the whole forest anyways
            print(f"\nThe chosen patch size of {patchx}x{patchy}um is larger than/equals the desired forest itself, which is {self.forest_length}x{self.forest_width}um.")
            print("Thus, the forest will not contain patches.\n")
            return self.get_random_forest(seed = seed, threshold = threshold, patches_list = None, blender_export = blender_export)
        

        else:
            boundary_width = 10     # width of the boundary region between patches
            boundary_regions = []
            xshift, yshift = 0, 0
            
            for x in range(self.forest_length):
                if x % patchx == 0 and x != 0:
                    boundary_v = Polygon([(x + xshift, -0.5*self.forest_width), (x + boundary_width + xshift, -0.5*self.forest_width), (x + boundary_width + xshift, 1.5*self.forest_width), (x + xshift, 1.5*self.forest_width), (x + xshift, -0.5*self.forest_width)])
                    boundary_regions.append(boundary_v)
                    xshift += boundary_width

            for y in range(self.forest_width):
                if y % patchy == 0 and y != 0:
                    boundary_h = Polygon([(-0.5*self.forest_length, y + yshift), (1.5*self.forest_length, y + yshift), (1.5*self.forest_length, y + boundary_width + yshift), (-0.5*self.forest_length, y + boundary_width + yshift), (-0.5*self.forest_length, y + yshift)])
                    boundary_regions.append(boundary_h)
                    yshift += boundary_width

            return self.get_random_forest(shiftx = xshift, shifty = yshift, seed = seed, threshold = threshold, patches_list = boundary_regions, blender_export = blender_export)
        

    def get_avg_vert_refr_idx(self, polygon_list, stepsize, n_air, delta):
        '''This function takes an array or forest or prisms as input and returns the average refractive index for every y-coordinate
        along the x-axis. The second input parameter defines the sampling distance in x-direction. Lastly, the refractive index of the 
        surrounding medium and the refractive index decrement of the prisms are to be specified.
        Returns the scan steps and average refractive index for each step the mean and standard deviation of all average refr. indices
        as well as the amount of scans which resulted in non-sense data
        '''
        super_polygon = MultiPolygon(polygon_list).buffer(0)
        minx, miny, maxx, maxy = super_polygon.bounds 
        x_offset = (maxx - minx) * 0.2
        y_offset = (maxy - miny) * 0.2
        steps = int((maxx - minx + x_offset) / stepsize)
        # print(steps)
        scan_steps, avg_refr_idx = [], []
        thrown_away = 0
        for s in trange(steps):
            # plt.figure(1) #enable all matplotlib commands to obtain a plot of each scan line

            # the principle is easy: First creating a scan line as shapely.geometry.linestring.LineString object
            v_scan_line = LineString([(minx + stepsize * s - 0.5*x_offset, miny - 0.5*y_offset), (minx + stepsize * s - 0.5*x_offset, maxy + 0.5*y_offset)])
            # then, the intersection of the scan line with the MultiPolygon is calculated
            intersections = v_scan_line.intersection(super_polygon)
            
            # unfortunately, the variable "intersections" is sometimes of shapely.geometry.linestring.LineString,
            # other times of shapely.geometry.linestring.MultiLineString type. Since both objects offer different methods,
            # a case distinction needs to be made

            # in any case, the assumption is: any intersection point of odd index must be the begining of a prism,
            # whereas any intersection point of even index must be the end of a prism
            # since the overall dimensions of the forest are known, this then allows for calculating the distance inside and outside
            # of prisms along the scan line at a given x-coordinate

            if intersections.geom_type == "LineString":
                if len(intersections.coords) == 0:  # no intersections found <-> no prisms at this x-coordinate
                    avg_refr_idx.append(n_air)
                    scan_steps.append(minx + stepsize * s - 0.5*x_offset)
                
                if len(intersections.coords) > 0: # there are prisms at this x-coordinate
                    if len(intersections.coords) % 2 == 0: # as mentioend above, the assumption of odd and even indices
                        prism_length = 0
                        for i in range(int(len(intersections.coords) / 2)):
                            _ = LineString([(intersections.coords[0][0], intersections.coords[0][1]), (intersections.coords[1][0], intersections.coords[1][1])])
                            prism_length += _.length
                            list(intersections.coords).pop(0), list(intersections.coords).pop(1)
                            vacuum_length = self.forest_width - prism_length
                            # plt.plot([intersections.coords[0][0], intersections.coords[1][0]], [intersections.coords[0][1], intersections.coords[1][1]])
                        avg_refr_idx.append((prism_length * (n_air - delta) + vacuum_length * n_air) / (prism_length + vacuum_length))
                        scan_steps.append(minx + stepsize * s - 0.5*x_offset)
                    
                    else: # if there is an overall odd amount of intersectiosn found, which should not be 
                        thrown_away += 1
                        

            elif intersections.geom_type == "MultiLineString": # see above
                prism_length = 0
                for line in intersections:
                    if len(line.coords) % 2 == 0:
                        for i in range(int(len(line.coords) / 2)):
                            _ = LineString([(line.coords[0][0], line.coords[0][1]), (line.coords[1][0], line.coords[1][1])])
                            prism_length += _.length
                            list(line.coords).pop(0), list(line.coords).pop(1)
                            vacuum_length = self.forest_width - prism_length
                            # plt.plot([line.coords[0][0], line.coords[1][0]], [line.coords[0][1], line.coords[1][1]])

                avg_refr_idx.append((prism_length * (n_air - delta) + vacuum_length * n_air) / (prism_length + vacuum_length))
                scan_steps.append(minx + stepsize * s - 0.5*x_offset)      

            else:
                print(f"Something weird happened. Geometry of type '{intersections.geom_type}' identified.")


        # sorting out all average indices stemming from scan areas outside of the forest, then calculate the mean and std
        avg_refr_idx_clipped = list(np.copy(avg_refr_idx))

        for i in range(len(avg_refr_idx_clipped)):
            if avg_refr_idx_clipped[i] == float(n_air):
                avg_refr_idx_clipped.pop(i)
            else:
                break

        for i in range(len(avg_refr_idx_clipped))[::-1]:
            if avg_refr_idx_clipped[i] == float(n_air):
                avg_refr_idx_clipped.pop(i)
            else:
                break
        
        mean = np.mean(avg_refr_idx_clipped)
        std = np.std(avg_refr_idx_clipped)

        return scan_steps, avg_refr_idx, mean, std, thrown_away

    def get_avg_hori_refr_idx(self, polygon_list, stepsize, n_air, delta):
        '''See above: get_avg_vert_refr_idx(). The same principal, only now scan lines are in x-direction and scan along the y-axis
        '''
        super_polygon = MultiPolygon(polygon_list).buffer(0)
        minx, miny, maxx, maxy = super_polygon.bounds 
        x_offset = (maxx - minx) * 0.2
        y_offset = (maxy - miny) * 0.2
        steps = int((maxy - miny + y_offset) / stepsize)
        # print(steps)
        scan_steps, avg_refr_idx = [], []
        thrown_away = 0
        for s in trange(steps):
            # plt.figure(1)
          
            h_scan_line = LineString([(minx - 0.5*x_offset, miny + stepsize * s - 0.5*y_offset), (maxx + 0.5*x_offset, miny + stepsize * s - 0.5*y_offset)])           
            intersections = h_scan_line.intersection(super_polygon)
            
            if intersections.geom_type == "LineString":
                if len(intersections.coords) == 0:
                    avg_refr_idx.append(n_air)
                    scan_steps.append(miny + stepsize * s - 0.5*y_offset)
                
                if len(intersections.coords) > 0:
                    if len(intersections.coords) % 2 == 0:
                        prism_length = 0
                        for i in range(int(len(intersections.coords) / 2)):
                            _ = LineString([(intersections.coords[0][0], intersections.coords[0][1]), (intersections.coords[1][0], intersections.coords[1][1])])
                            prism_length += _.length
                            list(intersections.coords).pop(0), list(intersections.coords).pop(1)
                            vacuum_length = self.forest_length - prism_length
                            # plt.plot([intersections.coords[0][0], intersections.coords[1][0]], [intersections.coords[0][1], intersections.coords[1][1]])
                        avg_refr_idx.append((prism_length * (n_air - delta) + vacuum_length * n_air) / (prism_length + vacuum_length))
                        scan_steps.append(miny + stepsize * s - 0.5*y_offset)
                    
                    else:
                        thrown_away += 1
                        

            elif intersections.geom_type == "MultiLineString":
                prism_length = 0
                for line in intersections:
                    if len(line.coords) % 2 == 0:
                        for i in range(int(len(line.coords) / 2)):
                            _ = LineString([(line.coords[0][0], line.coords[0][1]), (line.coords[1][0], line.coords[1][1])])
                            prism_length += _.length
                            list(line.coords).pop(0), list(line.coords).pop(1)
                            vacuum_length = self.forest_length - prism_length
                            # plt.plot([line.coords[0][0], line.coords[1][0]], [line.coords[0][1], line.coords[1][1]])

                avg_refr_idx.append((prism_length * (n_air - delta) + vacuum_length * n_air) / (prism_length + vacuum_length))
                scan_steps.append(miny + stepsize * s - 0.5*y_offset)              

            else:
                print(f"Something weird happened. Geometry of type '{intersections.geom_type}' identified.")
                
        avg_refr_idx_clipped = list(np.copy(avg_refr_idx))

        for i in range(len(avg_refr_idx_clipped)):
            if avg_refr_idx_clipped[i] == float(n_air):
                avg_refr_idx_clipped.pop(i)
            else:
                break

        for i in range(len(avg_refr_idx_clipped))[::-1]:
            if avg_refr_idx_clipped[i] == float(n_air):
                avg_refr_idx_clipped.pop(i)
            else:
                break

        mean = np.mean(avg_refr_idx_clipped)
        std = np.std(avg_refr_idx_clipped)

        return scan_steps, avg_refr_idx, mean, std, thrown_away


    def get_intersection_list(self, ray, super_polygon, ray_tracing_class):
        '''Unfinished, also does not belong in this class since it is concerned with ray tracing '''
        out_list = []
        try:
            intersec = ray.intersection(super_polygon)      
            if intersec.geom_type == "LineString":
                if len(intersec.coords) > 0:
                    for ii in intersec.coords:
                        out_list.append(ii)
                    return out_list

                else:
                    print("No intersections found")
                    return None


            if intersec.geom_type == "MultiLineString":
                if len(intersec) > 0:
                    for ls in intersec:
                        for ii in ls.coords:
                            out_list.append(ii)
                        return out_list

                else:
                    print("No intersections found")
                    return None  

        except:
            return False


class AngleConverter(PrismDistribution, Snell):

    '''Angles calculated in Snell's law are always relative to the normal of the corresponding surface, at which refraction takes place.
    This is inpractical, though.
    For easier handling, this converter uses the top angle (alpha, in degrees) of the prisms and converts between the angle relative to the optical axis 
    (which would be the z-axis in the SigRay setup) and angles from Snell's law.
    '''

    def __init__(self, n_air, delta, n_prisms, separation, width, alpha, forest_length, forest_width):
        PrismDistribution.__init__(self, n_prisms, separation, width, alpha, forest_length, forest_width)
        self.n_air = n_air
        self.d = delta
        n_pol = n_air - self.d
        self.n_pol = n_pol

    def opt_ax_to_snell(self, theta_inc):   # convert from angles relative to the optical axis to angles relative to the normal of the surface where refraction is to take place
        return theta_inc + (self.alpha / 2)

    def snell_to_opt_ax(self, theta_out):   # convert from angles relative to the normal of the surface where refraction is to take place to angles relative to the optical axis 
        return (self.alpha / 2) - theta_out

    def get_snell_angles(self, theta_inc_snell, theta_1): # returns all four Snell angles for a ray travelling through a prism
        theta_2 = Snell.air_to_polymer(self, theta_1)
        theta_3 = self.alpha - theta_2
        theta_4 = Snell.polymer_to_air(self, theta_3)

        theta_rel = theta_4 - self.alpha
        delta_theta_rel = theta_inc_snell + theta_rel

        delta_theta_abs = self.snell_to_opt_ax(theta_4)
        return theta_1, theta_2, theta_3, theta_4, delta_theta_rel, delta_theta_abs


class RayTracing: #(AngleConverter):
    '''This class is concerned with ray tracing. As of Nov. 08, 2021 there are a few limitations:
    1. Ray Tracing was only implemented for horizontal arrays of prisms
    2. The case of a ray hitting the bottom of a prism is not defined (because it is very unlikely with X-rays)
    3. Fundamental limitation of the shapely library: restriction to 2D geometrical problems
    '''

    def __init__(self, theta_inc):
        # AngleConverter.__init__(self)
        self.theta_inc = theta_inc      # relative to optical axis

    def get_source_ray(self, theta_inc, b, source_coords, length):
        #returns the starting ray (source -> first prism), shapely.geometry.linestring.LineString object
        x, y = source_coords
        m = np.tan(math.radians(theta_inc))
        end_coords = (x + length, m * (x + length) + b)
        return LineString([source_coords, end_coords])

    def get_ray(self, angle, intersec_coords, length):
        # returns ray based on it's origin (intersection corrdinates), slope (angle, rel. to opt. ax.) and length,
        # shapely.geometry.linestring.LineString object
        x, y = intersec_coords
        m = np.tan(math.radians(angle))
        new_offset = y - m * x
        end_coords = (x + length, m * (x + length) + new_offset)
        return LineString([intersec_coords, end_coords])

    def get_plot_from_LineString(self, linestring):
        x_list = [linestring.coords[0][0], linestring.coords[1][0]]
        y_list = [linestring.coords[0][1], linestring.coords[1][1]]
        return x_list, y_list

    def to_dataframe(self, df, prism_centroid, width, separation, alpha, source, theta1, theta2, theta3, theta4, thetarel, deltatheta, intlx, intly, intrx, intry):
        # creates dataframe based on inout data
        cols =      [
                    "Prism Centroid X [um]", "Prism Centroid Y [um]", 
                    "Width [um]", "Separation [um]", "Alpha [deg]", 
                    "Source X [um]", "Source Y [um]", 
                    "Theta 1 [deg]", "Theta 2 [deg]", "Theta 3 [deg]", "Theta 4 [deg]", 
                    "Delta Theta rel. [deg]", "Delta Theta abs. [deg]", 
                    "Intersec In X [um]", "Intersec In Y [um]", 
                    "Intersec Out X [um]", "Intersec Out Y [um]"
                    ]
            
        prism_centroid_x, prism_centroid_y = prism_centroid
        source_x, source_y = source
        data = [prism_centroid_x, prism_centroid_y, width, separation, alpha, source_x, source_y, theta1, theta2, theta3, theta4, thetarel, deltatheta, intlx, intly, intrx, intry]


        temp_df = pd.DataFrame([data], columns = cols)
        df = df.append(temp_df)
        
        return df

class Detector:
    '''Initialized using the pixel size of the detector used (Sigray lab, Si det.: 55um) and its distance.
    The angle gamma is the angle relative to the opt. ax. the ray leaves the last prism with'''
    def __init__(self, gamma, distance, pixelsize):
        self.gamma = gamma
        self.distance = distance
        self.pixelsize = pixelsize
    
    def get_spot_difference(self):
        # returns the distance between an unrefracted, direct beam parallel to the opt. ax. and a refracted beam.
        # makes only sense when used for incident angles rel. to opt. ax. of 0 deg
        spots = np.tan(math.radians(self.gamma)) * self.distance
        pixels = spots / self.pixelsize
        return spots, pixels


class MirageDeflector:
    '''Based on the idea of H. Chapman, an array of prisms without separation between two adjacent prisms can be dealt with as a
    mirage deflector. See supplementary files.
    Only the last function needs to be called. Returns the total deflection expected.'''
    def __init__(self, theta_inc, delta, offset, alpha, width, separation, n_prisms):
        self.theta0 = theta_inc
        self.t = width
        self.delta = delta
        self.y0 = offset
        self.alphaHC = alpha / 2   
        self.separation = separation
        self.n_prisms = n_prisms


    def linear_profile_gradient(self):
        return 2 * self.delta * math.degrees(np.tan(math.radians(self.alphaHC))) / self.t

    def deflector_length(self):
        return self.n_prisms * (self.t + self.separation) - self.separation

    def get_capital_delta(self):
        # retuns total deflection expected
        g = self.linear_profile_gradient()
        a = g * self.deflector_length() / (1 - g * self.y0)
        b = self.theta0 * (g * self.deflector_length())**2 / (2 * (1 - g * self.y0)**2)
        return a + b

class PyhsicalProperties:
    '''This class is used mainly to deal with Lambert Beer's law and optical and geometric path lengths.
    The input dataframe is a dataframe resulting from a completed ray tracing simulation, in the style of
    the return of RayTracing.to_dataframe(). delta and beta are the real and imagrinary part of the complex refractive index'''
    def __init__(self, dataframe, n_air, delta, beta):
        self.dataframe = dataframe
        self.n_air = n_air
        self.delta = delta
        self.n_pol = n_air - delta
        self.beta = beta

    def get_gpl_and_opl(self, all_path_lengths = False):
        '''Based on the obtained intersection coordinates,
        both the Geometric Path Length (GPL) and Optical Path Lengt (OPL)
        are calculated.
        Using the boolean input parameter, only the GPL throughout all prisms
        can be returned. Such an array could then be used e.g. as input array
        for the following function down below.
        '''

        geom_path_length_tot = ()
        geom_path_length_prism = 0
        geom_path_length_air = 0
        
        for i, (inx, outx, iny, outy) in enumerate(zip(self.dataframe["Intersec In X [um]"], self.dataframe["Intersec Out X [um]"], self.dataframe["Intersec In Y [um]"], self.dataframe["Intersec Out Y [um]"])):
            geom_path_length_prism += np.sqrt((outx-inx)**2 + (outy-iny)**2)
            geom_path_length_tot += ((inx, iny),) 
            geom_path_length_tot += ((outx, outy),) 


        geom_path_length_tot = MultiLineString([geom_path_length_tot]).length
        geom_path_length_air = geom_path_length_tot - geom_path_length_prism

        opt_path_length_prism = geom_path_length_prism * self.n_pol
        opt_path_length_air = geom_path_length_air * self.n_air
        opt_path_length_tot = opt_path_length_prism + opt_path_length_air

        # print(f"\n\t Prism:\tGeom. Path Length: {geom_path_length_prism}um,\t Optical Path Length: {opt_path_length_prism}um")
        # print(f"\t Air:\tGeom. Path Length: {geom_path_length_air}um,\t Optical Path Length: {opt_path_length_air}um")
        # print(f"\t ___________________________________________________________________________________________________")
        # print(f"\t Total:\tGeom. Path Length: {geom_path_length_tot}um,\t Optical Path Length: {opt_path_length_tot}um\n")

        # deflection = data.append(list(self.dataframe["Delta Theta abs. [deg]"])[-1]) # this is not necessary here, but might be useful
        if all_path_lengths == True:
            return geom_path_length_tot, opt_path_length_tot, geom_path_length_air, opt_path_length_air, geom_path_length_prism, opt_path_length_prism
        
        else:
            return geom_path_length_tot, opt_path_length_tot


    def get_intensity(self, step = 1, intensity_0 = 1, threshold = 10):
        '''Automatically calls the function above to get the GPL inside the prisms.
        The stepsize for the array resulting from the GPL is by default set to 1.
        The starting intensity is set to 1 by default, too.
        Assumes no absorption of photons while the beam is in air in between prisms.
        Also, an attempt to calculate the half-life of the incident Intensity is done.
        This function returns an array of resulting intensities, the half-life (tubple of floats or None)
        and a boolean, indicating if the half-life was reached inside the specified input array 
        ( == if the half-life was reached during the geometric path lenght in the prisms).
        '''
        
        _1, _2, _3, _4, geom_path_length_prism, _5 = self.get_gpl_and_opl(all_path_lengths = True)
        distances = np.arange(0, geom_path_length_prism + step, step)
        
        def f(x):
            return intensity_0 * np.exp(-1 * self.beta * x)

        def f_fit(x, a, c):
            return a * np.exp(-1 * self.beta * x) + c

        out_list = [f(x) for x in distances]       
        popt, pcov = curve_fit(f_fit, distances, out_list)
        a, c = popt
        distances_high_res = np.arange(distances[0], distances[-1], 0.01 * (abs(distances[0] - distances[-1]) / len(distances)))
        out_list_high_res = [f_fit(x, a, c) for x in distances_high_res]
      
        for i, (ii, jj) in enumerate(zip(distances_high_res, out_list_high_res)):
            try:
                if out_list_high_res[i] >= out_list[0] / np.exp(1) > out_list_high_res[i + 1]:
                    half_life = (ii, jj)
                    found_inside_array = True
                    break

            except(IndexError):
                print("\n\tThe half-life of the Intensity's exponential decay was not reached yet.")
                print("\tTrying to determine where it would approximately be...\n")

                x = distances[-1]
                check = False
                while check == False:
                    next_x = x + 0.5 * (abs(distances[0] - distances[-1]) / len(distances))
                    y_i = f_fit(x, a, c)
                    y_i_plus_1 = f_fit(next_x, a, c)

                    if y_i >= out_list[0] / np.exp(1) > y_i_plus_1:
                        half_life = (x, y_i)
                        found_inside_array = False
                        check == True
                        
                    if x >= threshold * (abs(distances[0] - distances[-1])):
                        print(f"\n\tThe calculation was stopped, since after {threshold} times the range of the input array the half-life would still not be reached.")
                        found_inside_array = False
                        half_life = None, None
                        break

                    else:
                        x = next_x


        return distances, out_list, half_life, found_inside_array
