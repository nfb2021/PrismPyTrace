import bpy
import numpy as np
import math
import mathutils
import time
import os


class Prism:

    """     ^"""
    """    / \\"""
    """   / ^ \\"""
    """  /  |  \\"""
    """ /'alpha'\\     <-- lenght of this side is calculated based on 'width' and 'alpha'"""
    """/         \\"""
    """----------- """
    """      ^"""
    """      |"""
    """This side is defined via 'width',"""
    """parallel to z-axis of Sigray defined"""
    """The angle opposite to this side is 'alpha'"""
    """'height' defines the distance between the two triangular sides of the prism"""


    def __init__(self, width, height, alpha):
        self.width = width
        self.height = height
        self.alpha = math.radians(alpha)


    def clear_scene(self):
        """This function clears the whole scene and all objects contained in it"""
        bpy.ops.object.select_all(action='SELECT')
        bpy.ops.object.delete(use_global=False)

    
    def define_prism(self, loc = (0, 0, 0), angle = None, base_width = None):
        """The default location assigned is (0, 0, 0). Using the 'update_coordinates'-function allows for reassignment of coordinates"""
        x, y, z = loc
        name = "prism"
        meshes = bpy.data.meshes
        if angle == None:
            angle = self.alpha
        else:
            angle = math.radians(angle)
        if base_width == None:
            base_width = self.width
        else:
            base_width = base_width
        
        

        points = [  [x, y, z], [x + base_width, y, z], [x + (base_width / 2), y + (base_width / (2 * np.tan(angle / 2))), z], 
                    [x, y, z + self.height], [x + base_width, y, z + self.height], [x + (base_width / 2), y + (base_width / (2 * np.tan(angle / 2))), z + self.height]   ]
        faces = [ [4,5,2],[1,0,3],[2,5,3],[4,3,5],[1,2,0],[1,4,2],[4,1,3],[0,2,3] ]

        shape_vertices = []
        for p in points:
            print(p)
            shape_vertices.append ( mathutils.Vector((p[0],p[1],p[2])) )

        new_mesh = bpy.data.meshes.new ( name + "_mesh" )
        new_mesh.from_pydata ( shape_vertices, [], faces )
        new_mesh.update()

        new_obj = bpy.data.objects.new ( name, new_mesh )
        return new_obj

    def link_prism(self, object):
        """Any created object in Blender needs to be linked to the scene, in order to be displayed"""
        bpy.context.collection.objects.link(object)

    def update_coordinates(self, new_location): 
        """This function allows for reassignment of coordinates"""
        return self.define_prism(loc = new_location)

    def update_alpha(self, new_alpha): 
        """This function allows for reassignment of the angle alpha"""
        return self.define_prism(angle = new_alpha)

    def update_width(self, new_width): 
        """This function allows for reassignment of the width of the prism"""
        return self.define_prism(base_width = new_width)
        
    def make_array(self, x, y, no_of_prisms, separation):
        for p in range(no_of_prisms):
            if p == 0:
                self.link_prism(self.update_coordinates((x, y, 0)))

            else:
                self.link_prism(self.update_coordinates( (p * (self.width + separation) + x, y, 0)))




