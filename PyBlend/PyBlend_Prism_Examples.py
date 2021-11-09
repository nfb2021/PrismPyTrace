# >>> This file contains the code for the prism forest to be implemented in Blender <<<
# >>> In order for this to function, you need the library 'pyblend_prism.py'. You can get this file from '...' <<<
# >>> Then,  copy the file to 'C:\Program Files\Blender Foundation\Blender x.x\x.x\scripts\modules\pyblend_prism.py' <<<
# >>> Once this is done, open Blender and click on 'General' to open a new work space <<<
# >>> Click on one of the two smaller windows on the right size of the blender main window (Should contain 'scene collection' and 'cube') <<<
# >>> Using the cursor you can adapt the size of each window accoridng to your liking <<<
# >>> Then, after selecting one of the two windows, press 'Shift' +  'F11' to open the editor <<<
# >>> Click on 'open' and select the 'Blender_Forest.py' file <<<
# >>> All code required is written below. Just run the script in Blender using 'Alt' + 'P' <<<


from pyblend_prism import Prism
import os

# get single prisms with different parameters
prism = Prism(20, 100, 60) 
prism.clear_scene()
prism_obj = prism.define_prism(loc = (20, 15, 0), angle = 30, base_width = 22)
prism.link_prism(prism_obj)

prism.link_prism(prism.update_alpha(120))
prism.link_prism(prism.update_coordinates((-10, 20, 0)))
prism.link_prism(prism.update_width(10))

# get 1D array of prisms
prism = Prism(20, 100, 60) 
prism.make_array(0, -30, 100, 5) # make_array(x_start, y_start, no_of_prisms, separation)

prism = Prism(40, 150, 120) 
prism.make_array(0, -50, 55, 5)        
