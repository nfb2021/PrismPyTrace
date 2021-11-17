# >>> This file contains the code for the prism forest to be implemented in Blender <<<
# >>> In order for this to function, you need the library 'pyblend_prism.py'. You can get this file from 'https://github.com/nfb2021/PrismPyTrace/tree/main/PyBlend' <<<
# >>> Then,  copy the file to 'C:\Program Files\Blender Foundation\Blender x.x\x.x\scripts\modules\pyblend_prism.py' <<<
# >>> Once this is done, open Blender and click on 'General' to open a new work space <<<
# >>> Click on one of the two smaller windows on the right size of the blender main window (Should contain 'scene collection' and 'cube') <<<
# >>> Using the cursor you can adapt the size of each window accoridng to your liking <<<
# >>> Then, after selecting one of the two windows, press 'Shift' +  'F11' to open the editor <<<
# >>> Click on 'open' and select the 'Blender_Forest.py' file <<<
# >>> All code required is written below. Just run the script in Blender using 'Alt' + 'P' <<<


from pyblend_prism import Prism
import os
prism = Prism(20, 100, 60) # initialize class with some standard values
prism.clear_scene()
with open(r'C:\Users\badernic\Documents\Ray_Tracing_Pulse_Compressor\PrismPyTrace\1636465631.6468964_Forest.txt', 'r') as f:
	lines = f.readlines()
	for l, line in enumerate(lines):
		if l == 0:
			continue
		x = float(line.split(';')[0])
		y = float(line.split(';')[1])
		z = float(line.split(';')[2])
		width = float(line.split(';')[3])
		alpha = float(line.split(';')[4])
		new_loc = (x, y, 0)
		prism_obj = prism.define_prism(loc = new_loc, angle = alpha, base_width = width)

		prism.link_prism(prism_obj)

