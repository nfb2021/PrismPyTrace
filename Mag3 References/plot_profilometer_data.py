# necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")

###################################################################################################################################
# Variables used for plotting

plt.rc('font', size = 15)
fontsize = 15
fontsize_legend = 11
fontsize_title = 15
markersize = 18
markersize_filled = 15
linewidth = 1
capsize = 10
dpi = 200
figsize = (10, 5)
colors = ['black', 'mediumblue', 'dodgerblue', 'slategray', 'teal', 'darkgreen', 'yellowgreen', 'saddlebrown', 'purple', "deeppink"]
###################################################################################################################################

original_files = [x for x in os.listdir() if x.endswith(".csv") and "at" not in x and "height" not in x]
resulting_files = [x for x in os.listdir() if x.endswith(".csv") and "height_results" in x]
sample = list(dict.fromkeys([x.split("_")[0] for x in original_files]))

if len(original_files) == len(resulting_files) and len(sample) == 1:
    print("\nAll provided files of the profilometer were evaluated.\n")
    sample = sample[0]
else:
    if len(sample) > 1:
        print("There seem to be files of more than one sample, being: ")
        for i in sample:
            print(f"\t{i}")

        print("Please sort this. Aborting program.")
        print("\n")
        exit()
    else:
        of = [x.split(".csv")[0] for x in original_files]
        rf = [x.split("_height_results.csv")[0] for x in resulting_files]
        missing = list(set(of).symmetric_difference(set(rf)))
        mf = [f"{x}.csv" for x in missing]
        print(f"\nIt seems as if following files were not evaluated:")
        for i in mf:
            print(f"\t{i}")
        print("Aborting program.")
        print("\n")
        exit()

image_heights = np.zeros(400).reshape((20,20))
image_std = np.zeros(400).reshape((20,20))
image_heights_b = np.zeros(400).reshape((20,20))
image_std_b = np.zeros(400).reshape((20,20))
image_heights_c = np.zeros(400).reshape((20,20))
image_std_c = np.zeros(400).reshape((20,20))


for item in resulting_files:
    data = pd.read_csv(item, skiprows = 0, header = 0, sep = ";")

    xx = list(data["x [mm]"])
    yy = list(data["y [mm]"])
    height = list(data["height mean [mm]"])
    std = list(data["height std [mm]"])
    height_b = list(data["height periphery mean [mm]"])
    std_b = list(data["height periphery std [mm]"])
    height_c = list(data["height center mean [mm]"])
    std_c = list(data["height center std [mm]"])

    for (X, Y, h, s, hb, sb, hc, sc) in zip(xx, yy, height, height_b, height_c, std, std_b, std_c):
        X_shifted = int((X + 50) / 5)
        Y_shifted = int((Y + 50) / 5)
        image_heights[Y_shifted][X_shifted] = h
        image_std[Y_shifted][X_shifted] = s
        image_heights_b[Y_shifted][X_shifted] = hb
        image_std_b[Y_shifted][X_shifted] = sb
        image_heights_c[Y_shifted][X_shifted] = hc
        image_std_c[Y_shifted][X_shifted] = sc


for i, (hh, ss) in enumerate(zip([image_heights, image_heights_b, image_heights_c], [image_std, image_std_b, image_std_c])):
    im_data_in, std_data_in = hh, ss

    plt.figure(1, figsize = (6, 6))
    masked_array = np.ma.array(im_data_in, mask=(im_data_in == 0))
    cmap = cm.CMRmap
    cmap.set_bad('gainsboro',1.)
    im = plt.imshow(masked_array, interpolation='nearest', cmap=cmap, origin = "lower")
    plt.colorbar(im, orientation='vertical', label = r"Height [$\AA$]")  
    plt.xticks([0, 5, 10, 15], [-50, -25, 0, 25])
    plt.yticks([0, 5, 10, 15], [-50, -25, 0, 25])
    plt.xlabel("X Axis [mm]", fontsize = fontsize)
    plt.ylabel("Y Axis [mm]", fontsize = fontsize)
    plt.minorticks_on()
    plt.gcf().subplots_adjust(bottom=0.15)
    if i == 0:
        plt.title(f"{sample}:\nHeight Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_height.png", dpi = dpi)
    if i == 1:
        plt.title(f"{sample}:\nHeight Periphery Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_height_periphery.png", dpi = dpi)
    if i == 2:
        plt.title(f"{sample}:\nHeight Center Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_height_center.png", dpi = dpi)



    plt.figure(2, figsize = (6, 6))
    masked_array = np.ma.array(std_data_in, mask=(std_data_in == 0))
    cmap = cm.CMRmap
    cmap.set_bad('gainsboro',1.)
    im = plt.imshow(masked_array, interpolation='nearest', cmap=cmap, origin = "lower")
    plt.colorbar(im, orientation='vertical', label = r"Standard Deviation [$\AA$]")  
    plt.xticks([0, 5, 10, 15], [-50, -25, 0, 25])
    plt.yticks([0, 5, 10, 15], [-50, -25, 0, 25])
    plt.xlabel("X Axis [mm]", fontsize = fontsize)
    plt.ylabel("Y Axis [mm]", fontsize = fontsize)
    # plt.tick_params(right = True, top = True, direction = 'out')
    plt.minorticks_on()
    plt.gcf().subplots_adjust(bottom=0.15)
    if i == 0:
        plt.title(f"{sample}:\nStandard Deviation Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_standard_deviation.png", dpi = dpi)
    if i == 1:
        plt.title(f"{sample}:\nStandard Deviation Periphery Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_standard_deviation_periphery.png", dpi = dpi)
    if i == 2:
        plt.title(f"{sample}:\nStandard Deviation Center Map", fontsize = fontsize_title)
        plt.savefig(f"{sample}_standard_deviation_center.png", dpi = dpi)
    

    plt.close("all")

# plt.show()