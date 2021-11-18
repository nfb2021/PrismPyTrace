import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR) # deals with a bug in matplotlib
from tqdm import trange
from bokeh.plotting import figure#, output_file, save, show
from scipy.signal import savgol_filter
from PyTrace import Prism, PrismDistribution, AngleConverter
# from PyTrace import params as p
plt.rc('font', size = 15)

###################################################################################################################################
# Variables used for the simulation
d = 0.8
n_air = 1
alpha = 60
width = 20
n_prisms = 10000
separation = width * 0.0
forest_length, forest_width = 1000, 250

###################################################################################################################################
# Variables used for plotting
fontsize = 15
fontsize_legend = 10
fontsize_title = 13
markersize = 18
markersize_filled = 15
linewidth = 1
capsize = 10
dpi = 100
figsize = (10, 5)

###################################################################################################################################
###################################################################################################################################


def main():
    prism = Prism(width, alpha)
    prism_distr = PrismDistribution(n_prisms,separation,width,alpha,forest_length,forest_width)

    forest = prism_distr.get_random_patchy_forest(seed = 2021, threshold = 1000, nanoscribe_patch_size = (250, 1000), blender_export = True)
    real_n_prisms = len(forest)
    
    stepsize = prism_distr.width*0.01

    hsteps, hrefr_idx, hmean, hstd, hthrown_away = prism_distr.get_avg_hori_refr_idx(forest, stepsize ,n_air,d)
    vsteps, vrefr_idx, vmean, vstd, vthrown_away = prism_distr.get_avg_vert_refr_idx(forest, stepsize ,n_air,d)
    
    hmean, hstd, vmean, vstd = float("{:.2f}".format(hmean)), float("{:.2f}".format(hstd)), float("{:.2f}".format(vmean)), float("{:.2f}".format(vstd))
    
    fig = plt.figure(1, figsize = figsize)
    ax = fig.add_subplot()
    plt.title(fr"#Prisms: {real_n_prisms}, $n_1$:{n_air},  $\delta$:{d},  $\alpha$:{alpha}$\degree$,  width:{width}$\mu$m", fontsize = fontsize_title)
    for prism in forest:
        ax.plot(*prism.exterior.xy, color = "grey") 

    ax2=ax.twiny()
    ax2.plot(hrefr_idx, hsteps, color="blue", linewidth = linewidth, label = "Av. Refr. Idx.") #, linestyle = "None", marker = "x", markersize = 1)
    ax2.axvline(x = hmean, ymin = 0, ymax = forest_length, color = "darkred", linewidth = linewidth, label = f"Mean = {hmean}")
    ax2.axvline(x = hmean + hstd, ymin = 0, ymax = forest_length, color = "red", linewidth = linewidth, label = f"Std = {hstd}")
    ax2.axvline(x = hmean - hstd, ymin = 0, ymax = forest_length, color = "red", linewidth = linewidth)

    ax.set_xlabel(r"Distance [$\mu$m]", fontsize = fontsize)
    ax.set_ylabel(r"Distance [$\mu$m]", fontsize = fontsize, color = "grey")   
    ax2.set_ylabel(r"Refractive Index", fontsize = fontsize, color = "blue")  

    ax.tick_params(axis='x', colors='grey')
    ax2.tick_params(axis='x', colors='blue')
    plt.legend(fontsize = fontsize_legend)




    fig = plt.figure(2, figsize = figsize)
    ax = fig.add_subplot()
    plt.title(fr"#Prisms: {real_n_prisms}, $n_1$:{n_air},  $\delta$:{d},  $\alpha$:{alpha}$\degree$,  width:{width}$\mu$m", fontsize = fontsize_title)
    for prism in forest:
        ax.plot(*prism.exterior.xy, color = "grey")    # set x-axis label

    ax2=ax.twinx()
    ax2.plot(vsteps, vrefr_idx, color="blue", linewidth = linewidth, label = "Av. Refr. Idx.") #, linestyle = "None", marker = "x", markersize = 1)
    ax2.axhline(y = vmean, xmin = 0, xmax = forest_width, color = "darkred", linewidth = linewidth, label = f"Mean = {vmean}")
    ax2.axhline(y = vmean + vstd, xmin = 0, xmax = forest_width, color = "red", linewidth = linewidth, label = f"Std = {vstd}")
    ax2.axhline(y = vmean - vstd, xmin = 0, xmax = forest_width, color = "red", linewidth = linewidth)

    ax.set_xlabel(r"Distance [$\mu$m]", fontsize = fontsize)
    ax.set_ylabel(r"Distance [$\mu$m]", fontsize = fontsize, color = "grey")   
    ax2.set_ylabel(r"Refractive Index", fontsize = fontsize, color = "blue")  

    ax.tick_params(axis='y', colors='grey')
    ax2.tick_params(axis='y', colors='blue')
    plt.axes().set_aspect('equal')
    plt.legend(fontsize = fontsize_legend)

   

    plt.figure(1000, figsize = figsize)#, figsize)
    plt.title(fr"#Prisms: {real_n_prisms}, $n_1$:{n_air},  $\delta$:{d},  $\alpha$:{alpha}$\degree$,  width:{width}$\mu$m", fontsize = fontsize_title)
    for prism in forest:
        plt.plot(*prism.exterior.xy, color = "grey")    # set x-axis label
    plt.xlabel(r"Distance [$\mu$m]", fontsize = fontsize)
    plt.ylabel(r"Distance [$\mu$m]", fontsize = fontsize)
    
    plt.figure(2000, figsize = figsize)#, figsize)
    plt.title(r"Verticllay: $\barn$" + fr" = {hmean} $\pm$ {hstd}", fontsize = 0.7*fontsize_title)
    plt.plot(hrefr_idx, hsteps, color="blue", linewidth = linewidth)#, label = "Av. Refr. Idx.")
    plt.axvline(x = float(hmean), ymin = 0, ymax = forest_length, color = "darkred", linewidth = linewidth, label = f"Mean = {hmean}")
    plt.axvline(x = float(hmean) + float(hstd), ymin = 0, ymax = forest_length, color = "red", linewidth = linewidth, label = f"Std = {hstd}")
    plt.axvline(x = float(hmean) - float(hstd), ymin = 0, ymax = forest_length, color = "red", linewidth = linewidth)
    plt.ylabel(r"Distance [$\mu$m]", fontsize = fontsize)
    plt.xlabel(r"Refractive Index", fontsize = fontsize) 

    plt.figure(3000, figsize = figsize)#, figsize)
    plt.title(r"Horizontally: $\barn$" + fr" = {vmean} $\pm$ {vstd}", fontsize = 0.7*fontsize_title)
    plt.plot(vsteps, vrefr_idx, color="blue", linewidth = linewidth)#, label = "Av. Refr. Idx.")
    plt.axhline(y = float(vmean), xmin = 0, xmax = forest_width, color = "darkred", linewidth = linewidth, label = f"Mean = {vmean}")
    plt.axhline(y = float(vmean) + float(vstd), xmin = 0, xmax = forest_width, color = "red", linewidth = linewidth, label = f"Std = {vstd}")
    plt.axhline(y = float(vmean) - float(vstd), xmin = 0, xmax = forest_width, color = "red", linewidth = linewidth)

    plt.xlabel(r"Distance [$\mu$m]", fontsize = fontsize)
    plt.ylabel(r"Refractive Index", fontsize = fontsize)  

    plt.legend(fontsize = fontsize_legend)
    plt.show()


if __name__ == '__main__':
    main()