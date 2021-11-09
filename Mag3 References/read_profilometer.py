# necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")

###################################################################################################################################
# here, parameters to be changed by the user are found

peak_detection_level = 12           # defines the threshold used to detect peaks (compare to blue line in Fig. 2), in units of standard deviations

threshold_mean_whole_square = 2         # defines the threshold used for calcualting the mean and standard deviation, in units of standard deviations
threshold_mean_periphery = 2
threshold_mean_center = 2

boundaries_whole_square = [0.1, 0.9]    # defines the regions used to calculate the corresponding mean, in percent of a whole segment (= "single measurement")
boundaries_periphery = [0.05, 0.1, 0.9, 0.95]
boundaries_center = [0.45, 0.55]

step_size = 5                           # size of a single square, in mm
rows_to_skip = 23
###################################################################################################################################
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
dpi = 100
figsize = (10, 5)
colors = ['black', 'mediumblue', 'dodgerblue', 'slategray', 'teal', 'darkgreen', 'yellowgreen', 'saddlebrown', 'purple', "deeppink"]
###################################################################################################################################

# this is the main function, where the evaluation happens. Only change this if you know what you are doing
def main(filename, manual_list = False):

    # loads data from file
    data = pd.read_csv(filename, skiprows = rows_to_skip, header = None, usecols = [0, 1], sep = ",")

    sample = filename.split("_")[0]
    name = filename.split(".csv")[0]
    x_start = float(filename.split("_")[1].split("x")[-1])
    y_start = float(filename.split("_")[2].split("y")[-1])
    y_end = float(filename.split("_")[-1].split("y")[-1].split(".csv")[0])

    if abs(y_start) > abs(y_end):
        data[0] = np.flip(np.array([-1*x for x in data[0]]))

    y_steps_mm = np.arange(y_start, y_end + step_size, step_size)
    y_steps = [x*1e3 for x in y_steps_mm]

    print(f"\n\tAssumed y-steps ({len(y_steps)} steps):\t\t\t{y_steps}")

    df_x_list, df_y_list = [], []
    df_mean_a_list, df_std_a_list, df_mean_b_list, df_std_b_list, df_mean_c_list, df_std_c_list = [], [], [], [], [], []
    # Note: variables ending with   _a refer to mean caclulations using     the whole square,
    #                               _b                                      the periphety,
    #                               _c                                      the center

    df_out = pd.DataFrame()

    xx, yy = list(data[0]), list(data[1])
    y_diff = np.diff(yy)
    xx1 = list(np.copy(xx))
    xx1.pop()

    # automatic peak detection
    peaks, _ = find_peaks(abs(y_diff), height = np.mean(abs(y_diff)) + peak_detection_level * np.std(abs(y_diff)), distance = 250)
    for i in peaks:
        plt.plot(xx1[i], abs(y_diff[i]), 'rx')

    peaks_x = [xx1[x] for x in peaks]
    peaks_y = [abs(y_diff[x]) for x in peaks]

    width = step_size * 1e3  # the square in microns
    cut_sites, cut_sites_idx = [xx[0]], []
    for p in range(len(peaks_x) - 1):
        if peaks_x[p+1] - peaks_x[p] < width/2:
            cut_sites.append(int((peaks_x[p+1] + peaks_x[p]) / 2))
    cut_sites.append(xx[-1]) 

         
    for c in cut_sites:
        cut_sites_idx.append(xx.index(min(xx, key=lambda x:abs(x-c))))


    dict_x, dict_y = {}, {}
    if manual_list == False:
        print(f"\tIdentified y-steps ({len(cut_sites) - 1} steps):\t\t\t{[x for x in cut_sites if x != cut_sites[-1]]}\n")
        for c in range(len(cut_sites_idx) - 1):
            dict_x[c] = xx[cut_sites_idx[c] : cut_sites_idx[c+1]]
            dict_y[c] = yy[cut_sites_idx[c] : cut_sites_idx[c+1]]

    if manual_list != False:
        print(f"\tManually specified y-steps ({len(manual_list)} steps):\t\t{manual_list}\n")
        # if y_end > y_start:
        manual_list.append(xx[-1])
        cut_sites_idx = [xx.index(min(xx, key=lambda x:abs(x-c))) for c in manual_list]

        for c in range(len(cut_sites_idx) - 1):
            dict_x[c] = xx[cut_sites_idx[c] : cut_sites_idx[c+1]]
            dict_y[c] = yy[cut_sites_idx[c] : cut_sites_idx[c+1]]

    df_x = pd.DataFrame.from_dict(dict_x, orient='index').transpose()
    df_y = pd.DataFrame.from_dict(dict_y, orient='index').transpose()

        
    def fit_func_lin(xlist, a, b):
        return [(a*x + b) for x in xlist]


    # plt.ion()
    plt.figure(1, figsize = figsize)
    plt.suptitle(f"Disected Measurement, {len(cut_sites_idx) - 1} segments identified", fontsize = fontsize_title)
    ax1 = plt.subplot(211)
    ax1.plot(data[0], data[1], linewidth = linewidth, color = "black")
    for c in range(len(cut_sites_idx)):
        ax1.vlines(x = xx[cut_sites_idx[c]], color = "red", ymin = min(yy), ymax = max(yy))
    ax1.set_ylabel(fr"Height [$\AA$]", fontsize = fontsize)

    ax2 = plt.subplot(212, sharex = ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_xlabel(fr"Distance in y-Direction [$\mu$$m$]", fontsize = fontsize)
    ax2.set_ylabel(fr"Height [$\AA$]", fontsize = fontsize)
    y_start_copy = y_start
    new_c = 0
    for c in range(len(cut_sites_idx) - 1):
        X, Y = df_x[c], df_y[c]
        if c >= len(colors):
            ax2.plot(X, Y, linewidth = linewidth, color = colors[new_c], label = f"({x_start}, {y_start_copy})")
            new_c += 1
        else:
            ax2.plot(X, Y, linewidth = linewidth, color = colors[c], label = f"({x_start}, {y_start_copy})")
        y_start_copy += step_size
        

    plt.savefig(f"{name}_segments.png", dpi = dpi)


    plt.subplots_adjust(hspace=.0)
    plt.legend(loc = "lower_a center", ncol = 5, fontsize = 7)

    new_i = 0
    for i in range(len(cut_sites_idx) - 1):
        xdata, ydata = list(df_x[i].dropna()), list(df_y[i].dropna())
        xlist = [xdata[0], xdata[-1]]
        ylist = [ydata[0], ydata[-1]]
        popt, pcov = curve_fit(fit_func_lin, xlist, ylist)
        a, b = popt

        
        y_fit = fit_func_lin(xdata, a, b)
        y_corr = np.subtract(ydata, y_fit)

        cols = ["x original", "y original", "y fit", "y corrected"]
        df = pd.DataFrame()

        df["x original"] = xdata
        df["y original"] = ydata
        df["y fit"] = y_fit
        df["y corrected"] = y_corr
        df.to_csv(f"{sample}_at_{x_start}_{y_start}.csv", sep = ";", header = True, index = False)


        plt.figure(1000 + i, figsize = figsize)
        plt.suptitle(f"Coordinates: ({x_start}, {y_start})\n", fontsize = fontsize_title)
        ax1 = plt.subplot(121)
        ax1.set_title(f"Original", fontsize = fontsize_title)
        if i >= len(colors):
            ax1.plot(xdata, ydata, color = colors[new_i], label = "Original data")
        else:
            ax1.plot(xdata, ydata, color = colors[i], label = "Original data")
        aformat, bformat = float("{:.2f}".format(a)), float("{:.2f}".format(b))
        ax1.plot(xdata, y_fit, color = "red", label = f"f(x) = {aformat}x + {bformat}")
        ax1.plot(xlist, ylist, "rx")
        ax1.set_xlabel(fr"Distance in y-Direction [$\mu$$m$]", fontsize = fontsize)
        ax1.set_ylabel(fr"Height [$\AA$]", fontsize = fontsize)
        plt.legend(fontsize = fontsize_legend)

        ax2 = plt.subplot(122)
        ax2.set_title("Corrected", fontsize = fontsize_title)
        if i >= len(colors):
            ax2.plot(xdata, y_corr, color = colors[new_i], linewidth = linewidth, marker = ".", label = "Corrected data")
            new_i += 1
        else:
            ax2.plot(xdata, y_corr, color = colors[i], linewidth = linewidth, marker = ".", label = "Corrected data")
        
        ax2.set_xlabel(fr"Distance in y-Direction [$\mu$$m$]", fontsize = fontsize)

        x_filtered_a = xdata[round(boundaries_whole_square[0]*len(xdata)):round(boundaries_whole_square[1]*len(xdata))]
        y_filtered_a = y_corr[round(boundaries_whole_square[0]*len(y_corr)):round(boundaries_whole_square[1]*len(y_corr))]
        x_filtered_a, y_filtered_a = list(x_filtered_a), list(y_filtered_a)
        lower_a, upper_a = np.mean(y_filtered_a) - threshold_mean_whole_square * np.std(y_filtered_a), np.mean(y_filtered_a) + threshold_mean_whole_square * np.std(y_filtered_a)
        x_temp, y_temp = [], []
        for j, (a, b) in enumerate(zip(x_filtered_a, y_filtered_a)):
            if lower_a <= b <= upper_a:
                x_temp.append(a)
                y_temp.append(b)
        x_filtered_a, y_filtered_a = x_temp, y_temp


        x_filtered_b_1 = xdata[round(boundaries_periphery[0]*len(xdata)):round(boundaries_periphery[1]*len(xdata))]
        y_filtered_b_1 = y_corr[round(boundaries_periphery[0]*len(y_corr)):round(boundaries_periphery[1]*len(y_corr))]
        x_filtered_b_2 = xdata[round(boundaries_periphery[2]*len(xdata)):round(boundaries_periphery[3]*len(xdata))]
        y_filtered_b_2 = y_corr[round(boundaries_periphery[2]*len(y_corr)):round(boundaries_periphery[3]*len(y_corr))]
        x_filtered_b = list(x_filtered_b_1) + list(x_filtered_b_2)
        y_filtered_b = list(y_filtered_b_1) + list(y_filtered_b_2)
        lower_b, upper_b = np.mean(y_filtered_b) - threshold_mean_periphery * np.std(y_filtered_b), np.mean(y_filtered_b) + threshold_mean_periphery * np.std(y_filtered_b)
        x_temp, y_temp = [], []
        for j, (a, b) in enumerate(zip(x_filtered_b, y_filtered_b)):
            if lower_b <= b <= upper_b:
                x_temp.append(a)
                y_temp.append(b)
        x_filtered_b, y_filtered_b = x_temp, y_temp


        x_filtered_c = xdata[round(boundaries_center[0]*len(xdata)):round(boundaries_center[1]*len(xdata))]
        y_filtered_c = y_corr[round(boundaries_center[0]*len(y_corr)):round(boundaries_center[1]*len(y_corr))]
        x_filtered_c, y_filtered_c = list(x_filtered_c), list(y_filtered_c)
        lower_c, upper_c = np.mean(y_filtered_c) - threshold_mean_center * np.std(y_filtered_c), np.mean(y_filtered_c) + threshold_mean_center * np.std(y_filtered_c)
        x_temp, y_temp = [], []
        for j, (a, b) in enumerate(zip(x_filtered_c, y_filtered_c)):
            if lower_c <= b <= upper_c:
                x_temp.append(a)
                y_temp.append(b)
        x_filtered_c, y_filtered_c = x_temp, y_temp


        final_mean_a, final_std_a = float("{:.2f}".format(np.mean(y_filtered_a))), float("{:.2f}".format(np.std(y_filtered_a)))
        ax2.plot(x_filtered_a, y_filtered_a, color = "gold", linestyle = "None", marker = "x", label = fr"mean: {final_mean_a}$\AA$ $\pm$ {final_std_a}$\AA$")
        
        final_mean_b, final_std_b = float("{:.2f}".format(np.mean(y_filtered_b))), float("{:.2f}".format(np.std(y_filtered_b)))
        ax2.plot(x_filtered_b, y_filtered_b, color = "violet", linestyle = "None", marker = "x", label = fr"mean: {final_mean_b}$\AA$ $\pm$ {final_std_b}$\AA$")

        final_mean_c, final_std_c = float("{:.2f}".format(np.mean(y_filtered_c))), float("{:.2f}".format(np.std(y_filtered_c)))
        ax2.plot(x_filtered_c, y_filtered_c, color = "darkturquoise", linestyle = "None", marker = ".", label = fr"mean: {final_mean_c}$\AA$ $\pm$ {final_std_c}$\AA$")

        
        plt.legend(fontsize = fontsize_legend)
        plt.tight_layout()
        plt.savefig(f"{sample}_at_{x_start}_{y_start}.png", dpi = dpi)

        df_x_list.append(x_start)
        df_y_list.append(y_start)

        df_mean_a_list.append(final_mean_a)
        df_std_a_list.append(final_std_a)

        df_mean_b_list.append(final_mean_b)
        df_std_b_list.append(final_std_b)

        df_mean_c_list.append(final_mean_c)
        df_std_c_list.append(final_std_c)

        y_start += step_size

    cutsites = []
    for c in range(len(cut_sites_idx)):
        cutsites.append(xx[cut_sites_idx[c]])


    
    df_out["x [mm]"] = df_x_list
    df_out["y [mm]"] = df_y_list
    df_out["height mean [mm]"] = df_mean_a_list
    df_out["height std [mm]"] = df_std_a_list
    df_out["height periphery mean [mm]"] = df_mean_b_list
    df_out["height periphery std [mm]"] = df_std_b_list
    df_out["height center mean [mm]"] = df_mean_c_list
    df_out["height center std [mm]"] = df_std_c_list

    df_out.to_csv(f"{name}_height_results.csv", sep = ";", header = True, index = False)

    # plt.close("all")
    plt.show()
    return cutsites

############################################################################################################################################
############################################################################################################################################
# the following part runs the script


if __name__ == '__main__':
    files = [x for x in os.listdir() if x.endswith(".csv") and "at" not in x and "height" not in x]
    print("\n")
    print(files)

    ############################################################################################################################################
    
    filename = "Z-211020B_x+15_y0_to_ y+40.csv"         # put the name of the file you want to evaluate here
    
    ############################################################################################################################################
    
    # for filename in files:
    print(f"\n{filename}")


    
    ############################################################################################################################################
    # add if-statements here

    if filename == "Z-211020B_x+40_y0_to_ y+20.csv":
        manual_list = [0, 4973, 9967, 14944, 19927]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+45_y0_to_ y+10.csv":
        manual_list = [0, 4985]  #even if filename suggests it, there is a step missing in the data
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+5_y0_to_ y+40.csv":
        manual_list = [0.0, 4994, 9980, 14972, 19968, 24970, 29967, 34973, 39958]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+10_y-45_to_ y-5.csv":
        manual_list = [-44998.5000474575, -39859, -34873, -29882, -24881, -19884, -14883, -9874, -4889]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+15_y-45_to_ y-5.csv":
        manual_list = [-44998.5000474575, -39995, -34947, -29955, -24958, -19959, -14959, -9950, -4964]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+25_y-40_to_ y-5.csv":
        manual_list = [-39998.6667103093, -34982, -29999, -25009, -20015, -15018, -10021, -5014]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+40_y-25_to_ y-5.csv":
        manual_list = [-24998.8889349743, -19951, -14970, -9984, -4993]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+5_y-45_to_ y-5.csv":
        manual_list = [-44998.5000474575, -39933, -34943, -29952, -24954, -19954, -14954, -9945, -4959]
        cutsites = main(filename, manual_list = manual_list)
    
    elif filename == "Z-211020B_x-20_y-45_to_ y-5.csv":
        manual_list = [-44998.5000474575, -39929, -34942, -29951, -24952, -19951, -14950, -9940, -4954]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x-35_y-35_to_ y-5.csv":
        manual_list = [-34998.8333697648, -29953, -24969, -19979, -14983, -9986, -4990]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x-45_y-25_to_ y-5.csv":
        manual_list = [-24998.8889349743, -19975, -14992, -10001, -5005]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x-5_y-45_to_ y-5.csv":
        manual_list = [-44998.5000474575, -39900, -34900, -29900, -24902, -19902, -14900, -9893, -4908]
        cutsites = main(filename, manual_list = manual_list)
    
    elif filename == "Z-211020B_x+15_y0_to_ y+40.csv":
        manual_list = [0.0, 4989, 9976, 14968, 19963, 24963, 29963, 34970, 39953]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+25_y0_to_ y+35.csv":
        manual_list =  [0.0, 4983, 9966, 14956, 19949, 24944, 29937, 34900]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+30_y0_to_ y+30.csv":
        manual_list =  [0.0, 4983, 9969, 14961, 19957, 24900, 29900]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x+35_y0_to_ y+25.csv":
        manual_list =  [0.0, 4978, 9960, 14950, 19943, 24900]
        cutsites = main(filename, manual_list = manual_list)

    elif filename == "Z-211020B_x-40_y0_to_ y+25.csv":
        manual_list =  [0.0, 5008, 9992, 14983, 19900, 24900]
        cutsites = main(filename, manual_list = manual_list)
    

    else:
        cutsites = main(filename)