import numpy as np
import math
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR) # deals with a bug in matplotlib
import pandas as pd
import os
from tqdm import trange
from shapely.geometry import LineString, MultiLineString
from scipy.signal import savgol_filter
from PyTrace import Prism, PrismDistribution, AngleConverter, RayTracing, Detector, PyhsicalProperties
plt.rc('font', size = 15)

###################################################################################################################################
# Variables used for the simulation

deltas = [3.57573435E-06, 3.41988334E-06]
betas = [5.31256283E-08, 5.29317106E-08]
energies = [23188.7852, 23211.2871]
			
deltas = [3.56042119E-06, 3.11976487E-06]
betas = [5.30091882E-08, 3.40426254E-07]
energies = [23202.2852, 23225.2969]

deltas = [3.73858029E-06, 3.73551097E-06]
betas = [3.18214262E-07, 3.17150466E-07]
energies = [23670.8516, 23693.1855]

deltas_all = [[3.57573435E-06, 3.41988334E-06], [3.56042119E-06, 3.11976487E-06], [3.73858029E-06, 3.73551097E-06]]
betas_all = [[5.31256283E-08, 5.29317106E-08], [5.30091882E-08, 3.40426254E-07], [3.18214262E-07, 3.17150466E-07]]
energies_all = [[23188.7852, 23211.2871], [23202.2852, 23225.2969], [23670.8516, 23693.1855]]

n_air = 1
alpha = 60
width = 20
separation = width * 0.0
n_prisms = 100
prism_start = (0, 0)
forest_length, forest_width = 200, 200
# b = 1
# theta_inc = 0.4        # for 100 prisms, b = 15 & theta_inc = -0.4 deg nearly the greates angle possible for this geometry
theta_inc_list = list(np.arange(-0.2, 0.2, step = 0.1))
b_list = list(np.arange(1, 15, (15-1)/len(theta_inc_list)))

# test runs with fewer values
# theta_inc_list = [-0.4, 0.4]
# b_list = [1, 15]

length = 1000
detector_distance = 3.5
pixelsize = 55e-6

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
figsize = (7, 5)

###################################################################################################################################
###################################################################################################################################

def main():
    prism = Prism(width, alpha)
    prism_distr = PrismDistribution(n_prisms, separation, width, alpha, forest_length, forest_width)
    prism_arr = prism_distr.get_h_array(prism_start)
   
    df = pd.DataFrame()
    cols = ["Theta_inc", "y0", 
            "Energy 1 [eV]", "delta 1 [ev]", "Deflection 1 [deg]", "OPL 1 [um]", "GPL 1 [um]", 
            "Energy 2 [eV]", "delta 2 [ev]", "Deflection 2 [deg]", "OPL 2 [um]", "GPL 2 [um]"]
    

    for u, (tt, bbb) in enumerate(zip(theta_inc_list, reversed(b_list))):
        theta_inc = tt
        b = bbb
        print(f"\n\t{u+1} of {len(theta_inc_list)}:")
        print(f"\tTheta_inc: {theta_inc}, b: {b}\n")
        ray_source = (-5, b)
        source_x, source_y = ray_source

        plt.figure(1 + u, figsize = figsize)
        for p, pp in enumerate(prism_arr):
            plt.plot(*pp.exterior.xy, color = "grey")
        plt.xlabel(r"Distance [$\mu$m]", fontsize = fontsize)
        plt.ylabel(r"Distance [$\mu$m]", fontsize = fontsize)
        plt.tick_params(right = True, top = True, direction = 'in')
        plt.tight_layout()

        data = [theta_inc, b]

        for i, (dd, bb) in enumerate(zip(deltas, betas)):
            plt.figure(1 + u, figsize = figsize)
            dd_formatted = "{:.3e}".format(dd)
            energy_formatted = "{:.3e}".format(energies[i])
            if i == 0:
                color = "red"
                label = r"$\delta_{1}$ = " + f"{dd_formatted}, $E_{1}$ = {energy_formatted}eV"
            if i == 1:
                color = "blue"
                label = r"$\delta_{2}$ = " + f"{dd_formatted}, $E_{2}$ = {energy_formatted}eV"

            converter = AngleConverter(n_air, dd, n_prisms, separation, width, alpha, forest_length, forest_width)
            rt = RayTracing(theta_inc)
            theta_inc_formatted = float("{:.2f}".format(rt.theta_inc))
            plt.title(r"Pulse of ~23.2keV $\pm$ 11.5eV on Rh" + f"\n" + fr"$n_1$:{converter.n_air}, $\alpha$:{prism.alpha}$\degree$,  width:{prism.width}$\mu$m,  $\Theta_{{inc}}$:{theta_inc_formatted}$\degree$,   source:{source_x, source_y}$\mu$m", fontsize = fontsize_title)

            data.append(energies[i])
            data.append(dd)

            for p in trange(len(prism_arr)):
                pp = prism_arr[p]
                pp_centroid = prism.get_centroid(pp)

                if p == 0:
                    sray = rt.get_source_ray(rt.theta_inc, b, ray_source, length)
                    sray_unrefr = sray
                    intersec_l = sray.intersection(pp) 

                    if len(intersec_l.coords) > 0: 
                        raylx, rayly = rt.get_plot_from_LineString(LineString([ray_source, (intersec_l.coords[0][0], intersec_l.coords[0][1])]))
                        plt.plot(raylx, rayly, color = color)
                        # plt.plot(intersec_l.coords[0][0], intersec_l.coords[0][1], "kx")

                        theta_inc_snell = converter.opt_ax_to_snell(rt.theta_inc)
                        theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta = converter.get_snell_angles(theta_inc_snell, theta_inc_snell)
                        theta_2_ax = -1 * converter.snell_to_opt_ax(theta_2)        # -1: convention
                        theta_4_ax = delta_theta
                        data_out_df = pd.DataFrame()

                        ray = rt.get_ray(theta_2_ax, (intersec_l.coords[0][0], intersec_l.coords[0][1]), length)
                        intersec_r = ray.intersection(pp)   
                        raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_l.coords[0][0], intersec_l.coords[0][1]), (intersec_r.coords[1][0], intersec_r.coords[1][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                        plt.plot(raypx, raypy, color = color)
                        # plt.plot(intersec_r.coords[1][0], intersec_r.coords[1][1], "kx")

                        data_out_df = rt.to_dataframe(data_out_df, pp_centroid, 
                        prism.width, prism_distr.separation, prism.alpha, ray_source,
                        theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta, 
                        intersec_l.coords[0][0], intersec_l.coords[0][1],
                        intersec_r.coords[1][0], intersec_r.coords[1][1])

                    else:
                        srayx, srayy = rt.get_plot_from_LineString(sray)
                        plt.plot(srayx, srayy, color = color)


                if 0 < p <= len(prism_arr) - 1:
                    if len(intersec_r.coords) > 1:
                        ray = rt.get_ray(theta_4_ax, (intersec_r.coords[1][0], intersec_r.coords[1][1]), length)
                        intersec_l = ray.intersection(pp)
                        if len(intersec_l.coords) > 0:
                            raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_r.coords[1][0], intersec_r.coords[1][1]), (intersec_l.coords[0][0], intersec_l.coords[0][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                            plt.plot(raypx, raypy, color = color)
                            # plt.plot(intersec_l.coords[0][0], intersec_l.coords[0][1], "kx")

                            theta_1_snell = prism.alpha - theta_4
                            theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta = converter.get_snell_angles(theta_inc_snell, theta_1_snell)
                            theta_2_ax = -1 * converter.snell_to_opt_ax(theta_2)        # -1: convention
                            theta_4_ax = delta_theta

                            ray = rt.get_ray(theta_2_ax, (intersec_l.coords[0][0], intersec_l.coords[0][1]), length)
                            intersec_r = ray.intersection(pp)  
                            if type(intersec_r) == MultiLineString: # no idea why this appears
                                intersec_r = intersec_r[1]          # this is an assumption
                            raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_l.coords[0][0], intersec_l.coords[0][1]), (intersec_r.coords[1][0], intersec_r.coords[1][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                            plt.plot(raypx, raypy, color = color)
                            # plt.plot(intersec_r.coords[1][0], intersec_r.coords[1][1], "kx")

                            data_out_df = rt.to_dataframe(data_out_df, pp_centroid, 
                            prism.width, prism_distr.separation, prism.alpha, ray_source,
                            theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta, 
                            intersec_l.coords[0][0], intersec_l.coords[0][1],
                            intersec_r.coords[1][0], intersec_r.coords[1][1])

                    else:
                        rayx, rayy = rt.get_plot_from_LineString(ray)
                        plt.plot(rayx, rayy, color = color)


                if p == len(prism_arr) - 1:
                    length_direct = 1.2 * (len(prism_arr) * (prism.width + prism_distr.separation) - prism_distr.separation)

                    unrefr_ray = rt.get_ray(rt.theta_inc,ray_source, length_direct)#, b, ray_source, length
                    unrefrx, unrefry = rt.get_plot_from_LineString(unrefr_ray)
                    length_rest = (length_direct - intersec_r.coords[1][0]) / (np.cos(math.radians(theta_4_ax)))
                    ray = rt.get_ray(theta_4_ax, (intersec_r.coords[1][0], intersec_r.coords[1][1]), length_rest)
                    raypx, raypy = rt.get_plot_from_LineString(ray) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)

                    plt.plot(raypx, raypy, color = color, label = label)
                    
                    eff_detector_distance =detector_distance - (len(prism_arr) * (prism.width + prism_distr.separation) - prism_distr.separation) * 1e-6
                    detector = Detector(theta_4_ax, eff_detector_distance,pixelsize)
                    diff_unrefr_refr_ray, result = detector.get_spot_difference()           # in meter
                    diff_unrefr_refr_ray = diff_unrefr_refr_ray * 1e6               # in microns


                    eff_detector_distance = detector_distance - (len(prism_arr) * (prism.width + prism_distr.separation) - prism_distr.separation) * 1e-6
                    detector = Detector(theta_4_ax, eff_detector_distance, pixelsize)
                    diff_unrefr_refr_ray, result = detector.get_spot_difference()           # in meter
                    diff_unrefr_refr_ray = diff_unrefr_refr_ray * 1e6               # in microns
                    # print(f"\n\tAt a detector distance of {detector_distance}m, spots the unrefracted and refracted beam are {diff_unrefr_refr_ray}um apart.")
                    # print(f"\tGiven a pixel size of {pixelsize}um, this corresponds to {result} pixels.\n")

                    if dd == deltas[-1]:
                        plt.plot(unrefrx, unrefry, color = "gold", linestyle = "dotted", label = "unrefracted beam")
                        plt.legend(fontsize =fontsize_legend, loc = "lower right")


                csv_name = f"Ray_Tracing_{prism_distr.n_prisms}_{prism.alpha}_{dd}_{prism.width}_{prism_distr.separation}_{rt.theta_inc}_{b}.csv"
                if csv_name in os.listdir():
                    os.remove(csv_name)
                data_out_df.to_csv(csv_name, sep = ";")

            phypo = PyhsicalProperties(data_out_df, n_air, dd, bb)
            distances, intensity, half_life, found_inside_array = phypo.get_intensity()
        
            geom_path_length_tot, opt_path_length_tot, geom_path_length_air, opt_path_length_air, geom_path_length_prism, opt_path_length_prism = phypo.get_gpl_and_opl(all_path_lengths = True)  

            data.append(list(data_out_df["Delta Theta abs. [deg]"])[-1])
            data.append(opt_path_length_tot)
            data.append(geom_path_length_tot)

            print(f"\n\tPrism:\tGeom. Path Length: {geom_path_length_prism}um,\t Optical Path Length: {opt_path_length_prism}um")
            print(f"\tAir:\tGeom. Path Length: {geom_path_length_air}um,\t\t Optical Path Length: {opt_path_length_air}um")
            print(f"\t___________________________________________________________________________________________________")
            print(f"\tTotal:\tGeom. Path Length: {geom_path_length_tot}um,\t\t Optical Path Length: {opt_path_length_tot}um\n")
            intensity_formatted = "{:.3e}".format((1 - (intensity[-1] / intensity[0])) * 100)
            print(f"\tThe fraction of photons absorpted by the prisms is {intensity_formatted}%\n")

        temp_df = pd.DataFrame([data], columns = cols)
        df = df.append(temp_df)

        plt.tight_layout()
        plt.savefig(f"Rh_pulse_{theta_inc}_{b}.png", dpi = dpi)
        plt.close("all")

    plt.show()

    csv_name = f"Rh_pulse_{n_prisms}_prisms_{energies[0]}_{energies[1]}.csv"
    if csv_name in os.listdir():
        os.remove(csv_name)
    df.to_csv(csv_name, sep = ";")


def comparison_plots():
    cols = ["Theta_inc", "y0", 
            "Energy 1 [eV]", "delta 1 [ev]", "Deflection 1 [deg]", "OPL 1 [um]", "GPL 1 [um]", 
            "Energy 2 [eV]", "delta 2 [ev]", "Deflection 2 [deg]", "OPL 2 [um]", "GPL 2 [um]"]

    colors = ["blue", "purple", "red"]
    
    for c in range(3):
        dd, bb, ee = deltas_all[c], betas_all[c], energies_all[c]


        data_in_df = pd.read_csv(f"Rh_pulse_{n_prisms}_prisms_{ee[0]}_{ee[1]}.csv", header = 1, sep = ";", usecols = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], names = cols)


        theta_inc = list(data_in_df["Theta_inc"])
        opl1 = list(data_in_df["OPL 1 [um]"])
        opl2 = list(data_in_df["OPL 2 [um]"])
        gpl1 = list(data_in_df["GPL 1 [um]"])
        gpl2 = list(data_in_df["GPL 2 [um]"])
        defl1 = list(data_in_df["Deflection 1 [deg]"])
        defl2 = list(data_in_df["Deflection 2 [deg]"])

        
        opl_diff = []
        for i, (a, b) in enumerate(zip(opl1, opl2)):
            opl_diff.append(abs(a - b))

        plt.figure(1003, figsize = figsize)
        plt.title(fr"E$_1$={ee[0]}eV, E$_2$={ee[1]}eV", fontsize = fontsize_title, y = 1.08)
        plt.plot(theta_inc, opl_diff, linewidth = linewidth, color = colors[c])
        plt.xlabel(r"$\Theta_{inc}$ [$\degree$]", fontsize = fontsize)
        plt.ylabel(r"$\Delta$ Optical Path Length [$\mu$$m$]", fontsize = fontsize)
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.gcf().subplots_adjust(bottom=0.15)
        # plt.legend(fontsize = fontsize_legend)
        plt.savefig(f"Rh_pulse_{n_prisms}_prisms_{ee[0]}_{ee[1]}_dOPL.png", dpi = dpi)

        gpl_diff = []
        for i, (a, b) in enumerate(zip(gpl1, gpl2)):
            gpl_diff.append(abs(a - b))

        plt.figure(1004, figsize = figsize)
        plt.title(fr"E$_1$={ee[0]}eV, E$_2$={ee[1]}eV", fontsize = fontsize_title, y = 1.08)
        plt.plot(theta_inc, gpl_diff, linewidth = linewidth, color = colors[c])
        plt.xlabel(r"$\Theta_{inc}$ [$\degree$]", fontsize = fontsize)
        plt.ylabel(r"$\Delta$ Geometric Path Length [$\mu$$m$]", fontsize = fontsize)
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.savefig(f"Rh_pulse_{n_prisms}_prisms_{ee[0]}_{ee[1]}_dGPL.png", dpi = dpi)
        # plt.legend(fontsize = fontsize_legend)

        defl_diff = []
        for i, (a, b) in enumerate(zip(defl1, defl2)):
            defl_diff.append(abs(a - b))

        plt.figure(1005, figsize = figsize)
        plt.title(fr"E$_1$={ee[0]}eV, E$_2$={ee[1]}eV", fontsize = fontsize_title, y = 1.08)
        plt.plot(theta_inc, defl_diff, linewidth = linewidth, color = colors[c])
        plt.xlabel(r"$\Theta_{inc}$ [$\degree$]", fontsize = fontsize)
        plt.ylabel(r"$\Delta$ Deflection [$\degree$]", fontsize = fontsize)
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        plt.gcf().subplots_adjust(bottom=0.15)
        # plt.legend(fontsize = fontsize_legend)
        plt.savefig(f"Rh_pulse_{n_prisms}_prisms_{ee[0]}_{ee[1]}_dDefl.png", dpi = dpi)

        plt.show()
        plt.close("all")

    
    





if __name__ == '__main__':
    main()

    comparison_plots()
    





