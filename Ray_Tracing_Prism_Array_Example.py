import numpy as np
import math
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR) # deals with a bug in matplotlib
import pandas as pd
import os
from tqdm import trange
from shapely.geometry import LineString, MultiLineString
from PyTrace import Prism, Snell, PrismDistribution, AngleConverter, RayTracing, Detector, MirageDeflector
plt.rc('font', size = 15)

###################################################################################################################################
# Variables used for the simulation
d = 1.33e-6
n_air = 1
alpha = 60
width = 20
separation = width * 0.0
n_prisms = 100
prism_start = (0, 0)
forest_length, forest_width = 125, 125
b = 7
ray_source = (-5, b)
theta_inc = 0
length = 1000
detector_distance = 3.5
pixelsize = 55e-6

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
    source_x, source_y =ray_source
    plt.figure(1, figsize =figsize)
    plt.title(fr"$n_1$:{converter.n_air},  $\delta$:{converter.d},  $\alpha$:{prism.alpha}$\degree$,  width:{prism.width}$\mu$m,  $\Theta_{{inc}}$:{rt.theta_inc}$\degree$,   source:{source_x, source_y}$\mu$m", fontsize = fontsize_title)
    plt.xlabel(r"Distance [$\mu$m]", fontsize =fontsize)
    plt.ylabel(r"Distance [$\mu$m]", fontsize =fontsize)
    plt.tick_params(right = True, top = True, direction = 'in')
    plt.tight_layout()
    for p, pp in enumerate(prism_arr):
        plt.plot(*pp.exterior.xy, color = "grey")

    for p in trange(len(prism_arr)):
        pp = prism_arr[p]
        pp_centroid = prism.get_centroid(pp)

        if p == 0:
            sray = rt.get_source_ray(rt.theta_inc,b,ray_source,length)
            sray_unrefr = sray
            intersec_l = sray.intersection(pp) 

            if len(intersec_l.coords) > 0: 
                raylx, rayly = rt.get_plot_from_LineString(LineString([ray_source, (intersec_l.coords[0][0], intersec_l.coords[0][1])]))
                plt.plot(raylx, rayly, color = "red")
                # plt.plot(intersec_l.coords[0][0], intersec_l.coords[0][1], "kx")

                theta_inc_snell = converter.opt_ax_to_snell(rt.theta_inc)
                theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta = converter.get_snell_angles(theta_inc_snell, theta_inc_snell)
                theta_2_ax = -1 * converter.snell_to_opt_ax(theta_2)        # -1: convention
                theta_4_ax = delta_theta
                data_out_df = pd.DataFrame()

                ray = rt.get_ray(theta_2_ax, (intersec_l.coords[0][0], intersec_l.coords[0][1]),length)
                intersec_r = ray.intersection(pp)   
                raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_l.coords[0][0], intersec_l.coords[0][1]), (intersec_r.coords[1][0], intersec_r.coords[1][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                plt.plot(raypx, raypy, color = "red")
                # plt.plot(intersec_r.coords[1][0], intersec_r.coords[1][1], "kx")

                data_out_df = rt.to_dataframe(data_out_df, pp_centroid, 
                prism.width, prism_distr.separation, prism.alpha,ray_source,
                theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta, 
                intersec_l.coords[0][0], intersec_l.coords[0][1],
                intersec_r.coords[1][0], intersec_r.coords[1][1])

            else:
                srayx, srayy = rt.get_plot_from_LineString(sray)
                plt.plot(srayx, srayy, color = "red")


        if 0 < p <= len(prism_arr) - 1:
            if len(intersec_r.coords) > 1:
                ray = rt.get_ray(theta_4_ax, (intersec_r.coords[1][0], intersec_r.coords[1][1]),length)
                intersec_l = ray.intersection(pp)
                if len(intersec_l.coords) > 0:
                    raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_r.coords[1][0], intersec_r.coords[1][1]), (intersec_l.coords[0][0], intersec_l.coords[0][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                    plt.plot(raypx, raypy, color = "red")
                    # plt.plot(intersec_l.coords[0][0], intersec_l.coords[0][1], "kx")

                    theta_1_snell = prism.alpha - theta_4
                    theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta = converter.get_snell_angles(theta_inc_snell, theta_1_snell)
                    theta_2_ax = -1 * converter.snell_to_opt_ax(theta_2)        # -1: convention
                    theta_4_ax = delta_theta

                    ray = rt.get_ray(theta_2_ax, (intersec_l.coords[0][0], intersec_l.coords[0][1]),length)
                    intersec_r = ray.intersection(pp)  
                    if type(intersec_r) == MultiLineString: # no idea why this appears
                        intersec_r = intersec_r[1]          # this is an assumption
                    raypx, raypy = rt.get_plot_from_LineString(LineString([(intersec_l.coords[0][0], intersec_l.coords[0][1]), (intersec_r.coords[1][0], intersec_r.coords[1][1])])) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)
                    plt.plot(raypx, raypy, color = "red")
                    # plt.plot(intersec_r.coords[1][0], intersec_r.coords[1][1], "kx")

                    data_out_df = rt.to_dataframe(data_out_df, pp_centroid, 
                    prism.width, prism_distr.separation, prism.alpha,ray_source,
                    theta_1, theta_2, theta_3, theta_4, theta_rel, delta_theta, 
                    intersec_l.coords[0][0], intersec_l.coords[0][1],
                    intersec_r.coords[1][0], intersec_r.coords[1][1])

            else:
                rayx, rayy = rt.get_plot_from_LineString(ray)
                plt.plot(rayx, rayy, color = "red")


        if p == len(prism_arr) - 1:
            length_direct = 1.2 * (len(prism_arr) * (prism.width + prism_distr.separation) - prism_distr.separation)

            unrefr_ray = rt.get_ray(rt.theta_inc,ray_source, length_direct)#, b, ray_source, length
            unrefrx, unrefry = rt.get_plot_from_LineString(unrefr_ray)
            length_rest = (length_direct - intersec_r.coords[1][0]) / (np.cos(math.radians(theta_4_ax)))
            ray = rt.get_ray(theta_4_ax, (intersec_r.coords[1][0], intersec_r.coords[1][1]), length_rest)
            raypx, raypy = rt.get_plot_from_LineString(ray) #intersec_r second cord pair bc fair pair is left intersection (== intersec_l)

            plt.plot(unrefrx, unrefry, color = "red", linestyle = "dotted", label = "unrefracted beam")
            plt.plot(raypx, raypy, color = "red", label = "refracted beam")
            plt.legend(fontsize =fontsize_legend, loc = "lower right")

            eff_detector_distance =detector_distance - (len(prism_arr) * (prism.width + prism_distr.separation) - prism_distr.separation) * 1e-6
            detector = Detector(theta_4_ax, eff_detector_distance,pixelsize)
            diff_unrefr_refr_ray, result = detector.get_spot_difference()           # in meter
            diff_unrefr_refr_ray = diff_unrefr_refr_ray * 1e6               # in microns




        csv_name = f"Ray_Tracing_{prism_distr.n_prisms}_{prism.alpha}_{prism.width}_{prism_distr.separation}_{rt.theta_inc}_{b}.csv"
        if csv_name in os.listdir():
            os.remove(csv_name)
        data_out_df.to_csv(csv_name, sep = ";", index = False)

    print(f"\n\tAt a detector distance of {detector_distance}m, spots the unrefracted and refracted beam are {diff_unrefr_refr_ray}um apart.")
    print(f"\tGiven a pixel size of {pixelsize}um, this corresponds to {result} pixels.\n")

    mirage_delta = mirdef.get_capital_delta()
    rt_delta = list(data_out_df["Delta Theta abs. [deg]"])[-1]
    print(f"\n\tAccording to the Ray Tracing Algorithm, the total deflection is {rt_delta}.")

    # Note: A Comparison to the Mirage Deflector only makes sense, if the separation between prisms is set to 0
    print(f"\tAccording to the Mirage Deflector Equation, the total deflection is {mirage_delta}.\n")


    # plt.axes().set_aspect('equal')
    plt.show()




if __name__ == '__main__':
    prism = Prism(width, alpha)
    prism_distr = PrismDistribution(n_prisms,separation,width,alpha,forest_length,forest_width)
    prism_arr = prism_distr.get_h_array(prism_start)

    converter = AngleConverter(n_air,d,n_prisms,separation,width,alpha,forest_length,forest_width)
    rt = RayTracing(theta_inc)
    mirdef = MirageDeflector(converter.opt_ax_to_snell(theta_inc),d,b,alpha,width,separation,n_prisms)


    main()

