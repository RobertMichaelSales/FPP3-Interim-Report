""" Created on Fri Oct 15 2021, Author: Robert Sales, File: tapping_coordinates.py """

import math, cmath
import numpy as np
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt

#==============================================================================

plt.rcParams.update({"text.usetex": False,
                     "font.family": "Times New Roman",
                     "font.size"  : 24,
                     "font.weight": "normal",
                     "figure.dpi" : 500
                     })

figsize = (12,2.5)
linewidth = '1.8'
color1 = '#D6083B'
color2 = '#0072CF'

#==============================================================================

#==============================================================================

def AnglesInEllipse(num, a, b):
    
    assert(num > 0)
    assert(a < b)
    angles = 2 * np.pi * np.arange(num) / num
    
    if a != b:
        e = (1.0 - a ** 2.0 / b ** 2.0) ** 0.5
        tot_size = sp.special.ellipeinc(2.0 * np.pi, e)
        arc_size = tot_size / num
        arcs = np.arange(num) * arc_size
        res = sp.optimize.root(
            lambda x: (sp.special.ellipeinc(x, e) - arcs), angles)
        angles = res.x 
        
    return angles

#==============================================================================

#==============================================================================

def Ellipse():
    
    r_major = 62.5 
    r_minor = 12.5    
    n = 60
    
    angles = AnglesInEllipse(n, r_minor, r_major)
    angles = [x for x in angles if np.rad2deg(x)<=90]
    angles.sort(reverse=True)
    
    coordinates = [(62.5 - r_major * np.sin(angle), r_minor * np.cos(angle)) for angle in angles]

    return coordinates

#==============================================================================

#==============================================================================

def FlatPlate():
    
    x_origin = 62.5
    y_origin = 12.5
    length = 150
    n = 32
    
    x_coordinates = list(np.linspace(x_origin, x_origin+length, n))
    x_coordinates.pop(0)
    
    coordinates = [(x,y_origin) for x in x_coordinates]
    
    return coordinates

#==============================================================================

#==============================================================================

def Circle():
    
    x_origin = 212.5
    radius = 12.5
    n = 7
    
    angles = list(np.linspace(math.pi/2, 0, n))
    angles.pop(0)

    coordinates = [(x_origin + radius * np.cos(angle), radius * np.sin(angle)) for angle in angles]
    
    return coordinates

#==============================================================================

#==============================================================================

def OutputTxt(coordinates):
    
    x_coordinates = (np.array([x[0] for x in coordinates], ndmin=2)).T
    y_coordinates = (np.array([x[1] for x in coordinates], ndmin=2)).T
    
    data = np.concatenate((x_coordinates,y_coordinates), axis=1)
    np.savetxt("tapping_coordinates.txt", data, delimiter='\t', fmt='%10.5f', newline='\n', encoding=None)
    
    return None

#==============================================================================

#==============================================================================

def PlotShape(coordinates, save_fig = False):
    
    x_coordinates_upper = [+x[0] for x in coordinates]
    y_coordinates_upper = [+x[1] for x in coordinates]
    x_coordinates_lower = [+x[0] for x in coordinates]
    y_coordinates_lower = [-x[1] for x in coordinates]
    
    x_section_1 = [+62.5, +62.5]
    y_section_1 = [+12.5, -12.5]
    
    x_section_2 = [+212.5, +212.5]
    y_section_2 = [+12.5, -12.5]
    
    plt.figure(figsize = figsize)  
    plt.plot(x_coordinates_upper, y_coordinates_upper, color = color2, marker = 'o', linestyle = 'solid')
    plt.plot(x_coordinates_lower, y_coordinates_lower, color = color2, marker = 'o', linestyle = 'solid')
    plt.plot(x_section_1, y_section_1, color = color2, marker = 'o', linestyle = 'solid')
    plt.plot(x_section_2, y_section_2, color = color2, marker = 'o', linestyle = 'solid')

    plt.fill_between(x_coordinates_upper,
                     y_coordinates_upper,
                     y_coordinates_lower,
                     alpha = 0.5,
                     facecolor = "gray",
                     antialiased = True
        )          

    plt.text(-5, -22, "LE", color = color2, weight = "bold")
    plt.text(220, -22, "TE", color = color2, weight = "bold")
    plt.title("Plate Geometry with Pressure Tapping Locations (in mm)")
    
    plt.grid()
    plt.axis('equal')
    
    if save_fig: plt.savefig('exp_1_plate_geometry.svg')
    plt.show()
    
    return None

#==============================================================================

#==============================================================================

coordinates_ellipse = Ellipse()
coordinates_flat_plate = FlatPlate()
coordinates_circle = Circle()

coordinates = coordinates_ellipse + coordinates_flat_plate + coordinates_circle

#OutputTxt(coordinates)
PlotShape(coordinates, save_fig = True)

