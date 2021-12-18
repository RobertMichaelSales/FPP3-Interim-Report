""" Created on Fri Oct 15 2021, Author: Robert Sales, File: plot_surface_pressure_distribution.py """

import os, csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.close('all')

#==============================================================================

plt.rcParams.update({"text.usetex": False,
                     "font.family": "Times New Roman",
                     "font.size"  : 24,
                     "font.weight": "normal",
                     "figure.dpi" : 500
                     })

matplotlib.rcParams['mathtext.fontset'] = 'stix'    

figsize = (12,5)
linewidth = '1.8'
markersize = 8
alpha = 1.0
color1 = '#D6083B'
color2 = '#0072CF'

#==============================================================================


#==============================================================================

def LoadCSV(filename):
    
    data = []
    filepath = os.path.join(os.getcwd().replace("Plots", "Data"), filename)
    
    with open(filepath, newline='') as file:
        reader = csv.reader(file, delimiter = '\t', quoting = csv.QUOTE_NONNUMERIC)
    
        for row in reader:
            data.append(row)
          
    return np.array(data, dtype = "float32")

#==============================================================================

#==============================================================================

def PressureDistribution(save_fig = False, reynold_mode = "h"):
    
    # Load in data
    filename_lookup = {"h":"high_reynolds.txt", "l":"low_reynolds.txt"}
    data = LoadCSV(filename_lookup[reynold_mode])
    
    # Load the coordinates
    coordinates = LoadCSV("tapping_coordinates.txt")
    
    # Separate readings of pressure and temperature
    p_data = data[:, 3:19]
    t_data = data[:,20:21]
    
    # Separate readings of stagnation/static pressures in the tunnel and inlet
    p_stag = p_data[:,0:1]
    p_stat = p_data[:,1:2]
    p_stat_inlet = [p_data[8,x] for x in range(2,8) if x !=4]

    # Calculate average wind tunnel conditions
    p_atm = float(data[-1,-1])
    p_stag_avg = float(sum(p_stag)/len(p_stag))
    p_stat_avg = float(sum(p_stat)/len(p_stat))
    t_stat_avg = float(sum(t_data)/len(t_data))
    
    # Calculate the fluid density (using thermofluid data)
    rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
        
    # Calculate the wind tunnel free-stream velocity
    velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5

    # Calculate the Reynold's number of the experiment
    Re = (rho * velocity * 0.225)/(17.9*10**-6)
    
    print(Re)
    
    # Initialise a pressure coefficient list
    p_coeff = []
    
    # Initialise the error/uncertainty lists
    y_p_error = []
    y_m_error = []
    n_readings = 200
    t_multiplier = 1.972 # for 200 readings
    
    p_bias_1 = 10 # static pressure measured
    p_bias_2 = 10 # total pressure inlet
    p_bias_3 = 10 # statoc pressure inlet
    
    s_x_1 = t_multiplier * ((3.7659)/(n_readings**0.5)) # static pressure measured
    s_x_2 = t_multiplier * ((2.0635)/(n_readings**0.5)) # total pressure inlet
    s_x_3 = t_multiplier * ((1.6572)/(n_readings**0.5)) # static pressure inlet

    # Loop through sucessive readings
    for row in range(p_data.shape[0]-1):
        
        # Calculate stagnation and static pressures
        p_stag_ref = p_data[row,0]
        p_stat_ref = p_data[row,1]
        
        # Loop through the DSA channels 
        for col in [x for x in range(2,16) if x != 4]:

            # Calculate and append the pressure coefficient
            p_stat = p_data[row,col]
            p_stat_coeff = (p_stat - p_stat_ref) / (p_stag_ref - p_stat_ref)
            p_coeff.append(p_stat_coeff)
            
            # Calculate and append the measurement errors
            b_cp_1 = ((1/(p_stag_ref - p_stat_ref))**2) * (p_bias_1**2)
            b_cp_2 = (((p_stat - p_stat_ref)/((p_stag_ref - p_stat_ref)**2))**2) * (p_bias_2**2)
            b_cp_3 = (((p_stat - p_stag_ref)/((p_stag_ref - p_stat_ref)**2))**2) * (p_bias_3**2)
            
            # Calculate and append the precision errors
            s_cp_1 = ((1/(p_stag_ref - p_stat_ref))**2) * (s_x_1**2)
            s_cp_2 = (((p_stat - p_stat_ref)/((p_stag_ref - p_stat_ref)**2))**2) * (s_x_2**2)
            s_cp_3 = (((p_stat - p_stag_ref)/((p_stag_ref - p_stat_ref)**2))**2) * (s_x_3**2)
            
            b_cp = (b_cp_1 + b_cp_2 + b_cp_3)**0.5
            s_cp = (s_cp_1 + s_cp_2 + s_cp_3)**0.5
            u_cp = (b_cp**2 + s_cp**2)**0.5
            
            y_p_error.append(p_stat_coeff + u_cp)
            y_m_error.append(p_stat_coeff - u_cp)
            
    # Split the pressure coefficient surfaces
    p_coeff_upper_surf = np.array(p_coeff[0:53])
    p_coeff_lower_surf = np.array(p_coeff[0:1] + p_coeff[53:] + p_coeff[52:53])
             
    # Split the errors
    y_p_error_upper_surf = np.array(y_p_error[0:53])
    y_p_error_lower_surf = np.array(y_p_error[0:1] + y_p_error[53:] + y_p_error[52:53])
    
    y_m_error_upper_surf = np.array(y_m_error[0:53])
    y_m_error_lower_surf = np.array(y_m_error[0:1] + y_m_error[53:] + y_m_error[52:53])
    
    # Initiate plotting
    fig, ax = plt.subplots(1, 1, figsize = figsize) 
    
    # Plot the plate section cutoffs
    x_section_1 = np.array([62.5/225, 62.5/225])
    y_section_1 = np.array([-0.80, +0.05])
    
    x_section_2 = np.array([212.5/225, 212.5/225])
    y_section_2 = np.array([-0.80, +0.05])

    ax.plot(x_section_1,
             y_section_1,
             color = "gray", 
             marker = 'D', 
             linestyle = 'dashed',
             markersize = markersize, 
             linewidth = linewidth,
             label = "_nolegend_"
             )
    
    ax.plot(x_section_2, 
             y_section_2, 
             color = "gray", 
             marker = 'D', 
             linestyle = 'dashed',
             markersize = markersize, 
             linewidth = linewidth,
             label = "_nolegend_"
             )   
    
    # Plot the selected data
    ax.plot(coordinates[:,0]/225,
             p_coeff_lower_surf,
             linewidth = linewidth,
             linestyle = 'solid',
             marker = 's',
             color = color2,
             markersize = markersize, 
             alpha = alpha)
    
    ax.plot(coordinates[:,0]/225,
             p_coeff_upper_surf,
             linewidth = linewidth,
             linestyle = 'solid',
             marker = 'o',
             color = color1,
             markersize = markersize, 
             alpha = alpha)
    
    ax.fill_between(coordinates[:,0]/225,
                     y_p_error_lower_surf,
                     y_m_error_lower_surf,
                     alpha = 0.2,
                     facecolor = color2,
                     #antialiased = True
        ) 
    
    ax.fill_between(coordinates[:,0]/225,
                     y_p_error_upper_surf,
                     y_m_error_upper_surf,
                     alpha = 0.2,
                     facecolor = color1,
                     #antialiased = True
        )          

    ax.text(0.5, -0.625, "FLAT REGION", color = "gray", weight = "bold")     
    
    # Define plot aesthetics
    plt.title("Static Pressure Distribution                   (at Re = {:.2e})".format(int(round(Re, -2))))
    plt.legend(["Upper Surface", "Lower Surface"], ncol = 2, loc = 'lower center', shadow = True)
    plt.xlabel(r"Nondimensionalised Axial Chord Position ($x /C_{x}$)")
    plt.ylabel(r"Pressure Coefficient, $C_{p}$")
    ax.set_xlim(-0.02,1.02)
    ax.set_ylim(0.38, -0.85)
    plt.grid()
    
    savename_lookup = {"h":"exp_1_surface_pressure_dist_high_re.svg",
                       "l":"exp_1_surface_pressure_dist_low_re.svg"} 
    
    if save_fig: plt.savefig(savename_lookup[reynold_mode])
    
    plt.show() 
    
#==============================================================================

PressureDistribution(save_fig = False, reynold_mode = "h")