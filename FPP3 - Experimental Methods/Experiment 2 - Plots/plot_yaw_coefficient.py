""" Created on Sun Oct 31 2021, Author: Robert Sales, File: plot_yaw_coefficient.py """

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

def YawCoefficient(p_l, p_r, p_c):
    
    numer = p_l - p_r
    denom = p_c - 0.5*(p_l + p_r)
    
    yaw_coefficient =  numer/denom
                          
    return yaw_coefficient
#==============================================================================

#==============================================================================

def PlotYawCoefficient(save_fig = False):
    
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       "PRB-L":3,
                       "PRB-C":4,
                       "PRB-R":6,
                       "SPP":  7
                       }
        
        marker = ['s', 'v', 'D', 'o']
        
        # Initiate plotting
        fig, ax = plt.subplots(1, 1, figsize = figsize)   
        legend_list = []
        
        for i, setting in enumerate(["low", "med", "high"]):
            
            # Configure the filename
            filename = "cal_fine_{}_speed.txt".format(setting)
            
            # Load in data from filename
            data = LoadCSV(filename)
            
            # Separate readings of pressure and temperature versus time (col 0)
            limit = 3
            angl_readings = data[limit:-limit, 2:3]    # Can be from 1:-1 but not ideal
            pres_readings = data[limit:-limit, 8:15]
            temp_readings = data[limit:-limit, 25:26]
            
            # Extract average conditions
            p_stag_avg = float(sum(pres_readings[:,0])/len(pres_readings[:,0]))
            p_stat_avg = float(sum(pres_readings[:,1])/len(pres_readings[:,1]))
            t_stat_avg = float(sum(temp_readings)/len(temp_readings))
            p_atm = float(data[0,0])

            # Calculate the fluid density (using thermofluid data)
            rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
                
            # Calculate the wind tunnel free-stream velocity
            velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5
            print(velocity)
        
            # Calculate the Reynold's number of the experiment
            Re = (rho * velocity * 0.225)/(17.9*10**-6)    
            legend_list.append("Re = {:.2e}".format(int(round(Re, -2))))
            
            # Create an array to store new coefficient values
            yaw_coefficient = np.zeros(angl_readings.shape)
            
            # Initialise the error/uncertainty lists
            y_p_error = np.zeros(angl_readings.shape)
            y_m_error = np.zeros(angl_readings.shape)
            n_readings = 200
            t_multiplier = 1.972 # for 200 readings
            
            p_bias_1 = 10 # p_l
            p_bias_2 = 10 # p_c
            p_bias_3 = 10 # p_r
            
            s_x_1 = t_multiplier * ((3.7659)/(n_readings**0.5)) # p_l
            s_x_2 = t_multiplier * ((3.7659)/(n_readings**0.5)) # p_c
            s_x_3 = t_multiplier * ((3.7659)/(n_readings**0.5)) # p_r

            # Loop through all rows and calculate coefficients
            for row in range(0,len(angl_readings)):
                
                p_l = pres_readings[row, (lookup_dict["PRB-L"]-1)]
                p_c = pres_readings[row, (lookup_dict["PRB-C"]-1)]
                p_r = pres_readings[row, (lookup_dict["PRB-R"]-1)]
                  
                yaw_coefficient[row] = YawCoefficient(p_l, p_r, p_c)
                
                # Calculate and append the measurement errors
                b_cg_1 = ((p_c - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_1**2)
                b_cg_2 = ((p_l - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_2**2)
                b_cg_3 = ((p_l - p_c)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_3**2)
                
                # Calculate and append the precision errors
                s_cg_1 = ((p_c - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_1**2)
                s_cg_2 = ((p_l - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_2**2)
                s_cg_3 = ((p_l - p_c)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_3**2)
                
                b_cg = (b_cg_1 + b_cg_2 + b_cg_3)**0.5
                s_cg = (s_cg_1 + s_cg_2 + s_cg_3)**0.5
                u_cg = (b_cg**2 + s_cg**2)**0.5
            
                y_p_error[row] = (yaw_coefficient[row] + u_cg)
                y_m_error[row] = (yaw_coefficient[row] - u_cg)
            
            alpha = [0.2, 0.4, 0.6]
            
            ax.fill_between(angl_readings[:,0],
                            y_p_error[:,0],
                            y_m_error[:,0],
                            alpha = alpha[i]
                            ) 
            
            # Plot the selected data
            ax.plot(angl_readings,
                    yaw_coefficient,
                    linewidth = linewidth,
                    linestyle = "solid",
                    markersize = markersize
                    )
   
        # Define plot aesthetics
        plt.title(r"Yaw Coefficient vs. Relative Flow Angle")
        plt.legend(legend_list, shadow = True)
        plt.xlabel(r"Relative Flow Angle ($^\circ$)")
        plt.ylabel(r"$C_\alpha$")
        ax.set_xlim(-32, 32)
        ax.set_ylim(-5,5)
        plt.grid()
        
        if save_fig: plt.savefig('exp_2_yaw_coefficient.svg')
        plt.show() 
        
#==============================================================================

#==============================================================================

def SaveYawCoefficient():
    
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       "PRB-L":3,
                       "PRB-C":4,
                       "PRB-R":6,
                       "SPP":  7
                       }
        
        for setting in ["high", "med", "low"]:
            
            # Configure the filename
            filename = "cal_fine_{}_speed.txt".format(setting)
            
            # Load in data from filename
            data = LoadCSV(filename)
            
            # Separate readings of pressure and temperature versus time (col 0)
            angl_readings = data[1:-1, 2:3]    # Can be from 1:-1 but not ideal
            pres_readings = data[1:-1, 8:15]
            temp_readings = data[1:-1, 25:26]
            
            # Extract average conditions
            p_stag_avg = float(sum(pres_readings[:,0])/len(pres_readings[:,0]))
            p_stat_avg = float(sum(pres_readings[:,1])/len(pres_readings[:,1]))
            t_stat_avg = float(sum(temp_readings)/len(temp_readings))
            p_atm = float(data[0,0])
                
            '''
            print("Average Temperature : {:.5f} ".format(t_stat_avg)+chr(176)+"C.")
            print("Atmospheric Pressure: {} Pa.".format(p_atm))
            '''
            
            # Calculate the fluid density (using thermofluid data)
            rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
                
            # Calculate the wind tunnel free-stream velocity
            velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5
        
            # Calculate the Reynold's number of the experiment
            Re = (rho * velocity * 0.225)/(17.9*10**-6)    
            print(Re)
            
            # Create arrays to store new coefficient values
            yaw_coefficient = np.zeros(angl_readings.shape)
                
            # Loop through all rows and calculate coefficients
            for row in range(0,len(angl_readings)):
                
                p_l = pres_readings[row, (lookup_dict["PRB-L"]-1)]
                p_c = pres_readings[row, (lookup_dict["PRB-C"]-1)]
                p_r = pres_readings[row, (lookup_dict["PRB-R"]-1)]
                  
                yaw_coefficient[row] = YawCoefficient(p_l, p_r, p_c)
            
            data = np.concatenate((angl_readings,yaw_coefficient), axis=1)
            np.savetxt("yaw_coeff_Re_{}.txt".format(setting), data, delimiter='\t', fmt='%10.5f', newline='\n', encoding=None)

#==============================================================================

PlotYawCoefficient(save_fig = False)
#SaveYawCoefficient()