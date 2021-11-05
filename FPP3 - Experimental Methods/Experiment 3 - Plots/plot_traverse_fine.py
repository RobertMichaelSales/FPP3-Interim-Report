""" Created on Sun Oct 31 2021, Author: Robert Sales, File: plot_dynamic_pressure_coefficient.py """

import os, csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

def Interpolate(mode, p_l, p_c, p_r, p_0, p_s):
    
    yaw_probe = GetYawProbe(mode, p_l, p_c, p_r)
    
    p_0_probe, p_0_coeff = Getp_0Probe(mode, p_l, p_c, p_r, yaw_probe)
    
    p_s_probe = Getp_sProbe(mode, p_l, p_c, p_r, yaw_probe, p_0_probe)
    
    return yaw_probe, p_0_probe, p_s_probe, p_0_coeff
#==============================================================================

#==============================================================================

def GetYawProbe(mode, p_l, p_c, p_r):
    
    # Configure the filename
    filename = "yaw_coeff_Re_{}.txt".format(mode)
        
    # Load in data from filename
    data = LoadCSV(filename)
    
    # Calculate the local yaw coefficient
    yaw_coeff = (p_l - p_r)/(p_c - 0.5*(p_l + p_r))
    
    # Set up the interpolation data
    yaw_angle_data = data[:,0]
    yaw_coeff_data = data[:,1]
    
    # Set up the interpolation function
    f_interp = interp1d(yaw_coeff_data, yaw_angle_data, assume_sorted = False)
    yaw_probe = f_interp(yaw_coeff)
    
    #print("yaw coeff: {}".format(yaw_coeff))
    #print("yaw_probe: {}".format(yaw_probe))
   
    return yaw_probe
#==============================================================================

#==============================================================================

def Getp_0Probe(mode, p_l, p_c, p_r, yaw_probe):
    
    # Configure the filename
    filename = "stagnation_pressure_coeff_Re_{}.txt".format(mode)
        
    # Load in data from filename
    data = LoadCSV(filename)
    
    # Set up the interpolation data
    yaw_angle_data = data[:,0]
    p_0_coeff_data = data[:,1]
    
    # Set up the interpolation function
    f_interp = interp1d(yaw_angle_data, p_0_coeff_data, assume_sorted = False)
    p_0_coeff = f_interp(yaw_probe)
    
    p_0_probe = p_c + (p_0_coeff*(p_c - 0.5*(p_l + p_r)))
    
    #print("p_0_coeff: {}".format(p_0_coeff))
    #print("p_0_probe: {}".format(p_0_probe))
    
    return p_0_probe, p_0_coeff
#==============================================================================

#==============================================================================

def Getp_sProbe(mode, p_l, p_c, p_r, yaw_probe, p_0_probe):
    
    # Configure the filename
    filename = "dynamic_pressure_coeff_Re_{}.txt".format(mode)
        
    # Load in data from filename
    data = LoadCSV(filename)
    
    # Set up the interpolation data
    yaw_angle_data = data[:,0]
    p_d_coeff_data = data[:,1]
    
    # Set up the interpolation function
    f_interp = interp1d(yaw_angle_data, p_d_coeff_data, assume_sorted = False)
    p_d_coeff = f_interp(yaw_probe)
    
    p_s_probe = p_0_probe - (p_d_coeff*(p_c - 0.5*(p_l + p_r)))
    
    #print("p_d_coeff: {}".format(p_d_coeff))
    #print("p_s_probe: {}".format(p_s_probe))
    
    return p_s_probe
#==============================================================================

#==============================================================================

def PlotYawProbe(save_fig = False):
    
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       "PRB-L":3,
                       "PRB-C":4,
                       "PRB-R":6
                       }
        
        # Initiate plotting
        fig, ax = plt.subplots(1, 1, figsize = figsize)   
        
        mode = "high"
        legend_list = []
        
        for i, setting in enumerate(["12", "25"]):
            
            # Configure the filename
            filename   = "wake_fine_{}ms.txt".format(setting)
        
            # Load in data from filename
            data = LoadCSV(filename)
            
            # Separate readings of pressure and temperature versus time (col 0)
            trav_readings = data[:, 3:4]
            pres_readings = data[:, 8:14]
            temp_readings = data[:, 25:26]
            
            # Extract average conditions
            p_stag_avg = float(sum(pres_readings[:,0])/len(pres_readings[:,0]))
            p_stat_avg = float(sum(pres_readings[:,1])/len(pres_readings[:,1]))
            t_stat_avg = float(sum(temp_readings)/len(temp_readings))
            p_atm = float(data[0,0])
            
            # Calculate the fluid density (using thermofluid data)
            rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
                
            # Calculate the wind tunnel free-stream velocity
            velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5
        
            # Calculate the Reynold's number of the experiment
            Re = (rho * velocity * 0.225)/(17.9*10**-6)
            legend_list.append("Re = {:.2e}".format(int(round(Re, -2))))
            
            # Create arrays to store yaw, stagnation and static pressure values
            yaw_probe_values = np.zeros(trav_readings.shape)
            #p_0_probe_values = np.zeros(trav_readings.shape)
            #p_s_probe_values = np.zeros(trav_readings.shape)
            
            # Initialise the error/uncertainty lists
            y_p_error = np.zeros(trav_readings.shape)
            y_m_error = np.zeros(trav_readings.shape)
            n_readings = 200
            t_multiplier = 1.972 # for 200 readings
            d_alpha = (40)/(2.15646 + 2.09674)
            
            p_bias_1 = 10 # p_l
            p_bias_2 = 10 # p_c
            p_bias_3 = 10 # p_r
            
            # In the free stream
            s_x_1_fs = t_multiplier * ((2.0043)/(n_readings**0.5)) # p_l
            s_x_2_fs = t_multiplier * ((1.2172)/(n_readings**0.5)) # p_c
            s_x_3_fs = t_multiplier * ((1.7571)/(n_readings**0.5)) # p_r
            
            '''      
            # In the wake region
            s_x_1_wk = t_multiplier * ((7.4288)/(n_readings**0.5)) # p_l
            s_x_2_wk = t_multiplier * ((10.152)/(n_readings**0.5)) # p_c
            s_x_3_wk = t_multiplier * ((13.478)/(n_readings**0.5)) # p_r
            '''
            
            for row in range(0,len(trav_readings)):
                
                p_l = pres_readings[row, (lookup_dict["PRB-L"]-1)]
                p_c = pres_readings[row, (lookup_dict["PRB-C"]-1)]
                p_r = pres_readings[row, (lookup_dict["PRB-R"]-1)]
                p_0 = pres_readings[row, (lookup_dict["TPR"]-1)]
                p_s = pres_readings[row, (lookup_dict["SPR"]-1)]
                
                yaw_probe, p_0_probe, p_s_probe, p_0_coeff = Interpolate(mode, p_l, p_c, p_r, p_0, p_s)
                
                yaw_probe_values[row] = yaw_probe
                #p_0_probe_values[row] = p_0_probe
                #p_s_probe_values[row] = p_s_probe
                
                # Calculate and append the measurement errors
                b_yaw_1 = ((p_c - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_1**2)
                b_yaw_2 = ((p_l - p_c)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_2**2)
                b_yaw_3 = ((p_l - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (p_bias_3**2)
                
                # Calculate and append the precision errors
                s_yaw_1 = ((p_c - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_1_fs**2)
                s_yaw_2 = ((p_l - p_c)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_2_fs**2)
                s_yaw_3 = ((p_l - p_r)/((p_c - 0.5*(p_l + p_r))**2))**2 * (s_x_3_fs**2)
                
                b_yaw = (b_yaw_1 + b_yaw_2 + b_yaw_3)**0.5
                b_ang = ((d_alpha**2) * (b_yaw**2))**0.5
                
                s_yaw = (s_yaw_1 + s_yaw_2 + s_yaw_3)**0.5
                s_ang = ((d_alpha**2) * (s_yaw**2))**0.5
                
                u_yaw = (b_ang**2 + s_ang**2)**0.5
            
                y_p_error[row] = (yaw_probe_values[row] + u_yaw)
                y_m_error[row] = (yaw_probe_values[row] - u_yaw)
            
            alpha = [0.2, 0.2]
            
            # Plot the selected data        
            ax.fill_between((((trav_readings[:,0])-81.5)/25),
                            y_p_error[:,0],
                            y_m_error[:,0],
                            alpha = alpha[i]
                            ) 
            
            # Plot the selected data
            ax.plot((((trav_readings[:,0])-81.5)/25),
                    yaw_probe_values,
                    linewidth = linewidth,
                    linestyle = "solid",
                    markersize = markersize
                    )
                    
        # Define plot aesthetics
        plt.title(r"Wake Traverse Relative Flow Angle")
        plt.legend(legend_list, shadow = True, loc = "upper right", ncol = 1)
        plt.xlabel(r"Normalised Traverse Distance ($y/t$)") 
        plt.ylabel(r"Relative Flow Angle ($^\circ$)")
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-25, +25)

        plt.grid()
        if save_fig: plt.savefig('exp_3_yaw_measurement.svg')
        plt.show() 
        
#==============================================================================

#==============================================================================

def PlotStgProbe(save_fig = False, loss_coefficient = True):
    
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       "PRB-L":3,
                       "PRB-C":4,
                       "PRB-R":6
                       }
        
        # Initiate plotting
        fig, ax = plt.subplots(1, 1, figsize = figsize)   
        
        mode = "high"
        legend_list = []
        
        for i, setting in enumerate(["12", "25"]):
            
            # Configure the filename
            filename   = "wake_fine_{}ms.txt".format(setting)
        
            # Load in data from filename
            data = LoadCSV(filename)
            
            # Separate readings of pressure and temperature versus time (col 0)
            trav_readings = data[:, 3:4]
            pres_readings = data[:, 8:14]
            temp_readings = data[:, 25:26]
            
            # Extract average conditions
            p_stag_avg = float(sum(pres_readings[:,0])/len(pres_readings[:,0]))
            p_stat_avg = float(sum(pres_readings[:,1])/len(pres_readings[:,1]))
            t_stat_avg = float(sum(temp_readings)/len(temp_readings))
            p_atm = float(data[0,0])
            
            # Calculate the fluid density (using thermofluid data)
            rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
                
            # Calculate the wind tunnel free-stream velocity
            velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5
        
            # Calculate the Reynold's number of the experiment
            Re = (rho * velocity * 0.225)/(17.9*10**-6)
            legend_list.append("Re = {:.2e}".format(int(round(Re, -2))))
            
            # Create arrays to store yaw, stagnation and static pressure values
            yaw_probe_values = np.zeros(trav_readings.shape)
            p_0_probe_values = np.zeros(trav_readings.shape)
            p_s_probe_values = np.zeros(trav_readings.shape)
            p_0_loss_coeff   = np.zeros(trav_readings.shape)
            mass_flow_values = np.zeros(trav_readings.shape)
            average_loss_num = np.zeros(trav_readings.shape)
            average_loss_den = np.zeros(trav_readings.shape)
            
            # Initialise the error/uncertainty lists
            y_p_error = np.zeros(trav_readings.shape)
            y_m_error = np.zeros(trav_readings.shape)
            n_readings = 200
            t_multiplier = 1.972 # for 200 readings
            d_p0 = 0
            
            p_bias_1 = 10 # p_0
            p_bias_2 = 10 # p_s
            p_bias_3 = 10 # p_l
            p_bias_4 = 10 # p_c
            p_bias_5 = 10 # p_r
            
            
           # In the free stream
            s_x_1_fs = t_multiplier * ((1.9559)/(n_readings**0.5)) # p_0, 
            s_x_2_fs = t_multiplier * ((1.6455)/(n_readings**0.5)) # p_s
            s_x_3_fs = t_multiplier * ((2.0043)/(n_readings**0.5)) # p_l
            s_x_4_fs = t_multiplier * ((1.2171)/(n_readings**0.5)) # p_c
            s_x_5_fs = t_multiplier * ((1.7571)/(n_readings**0.5)) # p_r

            '''      
            # In the wake region
            s_x_1_wk = t_multiplier * ((1.7674)/(n_readings**0.5)) # p_0
            s_x_2_fs = t_multiplier * ((1.4562)/(n_readings**0.5)) # p_s
            s_x_2_wk = t_multiplier * ((7.4288)/(n_readings**0.5)) # p_l
            s_x_3_wk = t_multiplier * ((10.152)/(n_readings**0.5)) # p_c
            s_x_4_fs = t_multiplier * ((13.478)/(n_readings**0.5)) # p_r

            '''
                    
            for row in range(0,len(trav_readings)):
                
                p_l = pres_readings[row, (lookup_dict["PRB-L"]-1)]
                p_c = pres_readings[row, (lookup_dict["PRB-C"]-1)]
                p_r = pres_readings[row, (lookup_dict["PRB-R"]-1)]
                p_0 = pres_readings[row, (lookup_dict["TPR"]-1)]
                p_s = pres_readings[row, (lookup_dict["SPR"]-1)]
                
                yaw_probe, p_0_probe, p_s_probe, p_0_coeff = Interpolate(mode, p_l, p_c, p_r, p_0, p_s)
                
                yaw_probe_values[row] = yaw_probe
                p_0_probe_values[row] = p_0_probe
                p_s_probe_values[row] = p_s_probe
                mass_flow_values[row] = 2*(p_0_probe - p_s_probe)
                
                # Calculate the loss coefficient
                correction = [0.09725651, 0.06702253]
                p_0_loss_coeff[row] = (p_0 - p_0_probe)/(p_0 - p_s) + correction[i]
                
                # Calculate mass averaged stagnation pressure loss
                average_loss_num[row] = mass_flow_values[row] * p_0_loss_coeff[row]
                average_loss_den[row] = mass_flow_values[row]
                
                # Calculate and append the measurement errors
                b_p0_1 = ((0)) * (p_bias_1**2)
                b_p0_3 = ((-0.5 * p_0_coeff)**2) * (p_bias_3**2)
                b_p0_4 = ((1 + p_0_coeff)**2) * (p_bias_4**2)
                b_p0_5 = ((-0.5 * p_0_coeff)**2) * (p_bias_5**2)
                
                # Calculate and append the precision errors
                s_p0_1 = ((0)) * (s_x_1_fs**2)
                s_p0_3 = ((-0.5 * p_0_coeff)**2) * (s_x_3_fs**2)
                s_p0_4 = ((1 + p_0_coeff)**2) * (s_x_4_fs**2)
                s_p0_5 = ((-0.5 * p_0_coeff)**2) * (s_x_5_fs**2)
                
                b_p0 = (b_p0_1 + b_p0_3 + b_p0_4 + b_p0_5)**0.5
                s_p0 = (s_p0_1 + s_p0_3 + s_p0_4 + s_p0_5)**0.5
                
                if loss_coefficient:
                    
                    b_yp_1 = (((p_0_probe - p_s)/(p_0 - p_s)**2)**2) * (p_bias_1**2)
                    b_yp_2 = (((-1)/(p_0 - p_s))**2) * (b_p0**2)
                    b_yp_3 = (((p_0 - p_0_probe)/(p_0 - p_s)**2)**2) * (p_bias_2**2)
                    
                    s_yp_1 = (((p_0_probe - p_s)/(p_0 - p_s)**2)**2) * (s_x_1_fs**2)
                    s_yp_2 = (((-1)/(p_0 - p_s))**2) * (s_p0**2)
                    s_yp_3 = (((p_0 - p_0_probe)/(p_0 - p_s)**2)**2) * (s_x_2_fs**2)
                    
                    b_yp = (b_yp_1 + b_yp_2 + b_yp_3)**0.5
                    s_yp = (s_yp_1 + s_yp_2 + s_yp_3)**0.5
                    
                    u_yp = (b_yp**2 + s_yp**2)**0.5
                    
                    y_p_error[row] = (p_0_loss_coeff[row] + u_yp)
                    y_m_error[row] = (p_0_loss_coeff[row] - u_yp)
                    
                else:
                    u_p0 = (b_p0**2 + s_p0**2)**0.5
            
                    y_p_error[row] = (p_0_probe_values[row] + u_p0)
                    y_m_error[row] = (p_0_probe_values[row] - u_p0)
            
            alpha = [0.2, 0.2]
            
            if loss_coefficient: 
                #Plot the selected data        
                ax.fill_between((((trav_readings[:,0])-81.5)/25),
                                y_p_error[:,0],
                                y_m_error[:,0],
                                alpha = alpha[i]
                                ) 
                
                # Plot the selected data
                ax.plot((((trav_readings[:,0])-81.5)/25),
                        p_0_loss_coeff,
                        linewidth = linewidth,
                        linestyle = "solid",
                        markersize = markersize
                        )
                    
                # Define plot aesthetics
                plt.title(r"Wake Traverse Stagnation Pressure Loss Coefficient")
                plt.legend(legend_list, shadow = True, loc = "upper right", ncol = 1)
                plt.xlabel(r"Normalised Traverse Distance ($y/t$)") 
                plt.ylabel(r"Stagnation Loss Coeff., $Y_p$")
                ax.set_xlim(-2.5, 2.5)
                ax.set_ylim(-0.20, 1.05)
            
            else:
                # Plot the selected data        
                ax.fill_between((((trav_readings[:,0])-81.5)/25),
                                y_p_error[:,0],
                                y_m_error[:,0],
                                alpha = alpha[i]
                                ) 
                
                # Plot the selected data
                ax.plot((((trav_readings[:,0])-81.5)/25),
                        p_0_probe_values,
                        linewidth = linewidth,
                        linestyle = "solid",
                        markersize = markersize
                        )
                    
                # Define plot aesthetics
                plt.title(r"Wake Traverse Guage Stagnation Pressure Distribution")
                plt.legend(legend_list, shadow = True, loc = "upper right", ncol = 1)
                plt.xlabel(r"Normalised Traverse Distance ($y/t$)") 
                plt.ylabel(r"Stagnation Pressure ($Pa$)")
                ax.set_xlim(-2.5, 2.5)
                #ax.set_ylim(-18, +18)

            # Determine mass averaged pressure loss coefficient            
            average_loss_coeff_num = np.trapz(average_loss_num[:,0], x = (((trav_readings[:,0])-81.5)/25))
            average_loss_coeff_den = np.trapz(average_loss_den[:,0], x = (((trav_readings[:,0])-81.5)/25))
            
            #print(min(p_0_loss_coeff))
            
            average_loss_coeff = average_loss_coeff_num / average_loss_coeff_den
            print("Average pressure loss coefficient at {}m/s = {}".format(setting, average_loss_coeff))
        
        plt.grid()
        
        if save_fig: plt.savefig('exp_3_stag_loss_coefficient.svg')
        plt.show() 
        
#==============================================================================

PlotYawProbe(save_fig = False)
PlotStgProbe(save_fig = False)



