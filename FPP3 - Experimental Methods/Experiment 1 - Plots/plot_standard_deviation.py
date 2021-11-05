""" Created on Fri Oct 15 2021, Author: Robert Sales, File: plot_standard_dev.py """

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

def CalculateStdDev(data):
     
    mean = sum(data) / len(data)
    deviations = [(measurement - mean)**2 for measurement in data]
    variance = sum(deviations)/(len(data)-1)
    standard_deviation = variance**0.5
                          
    return standard_deviation

#==============================================================================

#==============================================================================

def StandardDeviation(save_fig = False, mode = "p", interval = 1):
    
    # Load in data
    data = LoadCSV("standard_deviation.txt")
    
    # Separate readings of pressure and temperature versus time (col 0)
    p_data = data[:,3 :19]
    t_data = data[:,20:21]
    
    # Extract average conditions
    p_stag_avg = float(sum(p_data[:,0])/len(p_data[:,0]))
    p_stat_avg = float(sum(p_data[:,1])/len(p_data[:,1]))
    t_stat_avg = float(sum(t_data)/len(t_data))
    p_atm = float(data[-1,-1])
    
    # Calculate the fluid density (using thermofluid data)
    rho = (p_atm + p_stat_avg)/(287.0 * (t_stat_avg + 273.15))
        
    # Calculate the wind tunnel free-stream velocity
    velocity = ((2*(p_stag_avg - p_stat_avg))/rho)**0.5

    # Calculate the Reynold's number of the experiment
    Re = (rho * velocity * 0.225)/(17.9*10**-6)    
    
    # Create matrices to store std.dev calculations versus sample size (col 0)
    p_std_dev = np.concatenate((data[:,2:3],np.zeros(p_data.shape)),axis=1)
    t_std_dev = np.concatenate((data[:,2:3],np.zeros(t_data.shape)),axis=1)
        
    # Calculate std.devs for pressures
    if mode == "p":
        
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       #"ISP-1":3,
                       #"ISP-2":4,
                       #"ISP-3":6,
                       #"ISP-4":7,
                       #"ISP-5":8,
                       "00-LE":9,
                       #"05-TOP":10,
                       #"10-TOP":11,
                       #"15-TOP":12, 
                       #"25-TOP":13,
                       #"35-TOP":14,
                       #"49-TOP":15,
                       "52-TE" :16,
                       }
        
        # Loop through all tappings in the lookup dictionary
        for tapping in lookup_dict:
            
            # Gather tapping data
            p_tapping = p_data[:,lookup_dict[tapping]-1]
            
            # Loop through all the std.dev sample sizes
            for sample_size in range(1,p_data.shape[0]):
                
                # Obtain a sample from the tapping data
                p_sample = p_tapping[0:sample_size+1]
                # Calculate the sample std.dev and save
                p_std_dev[sample_size,lookup_dict[tapping]] = CalculateStdDev(p_sample)
    
    # Calculate std.devs for temperatures
    if mode == "t":
        
        # Loop through all the std.dev sample sizes
        for sample_size in range(1,t_data.shape[0]):
            
            # Obtain a sample from the tapping data
            t_sample = t_data[0:sample_size+1]
            # Calculate the sample std.dev and save
            t_std_dev[sample_size,1] = CalculateStdDev(t_sample)
    
    # Initiate plotting
    fig, ax = plt.subplots(1, 1, figsize = figsize) 
    
    # Plot pressure std.devs
    if mode == "p":
    
        # Loop through all tappings in the lookup dictionary
        for i, tapping in enumerate(lookup_dict):
            
            marker = ['s', 'v', 'D', 'o']

            # Select the data to plot
            data_to_plot = [p_std_dev[x, lookup_dict[tapping]] for x in range(1,p_std_dev.shape[0], interval)]         
            
            # Plot the selected data
            ax.plot(list(range(1, p_std_dev.shape[0], interval)),
                     data_to_plot,
                     linewidth = linewidth,
                     marker = marker[i],
                     linestyle = "dashed",
                     markersize = markersize)
            
        # Define plot aesthetics
        plt.title("Standard Deviation vs. Sample Size     (at Re = {:.2e})".format(int(round(Re, -2))))  
        plt.legend([r"$P_{0,ref}$", r"$P_{s,ref}$", r"LE Tapping", r"TE Tapping"], loc = "upper right", ncol = 2, shadow = True)
        plt.xlabel("DSA Sample Size")
        plt.ylabel(r"Standard Deviation (Pa)")
        ax.set_ylim([-5,305])
        ax.set_ylim([0,8.5])        
        plt.grid()
        
        if save_fig: plt.savefig('exp_1_standard_deviation_pressure.svg')
        plt.show()     
        
    # Plot temperature std.devs
    if mode  == "t":
        
        # Select the data to plot
        data_to_plot = [t_std_dev[x, 1] for x in range(1,t_std_dev.shape[0], interval)]         
        # Plot the selected data
        ax.plot(list(range(1, t_std_dev.shape[0], interval)),
                 data_to_plot,
                 linewidth = linewidth, 
                 marker = "s",
                 linestyle = "dashed",
                 markersize = markersize)

        # Define plot aesthetics
        plt.title("Standard Deviation vs. Sample Size     (at Re = {:.2e})".format(int(round(Re, -2))))  
        plt.legend([r"$T_{s,inlet}$"], loc = "lower right", shadow = True)
        plt.xlabel("DSA Sample Size")
        plt.ylabel(r"Standard Deviation ("+chr(176)+"C)")
        ax.set_xlim([-5,305])
        ax.set_ylim([0.0015,0.012])        
        plt.grid()
        
        if save_fig: plt.savefig('exp_1_standard_deviation_temp.svg')
        plt.show()

    #return p_std_dev
    
#==============================================================================

a = StandardDeviation(save_fig = True, mode = "t", interval = 6)
