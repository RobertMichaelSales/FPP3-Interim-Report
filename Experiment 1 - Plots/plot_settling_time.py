""" Created on Fri Oct 15 2021, Author: Robert Sales, File: plot_settling_time.py """

import os, csv
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================

plt.rcParams.update({"text.usetex": False,
                     "font.family": "Times New Roman",
                     "font.size"  : 14,
                     "font.weight": "normal",
                     "figure.dpi" : 500
                     })

figsize = (12,6)
linewidth = '1.0'

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

def SettlingTime(save_fig = False, mode = "p", t_min = 0, t_max = 1305):
    
    # Load in data
    data = LoadCSV("settling_time.txt")
    
    # Separate readings of pressure and temperature versus time (col 0)
    pres_readings = np.concatenate((data[:,0:1],data[:,3 :19]), axis=1)
    temp_readings = np.concatenate((data[:,0:1],data[:,20:21]), axis=1)
    
    # Extract average conditions
    temp_average = sum(temp_readings[:,1])/len(temp_readings[:,1])
    pres_atmosph = data[0,-1]
     
    # Initiate plotting
    plt.figure(figsize = figsize)
    
    # Plot pressures
    if mode == "p":
        
        # Define a lookup dictionary to identify each tapping channel
        lookup_dict = {"TPR"  :1,
                       "SPR"  :2,
                       "ISP-1":3,
                       "ISP-2":4,
                       "ISP-3":6,
                       "ISP-4":7,
                       "ISP-5":8,
                       "00-LE":9,
                       "05-TOP":10,
                       "10-TOP":11,
                       "15-TOP":12, # Takes longer to settle, looks less noisy 
                       "25-TOP":13,
                       "35-TOP":14,
                       "49-TOP":15,
                       "52-TE" :16,
                       }
        
        # Loop through all tappings in the lookup dictionary
        for tapping in lookup_dict:
            
            # Plot the data
            plt.plot(pres_readings[t_min:t_max,0], pres_readings[t_min:t_max,lookup_dict[tapping]],linewidth = linewidth)
        
        # Define plot aesthetics
        plt.legend(lookup_dict.keys(), bbox_to_anchor = (0.975,0.85), loc = "upper right", ncol = 5, shadow = True)
        plt.xlabel("Time (seconds)")
        plt.ylabel("Pressure (Pascals)")
        plt.grid()
        plt.show()

    # Plot temperature
    if mode == "t":
        
        plt.plot(temp_readings[t_min:t_max,0], temp_readings[t_min:t_max,1],linewidth = linewidth)
    
        # Define plot aesthetics
        plt.legend("Tunnel Inlet Temperature", loc = "lower right", shadow = True)
        plt.xlabel("Time (seconds)")
        plt.ylabel("Temperature ("+chr(176)+"C)")
        plt.grid()
        plt.show()
    
    # Save the plot as svg figure
    if save_fig: plt.savefig('settling_time.svg')

#==============================================================================

#'''
SettlingTime(save_fig = False, mode = "t", t_min = 9, t_max = 1305)
#'''