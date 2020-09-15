import matplotlib.pyplot as plt
import numpy as np
import math

from lmfit import Model
from matplotlib.backends.backend_pdf import PdfPages

def plotfile():
#------------------------------------------------------------------------------------------------------------------
# This function plots magnetization and spin current data
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------------------------------------------
   LZ=256    #length of the system
   alpha=0.05 # damping constant
   omega=0.5 # excitation frequency
   output_number=10 # number of timesteps where output was created
   timestep=100000 # number of steps between outputs

   #------------------------------------------------------------------------------------------------------------------
   # reading the input files of the simulation: Spin_STM-File, NN-table, DMI-vector
   #------------------------------------------------------------------------------------------------------------------
   job="mag_0.500000.dat"
   infile=open(job,"r")

   #------------------------------------------------------------------------------------------------------------------------------
   # defining output arrays
   #--------------------------------------------------------------------------------------------------------------------------
   magx_array=[]
   magy_array=[]
   magz_array=[]
   pos_array=[]

   #--------------------------------------------------------------------------------------------------------------------------------
   # reading and extracting the input file
   #-----------------------------------------------------------------------------------------------------------------------------------
   for index, line in enumerate(infile):
       if(index>LZ*(output_number-1)):
           current_line=line.split()
           pos_array.append(int(current_line[0]))
           magx_array.append(float(current_line[1]))
           magy_array.append(float(current_line[2]))
           magz_array.append(float(current_line[3]))
   #--------------------------------------------------------------------------------------------------------------------------------
   # fit of the data
   #-----------------------------------------------------------------------------------------------------------------------------------
   def wave(x, a, b,c,d):
       return a*np.cos(b*x + c)*np.exp(-x/d)

   gmodel = Model(wave)
   result = gmodel.fit(magy_array, x=pos_array, a=0.1, b=0.5, c=3.0, d=10)
   params = gmodel.make_params()
   print(result.fit_report())

   result2 = gmodel.fit(magz_array, x=pos_array, a=0.1, b=0.5, c=3.0, d=10)
   params2 = gmodel.make_params()
   print(result2.fit_report())

   #---------------------------------------------------------------------------------------------------------------------------
   #creating the plots
   #--------------------------------------------------------------------------------------------------------------------------
   pp = PdfPages('mag_0.5.pdf')
   plt.plot(pos_array,magy_array, label='m_y')
   plt.plot(pos_array,magz_array, label='m_z')
   plt.plot(pos_array, result.best_fit, label='Fit: m_y')
   plt.plot(pos_array, result2.best_fit, label='Fit: m_z')
   plt.xlabel('position (a)')
   plt.ylabel('magnetization comp.')
   plt.savefig(pp, format='pdf')
   pp.close()
   plt.legend()
   plt.show()


plotfile()
