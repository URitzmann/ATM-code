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

   timestep=100000 # number of steps

   #------------------------------------------------------------------------------------------------------------------
   # reading the input files of the simulation: Spin_STM-File, NN-table, DMI-vector
   #------------------------------------------------------------------------------------------------------------------
   job="mag.dat"
   infile=open(job,"r")

   #------------------------------------------------------------------------------------------------------------------------------
   # defining output arrays
   #--------------------------------------------------------------------------------------------------------------------------
   magx_array=[]
   magy_array=[]
   magz_array=[]
   time_array=[]

   #--------------------------------------------------------------------------------------------------------------------------------
   # reading and extracting the input file
   #-----------------------------------------------------------------------------------------------------------------------------------
   for index, line in enumerate(infile):
       current_line=line.split()
       time_array.append(int(current_line[0]))
       magx_array.append(float(current_line[1]))
       magy_array.append(float(current_line[2]))
       magz_array.append(float(current_line[3]))
   #--------------------------------------------------------------------------------------------------------------------------------
   # fit of the data
   #-----------------------------------------------------------------------------------------------------------------------------------
   def mag(x,a):
       return a+0*x

   gmodel = Model(mag)
   result = gmodel.fit(magx_array, x=time_array, a=0.5)
   params = gmodel.make_params()
   print(result.fit_report())

   #---------------------------------------------------------------------------------------------------------------------------
   #creating the plots
   #--------------------------------------------------------------------------------------------------------------------------
   pp = PdfPages('mag_0.5.pdf')
   plt.plot(time_array,magx_array, label='m_x')
   plt.plot(time_array, result.best_fit, label='Fit')
   plt.xlabel('time')
   plt.ylabel('magnetization comp.')
   plt.savefig(pp, format='pdf')
   pp.close()
   plt.legend()
   plt.show()


plotfile()
