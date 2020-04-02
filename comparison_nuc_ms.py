from ms_functions import *
from constants import *
from rk45_integration import *

import matplotlib.pyplot as plt
import numpy as np


ms=np.loadtxt("Generated_params.csv")

x=ms[:,2] #T 
y=ms[:,4] #L
r=ms[:,0] #r
m=ms[:,3] #M


n=list(map(lambda z, k: (z/(4*np.pi*sigma*(k**2)))**(1/4), y, r))
plt.plot(n,[i/L_sun for i in y],'bo', label='Calculated')
d = list(x)[:18] + n[18:]
plt.plot(d,[i/L_sun for i in y],'bo')
plt.xlim(30000,1000)
plt.xlabel('T (K)')
plt.ylabel(r'$L/L_{sun}$')
plt.xscale('log')
plt.yscale('log')
plt.title('HR Diagram')


nuc=np.loadtxt("Generated_params_Nuc_1e-9.csv")
x_1=nuc[:,2] #T 
y_1=nuc[:,4] #L
r_1=nuc[:,0] #r
m_1=nuc[:,3] #M

n_1=list(map(lambda g, l: (g/(4*np.pi*sigma*(l**2)))**(1/4), y_1, r_1))
#plt.plot(n_1,[i/L_sun for i in y_1],'k', label='Calculated')
d_1 = list(x_1)[:18] + n[18:]
plt.plot(d_1,[i/L_sun for i in y_1],'ko')
plt.xlim(30000,1000)
plt.xlabel('T (K)')
plt.ylabel(r'$L/L_{sun}$')
plt.xscale('log')
plt.yscale('log')
plt.title('HR Diagram PP 1e-9')
plt.savefig('Nuclear_Comparison_PP_1e-9.png')
plt.show()
