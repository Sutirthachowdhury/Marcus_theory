# Parameters
import math
def marcus(g): 
 Vda = (0.01/27.2114) #diabatic coupling
 pi  =  math.pi
 lam =  0.65/27.2114 #solvent reorganization energy
 beta = 1052.85 # temperature of the system

 #--- prefactor term---------------
 A  = (Vda**2) * (pi*beta/lam) **0.5 
 A  = A*41341.37 
 #-------------------------------

 #---- FREE ENERGY SCAN---------
 dG = -g 
 alpha = beta/(4*lam) 
 #---------------------------

 #-------- exponential part-------------
 Exp = math.exp(-alpha*(dG+lam)**2) 


 #---- total rate------
 K = A*Exp
 #----------------

 print K, dG 
 return K

mar = open("marxus.txt","w+") 
gc = 0.0
while gc<=0.045*2:
        k1 = marcus(gc)
        k2 = marcus(gc-1.0/27.2114)
        mar.write("%s\t%s\n"%(gc*27.2114,k1+k2))
        gc+=0.001
mar.close()

# See Marcus theory at  : https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/15%3A_Energy_and_Charge_Transfer/15.05%3A_Marcus_Theory_for_Electron_Transfer
# this is standard Marcus theory