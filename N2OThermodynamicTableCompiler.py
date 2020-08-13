"""
OOrtiz 
N2O Thermodynamic Data Gathering
v:0.0
"""
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy import *


class N20ThermoCompile():
    def __init__(self, TrangeMin, TrangeMax, Tincrement, PrangeMin, PrangeMax, Pincrement):
        self.TrangeMin = TrangeMin
        self.TrangeMax = TrangeMax
        self.Tincrement = Tincrement
        self.PrangeMin = PrangeMin
        self.PrangeMax = PrangeMax
        self.Pincrement = Pincrement
        self.conversions() #conversions occur at start
        self.ReadMe()

    def conversions(self):
        C2K = 273 #Celceius to Kelvin 
        self.TrangeMin += C2K
        self.TrangeMax += C2K

    def dataGather(self):
        print("Gathering Data...")
        #Range from Low to High temperature entered
        for temp in range(self.TrangeMin, self.TrangeMax + self.Tincrement, self.Tincrement):
            dataStore = [[],[], [],[],[],[],[]] #T , Saturation values, P_n, Rho_n, s_n, h_n, u_n
            Tcurr = temp 
            dataStore[0].append(Tcurr)
            #Liquid saturation values
            sv0 = PropsSI('P','T', Tcurr, 'Q', 0, 'N2O')#P_sat_l
            sv2 = PropsSI('D','T', Tcurr, 'Q', 0, 'N2O')#Rho_sat_l 
            sv4 = PropsSI('S','T', Tcurr, 'Q', 0, 'N2O')#s_sat_l
            sv6 = PropsSI('H','T', Tcurr, 'Q', 0, 'N2O')#h_sat_l
            sv8 = PropsSI('U','T', Tcurr, 'Q', 0, 'N2O')#u_sat_l
            #Vapor saturation values
            sv1 = PropsSI('P','T', Tcurr, 'Q', 1, 'N2O')#P_sat_v
            sv3 = PropsSI('D','T', Tcurr, 'Q', 1, 'N2O')#Rho_sat_v
            sv5 = PropsSI('S','T', Tcurr, 'Q', 1, 'N2O')#s_sat_v
            sv7 = PropsSI('H','T', Tcurr, 'Q', 1, 'N2O')#h_sat_v 
            sv9 = PropsSI('U','T', Tcurr, 'Q', 1, 'N2O')#u_sat_v
            #appends first index
            for j in range(0,10):
                satIndex = 'sv'+str(j)
                satIndex = eval(satIndex)
                dataStore[1].append(satIndex)
            for pressure in range(self.PrangeMin, self.PrangeMax + self.Pincrement, self.Pincrement):
                rho = PropsSI('D','T',Tcurr,'P',pressure,'N2O') #pressure
                s = PropsSI('S','T',Tcurr,'P',pressure,'N2O') #entropy
                h = PropsSI('H','T',Tcurr,'P',pressure,'N2O') #enthalpy
                u = PropsSI('U','T',Tcurr,'P',pressure,'N2O') #internal energy
                dataStore[2].append(pressure)
                dataStore[3].append(rho)
                dataStore[4].append(s)
                dataStore[5].append(h)
                dataStore[6].append(u)

            #creates a new txt file for N2O table
            f = open("C:/Users/m_i_d/Desktop/N2O Thermodynamic Property Tables/"+str(temp)+"K Table.txt","w")
            #write current temperature
            f.write("Temperature\n"+str(dataStore[0][0])+"\n")
            #writes title headings 
            f.write("P Sat l\tP Sat v\tRho Sat l\tRho Sat v\ts Sat l\ts Sat v\t\
                h Sat l\th Sat v\tu Sat l\tu Sat v\n")
            #wrires saturation values 
            for space in range(0, len(dataStore[1])):                
                f.write("{:.5f}".format(dataStore[1][space]))
                if space + 1 < len(dataStore[1]): #tabs in between each value
                    f.write("\t")
                else: #unless its the last value in the array
                    f.write("\n")
            f.write("Pressure\tDensity\tEntropy\tEnthalpy\tInternal Energy\n") #writes title headings
            for line in range(0, len(dataStore[2])):
                f.write("{:.2f}".format(dataStore[2][line])+"\t{:.5f}".format(dataStore[3][line])+"\t{:.5f}".format(dataStore[4][line])+\
                    "\t{:.5f}".format(dataStore[5][line])+"\t{:.5f}".format(dataStore[6][line])+"\n")
            f.close()

        print("Data Gathering Completed...")

    def ReadMe(self):
            rm = open("C:/Users/m_i_d/Desktop/N2O Thermodynamic Property Tables/@@Read Me@@.txt","w")
            rm.write("All thermodynamic information in the following tables has been \n gathered from CoolProp\
                @ https://CoolProp.org or https://github.com/CoolProp/CoolProp.")
            rm.close()        


#       (TrangeMin, TrangeMax, Tincrement, PrangeMin, PrangeMax, Pincrement)
NTC = N20ThermoCompile(-50, 36, 1, 100000, 9000000, 25000) #Temp in C, P in Pa 

NTC.dataGather()
