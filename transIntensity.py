# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:01:07 2019

@author: zhxxu
"""

import numpy as np
import scipy.constants as const

def f_unit(value,unit):
    # output unit: {f: Hz, omega: rad Hz, lambda: nm}
    if (unit == "nm"):
        f = const.c*1e9/value
    elif (unit == "Hz"):
        f = value
    elif (unit == "cm-1"):
        f = const.c*1e2 * value
    else:
        print("error of unit in f_unit")
        return -1
    return {"f":f, "omega":2*np.pi*f, "lambda":const.c*1e9/f}

def mu_unit(value,unit):
    if (unit == "D" or unit == "Debye"):
        mu = value
    elif (unit == "eA"):
        mu = value / 0.208193
    elif (unit == "ea0"):
        mu = value / 0.393430
    elif (unit == "Cm"):
        mu = value / 3.33564e-30
    else:
        print("error of unit in mu_unit")
        return -1
    return {"D":mu, "Debye":mu, "Cm":mu*3.33564e-30, "ea0":mu*0.393430, "eA":mu*0.208193}
    
def sigma_unit(value,unit1,unit2,unit3):
    # unit1: cm2, m2, unit2: 1, rad, unit3: s-1, Hz, cm-1
    sigma = value
    if (unit1 == "cm2"):
        sigma *= 1e-4
    elif (unit1 != "m2"):
        print("error of unit1 in sigma_unit")
        return -1
    if (unit2 == "rad"):
        sigma /= 2*np.pi
    elif (unit2 != "1"):
        print("error of unit2 in sigma_unit")
        return -1        
    if (unit3 == "cm-1"):
        sigma *= const.c*1e2
    elif (unit3 not in ["Hz", "s-1"]):
        print("error of unit3 in sigma_unit")
        return -1   
    # sigma in unit m2s-1       
    return {"m2s-1":sigma, "cm2s-1":sigma*1e4, "m2rads-1":sigma*2*np.pi, "cm2rads-1":sigma*2*np.pi*1e4,
            "m2cm-1":sigma/(const.c*1e2), "cm2cm-1":sigma*1e2/const.c}

def A_to_sigma(g1,g2,freq,A):
    # output unit: m^2 rad s^-1
    return sigma_unit(g2 * freq["lambda"]**2 * A / (4 * g1) * 1e-18, "m2", "rad", "s-1")

def A_to_f(g1,g2,freq,A):
    # output unit: 1
    return g2/g1 * (2*const.pi*const.epsilon_0*const.m_e*const.c**3) / (freq["omega"]**2 * const.e**2) * A

def A_to_mu(freq,A):
    # output unit: using the unit as key
    mu2 = 3*const.epsilon_0*const.h*const.c**3 / (2*freq["omega"]**3) * A
    mu = np.sqrt(mu2)  
    return mu_unit(mu,"Cm")

def sigma_to_A(g1,g2,freq,sigma):
    A = sigma["m2rads-1"]*1e18*4*g1/ (g2*freq["lambda"]**2)
    return A

def f_to_A(g1,g2,freq,f):
    A = f/(g2/g1 * (2*const.pi*const.epsilon_0*const.m_e*const.c**3) / (freq["omega"]**2 * const.e**2))
    return A

def mu_to_A(freq,mu):
    A = mu["Cm"]**2 / (3*const.epsilon_0*const.h*const.c**3 / (2*freq["omega"]**3))
    return A
 
class transition(object):
    def __init__(self, know_type=None, know_value=None, g1=None, g2=None, freq=None):
        if g1 == None:
            self.g1 = int(input("input g1:\n"))
        else:
            self.g1 = g1
        if g2 == None:
            self.g2 = int(input("input g2:\n")) 
        else:
            self.g2 = g2
        if type(freq) != dict:
            inp = input("input frequency [value, unit], supported unit: cm-1, Hz, nm, example: 140.1 nm:\n").split()
            self.freq = f_unit(float(inp[0]),inp[1])
        else:
            self.freq = freq
        if (know_type == None or know_type not in [1,2,3,4,"A","f","sigma","mu"]):
            self.know_type = int(input("input the type (number) of information given:\n"
                                  "1.A, Einstein coefficient\n"
                                  "2.f, Oscillator strength\n"
                                  "3.sigma, Cross section\n"
                                  "4.mu, Transition dipole moment\n"))
        else:
            self.know_type = know_type
            
        if self.know_type in [1,"A"]: 
            if (know_value == None):
                self.A = float(input("input A Einstein coefficient in unit s-1:\n"))
            else:
                self.A = know_value
                
        elif self.know_type in [2,"f"]:
            if (know_value == None):
                self.f = float(input("input f value, Oscillator strength in unitless:\n"))
            else:
                self.f = know_value
                
        elif self.know_type in [3,"sigma"]:
            if (know_value == None):
                inp = input("input sigma, Cross section, with unit: [cm2 or m2] [1 or rad] [s-1 or Hz or cm-1]:\n"
                            "example, 2.9e-5 m2 rad s-1\n").split()
                self.sigma = sigma_unit(float(inp[0]),inp[1],inp[2],inp[3])
            else:
                self.sigma = know_value
        
        elif self.know_type in [4,"mu"]:
            if (know_value == None):
                inp = input("input mu, Transition dipole moment in unit D/Debye/Cm/ea0/eA:\n").split()
                self.mu = mu_unit(float(inp[0]),inp[1])
            else:
                self.mu = know_value
        
        if self.know_type in [1,"A"]:
            self.A = self.A
        elif know_type in [2,"f"]:
            self.A = f_to_A(self.g1,self.g2,self.freq,self.f)
        elif self.know_type in [3,"sigma"]:
            self.A = sigma_to_A(self.g1,self.g2,self.freq,self.sigma)
        elif self.know_type in [4,"mu"]:
            self.A = mu_to_A(self.freq,self.mu)
        
        self.sigma = A_to_sigma(self.g1,self.g2,self.freq,self.A)
        self.f = A_to_f(self.g1,self.g2,self.freq,self.A)
        self.mu = A_to_mu(self.freq,self.A)
        
        print("Calculation finished with code 0")

print("to start, type \"a = transition()\" and follow the steps\n"
      "then access data with a.A, a.f, a.mu, and a.sigma\n")

"""
Test data:
A = 2.2e8
freq = f_unit(422.7,"nm")
g2 = 3
g1 = 1
sigma = A_to_sigma(g1,g2,freq,A) # 2.9e-5 m2rads-1
f = A_to_f(g1,g2,freq,A)         # 1.7
mu = A_to_mu(freq,A)             # 1.5 eA = 2.4e-29 Cm

A_s = sigma_to_A(g1,g2,freq,sigma)
A_f = f_to_A(g1,g2,freq,f)
A_m = mu_to_A(freq,mu)
"""