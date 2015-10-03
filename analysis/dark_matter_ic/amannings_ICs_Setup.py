import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt
from scipy.misc import derivative
from scipy import integrate
from scipy.optimize import newton, brentq, bisect


# Derivatives and equations for setting up NFW profile and distribution
"""Will give the derivatives that will be plugged into the distribution
function(B&T 4-140b). Also defines density profile for NFW profile.
Uses (psider) equations found in art of computational science piet hut and jun makino
Equations also found in Kazantzidis, Magorrian, Moore 2004 paper.
Values are in solar mass/pc units. cgs units made numbers too large and crashed derivative function"""

# Radius of system
R = 3528.45

#scale density
rhos = 0.02705

#scale radius
rs = 795.0

#virial radius
rvir = 1411.38

#decay radius
rdec = .5 * rvir

# Gravitational Constant
G = 4.302e-3

#--------------------------------------------------------------------------------------#
#                                                                                      #   
#                                         FUNCTIONS                                    #
#--------------------------------------------------------------------------------------#
"""
Functions to calculate density, potential and the energy distribution
"""



# -----------------------------------------------------  #
#                                                        # 
#                         DENSITY                        #
#                                                        #
# -----------------------------------------------------  # 
""" 
Defines density for the NFW profile.
Following three functions pertain to rho(r).
The derivative is taken with respect to r
to be used in the distribution
function. 
"""



def calculate_epsilon(alpha, beta, gamma, virial_radius, decay_radius, scale_radius):
    
    c = virial_radius/scale_radius
    
    d = (-gamma - beta*(c)**alpha) / (1 + c**alpha)

    e = d + (virial_radius/decay_radius)
    
    return e

e = calculate_epsilon(1, 3, 1, rvir, rdec, rs)

# NFW density function
#  #  #  #  #  #  #  #

def calculate_density(radius_array):
    """
    This function takes an array of radii as input
    Uses this info to calculate the density at each
    of these values. outputs them as an array called 
    dens_array. Else statement uses exponential cutoff for
    radii above virial radius (rvir).
    """

    #array of values density of profile
    dens_array = []



    if np.size(radius_array) == 1:
            
        if radius_array < rvir:
            return (rhos/((radius_array/rs)*((1 + radius_array/rs)**2)))
        else:
            return (rhos/((rvir/rs)*((1 + rvir/rs)**2))) *  np.exp(-(radius_array-rvir)/rdec) *((radius_array/rvir)**e)


    
    else:
        
        for radius in radius_array:
            if radius == 0: radius = 1.0 
            #don't want the value to equal zero because we 
            #would then divide by zero and get infinity, so change it to one. 
            #close enough to zero. then continues through the rest of the values
  
            if radius < rvir:
                denominator = (radius/rs)*((1.0 + radius/rs)**2) 
                dens = rhos/denominator         
                dens_array.append(dens)


            else:
                denominator = (rvir/rs)*((1.0 + rvir/rs)**2)
                dens = rhos/denominator
                decay = np.exp(-(radius-rvir)/rdec)
                densdecay = dens*((radius/rvir)**e)*decay
                dens_array.append(densdecay)


    return dens_array


def f1(radius):   
    radius[0] = 1
    if radius.any() < rvir:
        return (rhos/((radius/rs)*((1 + radius/rs)**2)))
    else:
        return (rhos/((rvir/rs)*((1 + rvir/rs)**2))) *  np.exp(-(radius-rvir)/rdec) *((radius/rvir)**e)    


# Derivative of Density function
# drho_dr
#    #    #    #    #    #    #
def der(func, radius_array):
    """
    Takes a function and an array of radii as arguments. 
    takes a derivative of that function, and finds the 
    value of the derivative at each of the values in the 
    radius_array. will take calculate_density as input
    """
    #array of derivative values
    derivative1 = []
    # initialize value
    val = 0

    derivative1 = derivative(func, radius_array)

    derivative1 = np.array(derivative1)
    return derivative1



# Second derivative
#  #  #  #  #  #  # 
def der2(func, radius_array):
    """
    does same thing as previous but is second order derivative
    """    
    #define drho_dr_2 as global to be used in finding the f(e) integrand
    derivative2 = []
    # initialize value
    val = 0

    derivative2 = derivative(func, radius_array, n=2)

    derivative2 = np.array(derivative2)
    return derivative2






# -----------------------------------------------------------------------  #
#                                                                          #
#                                                                          #
#                             PSI/POTENTIAL                                #
#                                                                          #
#  ---------------------------------------------------------------------   #




"""
the below code goes through psi, dPsi/dr, and d2Psi/dr2
"""


# defines the integrand found in the psi and dPsi/dr functions as a callable function

def g1(radius):
    if radius < rvir:
        return radius**2 * (rhos/((radius/rs)*((1 + radius/rs)**2)))
    else:
        return radius**2 * (rhos/((rvir/rs)*((1 + rvir/rs)**2))) *  np.exp(-(radius-rvir)/rdec) *((radius/rvir)**e)



# defines density function. difference from above is the radius is not squared. Used as input for the psi function
def g2(radius):
    if radius < rvir:
        return radius * (rhos/((radius/rs)*((1 + radius/rs)**2)))
    else:
        return radius * (rhos/((rvir/rs)*((1 + rvir/rs)**2))) *  np.exp(-(radius-rvir)/rdec) *((radius/rvir)**e)





# defines psi
#  #  #  #  #  #
def calculate_psi(radius_array, function1, function2):
    """
    This function calculated the value of psi. 
    it takes an array of radii and two different functions as 
    input because the integrands are different for the two
    integrals in this function
    """


    # constant
    cons = -4 * math.pi * G
    # array of integral values (1)
    integral_a = []
    # array of integral values (2)
    integral_b = []
    # array of the values of the two integrals added together
    psi_values = []
    # array to use for limit of second integral since it is from R to infinity or very large number (make same length as 'r' array)

    # error arrays
    error1 = []
    error2 = []
    
    
    integral_a = np.zeros(np.size(radius_array)) # zeros of length of radius array
    integral_b = np.zeros(np.size(radius_array)) # zeros of length of radius array
    psi_values = np.zeros(np.size(radius_array)) # zeros of length of radius array


    #in case there is only one radius value that we want to know the value for, will take and make into a longer array size of two, but return only the second value of the integral
    if np.size(radius_array) == 1: 
        radius_array = np.array([0.001 , radius_array])
        integral_a = np.zeros(np.size(radius_array)) # zeros of length of radius array
        integral_b = np.zeros(np.size(radius_array)) # zeros of length of radius array
        psi_values = np.zeros(np.size(radius_array))
    
    

    for i in np.arange(1, np.size(radius_array)):

        #does first integral found in the potential function and multiplies by 1/same r that is being used in integral
        
        integral_1 ,error1 = integrate.quad(function1, radius_array[0], radius_array[i], full_output=0)
        integral_1 = integral_1 * (1/radius_array[i])
        integral_a[i] = integral_1
        
        #second integral
        
        integral_2 ,error2 = integrate.quad(function2, radius_array[i], float('Inf'), full_output = 0)
        integral_b[i] = integral_2
        
        # multiplies the sum of the two integrals by -4piG
        integral = cons * (integral_a[i] + integral_b[i])
        psi_values[i] = integral


    # only returns the integral value and not the zero value
    if np.size(psi_values) == 2: 
        return psi_values[1]
    else:
        return psi_values


        
           

#  Defines dpsi/dr
#  #  #  #  #  #  #
def psider(function, radius_array):
    """
    Defines dPsi/dr defined as G/R**2 * M(R). 
    g(radius) above defines the integrand so that it can easily be integrated. 
    Integrates over the range of zero to R but takes integral at each r value. 
    Integral from previous r value to that r.
    so as to not repeat what has already been calculated
    Takes a function and an array of radii as input
    """
    
    #four pi
    fourpi = 4 * math.pi
    # array for integral values 
    psider_vals = []
    # separate array for error values
    error3 = []

    psider_vals = np.zeros(np.size(radius_array)) # zeros of length of radius array

    #in case there is only one radius value that we want to know the value for, will take and make into a longer array size of two, but return only the second value of the integral
    if np.size(radius_array) == 1: 
        radius_array = np.array([0.001 , radius_array])
        psider_vals = np.zeros(np.size(radius_array)) # zeros of length of radius array

   
    for i in np.arange(np.size(radius_array)):

        integral ,error3 = integrate.quad(function, 0, radius_array[i], full_output=0)
        integral = integral * (G/radius_array[i]**2) * fourpi 
        psider_vals[i] = integral 
        


    if np.size(psider_vals) == 2: 
        return psider_vals[1]
    else:
        return psider_vals

    


# Defines d2Psi/dr2
#  #  #  #  #  #  #
def psider2(radius_array):
    """      
    Second derivate of Psi.
    takes radius value as input, can take multiple values. 
    outputs second derivative in an array
    calls calculate_density and psider
    """
    #four pi
    fourpi = 4 * math.pi
    #constant
    cons = G * fourpi
    #array of second derivative values Psi
    #psider2_vals = []
 

    print 'rarray',radius_array
    if np.size(radius_array) == 1: 
        radius_array = np.array([0.001 , radius_array])
        density_values = calculate_density(radius_array)
        psi_derivative = psider(g1, radius_array)
        val = (cons * density_values[1]) - ((2/radius_array[1]) * psi_derivative)
        psider2_vals = val  



    else:
       
        print "inside psider2"
        density_values = calculate_density(radius_array)
        print np.size(density_values)
        print "dens", density_values
        psi_derivative = psider(g1, radius_array)
        print "psider", psi_derivative
        print np.size(psi_derivative)
        psider2_vals = np.zeros(np.size(radius_array)) # zeros of length of radius array
            

        for i in np.arange(np.size(radius_array)):
            val = (cons * density_values[i]) - ((2/radius_array[i]) * psi_derivative[i])
            psider2_vals[i] = val


 
    return psider2_vals







# -------------------------------------------------------------------------------------------------------  #
#                                                                                                          #
#                                                                                                          #
#                                       DISTRIBUTION FUNCTION                                              #
#                                                                                                          #
#                                                                                                          #
#  ------------------------------------------------------------------------------------------------------  #



"""
Here f(e) will be calculated. Will first have to calculate
the value of E epsilon, and do first and second order chain rules 
to find drho/dpsi & d2rho/dpsi2 and find some of the parts of the equation
separately then combine in the calculate_fofe function
"""


# Second order chain rule d2rho/dpsi2
#  #  #  #  #  #  #  #  #  #  #  #  #  

def second_order_chain_rule(radius):
    """
    This function defines part of the integrand for f(e).
    defines the second term d2rho/dpsi2 using the second order chain rule
    take radius as input. calls the outside functions and evalutes them inside
    """

    term1 = der2(calculate_density, radius) / (psider(g1, radius)**2) #first term
    term2 = der(calculate_density, radius) / psider2(radius) #second term
    chain_rule_2 = term1 + term2 
            


    return chain_rule_2




# first order chain rule drho/dpsi
#  #  #  #  #  #  #  #  #  #  #  #

def first_order_chain_rule(array1, array2):
    """
    takes arrays of derivative values as input
    outputs drho/dpsi to be used in f(e) (added term)
    used only evaluated at psi = 0 (so only one value)
    array1 = drho_dr
    array2 = dpsi_dr    
    """
    # array for drho/dpsi values

    chain1 = []
    

    chain1 = np.zeros(1) # zeros of length of radius array


    value = array1[0]/array2[0]
    chain1[0] = value

    print "chain1", chain1
    return chain1
    



def root(radius_array, energy_val):
    
    """
    this function takes an array of radii as input and an energy value
    defines the function that will be used in the root finder
    in order to find the value of r where E(epsilon) is the largest, 
    so it can be used as limit of integration and plugged into integrand
    root finder is executed in the next function for fofe
    """
    

    return energy_val + calculate_psi(radius_array, g1, g2)




def calculate_fofe(energy_array, chain_rule_1_array):
    
    """
    this function calculates the energy distribution
    uses simpson method for numerical integration given an array of values
    energy_array = energy_values
    chain_rule_1_array = drho_dpsi
    calculates the rmax for the given value of E. 
    uses this rmax to calculate the values within the integral then the integral itself.  
    integral is taken with respect to r so there is an extra dpsi/dr term seen in the "inte" term
    """


    #constant
    cons = 1/(math.sqrt(8)*(math.pi**2))
    error4 = []
    fofe = np.zeros(np.size(energy_array))
    i = 0
  
    
    for E in energy_array:
        
        try:
         
            #use root finder to find value of rmax for given E value
            maximum_r = brentq(root, 0.001, R, args=(E,))
        

            #define added term that is supposedly equal to 0    
            added_term = (1/math.sqrt(E)) * chain_rule_1_array[0]
         
            def integrand(radius):
                if (E + calculate_psi(radius, g1, g2) == 0):
                    return math.sqrt(E) * (2/math.sqrt(2)) * second_order_chain_rule(radius)

                else:
                    return (1 / math.sqrt(E - calculate_psi(radius, g1, g2))) * second_order_chain_rule(radius) * psider(g1, radius)

   
       
          
            integral ,error4 = integrate.quad(integrand, 0, maximum_r)
            integral = integral + added_term
            integral = integral * cons
            fofe[i] = integral
            
            i = i + 1

   

        except ValueError:
        
            print "Energy value has no zero"
     

           
    return fofe



 


#----------------------------------------------------------------------------------------#
#                                 COMPUTING STUFF                                        #
#                                                                                        #
# now calling functions and storing the values in variables                              #
# also plots the density function, derivatives, and psi functions                        #
#----------------------------------------------------------------------------------------#


# defining variable r in range from 0 to R
# 3rd number determines the number of values in the array

r= np.linspace(0.0, R, 10)
r[0] = r[0] + 1 # changes first number to small number so that there is no division by 0


energy_values = np.linspace(0, 600, 10)    


density = calculate_density(r)
print "density", density

drho_dr = der(f1, r)
print "drho/dr", drho_dr

drho_dr_2 = der2(f1, r)
print "drho/dr_2", drho_dr_2


psi = calculate_psi(r, g1, g2)
print "psi", psi

dpsi_dr = psider(g1, r)
print "dpsi/dr", dpsi_dr

dpsi_dr_2 = psider2(r)
print "psider2", dpsi_dr_2

drho_dpsi = first_order_chain_rule(drho_dr, dpsi_dr)
#print "drho/dpsi", drho_dpsi


forfe = calculate_fofe(energy_values, drho_dpsi)
print "fofe", forfe


"""
# dictionary that stores all the information pertaining to a specific radius
information_on_radius = {}


for i in range(len(r)):
    information_on_radius[r[i]] = [density[i], drho_dr[i], drho_dr_2[i], psi[i], dpsi_dr[0], dpsi_dr_2[i], d2rho_dpsi2[i]]






"""







plt.plot(r, density, 'r')
plt.semilogy()
plt.ylabel("log Density $({SM}/{pc^3})$")
plt.xlabel("Radius (pc)")
plt.title("Density Plot")
plt.show()


plt.plot(r, abs(drho_dr), 'g')  
plt.semilogy()
plt.ylabel('log Derivative $({$M_{\odot}}/{pc^4})$')
plt.xlabel("Radius (pc)")
plt.title("Derivative of Density (absolute value)")
plt.show()



plt.plot(r, drho_dr_2, 'b')  
plt.semilogy()
plt.ylabel("log 2nd Derivative $({M_{\odot}}/{pc^5})$")
plt.xlabel("Radius (pc)")
plt.title("2nd Derivative of Density")
plt.show()


plt.plot(r, psi, 'g')
plt.ylabel("Gravitational Potential $({M_{\odot}*pc^2}/{s^2})$")
plt.xlabel("Radius (pc)")
plt.title("Psi")
plt.show()


plt.plot(r, dpsi_dr, 'b')
plt.ylabel("dPsi/dr")
plt.xlabel("Radius (pc)")
plt.title("dPsi/dr")
plt.show()


plt.plot(r, dpsi_dr_2, 'r')
plt.semilogy()
plt.ylabel(" log 2nd Derivative of Psi")
plt.xlabel("Radius (pc)")
plt.title("2nd Derivative of Psi")
plt.show()


plt.plot(r, forfe, 'g')
#plt.semilogy()
#plt.semilogx()
plt.ylabel("Probability")
plt.xlabel("Radius (pc)")
plt.title("Energy distribution f(e) (absolute value)")
plt.show()

"""
plt.plot(r, E_array, 'b')
plt.semilogy()
plt.ylabel("Energy (eV)")
plt.xlabel("Radius (pc)")
plt.title("Energy")
plt.show()
"""













#Alexandra Mannings
