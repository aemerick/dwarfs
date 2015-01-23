import numpy as np


def spherical_NFW(rProf, sim_TCloud, sim_TAmbient, sim_bParam,
                  sim_RL, sim_rhoRL, sim_rhoCenter,
                  sim_rho1rm, sim_rho2rm, sim_gamma = 1.666667, sim_mu = 1.0):

    rho = np.zeros(np.size(rProf))
    pressure = np.zeros(np.size(rProf))
    temperater = np.zeros(np.size(rProf))


    # define some constants
    sim_pmass = 1.6726231E-24 # g
    sim_boltz = 1.380658E-16
    sim_pi    = np.pi
    sim_newton = 6.67259E-8
    

    kfactor = sim_boltz / (sim_mu * sim_pmass)

    cs1 = (sim_gamma * kfactor * sim_TCloud)**0.5
    cs2 = (sim_gamma * kfactor * sim_TAmbient)**0.5

    cPhi = 4.0 * sim_pi * sim_newton * sim_rhoCenter * sim_bParam**3


    cRho2 = sim_rhoRL * np.exp((-cPhi / (cs2*cs2*sim_RL))*np.log(1.0 + sim_RL / sim_bParam))


    # uses the bisection search method to calculate rm 
    RM = find_rm(cRho2, cs2, sim_rho2rm, cPhi, sim_bParam, 6.171E21)

    cRho1 = sim_rho1rm * np.exp((-cPhi/(cs1*cs1*RM))*np.log(1.0 + (RM/sim_bParam)))
    
    rho1 = np.zeros(np.size(rProf)); rho2 = np.zeros(np.size(rProf))
    ## calculate rho1 and rho2
    rho1[rProf == 0] = cRho1 * np.exp((cPhi/(cs1*cs1)) / sim_bParam)
    rho2[rProf == 0] = cRho2 * np.exp((cPhi/(cs2*cs2)) / sim_bParam)
    
    rho1[rProf > 0] = cRho1 * np.exp( cPhi / (cs1*cs1*rProf[rProf>0])*\
                              np.log(1.0 + rProf[rProf>0] / sim_bParam))
    rho2[rProf > 0] = cRho2 * np.exp( cPhi / (cs2*cs2*rProf[rProf>0])*\
                              np.log(1.0 + rProf[rProf>0] / sim_bParam))
   
    rho[rProf < RM] = rho1[rProf < RM]
    rho[rProf > RM] = rho2[rProf > RM]
    
    ##
    pressure[rProf < RM] = rho1[rProf<RM] * kfactor * sim_TCloud
    pressure[rProf > RM] = rho2[rProf>RM] * kfactor * sim_TAmbient 
    
   
    temperature = pressure / (rho * kfactor)    


    return rho, pressure, temperature
    
                              
def find_rm(cp2, cs2, rho2rm, cPhi, bparam, length):
    
    eps = 1.0E-7
    nmax = 2000
    rhi = 2.0*length
    rlo = 0.0
    ic = 0

    rhomid = 1.0E-32

    while (abs(rhomid - rho2rm)/rho2rm > eps) and (ic <= nmax):

        rmid   = 0.5 * (rhi + rlo) 
        rhomid = cp2 * np.exp( (cPhi/(cs2*cs2*rmid))*np.log(1.0 + rmid/bparam))

        if (rhomid > rho2rm):
           rlo = rmid
        else:
           rhi = rmid
    
        ic = ic + 1

    return rmid
