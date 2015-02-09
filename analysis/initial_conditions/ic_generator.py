import numpy as np
import cgs as cgs
# this is mosty an exact copy of what is coded up in fortran in FLASH
# this is done mostly ingorning the python-ness of python in order to 
# make sure the logic flows perfectly as coded


def spherical_NFW(rProf, sim_TCloud, sim_TAmbient, sim_bParam,
                  sim_RL, sim_rhoRL, sim_rhoCenter,
                  sim_rho1rm, sim_rho2rm, sim_crho1, sim_gamma = 1.666667, sim_mu = cgs.mu,M_DM=None):

    rho = np.zeros(np.size(rProf))
    pressure = np.zeros(np.size(rProf))
    temperater = np.zeros(np.size(rProf))


    # define some constants
    sim_pmass = cgs.mp
    sim_boltz = cgs.kb
    sim_pi    = np.pi
    sim_newton = cgs.G
    

    kfactor = sim_boltz / (sim_mu * sim_pmass)

    cs1 = (sim_gamma * kfactor * sim_TCloud)**0.5
    cs2 = (sim_gamma * kfactor * sim_TAmbient)**0.5

    if M_DM == None:
        cPhi = 4.0 * sim_pi * sim_newton * sim_rhoCenter * sim_bParam**3
    else:
        cPhi = M_DM * sim_newton

#    cRho2 = sim_rhoRL * np.exp((-cPhi / (cs2*cs2*sim_RL))*np.log(1.0 + sim_RL / sim_bParam))

#    print 'before', cRho2
 #   factor = 1.0E4
 #   cRho2 = sim_rho2rm * (1.0 + cPhi/(cs2*cs2)*\
 #                 (np.log(1.0/(factor*sim_bParam))-np.log(1.0/sim_bParam))/(factor*sim_bParam))

 #   print 'after', cRho2
    # uses the bisection search method to calculate rm 
#    RM = find_rm_2(cRho2, cs2, sim_rho2rm, cPhi, sim_bParam, 2.0*cgs.kpc)

#    cRho1 = sim_rho1rm * np.exp((-cPhi/(cs1*cs1*RM))*np.log(1.0 + (RM/sim_bParam)))
    

     
    cRho1 = sim_crho1 * np.exp(-(cPhi/(cs1*cs1))/sim_bParam)

    RM = find_rm_2(cRho1, cs1, cs2,sim_rho2rm, cPhi, sim_bParam,2.0*cgs.kpc)



    rho1 = np.zeros(np.size(rProf)); rho2 = np.zeros(np.size(rProf))
    ## calculate rho1 and rho2
    rho1[rProf == 0] = cRho1 * np.exp((cPhi/(cs1*cs1)) / sim_bParam)
#    rho2[rProf == 0] = cRho2 * np.exp((cPhi/(cs2*cs2)) / sim_bParam)
 
    rho2[rProf >= 0] = sim_rho2rm
   
    rho1[rProf > 0] = cRho1 * np.exp( cPhi / (cs1*cs1*rProf[rProf>0])*\
                              np.log(1.0 + rProf[rProf>0] / sim_bParam))
#    rho2[rProf > 0] = cRho2 * np.exp( cPhi / (cs2*cs2*rProf[rProf>0])*\
#                              np.log(1.0 + rProf[rProf>0] / sim_bParam))
   
    rho[rProf < RM] = rho1[rProf < RM]
    rho[rProf > RM] = rho2[rProf > RM]
    
    ##
    pressure[rProf < RM] = rho1[rProf<RM] * kfactor * sim_TCloud
    pressure[rProf > RM] = rho2[rProf>RM] * kfactor * sim_TAmbient 
    
   
    temperature = pressure / (rho * kfactor)    


    return RM,rho, pressure, temperature
    
                              
def find_rm_2(cp1, cs1, cs2, rho2rm, cPhi, bparam, length):
    
    eps = 1.0E-7
    nmax = 2000
    rhi = length
    rlo = 0.0
    ic = 0

    rhomid = 1.0E-32

    while (abs(rhomid - rho2rm)/rho2rm > eps) and (ic <= nmax):

        rmid   = 0.5 * (rhi + rlo) 
        #rhomid = cp1 * np.exp( (cPhi/(cs1*cs1*rmid))*np.log(1.0 + rmid/bparam))

        rhomid = cp1 * np.exp((cPhi/(cs1*cs1*rmid))*np.log(1.0+rmid/bparam)) *(cs1*cs1)/(cs2*cs2)


        if (rhomid > rho2rm):
           rlo = rmid
        else:
           rhi = rmid
    
        ic = ic + 1

    return rmid

def find_rm_gatto(rho_center, cPhi, bparam, Pcorona, kTdwarf):

    eps  = 1E-7
    nmax = 2000

    rhi = 3.0856E22 
    rlo = 0.0
    ic  = 0



    Pdwarf = 100.0*Pcorona

    print "corona dwarf", Pcorona, Pdwarf
  
    while ((( np.abs(Pdwarf - Pcorona)/Pcorona > eps) and (ic <= nmax))):

        rmid   = 0.5 * (rhi + rlo) 

        rhomid = rho_center * np.exp(- (cPhi)*(1.0-np.log(1.0 + rmid/bparam)/(rmid/bparam)))
        Pdwarf = rhomid * kTdwarf

        if ((Pdwarf-Pcorona) > 0.0 ):
           rlo = rmid
        else:
           rhi = rmid
    
        ic = ic + 1



    return rmid



def find_rm(cp2, cs2, rho2rm, cPhi, bparam, length):

    eps = 1.0E-7
    nmax = 2000
    rhi = length
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




