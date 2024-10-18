import os, errno
import numpy as np


class utils(object):

    def snp2S(data, tag=None):
    
        '''
        Convert the s2p (Touchstone) in Scattering parameters (S-matrix)
        

        SnP (Touchstone) File Format
        
        S2P: each record contains 1 stimulus value and 4 S-parameters (total of 9 values)
        Stim (freq)  Real (S11)  Imag(S11)  Real(S21)  Imag(S21)  Real(S12)  Imag(S12)  Real(S22)  Imag(S22)
        
        '''

        data  = np.array(data)
        ncols = 9
        
        if data.shape[1]<ncols:
            data=np.hstack((data, np.zeros((data.shape[0], ncols-data.shape[1]))))
        
        idxs   = [i for i in range(data.shape[1]-1) if i%2!=0]
        n      = int(len(idxs)/2)  # Replace with any even number
        labels = [f"S{i}{j}" for i in range(1, n+1) for j in range(1, n+1)]

        #print(len(idxs)/2)
        #print(data.shape)
        
        SM={labels[k]:data[:,i]+1j*data[:,i+1] for k, i in enumerate(idxs)}

        '''
        SM={'freq': data[:,0],
            'S11' : data[:,1]+1j*data[:,2],
            'S12' : data[:,3]+1j*data[:,4],
            'S21' : data[:,5]+1j*data[:,6],
            'S22' : data[:,7]+1j*data[:,8]}
        '''
        
        SM.update({'freq': data[:,0]})
        
        return SM.get(tag, SM)
    

    def snp2Y(data, tag=None, Z0=50):
    
        '''
        Convert the s2p (Touchstone) in admittance parameters (Y-matrix)
        
        '''
        
        S  = utils.snp2S(data)
        dS = (1+S['S11'])*(1+S['S22'])-S['S12']*S['S21']

        Y={'freq': S['freq'],
           'Y11' : ((1-S['S11'])*(1+S['S22'])+S['S12']*S['S21'])/dS*1/Z0,
           'Y12' : 2*S['S12']/dS*Z0,
           'Y21' : 2*S['S21']/dS*Z0  ,
           'Y22' : ((1+S['S11'])*(1-S['S22'])+S['S12']*S['S21'])/dS*1/Z0
           }

        return Y.get(tag, Y)


    def snp2Z(data, tag=None, Z0=50):
    
        '''
        Convert the s2p (Touchstone) in Impedance parameters (Z-matrix)
        
        '''
        
        S  = utils.snp2S(data)
        dS = (1-S['S11'])*(1-S['S22'])-S['S12']*S['S21']

        Z={'freq': S['freq'],
           'Z11' : ((1+S['S11'])*(1-S['S22'])+S['S12']*S['S21'])/dS*Z0,
           'Z12' : -2*S['S12']/dS*Z0,
           'Z21' : -2*S['S21']/dS*Z0  ,
           'Z22' : ((1-S['S11'])*(1+S['S22'])+S['S12']*S['S21'])/dS*Z0
           }

        return Z.get(tag, Z)

    

    def snp2VSWR(data, tag=None):
        
        S  = utils.snp2S(data)

        V={'freq' : S['freq'],
           'VSWR1': (1+np.abs(S['S11']))/(1-np.abs(S['S11'])),
           'VSWR2': (1+np.abs(S['S22']))/(1-np.abs(S['S22']))
           }
                
        return V.get(tag, V)


    def snp2Zin(data, tag=None, Z0=50):
        
        S  = utils.snp2S(data)

        Z={'freq': S['freq'],
           'Zin1': (1+S['S11'])/(1-S['S11'])*Z0,
           'Zin2': (1+S['S22'])/(1-S['S22'])*Z0,
           } 
                
        return Z.get(tag, Z)


    def snp2LC(data, tag=None, Z0=50):

        '''
        Capacitance1 is the capacitance of a one-port or two-port circuit assuming a series RC
        https://www.sonnetsoftware.com/support/help-18/Sonnet_Suites/Sonnet%20Suites%20Documentation.html?Capacitance1.html

        Capacitance2 is the series capacitance between any pair of ports assuming a Pi-model as shown below
        https://www.sonnetsoftware.com/support/help-18/Sonnet_Suites/Sonnet%20Suites%20Documentation.html?Inductance1.html

        Inductance1 is the inductance of a one-port or two-port circuit assuming a series RL
        https://www.sonnetsoftware.com/support/help-18/Sonnet_Suites/Sonnet%20Suites%20Documentation.html?Inductance1.html

        Inductance2 is the series inductance between any pair of ports assuming a PI-model as shown below
        https://www.sonnetsoftware.com/support/help-18/Sonnet_Suites/Inductance2.html

        '''

        
        Y = utils.snp2Y(data, Z0=Z0)

        LC = {'freq': Y['freq'],
              'C1'  : -1/(2*np.pi*Y['freq']*np.imag(1/Y['Y11'])),
              'C2'  : -1/(2*np.pi*Y['freq']*np.imag(1/Y['Y22'])),
              'L1'  :  1/(2*np.pi*Y['freq'])*np.imag(1/Y['Y11']), 
              'L2'  : -1/(2*np.pi*Y['freq'])*np.imag(1/Y['Y21'])
              } 
        
        
        
        return LC.get(tag, LC)
    
    def YtoC1(f, Y11):

        '''
        Capacitance1 is the capacitance of a one-port or two-port circuit assuming a series RC
        https://www.sonnetsoftware.com/support/help-17/Sonnet_Suites/Sonnet%20Suites%20Documentation.html?Capacitance1.html

        '''
        return -1/(2*np.pi*f*np.imag(1/Y11))



    
    def YtoL1(f, Y11):

        '''
        Inductance1 is the inductance of a one-port or two-port circuit assuming a series RL
        https://www.sonnetsoftware.com/support/help-17/Sonnet_Suites/Sonnet%20Suites%20Documentation.html?Inductance1.html

        '''
        return 1/(2*np.pi*f)*np.imag(1/Y11)


    def YtoL2(f, Y22):

        '''
        Inductance2 is the series inductance between any pair of ports assuming a PI-model as shown below
        https://www.sonnetsoftware.com/support/help-17/Sonnet_Suites/Inductance2.html

        '''
        return -1/(2*np.pi*f)*np.imag(1/Y22)

    
    def snp2ABCD(data, tag=None, Z0=50):
        
        S  = utils.snp2S(data)
        dS = 2*S['S21']

        ABCD={'freq': S['freq'],
              'A'   : ((1+S['S11'])*(1-S['S22'])+S['S12']*S['S21'])/dS,
              'B'   : ((1+S['S11'])*(1+S['S22'])-S['S12']*S['S21'])/dS*Z0,
              'C'   : ((1-S['S11'])*(1-S['S22'])-S['S12']*S['S21'])/dS/Z0,
              'D'   : ((1-S['S11'])*(1+S['S22'])+S['S12']*S['S21'])/dS
              } 
                
        return ABCD.get(tag, ABCD)


    def toS21(I,Q, isdB=False):
        '''
        From I, Q to S21

        '''
        return 20*np.log10(np.sqrt(I**2+Q**2)) if isdB else np.sqrt(I**2+Q**2)


    def pol2cart(radii, angles):
        '''
        From polar to rectangular (or cartesian) 

        '''
        return radii * np.exp(1j*angles)

    def cart2pol(x):
        '''
        From rectangular (or cartesian) to polar 

        '''

        return abs(x), np.angle(x)

    

    def sid2TL(data, tag=None):

        
        '''
        Sonnet sid response format
        
        P <number> 56.6767709 0.0 189.369804 0.0 0.0 4.55225e-4
        P <number> <Eff real> <Eff imag> <Z0 real> <Z0 real> 0.0 4.55225e-4
        
        ! P1 F=2.0 Eeff=(56.6771 + j0) Z0=(189.371 + j0)
        ! P2 F=2.0 Eeff=(56.6789 + j0) Z0=(189.368 + j0)
        '''
        data=np.array(data)

        TL=dict()
        
        for p in np.unique(data[:,0]):

            plabel='P'+str(int(p))
            
            TL.setdefault(plabel,{})
            tmparr=data[np.where(data[:,0]==p)[0],:][:,1:]

            TL[plabel].update({'freq': tmparr[:,0]})
            TL[plabel].update({'eeff': np.sqrt(tmparr[:,1]**2+tmparr[:,2]**2)})
            TL[plabel].update({'Z0'  : np.sqrt(tmparr[:,3]**2+tmparr[:,4]**2)})
            
        
        return TL.get(tag, TL)


    

    def dBm2watt(dBm, Z=50, R=50):
        '''
        dBm to watt

        '''
        return 10**(dBm/10)/1000
    
    def dBm2voltage(dBm, Z=50, R=50):
        '''
        dBm to voltage

        '''
        return np.sqrt(Z/1000)*10**(dBm/20)

    def dBm2current(dBm, Z=50, R=50):
        '''
        dBm to current
        
        '''
        return np.sqrt(Z/1000)*10**(dBm/20)/R


    def current2dBm(I, Z=50,R=50):
        '''
        current to dBm
        
        '''
        return 20*np.log10((I*R)/np.sqrt(Z/1000))

    
    def voltage2dBm(V, Z=50, R=50):
        '''
        voltage to dBm
        
        '''
        return 20*np.log10(V/np.sqrt(Z/1000))


    def watt2dBm(W, Z=50, R=50):
        '''
        watt to dBm
        
        '''
        return 10*np.log10(1000*W)
    

    def toMag(i, r, isdB=False):
        '''
        From imaginary to mangitude 

        '''
        return 20*np.log10(np.sqrt(i**2+r**2)) if isdB else np.sqrt(i**2+r**2)


    def toS21(I, Q, isdB=False):
        '''
        From I, Q to S21

        '''
        return rfutils.toMag(I,Q, isdB) 

