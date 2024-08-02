
'''
 @file    sreader.py
 @author  Andrea Giachero <andrea.giachero@mib.infn.it> 
 @date    28 January 2023
 @brief   Collection of module for interacting with Sonnet file

'''

import os, errno
import numpy as np

from twpazer import utils


class sidreader(object):
    '''
    Sonnet sid file reader
    
    '''
    def __init__(self, filename=None, data=None, fmin=None, fmax=None):
        
        self.__filename = filename
        self.__data     = list() if data is None else data
        self.__read() 

        self.__check_frequency(fmin, fmax)
        self.__arrange()
        
        return

    def __iter__(self):
        return iter(self.__data)
    
    def __read(self):

        if self.__filename is None:
            return

        if not os.path.exists(self.__filename):
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), self.__filename)
        
        with open(self.__filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('!')))

            for line in lines:
                if 'SIDATA' in line:
                    self.__data.append(dict())
                if '=' in line:
                    self.__data[-1].update({line.strip().split('=')[0].strip(): float(line.strip().split('=')[1])})   
                if line[:1].isdigit():
                    if 'data' not in self.__data[-1]:
                        self.__data[-1].setdefault('data',[])
                    self.__data[-1]['data'].append(list(map(float, line.strip().split())))

                elif line[:1]=='P' and line[1:2].isspace() and 'UNDEF' not in line:
                    tag = line[:1]+line[2:3]
                    if tag not in self.__data[-1]:
                        self.__data[-1].setdefault(tag,[])
                    self.__data[-1][tag].append(list(map(float, line[4:].strip().split())))
                    self.__data[-1][tag][-1].insert(0, self.__data[-1]['data'][-1][0])


        return

    def __arrange(self):
        for i in range(len(self.__data)):
            for tag in (tag for tag in ['data', 'P1', 'P2'] if tag in self.__data[i]):
                self.__data[i][tag] = np.array(self.__data[i][tag])
                self.__data[i][tag] = self.__data[i][tag][self.__data[i][tag][:, 0].argsort()]

                if self.__fmin and self.__fmin > np.min(self.__data[i][tag][:,0]):
                    idx = np.argmin(np.abs(self.__fmin-self.__data[i][tag][:,0]))
                    self.__data[i][tag]=self.__data[i][tag][idx:,:]

                if self.__fmax and self.__fmax < np.max(self.__data[i][tag][:,0]):
                    idx = np.argmin(np.abs(self.__fmax-self.__data[i][tag][:,0]))
                    self.__data[i][tag]=self.__data[i][tag][:idx,:]
        
        return

    def __check_frequency(self, fmin, fmax):

        freqs = [np.abs(f) if f is not None else f for f in [fmin, fmax]]
        self.__fmin, self.__fmax = sorted(freqs) if None not in freqs else freqs
        
        return

    #def get_range(self):
    #    return self.__fmin, self.__fmax        

    
    def unique(self, opt=None):

        if not self.__data:
            return None

        def get_k(data):
            return sorted([dd for dd in data.keys() if dd not in ['data', 'P1', 'P2', 'pars']]) 

        def get_v(data):
            return [data[k] for k in get_k(data)]
            
        newdata=list()
        newdata.clear()

        import copy
        for dd in self.__data:
            d = copy.deepcopy(dd)
            if not newdata:
                newdata.append(d)
            else:
                vs = [get_v(n) for n in newdata]
                if get_v(d) not in vs:
                    newdata.append(d)
                else:
                    if opt=='merge':
                        idx = vs.index(get_v(d))
                        for tag in  ['data', 'P1', 'P2']:
                            if tag in newdata[idx] and tag in d:
                                newdata[idx][tag]=np.vstack((newdata[idx][tag], d[tag]))
                            elif tag not in newdata[idx] and tag in d:
                                newdata[idx].update({tag: d[tag]}) if tag in d else None
                            else:
                                pass
        
                
        return newdata
    
    def list_parameters(self):
        return [l for l in self.__data[0] if l not in ['data', 'P1', 'P2', 'pars']]

    def get_parameters(self, tag):
        return [d[tag] for d in self.__data if tag in d] 
    
    def get(self, **kargs):
        return self.__data if not kargs else [d for d in self.__data if all([k in d and kargs[k]== d[k] for k in kargs])]
        
    def compute(self, target='S'):
        {'S'   : self.__computeS,
         'Y'   : self.__computeY,
         'Z'   : self.__computeZ,
         'Zin' : self.__computeZin,
         'VSWR': self.__computeVSWR,
         'TL'  : self.__computeTL,
         'LC'  : self.__computeLC}.get(target, 'S')()
        
        return

    def __computeX(self, tags=['freq', 'S11', 'S21', 'S12', 'S22'], func=utils.snp2S, isimag=True):

        for tag in tags:
            for i in range(len(self.__data)):
                if 'pars' not in self.__data[i]: 
                    self.__data[i].setdefault('pars', {})
                    
                self.__data[i]['pars'].update({tag: func(self.__data[i]['data'], tag)})
                if tag != 'freq' and isimag:
                    self.__data[i]['pars'].update({tag+'dB' : utils.toS21(self.__data[i]['pars'][tag].imag, self.__data[i]['pars'][tag].real, isdB=True)})
                    self.__data[i]['pars'].update({tag+'mag': utils.toS21(self.__data[i]['pars'][tag].imag, self.__data[i]['pars'][tag].real, isdB=False)})
                    
        return

    def __computeS(self):
        return self.__computeX()
    
    def __computeY(self):
        return self.__computeX(tags=['freq', 'Y11', 'Y21', 'Y12', 'Y22'], func=utils.snp2Y)
    
    def __computeZ(self):
        return self.__computeX(tags=['freq', 'Z11', 'Z21', 'Z12', 'Z22'], func=utils.snp2Z)

    def __computeZin(self):
        return self.__computeX(tags=['freq', 'Zin1', 'Zin2'], func=utils.snp2Zin)
        
    def __computeVSWR(self):
        return self.__computeX(tags=['freq', 'VSWR1', 'VSWR2'], func=utils.snp2VSWR)

    def __computeTL(self):

        for i in range(len(self.__data)):
            for tag in (tag for tag in ['P1', 'P2'] if tag in self.__data[i]):
                if 'pars' not in self.__data[i]: 
                    self.__data[i].setdefault('pars', {})

                self.__data[i]['pars'].update({'Z0'+tag  : np.sqrt(self.__data[i][tag][:,3]**2+self.__data[i][tag][:,4]**2)})
                self.__data[i]['pars'].update({'Eeff'+tag: np.sqrt(self.__data[i][tag][:,1]**2+self.__data[i][tag][:,2]**2)})

        
        return
    

    def __computeLC(self):
        return self.__computeX(tags=['freq', 'C1', 'C2', 'L1', 'L2',], func=utils.snp2LC, isimag=False)

