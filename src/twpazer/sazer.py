'''
 @file    sazer.py
 @author  Andrea Giachero <andrea.giachero@mib.infn.it> 
 @date    21 November 2023
 @brief   Collection of base classes to analyze Sonnet data. 
          sazer: Sonnet Analyzer    

'''

import matplotlib.cm as cm
import matplotlib.colors as mcolors

import numpy as np
import pylab as plt
import os, errno, sys

from dartwarslab.base import sreader



whoami = lambda: sys._getframe(1).f_code.co_name
Z0val  = lambda L, C : np.sqrt(L/C)
getpar = lambda l: [v for v in l if type(v) != str] 


class Z0zer(object):
    '''
    Class to analyze the KIT transmission line simulations

    '''
    def __init__(self, filenameL, filenaleC, **kwargs):
        self.__pars={'npol'     : 3,
                     'ncell'    : 1,
                     'nfit'     : 200,
                     'filenameL': None,
                     'filenameC': None,
                     'savepath' : None,
                     'Lnorm'    : 1,
                     'Cnorm'    : 1,
                     'Znorm'    : 1,
                     'Lunit'    : 'H',
                     'Cunit'    : 'C',
                     'Zunit'    : '$\\Omega$',
                     'barunit'  : '$\mu$m',
                     'barvar'   : '$\ell$',
                     'barname'  : 'Finger Length',
                     'Ldigits'  : 3,
                     'Cdigits'  : 1,
                     'Zdigits'  : 1,
                     'xdigits'  : 0,
                     'issave'   : False,
                     'isfit'    : True,
                     'fittypeC' : 'poly',
                     'fittypeL' : 'poly'} 

        for filename, tag in zip([filenameL, filenaleC],['filenameL', 'filenameC']):
            if not os.path.exists(filename):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
            self.__pars.update({tag:filename})
            
        # Input parameters
        self.__pars.update(kwargs)

        # Results data
        self.__data=dict()
        
        return


    def update(self, **kwargs):
        self.__pars.update(kwargs)
        return

    def pars(self):
        return self.__pars


    def results(self):
        return self.__data

            
    def __computeZ0(self):


        print ('({name}) Compute Z0'.format(name = whoami()))
        
        
        for l in getpar(self.__data):    
            Z = Z0val(self.__data[l]['L']['L'], self.__data[l]['C']['C'])
            self.__data[l].setdefault('Z',{})
            self.__data[l]['Z'].update({'freq': self.__data[l]['L']['freq']})
            self.__data[l]['Z'].update({'Z'   : Z})
            
            f     = self.__data[l]['Z']['freq']
            p     = np.poly1d(np.polyfit(f , Z, self.__pars['npol']))
            f_fit = np.linspace(np.min(f), np.max(f), self.__pars['nfit']);

            self.__data[l]['Z'].update({'Zfit'    : p(f_fit)})
            self.__data[l]['Z'].update({'freq_fit': f_fit})
            self.__data[l]['Z'].update({'p': p.coef})
            self.__data[l]['Z'].update({'Z0': p(0)})
        
        return

    def __computeLC(self):
        
        filename={'C': self.__pars['filenameC'],
                  'L': self.__pars['filenameL']}

        for tag in ['L', 'C']:

            print ('({name}) Compute {tag}'.format(name = whoami(), tag=tag))
            
            sid  = sreader(data = sreader(filename.get(tag,'C')).unique('merge'),
                           fmin = self.__pars['fmin'] if 'fmin' in self.__pars else None,
                           fmax = self.__pars['fmax'] if 'fmax' in self.__pars else None)
            
            sid.compute('LC')
            
            for s in sid:
                f = s['pars']['freq']
                X = s['pars'][tag+'1']
                l = s[self.__pars['target']]
                
                self.__computeX(f, X, l, tag)
                
        return

    def __computeX(self, f, X, l, tag):
        
        #idxmax = np.argmin(abs(f-self.__pars['fmax'])) if 'fmax' in self.__pars else None 
        #idxmin = np.argmin(abs(f-self.__pars['fmin'])) if 'fmin' in self.__pars else None 

        p      = np.poly1d(np.polyfit(f, X/self.__pars['ncell'], self.__pars['npol']))
        f_fit  = np.linspace(np.min(f),np.max(f), self.__pars['nfit']);
        
        if l not in self.__data:
            self.__data.setdefault(l, {})

        self.__data[l].setdefault(tag,{})
        self.__data[l][tag].update({tag       : X/self.__pars['ncell'],
                                    tag+'fit' : p(f_fit),
                                    tag+'0'   : p(0),
                                    'freq'    : f,
                                    'freq_fit': f_fit,
                                    'p'       : p.coef})
        
        return



    def __fitLC(self):


        for tag in ['L', 'C']:
            print ('({name}) Fitting {tag}0'.format(name = whoami(), tag=tag))
            self.__fitX(tag)
            
        return

    def __fitX(self, tag):

        ls = getpar(self.__data)
        if 'exclude' in self.__pars:
            ls = [l for l in ls if l not in self.__pars['exclude']]
        
        ls_fit = np.linspace(np.min(ls),np.max(ls), self.__pars['nfit']);
        X0=[self.__data[l][tag][tag+'0'] for l in ls]
                

        if 'fittype'+tag not in self.__pars or self.__pars['fittype'+tag]=='poly': 
            print('({name}) Performing linear fit'.format(name = whoami()))
            
            p      = np.poly1d(np.polyfit(ls, X0, 2))
            X0_fit = p(ls_fit)
            X00    = p(0) 
            
            print('({name}) {tag}0 from polynomial fit : {tag} = {X:.2f} {u}'.format(name = whoami(),
                                                                                     tag  = tag,
                                                                                     X    = X00*self.__pars[tag+'norm'],
                                                                                     u    = self.__pars[tag+'unit']))

        if 'Z0fit' not in self.__data: 
            self.__data.setdefault('Z0fit', {})


        self.__data['Z0fit'].update({self.__pars['target']      : np.array(ls),
                                     self.__pars['target']+'fit': np.array(ls_fit),
                                     tag+'0'    : np.array(X0),
                                     tag+'00'   : np.array(X00),
                                     tag+'0fit' : np.array(X0_fit),
                                     'd'+tag+'0': [np.abs(np.max(X0)-np.min(X0)), np.abs(np.max(X0)-np.min(X0))/np.max(X0)],
                                     'p'+tag    : p.coef})
        
        return
        

    def __fitZ0(self):
        
        self.__data['Z0fit'].update({'Z0'   : Z0val(self.__data['Z0fit']['L0'],    self.__data['Z0fit']['C0']),
                                     'Z0fit': Z0val(self.__data['Z0fit']['L0fit'], self.__data['Z0fit']['C0fit'])})
        
        return

    
    def compute(self):
        self.__computeLC()
        self.__computeZ0()

        self.__fitLC()
        self.__fitZ0()

        return 
    



    def plot(self, tag='Z'):
        return {'Z' : self.__plotX,
                'L' : self.__plotX,
                'C' : self.__plotX,
                'Zall' : self.__plot_allX,
                'Lall' : self.__plot_allX,
                'Call' : self.__plot_allX,}.get(tag, 'Z')(tag)
    

    def __plot_allX(self, tag='Z0all'):

        lut={'Zall': 'Z',
             'Lall': 'L',
             'Call': 'C'}

        norm = [p for p in [p for p in self.__pars if p in ['Lnorm', 'Cnorm', 'Znorm']] if lut[tag] in p][0]
        unit = [p for p in [p for p in self.__pars if p in ['Lunit', 'Cunit', 'Zunit']] if lut[tag] in p][0]


        plt.figure()
        ws=[w for w in self.__data if type(w) != str]
        sm, cnorm, cmap = create_colormap(ws)
        for w in ws:
            plt.plot(self.__data[w][lut[tag]]['freq']/1e6,
                     self.__data[w][lut[tag]][lut[tag]]*self.__pars[norm], color=cmap(cnorm(w)),
                     marker='o', ls='none', zorder=100)
            plt.plot(self.__data[w][lut[tag]]['freq_fit']/1e6,
                     self.__data[w][lut[tag]][lut[tag]+'fit']*self.__pars[norm], color=cmap(cnorm(w))) 

        plt.title(self.__create_title())
        plt.ylabel('${s}$ [{u}]'.format(s=lut[tag], u=self.__pars[unit]))
        plt.xlabel('Frequency $f$ [MHz]')

        #create_cbar(sm, 1, r'Finger Length $\ell$ [$\mu$m]')
        create_cbar(sm, 1, r'{barname} {barvar} [{barunit}]'.format(barname = self.__pars['barname'],
                                                                    barvar  = self.__pars['barvar'],
                                                                    barunit = self.__pars['barunit']))
        
        self.__save_plot(self.__create_savename('scan_{X}_vs_f_vs_{w}_{post}'.format(X=lut[tag],
                                                                                     w=self.__pars['target'],
                                                                                     post='')))
        return
    
    
    def __plotX(self, tag='L'):
        
        norm = [p for p in [p for p in self.__pars if p in ['Lnorm', 'Cnorm', 'Znorm']] if tag in p][0]
        unit = [p for p in [p for p in self.__pars if p in ['Lunit', 'Cunit', 'Zunit']] if tag in p][0]
        
        plt.figure()
        if self.__pars['isfit']:
            plt.plot(self.__data['Z0fit'][self.__pars['target']],       self.__data['Z0fit'][tag+'0']*self.__pars[norm], marker='o', ls='none', zorder=100)
            plt.plot(self.__data['Z0fit'][self.__pars['target']+'fit'], self.__data['Z0fit'][tag+'0fit']*self.__pars[norm])
        else:
            plt.plot(self.__data['Z0fit'][self.__pars['target']], self.__data['Z0fit'][tag+'0']*self.__pars[norm], marker='o')
            
        plt.title(self.__create_title())
        plt.ylabel('${s} [${u}$]$'.format(s=tag, u=self.__pars[unit]))
        plt.xlabel(r'{barname} {barvar} [{barunit}]'.format(barname = self.__pars['barname'],
                                                            barvar  = self.__pars['barvar'],
                                                            barunit = self.__pars['barunit']))

        plt.ticklabel_format(useOffset=False, style='plain', axis='y')
        #plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
        
        self.__save_plot(self.__create_savename('trend_{X}_vs_{w}_'.format(X=tag,
                                                                           w=self.__pars['target'])))
        
        
        return


    def __create_title(self, titlestr=''):

        varconverter={'Lk':'L_k', 't' :'t', 'w' :'w', 's' :'s','d' :'d', 'eps': '\\varepsilon_r'}
        
        for l, u in zip(['Lk', 't', 'w', 's',  'd', 'eps'], ['pH/sq', 'nm', '$\\mu$m', '$\\mu$m',  'nm', '']):
            if l in self.__pars:
                titlestr+=r'${l}={v}$ {u} , '.format(l=varconverter[l], v=self.__pars[l], u=u)

        return titlestr[:-3]

    def __create_savename(self, savename=''):
        return savename+'_'.join([k+str(self.__pars[k]) for k in ['Lk', 't', 'w', 's', 'd', 'eps'] if k in self.__pars])

    def __save_plot(self, savename):
        if self.__pars['issave'] and 'saveplot' in self.__pars:
            for e in ['svg', 'pdf']:
                s=os.path.join(self.__pars['saveplot'], savename+'.{e}'.format(e=e))
                plt.gcf().savefig(s, bbox_inches='tight', transparent=True)
                print ('({name}) Saving plot {savename}'.format(name = whoami(), savename=s))
        
        return




class CellZer(object):
    '''

    Class to analyze the resuls from the Z0zer class

    '''

    def __init__(self, l, L, C, Z, **kwargs):

        self.__data=dict()
        self.__pars=dict()
        for t, v in zip(['l','L', 'C', 'Z',], [l, L, C, Z]):
            self.__data.update({t:v})


        self.__labels = ['Z0target', 'l'   , 'Z0'   ,'L'   , 'C'   ]
        self.__units  = ['[Ohm]'   , '[um]', '[Ohm]','[pH]', '[fF]']
        self.__norms  = [ 1        ,  1    ,  1,      1e12 ,  1e15 ] 

            
        self.__pars.update(kwargs)
        self.compute()
        self.__targets=dict()

        return

    def update(self, **kwargs):
        self.__pars.update(kwargs)
        return
    
    def pars(self):
        return self.__pars

    def targets(self):
        return self.__targets

    def results(self):
        return self.__data
    
    def add_Z0_target(self, Z0target):

        idx = np.argmin(np.abs(Z0target-self.__data['Z']))
        t   = [self.__data['l'][idx], self.__data['L'][idx], self.__data['C'][idx]]

        t.insert(1, np.sqrt(t[1]/t[2]))
                
        self.__targets.setdefault(Z0target, t)

        return

    def add_finger_length_target(self, ltarget):

        idx = np.argmin(np.abs(ltarget-self.__data['l']))
        Z0 = np.sqrt(self.__data['L'][idx]/self.__data['C'][idx])
        self.add_Z0_target(Z0)
        
        return 
    
    def compute(self):
        if 'w' in self.__pars and 'nsq' not in self.__labels:

            for l, v in zip([self.__labels, self.__units, self.__norms],['nsq', '[squares]', 1]):
                l.append(v)

            for t in self.__targets:
                self.__targets[t].append(self.__targets[t][0]/self.__pars['w'])

        if all(x in self.__pars for x in ['w', 'Lk']) and not all(x in self.__pars for x in ['nsq', 'Lf']): 
            for l, v in zip([self.__labels, self.__units, self.__norms],['Lk', '[pH]', 1e12]):
                l.append(v)
                
            for t in self.__targets:
                self.__targets[t].append(self.__pars['Lk'])

            for l, v in zip([self.__labels, self.__units, self.__norms],['Lf', '[nH]', 1e9]):
                l.append(v)
                
            for t in self.__targets:
                self.__targets[t].append(self.__pars['Lk']*self.__targets[t][-2])

            for l, v in zip([self.__labels, self.__units, self.__norms],['Cf', '[fF]', 1e15]):
                l.append(v)
            
            for t in self.__targets:
                self.__targets[t].append(self.__targets[t][3]/2)
                
        return 
    
    
    def dump(self, filename=None):
        
        formatter=''.join(['%-12s']*len(self.__labels))
        
        print(formatter%tuple(self.__labels))
        print(formatter%tuple(self.__units))
        
        for t in self.__targets:
            restarget=[a*b for a,b in zip([t]+self.__targets[t], self.__norms)]
            print (formatter%tuple(["%.2f" % i for i in restarget]))
            
        if filename:
            with open(filename, "w") as f:
                f.write(formatter%tuple(self.__labels)+'\n')
                f.write(formatter%tuple(self.__units)+'\n')
                
                for t in self.__targets:
                    restarget=[a*b for a,b in zip([t]+self.__targets[t], self.__norms)]
                    f.write(formatter%tuple(["%.2f" % i for i in restarget])+'\n')
                
        return 
    

    

def create_colormap(vals):
        
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    from matplotlib import gridspec
    
    cmap = cm.turbo
    norm = mcolors.Normalize(vmin=np.min(vals), vmax=np.max(vals))
    
    sm = cm.ScalarMappable(norm, cmap)
    sm.set_array([round(v,2) for v in vals])
    
    return sm, norm, cmap


def create_cbar(sm, n, label='Heat Map'):
    import matplotlib.ticker
    import warnings
    warnings.filterwarnings("ignore")

    
    xlabel_format = '{{:.{n}f}}'.format(n=n)
    
    cbar=plt.gcf().colorbar(sm, ax = plt.gca())
    cbar.set_label(label)

    fmt = matplotlib.ticker.StrMethodFormatter("{x}")
    cbar.ax.yaxis.set_major_formatter(fmt)    
    cbar.set_ticklabels([xlabel_format.format(t) for t in cbar.ax.get_yticks()])

    warnings.filterwarnings("default")
    
    return
