#!/usr/bin/env python


'''
 @file    Z0analyzer.py
 @author  Andrea Giachero <andrea.giachero@mib.infn.it> 
 @date    13 December 2023
 @brief   Script to analyze the transmission line Sonnet simulation ad to extrapolate Z0
          Usage: ./create_KIT.py --cfg path/where/cfgfile/is/cfgfile.json 
'''

from twpazer.sazer import Z0zer, CellZer
from dartwarslab.base import dic2file, file2dic

#from dartwarslab.azer import Z0zer, CellZer
#from dartwarslab.base import dic2file, file2dic

import numpy as np

import sys, os
import argparse

whoami = lambda: sys._getframe(1).f_code.co_name


def Z0analyzer(cfg):

    #cfg = file2dic(cfgfile , writer='json').read()

    for c in cfg:
        print ('({name}) Processing simulation \"{sim}\"'.format(name = whoami(), sim=c))

        cfg[c]['simpath']=os.path.join(cfg[c]['pathfile'], c)
        if not os.path.exists(cfg[c]['simpath']):
            os.makedirs(cfg[c]['simpath'])
            print ('({name}) Path {path} created'.format(name = whoami(), path=cfg[c]['simpath']))
            
        for d in ['savepath', 'saveplot']:
            cfg[c][d]=os.path.join(cfg[c]['simpath'], d)
            if not os.path.exists(cfg[c][d]):
                os.makedirs(cfg[c][d])
                print ('({name}) Path {path} created'.format(name = whoami(), path=cfg[c][d]))
                
        Z0res=Z0zer(os.path.join(cfg[c]['pathfile'], cfg[c]['filenameL']),
                    os.path.join(cfg[c]['pathfile'], cfg[c]['filenameC']),
                    **{k: cfg[c][k] for k in cfg[c] if k not in ['filenameL','filenameC']})
        
        target=Z0res.pars()['target']
        
        import warnings
        warnings.filterwarnings("ignore")
        
        Z0res.compute()

        import pylab as plt
        plt.rcParams.update({"axes.grid" : True})

        Z0res.update(Ldigits=0, Cdigits=0, Zdigits=0)
        for tag in ['Lall' , 'Call', 'Zall']:
            Z0res.plot(tag)
        
        Z0res.update(Ldigits=2, Cdigits=0, Zdigits=0, isfit=True)
        for x in ['L', 'C', 'Z']:
            Z0res.plot(x)

            
        if 'Ztarget' in cfg[c] and cfg[c]['Ztarget'] and isinstance(cfg[c]['Ztarget'], list): 
            res=Z0res.results()


            print('target', target)

            cell=CellZer(target = res['Z0fit'][target+'fit'],
                         L      = res['Z0fit']['L0fit'],             
                         C      = res['Z0fit']['C0fit'],              
                         Z      = np.sqrt(res['Z0fit']['L0fit']/res['Z0fit']['C0fit']),
                         tname  = target,
                         tunit  = Z0res.pars()['tunit'])

            cell.update(Lk = Z0res.pars()['Lk']*1e-12, w = Z0res.pars()['w'])

            restarget={}
            for Z in cfg[c]['Ztarget']:
                cell.add_Z0_target(Z)
                
            
            if target in cfg[c] and cfg[c][target] and isinstance(cfg[c][target], list): 
                for l in cfg[c][target]:
                    cell.add_finger_length_target(l)

            cell.compute()
            cell.dump(os.path.join(os.path.join(cfg[c]['savepath'], Z0res._Z0zer__create_savename()+'.txt')))

            t = cell.targets()
            for Z, yshift in zip(cfg[c]['Ztarget'], [1, 0]):
                if Z in t:
                    l  = t[Z][0]
                    Z0 = t[Z][1]
                    plt.plot(l, Z0, marker='o', ls='none', color='k', zorder=100)
                    plt.axhline(y=Z0, color='black', lw=0.9, ls='--')
                    plt.axvline(x=l,  color='black', lw=0.9, ls='--')

                    xtext=plt.gca().get_xlim()[1]-0.20*(plt.gca().get_xlim()[1]-plt.gca().get_xlim()[0])
                    plt.text(xtext, Z0+4, '$Z_0 = {Z:.1f}\\,\\Omega$'.format(Z=Z0))
                    
                    ytext=plt.gca().get_ylim()[1]-0.25*(plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])
                    plt.text(l+0.50,
                             Z0+4,
                             '{tvar} = {l:0.2f} {tunit}'.format(tvar  = Z0res.pars()['tvar'],
                                                                tunit = Z0res.pars()['tunit'],
                                                                l     = l))
                
            
            for e in ['svg', 'pdf']:
                savename=os.path.join(cfg[c]['saveplot'], 'trend_Z_vs_'+target+'_'+Z0res._Z0zer__create_savename()+'_target.'+e)
                print ('({name}) Saving plot {savename}'.format(name = whoami(), savename=savename))
                plt.gcf().savefig(savename, bbox_inches='tight', transparent=True)
                
            savename=os.path.join(os.path.join(cfg[c]['savepath'], Z0res._Z0zer__create_savename()+'.h5'))
            print ('({name}) Saving results {savename}'.format(name = whoami(), savename=savename))
            dic2file(savename).save({'fit' : Z0res.results()['Z0fit'], 
                                     'pars': Z0res.pars(),
                                     'l'   : {k: v for k, v in Z0res.results().items() if type(k) != str }})
            
            savename=os.path.join(os.path.join(cfg[c]['savepath'], Z0res._Z0zer__create_savename()+'_reduced.h5'))
            print ('({name}) Saving results {savename}'.format(name = whoami(), savename=savename))
            dic2file(savename).save({'L'   : Z0res.results()['Z0fit']['L0'], 
                                     'Lfit': Z0res.results()['Z0fit']['L0fit'],
                                     'C'   : Z0res.results()['Z0fit']['C0'], 
                                     'Cfit': Z0res.results()['Z0fit']['C0fit'],
                                     target+'fit': Z0res.results()['Z0fit'][target+'fit'],
                                     target      : Z0res.results()['Z0fit'][target]})

            
    return



def main():
    parser = argparse.ArgumentParser(description='Extrapolate Z0 from sonnet simulations')

    parser.add_argument("-c", "--cfg" , dest="cfgfile" , type=str  , help="Config file", required = True)
    #parser.add_argument("-s", "--save", dest="savepath", type=str  , help="Savepath"   , required = True)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.error('No argument given!')
        return 

    if args.cfgfile and not os.path.exists(args.cfgfile):
        parser.error('({name}) Filename {cfg} does not exist!'.format(name = whoami(),
                                                                      cfg  = args.cfgfile))
        return

    Z0analyzer(file2dic(args.cfgfile , writer='json').read())     
            
    return

if __name__ == '__main__':
    main()

    
