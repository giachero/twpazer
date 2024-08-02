import numpy as np
import scipy.signal    

def val2idx(val, vals):
    return np.argmin(np.abs(val-vals))

def get_compensation_factor(freq, data, fcenter, deltaf=1):

    frange = freq[val2idx(fcenter-deltaf, freq):val2idx(fcenter+deltaf, freq)]
    drange = data[val2idx(fcenter-deltaf, freq):val2idx(fcenter+deltaf, freq)]

    dataf = scipy.signal.savgol_filter(list(drange), 500, 3)

    newdata = drange-dataf
    idx     = val2idx(frange, fcenter)
    factor  = np.std(newdata[:idx])/np.std(newdata[idx:]) 
    
    return factor

def compensate_data(freq, data, fcenter, deltaf=1):
    factor  = get_compensation_factor(freq, data, fcenter, deltaf=1)
    dataf   = scipy.signal.savgol_filter(list(data), 500, 3)
    newdata = data-dataf

    idx = val2idx(freq, fcenter)
    newdata[idx:] = newdata[idx:]*factor 
    
    return freq, newdata+dataf, scipy.signal.savgol_filter(list(newdata+dataf), 500, 3)


def get_bandwidth(freq, data, threshold):


    idxmax=np.argmax(data)

    newdata=data-(np.max(data)-threshold)

    f1 = freq[val2idx(0, np.abs(newdata[:idxmax]))]
    f2 = freq[val2idx(0, np.abs(newdata[idxmax:]))+idxmax]
    
    return f1, f2


def estimate_ripple(f, data, fcenter, deltaf):

    idxmin = val2idx(fcenter-deltaf, f)
    idxmax = val2idx(fcenter+deltaf, f)

    newdata = data - scipy.signal.savgol_filter(list(data), 500, 3)
    
    ripplestd = np.std(newdata[idxmin:idxmax])
    rippleMm  = np.abs(np.max(newdata[idxmin:idxmax])-np.min(newdata[idxmin:idxmax]))
    
    return ripplestd, rippleMm
