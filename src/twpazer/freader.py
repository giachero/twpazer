'''
 @file    freader.py
 @author  Andrea Giachero <andrea.giachero@mib.infn.it> 
 @date    17 October 2024
 @brief   Collection of base classes to handle json and hdf5 files

 @updates 2024-10-17: first commit

'''


'''
From dictionary to hdf5 or json file

'''

import json, h5py
import numpy as np

class file2dic(object):
    def __init__(self, savename, writer='h5'):

        self.__savename = savename
        self.__writer   = writer
        
        return

    def read(self, writer='h5'):
        return {'h5'  : self.__hdf5_read,
                'json': self.__json_read}.get(self.__writer, 'h5')()

    def __hdf5_read(self):
        return h2dic(self.__savename).read()
    
    def __json_read(self):

        import json
        with open(self.__savename, "r") as f: 
            return json.load(f)
        
        return

'''
From hdf5 or json file to dictionary

'''
class dic2file(object):

    def __init__(self, savename, writer='h5'):

        self.__savename = savename
        self.__writer   = writer
        
        return

    def save(self, d):
        return {'h5'  : self.__hdf5_save,
                'json': self.__json_save}.get(self.__writer, 'h5')(d)

    
    def __hdf5_save(self, d):
        return dic2h(self.__savename).write(d)
    

    def __json_save(self, d):

        import json
        
        class NpEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                if isinstance(obj, np.floating):
                    return float(obj)
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)
            
            
        with open(self.__savename, 'w', encoding='utf-8') as fp:
            json.dump(d, fp, cls=NpEncoder, indent=4)
                      
        return
    
class dic2h(object):
    def __init__(self, filename):

        self.__filename=filename

        return

    def write(self, d):
        if not isinstance(d, dict):
            return

        with h5py.File(self.__filename, 'w') as f:
            self.__recursively_save_dict_contents_to_group(f,  d)
            
        return

    def __recursively_save_dict_contents_to_group(self, f, d):

        for key, item in d.items():
            
            key = str(key)
            if isinstance(item, list):
                item = np.array(item)

            if isinstance(item, (str, float, int, np.ndarray, np.int64, np.float64,  np.float32)):
                if isinstance (item, str):
                    dtype=h5py.special_dtype(vlen=str)
                elif isinstance (item, np.ndarray):
                    dtype=item.dtype
                else:
                    dtype=type(item)
                    
                f.create_dataset(key, data=item, dtype=dtype)

            elif isinstance(item, dict):
                f.create_group(key)
                self.__recursively_save_dict_contents_to_group(f[key], item)
                
        return
    

class h2dic(object):
    def __init__(self, filename):
        
        self.__filename=filename
        
        return
    
    def read(self):

        with h5py.File(self.__filename, 'r') as f:
            return self.__recursively_load_dict_contents_from_group(f)
        
        return

    def __recursively_load_dict_contents_from_group(self, f):
        ans = dict()
        for key, item in f.items():
            if isinstance(item, h5py._hl.dataset.Dataset):
                ans[key] = item[()].decode('utf-8') if isinstance (item[()], bytes) else item[()]
            elif isinstance(item, h5py._hl.group.Group):
                ans[key] = self.__recursively_load_dict_contents_from_group(f[key])
                
        return ans
