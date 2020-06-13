import numpy as np
import json
class NumpyEncoder(json.JSONEncoder):    
    """ Special json encoder for numpy types """    
    def default(self, obj):        
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,            
            np.int16, np.int32, np.int64, np.uint8,            
            np.uint16, np.uint32, np.uint64)):            
            return int(obj)        
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):            
            return float(obj)        
        elif isinstance(obj,(np.ndarray,)): #### This is the fix            
            return obj.tolist()        
        return json.JSONEncoder.default(self, obj)

def save_to_json(dict,filename):
    """For storing the dictionary of multiple spheres as json file"""
    finalpath = '../IOfile/json/'+filename+'.json'
    with open(finalpath,'w') as f:
        json.dump(dict, f, cls=NumpyEncoder,indent = 4)
    f = open(finalpath,'r')
    data = json.load(f)
    print(data)
    f.close()
