# numpy array to csv

import numpy as np
import math
import random

a = np.asarray([[1,2,3],[4,5,6],[7,8,9]])

num = math.floor((10**6)*np.random.rand())

num_str = str(num)

np.savetxt("Array_saved_as_txt_file_"+num_str+".csv",a,delimiter=",")
