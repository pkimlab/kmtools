# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 15:59:19 2014

@author: alexey
"""

# ipcluster start -n 4

from IPython.parallel import Client
rc = Client()



#%%
with rc[:].sync_imports():
    import numpy as np
    import pandas as pd
    import sqlalchemy as sa




#%%





#%%


