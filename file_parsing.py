# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:00:27 2022

@author: 2454683t
"""

import pandas as pd 
import glob

def findandimport(ext, name):
    files = glob.glob("**/" + name)
    
    dfs = []
    
    
    for file in files: 
        if file.endswith(ext):
            
            df = pd.read_table(file, sep = "\t", header = 0) 
            
            dfs.append((file, 
                        df))
    return(dfs)
    
plasmids = findandimport(ext = ".txt", name = "*plasmid*")




            
