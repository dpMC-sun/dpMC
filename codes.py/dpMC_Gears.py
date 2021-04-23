# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:34:23 2020
@institute: BMC∙LMU∙DE
@author: Sunberg
"""


import os
import fnmatch
import numpy as np
import pandas as pd

import random
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import keras.preprocessing.text as T
from keras.preprocessing.text import Tokenizer
from keras.utils import np_utils

from sklearn.metrics.pairwise import cosine_similarity

class Gears:
    def __init__(self,charge_ions_max = 2,seq_ln_min = 9,seq_ln_max = 15):
        self.charge_ions_max = charge_ions_max
        self.seq_ln_min = seq_ln_min
        self.seq_ln_max = seq_ln_max
        self.aa_vec = self.AAVecDict()
        self.AA_idx = dict(zip("ACDEFGHIKLMNPQRSTVWY",range(0,20)))
        self.AA_ed = "A C D E F G H I K L M N P Q R S T V W Y"


    def AAVecDict(self):
        aa_vec = {}
        s = "0ACDEFGHIKLMNPQRSTVWY"
        v = [0]*len(s)
        v[0] = 1
        for i in range(len(s)):
            aa_vec[s[i]] = list(v)
            v[i],v[(i+1) % 21] = 0,1
        return aa_vec

    def OnehotEncod(self, seq):
        one_hot_encod = np.zeros((self.seq_ln_max,21),dtype=np.float32)
        for i in range(0,len(seq)):
            one_hot_encod[i,0:21] = self.aa_vec[seq[i]]
        return one_hot_encod

    def mapping(self, file):
        with open(file,'r') as dpf:
            aas = []
            mcs = []
            
            for dpl in dpf:
                if dpl == "": break

                dpi = dpl.split("\t")
                seq = dpi[0]
                if len(seq) <= self.seq_ln_max:
                    mc = dpi[1]
                    mc_labels = np_utils.to_categorical(mc,2)
                    
                    ion = self.OnehotEncod(seq)
                    
                    aas.append(ion)
                    mcs.append(mc_labels)
                else:
                    pass

        return aas,mcs
    
    def gettingData(self, filenames):

        
        for filename in filenames:
            print(filename)
            ion,mcs = self.mapping(filename)

        ions = np.array(ion)
        mcs = np.array(mcs)

        print (ions.shape, mcs.shape)
        return ions,mcs


    def load_data(self,foldername):
        print("Loading:")

        filenames = []
        
        for file in os.listdir(foldername):
            if fnmatch.fnmatch(file, "*.dp"):
                filenames.append(os.path.join(foldername,file))
        return self.gettingData(filenames)

    def load_md(self, foldername):
        print("Loading:")

        filenames = ''

        for file in os.listdir(foldername):
            if fnmatch.fnmatch(file, "*_i.h5"):
                filenames = str(os.path.join(foldername, file))
                print(filenames)
        return filenames
    
    def mapping_pred(self, file):
        with open(file,'r') as dpf:
            aas = []

            for dpl in dpf:
                if dpl == "": break

                seq = dpl.split("\n")[0]
                if len(seq) <= self.seq_ln_max:
                    ion = self.OnehotEncod(seq)
                    
                    aas.append(ion)

                else:
                    pass
        return aas
    

if __name__ == '__main__':
    Gears()
