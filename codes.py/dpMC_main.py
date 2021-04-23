# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:25:23 2020
@institute: BMC∙LMU∙DE
@author: Sunberg
"""

import ipykernel
import os
import argparse
import re
import datetime

#from keras import *
from keras import backend as K
from keras.models import Input, Model
from keras.utils import multi_gpu_model
from keras.callbacks import ModelCheckpoint,EarlyStopping
from keras.layers import Conv1D,Bidirectional,LSTM,Activation,Dropout,MaxPooling1D
from keras.layers import Dense,TimeDistributed,concatenate,Masking,merge
import math



from keras.utils import CustomObjectScope
from keras.models import load_model
from keras.models import model_from_json

from dpMC_Gears import Gears
from dpMC_tcn import TCN


import fnmatch
from Bio import SeqIO
import pandas as pd
import numpy as np
import re

from keras import layers, models, optimizers
from keras import backend as K
from keras.utils import to_categorical
import matplotlib.pyplot as plt
import tensorflow as tf
from keras.layers.core import Flatten

from dpMC_preprocess_MQ import preprocess_MQ
from dpMC_preprocess_SP import preprocess_SP

###########################################################
dp = Gears()
###########################################################
def AAVecDict():
    aa_vec = {}
    s = "0ACDEFGHIKLMNPQRSTVWY"
    v = [0] * len(s)
    v[0] = 1
    for i in range(len(s)):
        aa_vec[s[i]] = list(v)
        v[i], v[(i + 1) % 21] = 0, 1
    return aa_vec


def OnehotEncod(seq):
    one_hot_encod = np.zeros((15, 21), dtype=np.float32)
    for i in range(0, len(seq)):
        one_hot_encod[i, 0:21] = AAVecDict()[seq[i]]
    return one_hot_encod


def mapping_pred(seqwinX):
    aas = []
    for dpl in seqwinX:
        if dpl == "": break
        seq = dpl
        if len(seq) <= 15:
            ion = OnehotEncod(seq)
            aas.append(ion)
        else:
            pass
    ions = np.array(aas)
    return ions

def dpMC_train(ions_train,mcs_train,md_train,wd):

    md = load_model(md_train)
    md.compile(optimizer='adam', loss='categorical_crossentropy',metrics=['accuracy'])
    md_path = wd + "/dpMC--{epoch:03d}-{loss:.5f}-{val_loss:.5f}.h5"  # the path which stores the trained models of each epoch

    checkpoint = ModelCheckpoint(md_path, monitor='val_loss', verbose=1, save_best_only=False,save_weights_only=False, mode='min')
    callbacks_list = [checkpoint]
    md.fit(ions_train, mcs_train, validation_split=0.05, epochs=20, shuffle=True,batch_size=64,verbose=1, callbacks=callbacks_list)

def dpMC_pred(fasta_file,md_pred,fasta_mc_dir):
    try:
        print('Loading and parsing Fasta file...')

        dict_fasta_pred = {rec.id.split('|')[1]: rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

        if len(dict_fasta_pred)>=1:
            pns = dict_fasta_pred.keys()
            seqwins = []
            ls_pos7 = []
            dpn7 = []

            n = 0
            for pn in pns:
                n += 1
                seq = str(dict_fasta_pred[pn])
                seqwin = [seq[(i.start(0) - 7):(i.end(0) + 7)] for i in re.finditer(r'K|R', seq) if i.end(0) > 7]
                ls_pos7m = [i.end(0) for i in re.finditer(r'K|R', seq) if i.end(0) > 7]
                ls_pos7.extend(ls_pos7m)
                dpn7.extend(pn for i in range(len(ls_pos7m)))
                seqwins.extend(seqwin)

            seqwinX = [seqwins[j] for j in range(len(seqwins)) if
                       'U' not in seqwins[j] and 'X' not in seqwins[j] and len(seqwins[j]) > 8]
            dpn7X = [dpn7[j] for j in range(len(seqwins)) if
                     'U' not in seqwins[j] and 'X' not in seqwins[j] and len(seqwins[j]) > 8]
            ls_pos7X = [ls_pos7[j] for j in range(len(seqwins)) if
                        'U' not in seqwins[j] and 'X' not in seqwins[j] and len(seqwins[j]) > 8]

            df_pos7X = pd.DataFrame({'pnm': dpn7X, 'pos': ls_pos7X})

            ls_pos1X = []
            dpn1X = []
            for pn in pns:
                seq = str(dict_fasta_pred[pn])
                ls_pos1m = [i.end(0) for i in re.finditer(r'K|R', seq) if i.end(0) <= 7]
                ls_pos1X.extend(ls_pos1m)
                dpn1X.extend(pn for i in range(len(ls_pos1m)))

            df_pos1X = pd.DataFrame({'pnm': dpn1X, 'pos': ls_pos1X})

            print('Predicting...')

            md = load_model(md_pred)

            ions_pred = mapping_pred(seqwinX)
            pdc_pred = np.round(md.predict(ions_pred),decimals=4)


            td0 = [i[1] for i in pdc_pred]
            td = [1 if i[1] > 0.5 else 0 for i in pdc_pred]

            df_pos7X['td0'] = td0
            df_pos7X['td'] = td
            datfs7 = df_pos7X

            df_pos1X['td'] = 1
            datfs = datfs7.append(df_pos1X,sort=True)
            datfs0 = datfs[datfs['td'] == 0]

            stats_lb = [';'.join([str(datfs.iloc[i, 0]), str(datfs.iloc[i, 1])]) for i in range(datfs.shape[0])]
            stats_val = [str(np.round(datfs.iloc[i, 3], decimals=4)) for i in range(datfs.shape[0])]
            stats = dict(zip(stats_lb, stats_val))

            cut0_seq = {i: np.array(datfs0[datfs0['pnm'] == i].iloc[:, 1]) for i in np.unique(datfs0['pnm'])}

            seqs_mc = []
            pnms_mc = []
            pos_start = []
            pos_end = []
            pos_KRss = []
            stats_KRss = []
            MCs = []

            for i in cut0_seq.keys():
                pos_KR = []
                pos_KRs = []
                stats_KRs = []
                cpos = sorted([0, (len(str(dict_fasta_pred[i])))] + sorted(cut0_seq[i]))
                win = [[cpos[j], cpos[j + 1]] for j in range(len(cpos) - 1)]
                pos_start.extend([(cpos[j] + 1) for j in range(len(cpos) - 1)])
                pos_end.extend([cpos[j + 1] for j in range(len(cpos) - 1)])
                seqs = [str(dict_fasta_pred[i])[0:(ii[1])] if ii[0] == 0 else str(dict_fasta_pred[i])[(ii[0]):(ii[1])] for
                        ii
                        in win]
                MC = [len(re.findall(r'K|R', j[0:(len(j) - 1)])) for j in seqs]
                pos_KR.extend([j.end(0) for j in re.finditer(r'K|R', str(dict_fasta_pred[i]))])

                n = 0
                id = 0
                for t in range(len(MC)):
                    pos_KRs.append(';'.join(str(s) for s in pos_KR[id:(MC[n] + 1 + id)]))
                    id = id + np.sum(MC[n]) + 1
                    n += 1

                for p in range(len(pos_KRs)):
                    stats_tmp = []
                    for ps in pos_KRs[p].split(';'):
                        if len(ps) != 0:
                            if int(ps) <= 7:
                                stats_tmp.append(str(1))
                            elif int(ps) == len(str(dict_fasta_pred[i])):
                                stats_tmp.append(str(0))
                            else:
                                try:
                                    stats_tmp.append(stats[';'.join([str(i), ps])])
                                except:
                                    stats_tmp.append(str(1))
                        else:
                            stats_tmp.append(str(0))

                    stats_KRs.append(';'.join(stats_tmp))

                pos_KRss.extend(pos_KRs)
                stats_KRss.extend(stats_KRs)
                MCs.extend(MC)
                seqs_mc.extend(seqs)
                pnms_mc.extend([i] * len(seqs))

            lens = [len(i) for i in seqs_mc]

            df_mcs = pd.DataFrame(
                {'ID': pnms_mc, 'Sequence': seqs_mc, 'PEPlength': lens, 'Pos_start': pos_start,
                 'Pos_end': pos_end, 'Pos_KR': pos_KRss, 'stats_KR': stats_KRss, 'MC': MCs})

            print('Writing results...')

            fasta_mc_path = fasta_mc_dir + 'dpMC_fasta_digested.txt'
            df_mcs.to_csv(fasta_mc_path, header=True, index=False, sep="\t", mode='w')

            K.clear_session()

        else:
            print('Please check your fasta file!')
    except:
        print('Please check your datasets!')


def main():

    while True:
        options = input('Please select the number of functions in dpMC：1. train; 2.digest; 3.exit\n')


        if options == '1':
            set_wd = input('Please set your working directory:')
            if os.path.exists(set_wd):
                os.chdir(set_wd)
            else:
                print('Please check your directory and try agian!')
                break

            time_now = datetime.datetime.now()
            md_mc_dir = os.getcwd() + '/dpMC/md/' + re.sub(':|\\.', '_', '_'.join(str(time_now).split(' '))) + '/'
            if os.path.exists(md_mc_dir):
                pass
            else:
                os.makedirs(md_mc_dir)

            fasta_ref = input('Please select the reference fasta for training[.fasta]:')
            md_train = input('Please select the pretrained model for fine-tuning[.h5]:')
            options2 = input('Please select the peptides file for training：1. MaxQuant; 2.Spectronaut; 3.exit\n')
            if options2 == '1':
                dats_mc_mq = input('Please select the peptides file of MaxQuant[peptides.txt]:')
                ions_train,mcs_train = preprocess_MQ(fasta_ref,dats_mc_mq)
                dpMC_train(ions_train,mcs_train,md_train,md_mc_dir)
            elif options2 == '2':
                dats_mc_sp = input('Please select the peptides file of Spectronaut[.xls]:')
                ions_train,mcs_train = preprocess_SP(fasta_ref,dats_mc_sp)
                dpMC_train(ions_train,mcs_train,md_train,md_mc_dir)
            elif options2 == '3':
                print('Thanks for using dpMC!')
                break

        elif options == '2':
            set_wd = input('Please set your working directory:')
            if os.path.exists(set_wd):
                os.chdir(set_wd)
            else:
                print('Please check your directory and try agian!')

            fasta_mc_dir = os.getcwd() + '/dpMC/Digested/'
            if os.path.exists(fasta_mc_dir):
                pass
            else:
                os.makedirs(fasta_mc_dir)

            fasta_pred = input('Please select the fasta file for digestion[.fasta]:')
            md_pred = input('Please select the model for digestion[.h5]:')

            dpMC_pred(fasta_pred, md_pred,fasta_mc_dir)

        elif options == '3':
            print('Thanks for using dpMC！')
            break
        else:
            print('Invalid input, please try again!')

if __name__ == '__main__':
    main()
