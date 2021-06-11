# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:25:23 2020
@institute: BMC∙LMU∙DE
@author: Sunberg
"""

##############################################################
from Bio import SeqIO
import pandas as pd
import numpy as np
from keras.utils import np_utils
import random
from random import shuffle
import re


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

def preprocess_SP(fasta_ref,dats_mc_sp):

    try:
        print('Loading and parsing Fasta file...')
        dict_fasta = {rec.id.split('|')[1]: rec.seq for rec in SeqIO.parse(fasta_ref, "fasta")}

        if len(dict_fasta)>=1:

            df_sp = pd.read_table(dats_mc_sp, encoding='unicode_escape', dtype=str)

            pgs = df_sp['ProteinGroups'].values.tolist()
            upgs = [s.split(';')[0] for s in pgs]
            seqs = df_sp['StrippedPeptide']
            dps = {'PG': upgs, 'SEQ': seqs}
            upgs_seqs = pd.DataFrame(dps).drop_duplicates()
            if upgs_seqs.shape[0] >= 1:
                upgs_seqs['len_seq'] = [len(upgs_seqs.iloc[i, 1]) for i in range(upgs_seqs.shape[0])]
                df1f = upgs_seqs[upgs_seqs['len_seq'] >= 7]

                jj = []
                for i in range(df1f.shape[0]):
                    if df1f.iloc[i, 0] in dict_fasta.keys():
                        try:
                            lc = re.search(df1f.iloc[i, 1], str(dict_fasta[df1f.iloc[i, 0].split(';')[0]])).span(0)[0] + 1
                        except:
                            jj.append(i)
                    else:
                        jj.append(i)

                df1f = df1f.drop(jj)

                seqqs = [df1f.iloc[i, 1] for i in range(df1f.shape[0])]
                PGs = [df1f.iloc[i, 0] for i in range(df1f.shape[0])]
                lastAA = [df1f.iloc[i, 1][len(df1f.iloc[i, 1]) - 1] for i in range(df1f.shape[0])]

                start_pos = [
                    re.search(df1f.iloc[i, 1], str(dict_fasta[df1f.iloc[i, 0].split(';')[0]])).span(0)[0] + 1
                    for i in range(df1f.shape[0])]

                end_pos = [re.search(df1f.iloc[i, 1], str(dict_fasta[df1f.iloc[i, 0].split(';')[0]])).span(0)[1]
                           for i in range(df1f.shape[0])]

                preAA = [str(dict_fasta[df1f.iloc[i, 0].split(';')[0]])[start_pos[i] - 2] for i in range(df1f.shape[0])]

                MC = [len(re.findall(r'K|R', df1f.iloc[i, 1])) - 1 for i in range(df1f.shape[0])]

                dft6_f = {'seq': df1f['SEQ'], 'preAA':preAA, 'LastAA': lastAA, 'start_pos': start_pos,
                          'end_pos': end_pos, 'name': df1f['PG'], 'MC': MC}
                df6_f = pd.DataFrame(dft6_f)
                df6_ff = df6_f.dropna(axis=0, how='any')

                df1_fKR = df6_ff.loc[df6_ff['LastAA'].isin(['K', 'R']) & df6_ff['preAA'].isin(['K', 'R'])]
                df1_fKR = df1_fKR.reindex(df1_fKR.seq.str.len().sort_values().index)
                df1_fKR.index = range(0, df1_fKR.shape[0])

                seqall = df1_fKR['seq'].tolist()

                seqFR = set(seqall)

                min_len = len(seqall[0])
                seqall.reverse()
                word_map = {}

                for s in seqall:
                    for cur_len in range(min_len, len(s)):
                        for i in range(len(s) - cur_len + 1):
                            sub_s = s[i: i + cur_len]
                            if sub_s not in word_map:
                                word_map[sub_s] = set()
                            word_map[sub_s].add(s)

                for s in seqall:
                    if s in word_map:
                        seqFR.remove(s)
                        for parent in word_map[s]:
                            if parent in seqFR:
                                seqFR.remove(parent)

                seqL = [i for i in df1_fKR['seq']]
                seq1tg = [seqL.index(i) for i in seqFR]

                dff_fKR11 = df1_fKR.iloc[seq1tg, :]

                ############################################
                df10_fKR1 = dff_fKR11[dff_fKR11['MC'] == 0]
                df11_fKR1 = dff_fKR11[dff_fKR11['MC'] == 1]
                df12_fKR1 = dff_fKR11[dff_fKR11['MC'] == 2]

                dff_fKR11_seqwin_last0 = [
                    str(dict_fasta[dff_fKR11.iloc[i, 5]])[(dff_fKR11.iloc[i, 4] - 8):(dff_fKR11.iloc[i, 4] + 7)]
                    for i in range(dff_fKR11.shape[0])]

                dff_fKR11_seqwin9_last0 = [dff_fKR11_seqwin_last0[i] for i in range(len(dff_fKR11_seqwin_last0)) if
                                           len(dff_fKR11_seqwin_last0[i]) >= 9]

                dff0_fKR11_seqwin9 = dff_fKR11_seqwin9_last0.copy()

                dff_fKR11_seqwin_pre0 = [
                    str(dict_fasta[dff_fKR11.iloc[i, 5]])[(dff_fKR11.iloc[i, 3] - 9):(dff_fKR11.iloc[i, 3] + 6)]
                    for i in range(dff_fKR11.shape[0])]

                dff_fKR11_seqwin9_pre0 = [dff_fKR11_seqwin_pre0[i] for i in range(len(dff_fKR11_seqwin_pre0)) if
                                          len(dff_fKR11_seqwin_pre0[i]) >= 9]

                dff0_fKR11_seqwin9.extend(dff_fKR11_seqwin9_pre0)
                ############################################
                df111_fKR1_seqwin = [str(dict_fasta[df11_fKR1.iloc[i, 5]])[
                                     ([m.start() for m in re.finditer('K|R', df11_fKR1.iloc[i, 0])][0] + df11_fKR1.iloc[
                                         i, 3] - 8):(
                                             [m.start() for m in re.finditer('K|R', df11_fKR1.iloc[i, 0])][0] + df11_fKR1.iloc[
                                         i, 3] + 7)]
                                     for i in range(df11_fKR1.shape[0])]

                df111_fKR1_seqwin9 = [df111_fKR1_seqwin[i] for i in range(len(df111_fKR1_seqwin)) if
                                      len(df111_fKR1_seqwin[i]) >= 9]

                ################################################
                df121_fKR1_seqwin = [str(dict_fasta[df12_fKR1.iloc[i, 5]])[
                                     ([m.start() for m in re.finditer('K|R', df12_fKR1.iloc[i, 0])][0] + df12_fKR1.iloc[
                                         i, 3] - 8):(
                                             [m.start() for m in re.finditer('K|R', df12_fKR1.iloc[i, 0])][0] + df12_fKR1.iloc[
                                         i, 3] + 7)]
                                     for i in range(df12_fKR1.shape[0])]

                df121_fKR1_seqwin9 = [df121_fKR1_seqwin[i] for i in range(len(df121_fKR1_seqwin)) if
                                      len(df121_fKR1_seqwin[i]) >= 9]

                df122_fKR1_seqwin = [str(dict_fasta[df12_fKR1.iloc[i, 5]])[
                                     ([m.start() for m in re.finditer('K|R', df12_fKR1.iloc[i, 0])][1] + df12_fKR1.iloc[
                                         i, 3] - 8):(
                                             [m.start() for m in re.finditer('K|R', df12_fKR1.iloc[i, 0])][1] + df12_fKR1.iloc[
                                         i, 3] + 7)]
                                     for i in range(df12_fKR1.shape[0])]

                df122_fKR1_seqwin9 = [df122_fKR1_seqwin[i] for i in range(len(df122_fKR1_seqwin)) if
                                      len(df122_fKR1_seqwin[i]) >= 9]

                shuffle(dff0_fKR11_seqwin9)

                df7f1 = df111_fKR1_seqwin9 + df121_fKR1_seqwin9 + df122_fKR1_seqwin9

                df7f0 = dff0_fKR11_seqwin9[0:len(df7f1)]

                dcMC0 = [df7f0[i] for i in range(len(df7f0)) if 'X' not in df7f0[i] and 'U' not in df7f0[i]]
                dcMC1 = [df7f1[i] for i in range(len(df7f1)) if 'X' not in df7f1[i] and 'U' not in df7f1[i]]

                dcMC0u = np.unique(dcMC0)
                dcMC1u = np.unique(dcMC1)

                dcMC01u_cp = [i for i in dcMC0u if i in dcMC1u]

                random.seed(0)
                dcMC0uf = [i for i in dcMC0u if i not in dcMC01u_cp]
                shuffle(dcMC0uf)

                dcMC1uf = [i for i in dcMC1u if i not in dcMC01u_cp]

                mp = min(len(dcMC0uf), len(dcMC1uf))
                dcMC0p = dcMC0uf[0:mp]

                dcMC1p = dcMC1uf[0:mp]

                mctg0 = [0] * mp
                mctg1 = [1] * mp

                dcMC01p = dcMC0p + dcMC1p
                mctg01 = mctg0 + mctg1

                datMC0 = {'sequence': dcMC01p, 'cls': mctg01}

                dfMC = pd.DataFrame(datMC0)

                dfMC.sample(frac=1, random_state=0)
                dfMCs = dfMC.sample(frac=1, random_state=0)

                ions_train = []
                mcs_train = []

                for i in range(dfMCs.shape[0]):
                    if dfMCs.iloc[i, :].empty: break
                    seq = dfMCs.iloc[i, 0]
                    if len(seq) <= 15:
                        mc = dfMCs.iloc[i, 1]
                        mc_labels = np_utils.to_categorical(mc, 2)
                        ion = OnehotEncod(seq)
                        ions_train.append(ion)
                        mcs_train.append(mc_labels)

                    else:
                        pass

                ions_train = np.array(ions_train)
                mcs_train = np.array(mcs_train)
                return ions_train, mcs_train

                print('start_training')

            else:
                print('Please check your Spectronaut file!')
        else:
            print('Please check your fasta file!')
    except:
        print('Please check your datasets!')

