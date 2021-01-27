# -*- coding: utf-8 -*-
"""
Created on Wed May  6 17:10:54 2020

@author: nguyentd
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
import matplotlib

# fs = 16
fs = 12
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fs }

sns.set(font_scale=2.2)
sns.set_style("whitegrid", {'ytick.left': True })
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex = False)
plt.rc('text', usetex=False)

current_palette = sns.color_palette("viridis_r",5)
palette = sns.color_palette()[0:5]

plt.close('all')

dpi=72
fig_w = 1600/dpi
fig_h = fig_w*2/3

plt.close('all')


# generate some data
x = np.arange(0, 10, 0.2)
y = np.sin(x)

# plot it
f, (a0, a1, a2 ) = plt.subplots(3, 3,  figsize=(fig_w, fig_h),
                                gridspec_kw={'height_ratios': [2.5,3, 1]},
                                tight_layout=True
                                )

# mng = plt.get_current_fig_manager()
# mng.frame.Maximize(True)

# figManager = plt.get_current_fig_manager()
# figManager.window.showMaximized()

#PSU
A4_PSU = pd.read_csv('analysis/data/raw/figure1/psu/A4_PSU_mu_0p001983_20200229.csv', sep=',');
A4_PSU= A4_PSU[(A4_PSU.PFPR==5) & (A4_PSU.TREATMENT_COVERAGE==0.4)]

medians = A4_PSU.groupby(['STARTING_PARTNER_DRUG_FREQ'])[['TIME_TO_1p_580Y','TIME_TO_10p_580Y','TIME_TO_25p_580Y','TIME_TO_10p_TF','TIME_TO_25p_TF']].median()*12

space = 0.3
text_offset_y1 = 0.1
text_offset_y2 = 0.23

text_offset=20
text_offset_x1=20

scale = 1
a_size = 12
marker_size = 100

for i in range(len(medians.index.values)):
    
    freq = medians.index.values[i]
    
    lines = [[(medians.TIME_TO_1p_580Y[freq], i*scale),(medians.TIME_TO_25p_580Y[freq],i*scale) ]]

    lines_tf = [[(medians.TIME_TO_10p_TF[freq], i*scale-space),(medians.TIME_TO_25p_TF[freq],i*scale-space) ]]

    lc = LineCollection(lines, colors=['k'])
    lc_tf = LineCollection(lines_tf, colors=['k'])
    
    a0[0].add_collection(lc)
    a0[0].add_collection(lc_tf)
    
    
    a0[0].scatter([medians.TIME_TO_1p_580Y[freq]],[i*scale], s=[marker_size], 
                  marker='o',c='k')
    a0[0].scatter([medians.TIME_TO_10p_580Y[freq],medians.TIME_TO_25p_580Y[freq]], 
                  [i*scale,i*scale] , s=[marker_size,marker_size], c='k',marker='s') 
    
    
    a0[0].scatter([medians.TIME_TO_10p_TF[freq],medians.TIME_TO_25p_TF[freq]], [i*scale-space,i*scale-space], 
                  s= [marker_size, marker_size], c='r',marker='s',
                  edgecolor='black', linewidth=1.5)

    if i==4:
        a0[0].annotate('10% TF', (medians.TIME_TO_10p_TF[freq]-16,i*scale-space-text_offset_y2),
                       size=a_size, c='r', weight='bold')
        a0[0].annotate('25% TF', (medians.TIME_TO_25p_TF[freq]-text_offset,i*scale-space-text_offset_y2),
                       size=a_size, c='r', weight='bold')
        a0[0].annotate('0.01', (medians.TIME_TO_1p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
                       size=a_size, weight='bold')
        a0[0].annotate('0.10', (medians.TIME_TO_10p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
                       size=a_size, weight='bold')
        a0[0].annotate('0.25 580Y', (medians.TIME_TO_25p_580Y[freq]+10,i*scale+text_offset_y1),
                       size=a_size, weight='bold')
    
# [medians.TIME_TO_1p_580Y[0],medians.TIME_TO_10p_580Y[0],medians.TIME_TO_25p_580Y[0]],[0,0,0]


start_time = 0
total_time= 30*12+1

a0[0].set_ylim(-0.5, scale*4.5)

a0[0].set_xlim(start_time, total_time)
a0[0].set_xticks(np.arange(start_time, total_time, 12*5.0))
a0[0].set_xticklabels([])

a0[0].set_yticks(np.arange(0, 5,1)-space/2)
a0[0].set_yticklabels(['0.00', '0.01', '0.10', '0.25', '0.5'])
a0[0].set_ylabel('Pre-existing partner\n drug resistance')

a0[0].grid(axis='y')
a0[0].spines['top'].set_visible(False)
a0[0].spines['right'].set_visible(False)
a0[0].spines['left'].set_visible(False)
a0[0].spines['bottom'].set_visible(False)

a0[0].set_title('PSU CIDD')

# plt.subplot_adjust(hspace=0.100)

#imperial Plot 1 

A4_IMPERIAL = pd.read_csv('analysis/data/raw/figure1/imperial/A4_IMPERIAL_magenta_1.3.0_2020022.csv', sep=',');
A4_IMPERIAL= A4_IMPERIAL[(A4_IMPERIAL.PFPR==5) & (A4_IMPERIAL.TREATMENT_COVERAGE==0.4)]

medians = A4_IMPERIAL.groupby(['STARTING_PARTNER_DRUG_FREQ'])[['TIME_TO_1p_580Y','TIME_TO_10p_580Y','TIME_TO_25p_580Y','TIME_TO_10p_TF','TIME_TO_25p_TF']].median()*12


for i in range(len(medians.index.values)):
    
    freq = medians.index.values[i]
    
    lines = [[(medians.TIME_TO_1p_580Y[freq], i*scale),(medians.TIME_TO_25p_580Y[freq],i*scale) ]]

    lines_tf = [[(medians.TIME_TO_10p_TF[freq], i*scale-space),(medians.TIME_TO_25p_TF[freq],i*scale-space) ]]

    lc = LineCollection(lines, colors=['k'])
    lc_tf = LineCollection(lines_tf, colors=['k'])
    
    a0[1].add_collection(lc)
    a0[1].add_collection(lc_tf)
    
    
    a0[1].scatter([medians.TIME_TO_1p_580Y[freq]],[i*scale], s=[marker_size], 
                  marker='o',c='k')
    a0[1].scatter([medians.TIME_TO_10p_580Y[freq],medians.TIME_TO_25p_580Y[freq]], 
                  [i*scale,i*scale] , s=[marker_size,marker_size], c='k',marker='s') 
    
    
    a0[1].scatter([medians.TIME_TO_10p_TF[freq],medians.TIME_TO_25p_TF[freq]], [i*scale-space,i*scale-space], 
                  s= [marker_size, marker_size], c='r',marker='s',
                  edgecolor='black', linewidth=1.5)
    # if i==4:
    #     a0[1].annotate('10%TF', (medians.TIME_TO_10p_TF[freq]-text_offset,i*scale-space-text_offset_y2),
    #                    size=a_size, c='r', weight='bold')
    #     a0[1].annotate('25% TF', (medians.TIME_TO_25p_TF[freq]-text_offset,i*scale-space-text_offset_y2),
    #                    size=a_size, c='r', weight='bold')
    #     a0[1].annotate('0.01', (medians.TIME_TO_1p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')
    #     a0[1].annotate('0.10', (medians.TIME_TO_10p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')
    #     a0[1].annotate('0.25', (medians.TIME_TO_25p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')

start_time = 0
total_time= 30*12+1

a0[1].set_ylim(-0.5, scale*4.5)

a0[1].set_xlim(start_time, total_time)
a0[1].set_xticks(np.arange(start_time, total_time, 12*5.0))
a0[1].set_xticklabels([])

a0[1].set_yticks(np.arange(0, 5,1)-space/2)
a0[1].set_yticklabels([])
# a0[1].set_ylabel('pre-existing res\n to partner drug')

a0[1].grid(axis='y')
a0[1].spines['top'].set_visible(False)
a0[1].spines['right'].set_visible(False)
a0[1].spines['left'].set_visible(False)
a0[1].spines['bottom'].set_visible(False)
a0[1].set_title('IMPERIAL')

#####
#MORU plot1

A4_MORU = pd.read_csv('analysis/data/raw/figure1/moru/A4_MORU_rev17.csv', sep=',');
A4_MORU= A4_MORU[(A4_MORU.PFPR==5) & (A4_MORU.TREATMENT_COVERAGE==0.4)]

medians = A4_MORU.groupby(['STARTING_PARTNER_DRUG_FREQ'])[['TIME_TO_1p_580Y','TIME_TO_10p_580Y','TIME_TO_25p_580Y','TIME_TO_10p_TF','TIME_TO_25p_TF']].median()*12


for i in range(len(medians.index.values)):
    
    freq = medians.index.values[i]
    
    lines = [[(medians.TIME_TO_1p_580Y[freq], i*scale),(medians.TIME_TO_25p_580Y[freq],i*scale) ]]

    lines_tf = [[(medians.TIME_TO_10p_TF[freq], i*scale-space),(medians.TIME_TO_25p_TF[freq],i*scale-space) ]]

    lc = LineCollection(lines, colors=['k'])
    lc_tf = LineCollection(lines_tf, colors=['k'])
    
    a0[2].add_collection(lc)
    a0[2].add_collection(lc_tf)
    
    
    a0[2].scatter([medians.TIME_TO_1p_580Y[freq]],[i*scale], s=[marker_size], 
                  marker='o',c='k')
    a0[2].scatter([medians.TIME_TO_10p_580Y[freq],medians.TIME_TO_25p_580Y[freq]], 
                  [i*scale,i*scale] , s=[marker_size,marker_size], c='k',marker='s') 
    
    
    a0[2].scatter([medians.TIME_TO_10p_TF[freq],medians.TIME_TO_25p_TF[freq]], [i*scale-space,i*scale-space], 
                  s= [marker_size, marker_size], c='r',marker='s',
                  edgecolor='black', linewidth=1.5)
    # if i==4:
    #     a0[1].annotate('10%TF', (medians.TIME_TO_10p_TF[freq]-text_offset,i*scale-space-text_offset_y2),
    #                    size=a_size, c='r', weight='bold')
    #     a0[1].annotate('25% TF', (medians.TIME_TO_25p_TF[freq]-text_offset,i*scale-space-text_offset_y2),
    #                    size=a_size, c='r', weight='bold')
    #     a0[1].annotate('0.01', (medians.TIME_TO_1p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')
    #     a0[1].annotate('0.10', (medians.TIME_TO_10p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')
    #     a0[1].annotate('0.25', (medians.TIME_TO_25p_580Y[freq]-text_offset_x1,i*scale+text_offset_y1),
    #                    size=a_size, weight='bold')

start_time = 0
total_time= 30*12+1

a0[2].set_ylim(-0.5, scale*4.5)

a0[2].set_xlim(start_time, total_time)
a0[2].set_xticks(np.arange(start_time, total_time, 12*5.0))
a0[2].set_xticklabels([])

a0[2].set_yticks(np.arange(0, 5,1)-space/2)
a0[2].set_yticklabels([])

a0[2].grid(axis='y')
a0[2].spines['top'].set_visible(False)
a0[2].spines['right'].set_visible(False)
a0[2].spines['left'].set_visible(False)
a0[2].spines['bottom'].set_visible(False)
a0[2].set_title('MORU')

#%%
#PSU plot 2
total_time = 649

# fig, axes = plt.subplots(1,3)
x=range(total_time)

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'psu'

for i in range(len(init_freqs)):
    freq = init_freqs[i]
    filename = 'analysis/data/raw/figure1/%s/A4_mu_0p001983_plas2_%s_comp_1p0_tc_0p4_pfpr_PFPR05_C580Y.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)

    a1[0].plot(quantile_init0p00[0.5],lw=2,color=current_palette[i])
    a1[0].fill_between(x, quantile_init0p00[0.75], quantile_init0p00[0.25],color=current_palette[i], alpha=0.1)
    

a1[0].set_ylim(0, 1)
start_time = 168
total_time=start_time + 30*12
a1[0].set_xlim(start_time, total_time)
a1[0].set_xticks(np.arange(start_time, total_time+1, 12*5.0))
# a1.set_xticklabels(np.arange(0, 61, 5))
a1[0].set_xticklabels([])

# a1.set_ylabel('580Y Frequency',fontsize=fs)
a1[0].set_ylabel('580Y Frequency')
# a1.set_title('[A4] TC=0.4  PFPR5 PSU',fontsize=18, fontweight='bold')

############

#Imperial plot 2
total_time = 481

# fig, axes = plt.subplots(1,3)
x=range(total_time)

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'imperial'

for i in range(len(init_freqs)):
    freq = init_freqs[i]
    filename = 'analysis/data/raw/figure1/%s/A4_init_res_freq_%s_tc_0p4_pfpr_PFPR05_C580Y_Imperial.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)

    a1[1].plot(quantile_init0p00[0.5],lw=2,color=current_palette[i])
    a1[1].fill_between(x, quantile_init0p00[0.75], quantile_init0p00[0.25],color=current_palette[i], alpha=0.1)
    

a1[1].set_ylim(0, 1)
start_time = 0
total_time=30*12
a1[1].set_xlim(start_time, total_time)
a1[1].set_xticks(np.arange(start_time, total_time, 12*5.0))
# a1.set_xticklabels(np.arange(0, 61, 5))
a1[1].set_xticklabels([])
a1[1].set_yticklabels([])

# a1.set_ylabel('580Y Frequency',fontsize=fs)
# a1[1].set_ylabel('580Y Frequency')
############

#MORU plot 2
total_time = 481

# fig, axes = plt.subplots(1,3)
x=range(total_time)

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'moru'

for i in range(len(init_freqs)):
    freq = init_freqs[i]
    filename = 'analysis/data/raw/figure1/%s/A4_init_res_freq_%s_tc_0p4_pfpr_PFPR05_C580Y_MORU.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)

    a1[2].plot(quantile_init0p00[0.5],lw=2,color=current_palette[i])
    a1[2].fill_between(x, quantile_init0p00[0.75], quantile_init0p00[0.25],color=current_palette[i], alpha=0.1)
    

a1[2].set_ylim(0, 1)
start_time = 0
total_time=30*12
a1[2].set_xlim(start_time, total_time)
a1[2].set_xticks(np.arange(start_time, total_time, 12*5.0))
# a1.set_xticklabels(np.arange(0, 61, 5))
a1[2].set_xticklabels([])
a1[2].set_yticklabels([])


#%%
#PSU plot 3

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'psu'

freq_criteria = 0.1


for i in range(len(init_freqs)):
    freq = init_freqs[i]
    
    filename = 'analysis/data/raw/figure1/%s/A4_mu_0p001983_plas2_%s_comp_1p0_tc_0p4_pfpr_PFPR05_C580Y.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)

    #find time point where median 580Y reach freq_criteria%
    tp = (quantile_init0p00 [0.5] >= freq_criteria).idxmax()
    
    # get the 580Y frequency of 100 runs at that time point
    t = df.loc[tp,:]
    
    #find the run that match with that time point
    run = np.where((t>=freq_criteria-0.005) & (t <=freq_criteria+0.005))[0][0]
    
    a2[0].plot(df[run], color = current_palette[i])


a2[0].set_ylim(0, 0.05)
start_time = 168;
total_time = start_time + 30*12;
a2[0].set_xlim(start_time, total_time)
a2[0].set_xticks(np.arange(start_time, total_time+1, 12*5.0))
a2[0].set_xticklabels(np.arange(0, round((total_time-start_time)/12+1), 5))


# a2.set_xlabel('Year',fontsize=fs)
a2[0].set_xlabel('Year')
# a2.set_ylabel('580Y Frequency',fontsize=fs)
a2[0].set_ylabel('580Y Frequency')

###########
#Imperial plot 3

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'imperial'


for i in range(len(init_freqs)):
    freq = init_freqs[i]
    filename = 'analysis/data/raw/figure1/%s/A4_init_res_freq_%s_tc_0p4_pfpr_PFPR05_C580Y_Imperial.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)
    
    tp = (quantile_init0p00[0.5] >= freq_criteria).idxmax()
    t = df.loc[tp,:]
    
    #find the run that match with that time point
    run = np.where((t>=freq_criteria-0.01) & (t <=freq_criteria+0.01))[0][0]
    
    a2[1].plot(df[run], color = current_palette[i])

# init_freqs_labels = ['0.00', '0.01', '0.10', '0.25','0.50']
# a2[1].legend(init_freqs_labels)

a2[1].set_ylim(0, 0.05)
start_time = 0;
total_time = 30*12
a2[1].set_xlim(start_time, total_time)
a2[1].set_xticks(np.arange(start_time, total_time+1, 12*5.0))
a2[1].set_xticklabels(np.arange(0, round((total_time-start_time)/12+1), 5))
a2[1].set_yticklabels([])


# a2.set_xlabel('Year',fontsize=fs)
a2[1].set_xlabel('Year')

#MORU plot 3

init_freqs = ['0p00', '0p01', '0p10', '0p25','0p50']
path = 'moru'


for i in range(len(init_freqs)):
    freq = init_freqs[i]
    
    filename = 'analysis/data/raw/figure1/%s/A4_init_res_freq_%s_tc_0p4_pfpr_PFPR05_C580Y_MORU.csv'%(path, freq)
    df =  pd.read_csv(filename, sep=',', header=None);
    quantile_init0p00 = df.quantile([.25, .5 , .75], axis=1).T.reset_index(drop=True)
    
    #find time point where median 580Y reach 1%
    tp = (quantile_init0p00[0.5] >= freq_criteria).idxmax()
    
    # get the 580Y frequency of 100 runs at that time point
    t = df.loc[tp,:]
    
    #find the run that match with that time point
    run = np.where((t>=freq_criteria-0.01) & (t <=freq_criteria+0.01))[0][0]
    
    a2[2].plot(df[run], color = current_palette[i])

a2[2].set_ylim(0, 0.05)
start_time = 0;
total_time = 30*12
a2[2].set_xlim(start_time, total_time)
a2[2].set_xticks(np.arange(start_time, total_time+1, 12*5.0))
a2[2].set_xticklabels(np.arange(0, round((total_time-start_time)/12+1), 5))
a2[2].set_yticklabels([])


# a2.set_xlabel('Year',fontsize=fs)
a2[2].set_xlabel('Year')

#legend

init_freqs_labels = ['0.00', '0.01', '0.10', '0.25','0.50']
a1[1].legend(init_freqs_labels)

plt.savefig("analysis/plots/fig1.png", bbox_inches = "tight", dpi = 300)
