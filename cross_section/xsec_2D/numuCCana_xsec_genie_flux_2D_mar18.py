#!/usr/bin/env python
# coding: utf-8

# In[1]:


from __future__ import division
import imp
import uproot
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
import collections
from array import array

inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/'
outputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted_improved/'+'xsec_flux_2D_mar18/' 
output_filedir = outputdir
input_filedir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/'
lib_function_dir = '/home/tmettler/Desktop/uBoone/do_plots/'

# helper functions
globale = imp.load_source('globale',lib_function_dir+'globale.py')
NuCC = imp.load_source('NuCC_function',lib_function_dir+'NuCC_function.py')
NuCC_w = imp.load_source('NuCCWeight_function',lib_function_dir+'NuCCWeight_function.py')

weight_list = [ 'expskin_FluxUnisim' ]
num_universes = 2
max_entries = 1000


# In[2]:


#!jupyter nbconvert --to script numuCCana_xsec_genie_other_2D_mar18.ipynb


# In[3]:


# initialte ROOT default canvas
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# # Flux cross section variation 100 unisims, 2D

# In[4]:


# Load input files
outputdir_png, outputdir_root, outputdir_pdf = NuCC.prepareOutput2(outputdir)
try:
    os.stat(output_filedir)
except:
    os.mkdir(output_filedir)
RootFile = ROOT.TFile(output_filedir+weight_list[0]+"_Flux_systematic.root","RECREATE");

#filename_overlay = 'NuCCana_overlay_V26_mar18.root'
filename_overlay = 'NuCCana_overlay_V26_mar18.root'
filename_data = 'NuCCana_data_V25.root'
filename_ext = 'NuCCana_ext_V25.root'
filename_dirt = 'NuCCana_dirt_V26_weight.root'
    
tree_name = 'numuCCAna'


# In[5]:


#Open all the trees of the four files (data, ext, dirt, overlay)

data, ext, dirt, overlay = NuCC.openTrees(inputdir, filename_data, filename_ext, filename_dirt, filename_overlay, tree_name)
NuCC.printNumberOfEntries(data,ext,dirt,overlay)

pot_overlay = NuCC.getPOT(inputdir,filename_overlay,tree_name)
pot_dirt =  NuCC.getPOT(inputdir,filename_dirt,tree_name)
#V25 files
pot_data =    7.644e+18  # best with tor875
data_trigger = 1838700.0 #2220362.0 #1854495.0 #4743794 # 1987072.0 # E1DCNT_wcut
ext_trigger =  18997529.0  #2120135 #5685315 # EXT

print 'POT: '
print 'Data:\t\t', pot_data
print 'Ext:\t\t', 0
print 'Overlay:\t', pot_overlay
print 'Dirt:\t\t', pot_dirt
print ''
sample = [data,ext,overlay,dirt]
scale = {data:1.0,ext:1.0,overlay:1.0,dirt:1.0}
name = {data:'data',ext:'ext',overlay:'overlay',dirt:'dirt'}

scale[data], scale[ext], scale[dirt], scale[overlay] = NuCC.calculateScale(data_trigger, ext_trigger, pot_data, pot_dirt, pot_overlay)

scale[dirt] = scale[dirt]
scale[overlay] = scale[overlay]
print 'Scalefactors: '
print 'Data:\t\t', scale[data]
print 'Ext:\t\t', scale[ext]
print 'Overlay:\t', scale[overlay]
print 'Dirt:\t\t', scale[dirt]


# In[6]:


if 1:
    filename_overlay = 'NuCCana_overlay_V26_flux.rootout4_small.root'
    #filename_overlay = 'NuCCana_overlay_V26_mar18_noflux.rootout4.root'
    filename_data = filename_data+'out4.root'
    filename_ext = filename_ext+'out4.root'
    filename_dirt = filename_dirt+'out4.root'

    tree_name = 't_out'

    data_out, ext_out, dirt_out, overlay_out = NuCC.openTreesOut(inputdir, filename_data, filename_ext, filename_dirt, filename_overlay, tree_name)
    NuCC.printNumberOfEntries(data_out,ext_out,dirt_out,overlay_out)

    sample_out = [data_out,ext_out,overlay_out,dirt_out]
    scale_out = {data_out:1.0,ext_out:1.0,overlay_out:1.0,dirt_out:1.0}
    name_out = {data_out:'data',ext_out:'ext',overlay_out:'overlay',dirt_out:'dirt'}

    scale_out[data_out], scale_out[ext_out], scale_out[dirt_out], scale_out[overlay_out] = NuCC.calculateScale(data_trigger, ext_trigger, pot_data, pot_dirt, pot_overlay)
    scale_out[dirt_out] = scale_out[dirt_out]
    scale_out[overlay_out] = scale_out[overlay_out]


# In[7]:


##### flux and number of tragets parameters###
flux = 1.16859e11/1.592e20 # flux per POT per cm2
print flux
flux = 7.3789785277e-10
print flux
roh_data = 1.3836 #data denisity g/cm3
roh_MC = 1.3954 # MC denisity g/cm3
mol = 39.95 # g for argon
N_A = 6.022140857e23 # molec/mol avogadro zahl
N_nucleons = 40.0
V_fid = ((254.8-10)-(-1.55+10))*((117.47-10)-(-115.53+10))*((1036.9-50)-(-0.1+10))
print 'Fiducial Volume: ', V_fid
##############################################

beam_flux = flux * pot_data
print 'Beam flux = {:.5e}'.format(beam_flux),' /cm2'
N_tot = roh_data*N_A*N_nucleons*V_fid/mol
print 'Number of target nuclei= {:.5e}'.format(N_tot),' /cm3'


# In[8]:


# Define signals

fidVol = '(Nu_Vx_sce>(-1.55+10) && Nu_Vx_sce<(254.8-10)) && (Nu_Vy_sce>(-115.53+10) && Nu_Vy_sce<(117.47-10)) &&(Nu_Vz_sce>(-0.1+10) && Nu_Vz_sce<(1036.9-50))'
MCfidVol = '(MCNu_Vx>(-1.55+10) && MCNu_Vx<(254.8-10)) && (MCNu_Vy>(-115.53+10) && MCNu_Vy<(117.47-10)) &&(MCNu_Vz>(-0.1+10) && MCNu_Vz<(1036.9-50))'
numu_signal = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG==13 && MCTrackPurity>0.5' # numu CC signal definition
numu_true = 'MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14' # numu CC signal definition
numu_nomu = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG!=13 && MCTrackPurity>0.5' # not an MC muon
numu_lowpur = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPurity<0.5' #low purity
numu_nc = 'fidVol && MCfidVol && MCNu_CCNC==1' # nutral current
numu_ov = 'fidVol && !MCfidVol' # out of fiducial
numu_other = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG!=14' # e.g anti nu or nue
#signal = 'MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14'
for x in sample:
    x.SetAlias('muon','(muon_candidate_key==track_key)')

energy_cut = ' && MCle_Energy>0.15'

numu_signal = numu_signal+energy_cut
numu_true = numu_true+energy_cut
numu_nomu = numu_nomu+energy_cut
numu_lowpur = numu_lowpur+energy_cut
numu_nc = numu_nc+energy_cut
numu_ov = numu_ov+energy_cut
numu_other = numu_other+energy_cut

num_fidVol = {}
for x in sample:
    x.SetAlias('fidVol',fidVol)
    x.SetAlias('MCfidVol',MCfidVol)
    x.SetAlias('numu_signal',numu_signal)
    x.SetAlias('numu_true',numu_true)
    x.SetAlias('numu_nomu',numu_nomu)
    x.SetAlias('numu_lowpur',numu_lowpur)
    x.SetAlias('numu_nc',numu_nc)
    x.SetAlias('numu_ov',numu_ov)
    x.SetAlias('numu_other',numu_other)
    num_fidVol[x] = x.GetEntries('fidVol && muon')*scale[x]
    
tot_num_fidVol = num_fidVol[ext]+num_fidVol[dirt]+num_fidVol[overlay]
overlay_signals = {'numu_signal','numu_nomu','numu_lowpur','numu_nc','numu_ov','numu_other'}


# In[9]:


for x in sample_out:
    x.SetAlias('muon','(track_key == key_muon)')

num_fidVol = {}
for x in sample_out:
    x.SetAlias('fidVol',fidVol)
    x.SetAlias('MCfidVol',MCfidVol)
    x.SetAlias('numu_signal',numu_signal)
    x.SetAlias('numu_true',numu_true)
    x.SetAlias('numu_nomu',numu_nomu)
    x.SetAlias('numu_lowpur',numu_lowpur)
    x.SetAlias('numu_nc',numu_nc)
    x.SetAlias('numu_ov',numu_ov)
    x.SetAlias('numu_other',numu_other)


# In[10]:


# Load the global variables for access of functions
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)
#NuCC.printGlobal()


# In[11]:


# initialte ROOT default canvas
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# In[12]:



track_start_border_x = '(TrackStart_x_sce <(-1.55+5) || TrackStart_x_sce > (254.8-5))'
track_end_border_x = '(TrackEnd_x_sce <(-1.55+5) || TrackEnd_x_sce > (254.8-5))'
track_start_border_y = '(TrackStart_y_sce <(-115.53+5) || TrackStart_y_sce > (117.47-5))'
track_end_border_y = '(TrackEnd_y_sce <(-115.53+5) || TrackEnd_y_sce > (117.47-5))'
track_start_border_z = '(TrackStart_z_sce <(0.1+5) || TrackStart_z_sce > (1036.9-5))'
track_end_border_z = '(TrackEnd_z_sce <(0.1+5) || TrackEnd_z_sce > (1039.9-5))'

track_end_uncontained = '(' + track_end_border_x + ' || ' + track_end_border_y + ' || ' + track_end_border_z+ ')'


data.SetAlias("track_end_uncontained",track_end_uncontained)
ext.SetAlias("track_end_uncontained",track_end_uncontained)
overlay.SetAlias("track_end_uncontained",track_end_uncontained)
dirt.SetAlias("track_end_uncontained",track_end_uncontained)
data_out.SetAlias("track_end_uncontained",track_end_uncontained)
ext_out.SetAlias("track_end_uncontained",track_end_uncontained)
overlay_out.SetAlias("track_end_uncontained",track_end_uncontained)
dirt_out.SetAlias("track_end_uncontained",track_end_uncontained)

data.SetAlias("crt_cut","(abs(crtt0_time+(crt_trig_corr_med)/1000-4)<0.9 || crtt0_time==-1)")
ext.SetAlias("crt_cut","(abs(crtt0_time+(crt_trig_corr_med)/1000-3.57+3.195-4)<0.9 || crtt0_time==-1)")
overlay.SetAlias("crt_cut","(abs(crtt0_time-3.95)<0.9 || crtt0_time==-1)")
dirt.SetAlias("crt_cut","(abs(crtt0_time-3.95)<0.9 || crtt0_time==-1)")
data_out.SetAlias("crt_cut","(abs(crtt0_time+(crt_trig_corr_med)/1000-3.95)<0.9 || crtt0_time==-1)")
ext_out.SetAlias("crt_cut","(abs(crtt0_time+(crt_trig_corr_med)/1000-3.57+3.195-3.95)<0.9 || crtt0_time==-1)")
overlay_out.SetAlias("crt_cut","(abs(crtt0_time-3.95)<0.9 || crtt0_time==-1)")
dirt_out.SetAlias("crt_cut","(abs(crtt0_time-3.95)<0.9 || crtt0_time==-1)")

crt_tom_cut = 'nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut'

data.SetAlias("crt_tom_cut",crt_tom_cut)
ext.SetAlias("crt_tom_cut",crt_tom_cut)
overlay.SetAlias("crt_tom_cut",crt_tom_cut)
dirt.SetAlias("crt_tom_cut",crt_tom_cut)
data_out.SetAlias("crt_tom_cut",crt_tom_cut)
ext_out.SetAlias("crt_tom_cut",crt_tom_cut)
overlay_out.SetAlias("crt_tom_cut",crt_tom_cut)
dirt_out.SetAlias("crt_tom_cut",crt_tom_cut)

weight_name = 'EventWeight*TunedCentralValue_Genie'


# In[13]:


from array import array
#mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
#mom_bins = [ 0.00, 0.225, 0.28, 0.33, 0.39, 0.52, 0.78, 1.21, 2.5]

mom_bins = {}
#mom_bins.Append()
mom_bins[0] = [ 0.00, 0.18, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[1] = [ 0.00, 0.18, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[2] = [ 0.00, 0.18, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[3] = [ 0.00, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[4] = [ 0.00, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[5] = [ 0.00, 0.30, 0.45, 0.77, 2.50 ]
mom_bins[6] = [ 0.00, 0.30, 0.45, 0.77, 1.28, 2.50 ]
mom_bins[7] = [ 0.00, 0.30, 0.45, 0.77, 1.28, 2.50 ]
mom_bins[8] = [ 0.00, 0.30, 0.45, 0.77, 1.28, 2.50 ]

theta_bins = [ -1.00, -0.50, 0.00, 0.28, 0.47, 0.63, 0.765, 0.865, 0.935, 1.00 ]

len_theta = len(theta_bins) - 1
len_mom = {}
sum_bins = 0
for i in range(len_theta):
    len_mom[i] = len(mom_bins[i]) - 1
    #print len_mom[i]
    sum_bins+=len_mom[i]


# In[14]:


#weight_list = [ 'AxFFCCQEshape_UBGenie','DecayAngMEC_UBGenie' ]
#weight_list = [ 'AxFFCCQEshape_UBGenie' , 'DecayAngMEC_UBGenie', 'NormCCCOH_UBGenie', 'NormNCCOH_UBGenie',\
#               'RPA_CCQE_Reduced_UBGenie','RPA_CCQE_UBGenie','ThetaDelta2NRad_UBGenie','Theta_Delta2Npi_UBGenie',\
#                'VecFFCCQEshape_UBGenie','XSecShape_CCMEC_UBGenie']
num_para = len(weight_list)
print 'Nuber of parameters: ',num_para


# In[15]:


#weight_name = 'EventWeight*TunedCentralValue_Genie'
title = 'true moum momentum efficiency'
nenner_cut = 'numu_true'
zahler_cut = 'fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1 && numu_signal'

h_1d_sel_cv = []
h_1d_gen_cv = []
h_1d_eff_cv = []

num_sel_cv = np.zeros(sum_bins+1)
num_gen_cv = np.zeros(sum_bins+1)
num_eff_cv = np.zeros(sum_bins+1)

for t_bin in range(len_theta):
    h_1d_sel_cv.append(ROOT.TH1F('h_1d_sel_cv['+str(t_bin)+']',"true momentum, selected (cv)",len_mom[t_bin],array('f',mom_bins[t_bin])))
    h_1d_gen_cv.append( ROOT.TH1F('h_1d_gen_cv['+str(t_bin)+']',"true momentum, generated (cv)",len_mom[t_bin],array('f',mom_bins[t_bin])) )
    h_1d_eff_cv.append( ROOT.TH1F('h_1d_eff_cv['+str(t_bin)+']',"true efficiency (cv)",len_mom[t_bin],array('f',mom_bins[t_bin])) )

bin_counter = 0
for t_bin in range(len_theta):
    this_cut = ' && cos(MCNu_leptonTheta)>'+str(theta_bins[t_bin])+' && cos(MCNu_leptonTheta)<'+str(theta_bins[t_bin+1])
    print this_cut
    globale.overlay_out.Draw('MCle_Energy>>h_1d_sel_cv['+str(t_bin)+']',weight_name+'*TunedCentralValue_Genie*('+zahler_cut+this_cut+')','',max_entries)
    globale.overlay_out.Draw('MCle_Energy>>h_1d_gen_cv['+str(t_bin)+']',weight_name+'*TunedCentralValue_Genie*('+nenner_cut+this_cut+')','',max_entries)
    h_1d_eff_cv[t_bin] = h_1d_sel_cv[t_bin].Clone()
    h_1d_eff_cv[t_bin].Divide(h_1d_gen_cv[t_bin])
    for m_bin in range(len_mom[t_bin]):
        num_sel_cv[bin_counter] = h_1d_sel_cv[t_bin].GetBinContent(m_bin+1)
        num_gen_cv[bin_counter] = h_1d_gen_cv[t_bin].GetBinContent(m_bin+1)
        num_eff_cv[bin_counter] = h_1d_eff_cv[t_bin].GetBinContent(m_bin+1)
        bin_counter+=1
    num_sel_cv[sum_bins] += h_1d_sel_cv[t_bin].GetBinContent(len_mom[t_bin]+1)
    num_gen_cv[sum_bins] += h_1d_gen_cv[t_bin].GetBinContent(len_mom[t_bin]+1)
    #num_eff[sum_bins] += h_1d_eff[t_bin].GetBinContent(len_mom[t_bin]+1)
    h_1d_sel_cv[t_bin].Write('h_1d_sel_cv['+str(t_bin)+']')
    h_1d_gen_cv[t_bin].Write('h_1d_gen_cv['+str(t_bin)+']')
    h_1d_eff_cv[t_bin].Write('h_1d_eff_cv['+str(t_bin)+']')

num_eff_cv[sum_bins] = num_sel_cv[sum_bins]/num_gen_cv[sum_bins]
    
np.save(output_filedir+'num_sel_cv',num_sel_cv)
np.save(output_filedir+'num_gen_cv',num_gen_cv)
np.save(output_filedir+'num_eff_cv',num_eff_cv)


# In[16]:


weight_name = 'EventWeight*TunedCentralValue_Genie'
title = 'true moum momentum efficiency'
nenner_cut = 'numu_true'
zahler_cut = 'fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1 && numu_signal'

h_1d_sel = []
h_1d_gen = []
h_1d_eff = []

num_sel = np.zeros((num_para,num_universes,sum_bins+1))
num_gen = np.zeros((num_para,num_universes,sum_bins+1))
num_eff = np.zeros((num_para,num_universes,sum_bins+1))


for para in range(num_para):
    h_1d_sel.append([])
    h_1d_gen.append([])
    h_1d_eff.append([])
    for uni in range(num_universes):
        h_1d_sel[para].append([])
        h_1d_gen[para].append([])
        h_1d_eff[para].append([])
        for t_bin in range(len_theta):
            h_1d_sel[para][uni].append(ROOT.TH1F('h_1d_sel['+str(para)+']['+str(uni)+']['+str(t_bin)+']',"true momentum, selected",len_mom[t_bin],array('f',mom_bins[t_bin])))
            h_1d_gen[para][uni].append( ROOT.TH1F('h_1d_gen['+str(para)+']['+str(uni)+']['+str(t_bin)+']',"true momentum, generated",len_mom[t_bin],array('f',mom_bins[t_bin])) )
            h_1d_eff[para][uni].append( ROOT.TH1F('h_1d_eff['+str(para)+']['+str(uni)+']['+str(t_bin)+']',"true efficiency",len_mom[t_bin],array('f',mom_bins[t_bin])) )

for para in range(num_para):
    for uni in range(num_universes):
        bin_counter = 0
        for t_bin in range(len_theta):
            this_cut = ' && cos(MCNu_leptonTheta)>'+str(theta_bins[t_bin])+' && cos(MCNu_leptonTheta)<'+str(theta_bins[t_bin+1])
            #print this_cut
            globale.overlay_out.Draw('MCle_Energy>>h_1d_sel['+str(para)+']['+str(uni)+']['+str(t_bin)+']',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+zahler_cut+this_cut+')','',max_entries)
            globale.overlay_out.Draw('MCle_Energy>>h_1d_gen['+str(para)+']['+str(uni)+']['+str(t_bin)+']',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+nenner_cut+this_cut+')','',max_entries)
            h_1d_eff[para][uni][t_bin] = h_1d_sel[para][uni][t_bin].Clone()
            h_1d_eff[para][uni][t_bin].Divide(h_1d_gen[para][uni][t_bin])
            for m_bin in range(len_mom[t_bin]):
                num_sel[para][uni][bin_counter] = h_1d_sel[para][uni][t_bin].GetBinContent(m_bin+1)
                num_gen[para][uni][bin_counter] = h_1d_gen[para][uni][t_bin].GetBinContent(m_bin+1)
                num_eff[para][uni][bin_counter] = h_1d_eff[para][uni][t_bin].GetBinContent(m_bin+1)
                bin_counter+=1
            num_sel[para][uni][sum_bins] += h_1d_sel[para][uni][t_bin].GetBinContent(len_mom[t_bin]+1)
            num_gen[para][uni][sum_bins] += h_1d_gen[para][uni][t_bin].GetBinContent(len_mom[t_bin]+1)
            #num_eff[sum_bins] += h_1d_eff[t_bin].GetBinContent(len_mom[t_bin]+1)
            h_1d_sel[para][uni][t_bin].Write('h_1d_sel['+str(para)+']['+str(uni)+']['+str(t_bin)+']')
            h_1d_gen[para][uni][t_bin].Write('h_1d_gen['+str(para)+']['+str(uni)+']['+str(t_bin)+']')
            h_1d_eff[para][uni][t_bin].Write('h_1d_eff['+str(para)+']['+str(uni)+']['+str(t_bin)+']')

        num_eff[para][uni][sum_bins] = num_sel[para][uni][sum_bins]/num_gen[para][uni][sum_bins]
    
np.save(output_filedir+'num_sel',num_sel)
np.save(output_filedir+'num_gen',num_gen)
np.save(output_filedir+'num_eff',num_eff)


# In[17]:


h_sel_all_cv = ROOT.TH1F('h_sel_all_cv','h_sel_all_cv',sum_bins+1,0,sum_bins+1)
h_gen_all_cv = ROOT.TH1F('h_gen_all_cv','h_gen_all_cv',sum_bins+1,0,sum_bins+1)
h_eff_all_cv = ROOT.TH1F('h_eff_all_cv','h_eff_all_cv',sum_bins+1,0,sum_bins+1)

for i in range(sum_bins+1):
    h_sel_all_cv.SetBinContent(i+1,num_sel_cv[i])
    h_gen_all_cv.SetBinContent(i+1,num_gen_cv[i])
    h_eff_all_cv.SetBinContent(i+1,num_eff_cv[i])
h_sel_all_cv.Write('h_sel_all_cv')
h_gen_all_cv.Write('h_gen_all_cv')
h_eff_all_cv.Write('h_eff_all_cv')


# In[18]:


h_sel_all = []
h_gen_all = []
h_eff_all = []

for para in range(num_para):
    h_sel_all.append([])
    h_gen_all.append([])
    h_eff_all.append([])
    for uni in range(num_universes):
        h_sel_all[para].append(ROOT.TH1F('h_sel_all['+str(para)+']['+str(uni)+']','h_sel_all',sum_bins+1,0,sum_bins+1))
        h_gen_all[para].append(ROOT.TH1F('h_gen_all['+str(para)+']['+str(uni)+']','h_gen_all',sum_bins+1,0,sum_bins+1))
        h_eff_all[para].append(ROOT.TH1F('h_eff_all['+str(para)+']['+str(uni)+']','h_eff_all',sum_bins+1,0,sum_bins+1))

        for i in range(sum_bins+1):
            h_sel_all[para][uni].SetBinContent(i+1,num_sel[para][uni][i])
            h_gen_all[para][uni].SetBinContent(i+1,num_gen[para][uni][i])
            h_eff_all[para][uni].SetBinContent(i+1,num_eff[para][uni][i])
        h_sel_all[para][uni].Write('h_sel_all['+str(para)+']['+str(uni)+']')
        h_gen_all[para][uni].Write('h_gen_all['+str(para)+']['+str(uni)+']')
        h_eff_all[para][uni].Write('h_eff_all['+str(para)+']['+str(uni)+']')


# In[19]:


#calculate cv true to reco matrix
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
bin_counter_true = 0
run = 1
if run==0:
    true_reco_matrix_cv = np.load(outputdir+'true_reco_matrix_cv.npy')

if run:
    true_reco_matrix_cv = np.zeros((sum_bins+1,sum_bins+1))
    for t_bin in range(len_theta):
        print 'next theta',t_bin
        #bin_counter_true+=1
        for m_bin in range(len_mom[t_bin]):
            bin_counter_reco = 0
            # define the true bins borders
            #plot for each true bin in mom and theta the distribution in reco for mom and theta, 42 times 42 bins
            mom_min = mom_bins[t_bin][m_bin]
            mom_max = mom_bins[t_bin][m_bin+1]
            theta_min = theta_bins[t_bin]
            theta_max = theta_bins[t_bin+1]
            #print mom_min, mom_max, theta_min, theta_max
            this_cut = ' && MCle_Energy>'+str(mom_min)+' && MCle_Energy<'+str(mom_max)+' && cos(MCNu_leptonTheta)>'+str(theta_min)+' && cos(MCNu_leptonTheta)<'+str(theta_max)
            for t2_bin in range(len_theta):
                #print 'new theta'
                #print len_mom[t2_bin], mom_bins[t2_bin]
                h_temp = ROOT.TH2F("h_temp",'h_temp',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                h_temp_r = ROOT.TH2F("h_temp_r",'h_temp_r',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                globale.overlay_out.Draw('cos(TrackTheta):TrackMomMCS_mom'+'>>h_temp',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& numu_signal && track_end_uncontained)','',max_entries)
                globale.overlay_out.Draw('cos(TrackTheta):TrackMomRange_mu'+'>>h_temp_r',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& numu_signal && !track_end_uncontained)','',max_entries)
                h_temp.Add(h_temp_r)
                #print h_temp.GetBinContent(1,1)
                #h_temp.Draw()
                #c1.Draw()
                for m2_bin in range(len_mom[t2_bin]):
                    #print bin_counter_true,bin_counter_reco
                    true_reco_matrix_cv[bin_counter_true,bin_counter_reco] = h_temp.GetBinContent(m2_bin+1,t2_bin+1)
                    #print bin_counter_true,bin_counter_reco, true_reco_matrix[bin_counter_true,bin_counter_reco]
                    bin_counter_reco+=1
                true_reco_matrix_cv[bin_counter_true,bin_counter_reco] += h_temp.GetBinContent(len_mom[t2_bin]+1,t2_bin+1)
                #raw_input("Press Enter to continue...")
                del h_temp
                del h_temp_r
                #bin_counter_reco+=1
            bin_counter_true+=1

    # fill true overflow bin
    bin_counter_reco = 0
    this_cut = ' && MCle_Energy>2.5'
    for t2_bin in range(len_theta):
        #print 'new theta'
        #print len_mom[t2_bin], mom_bins[t2_bin]
        h_temp = ROOT.TH2F("h_temp",'h_temp',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
        h_temp_r = ROOT.TH2F("h_temp_r",'h_temp_r',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
        globale.overlay_out.Draw('cos(TrackTheta):TrackMomMCS_mom'+'>>h_temp',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& numu_signal && track_end_uncontained)','',max_entries)
        globale.overlay_out.Draw('cos(TrackTheta):TrackMomRange_mu'+'>>h_temp_r',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& numu_signal && !track_end_uncontained)','',max_entries)
        h_temp.Add(h_temp_r)
        for m2_bin in range(len_mom[t2_bin]):
            #print bin_counter_true,bin_counter_reco
            true_reco_matrix_cv[sum_bins,bin_counter_reco] = h_temp.GetBinContent(m2_bin+1,t2_bin+1)
            #print bin_counter_true,bin_counter_reco, true_reco_matrix[bin_counter_true,bin_counter_reco]
            bin_counter_reco+=1
        true_reco_matrix_cv[sum_bins,bin_counter_reco] += h_temp.GetBinContent(len_mom[t2_bin]+1,t2_bin+1)
        #raw_input("Press Enter to continue...")
        del h_temp
        del h_temp_r


# In[20]:


# save cv true to reco matrix
h_true_reco_cv = ROOT.TH2F("h_true_reco_cv",'Migration matrix (cv)',sum_bins+1,0,sum_bins+1,sum_bins+1,0,sum_bins+1)
for i in range(sum_bins+1):
    for j in range(sum_bins+1):
        h_true_reco_cv.SetBinContent(i+1,j+1, true_reco_matrix_cv[i,j])
h_true_reco_cv.Write('h_true_reco_cv')
if run:
    np.save(outputdir+'true_reco_matrix_cv',true_reco_matrix_cv)


# In[ ]:


# calculate all genie true to reco matrix
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
run = 1
if run==0:
    true_reco_matrix = np.load(outputdir+'true_reco_matrix.npy')

if run:
    true_reco_matrix = np.zeros((num_para,num_universes,sum_bins+1,sum_bins+1))
    for para in range(num_para):
        for uni in range(num_universes):
            bin_counter_true = 0
            for t_bin in range(len_theta):
                #print 'next theta',t_bin
                #bin_counter_true+=1
                for m_bin in range(len_mom[t_bin]):
                    bin_counter_reco = 0
                    # define the true bins borders
                    #plot for each true bin in mom and theta the distribution in reco for mom and theta, 42 times 42 bins
                    mom_min = mom_bins[t_bin][m_bin]
                    mom_max = mom_bins[t_bin][m_bin+1]
                    theta_min = theta_bins[t_bin]
                    theta_max = theta_bins[t_bin+1]
                    #print mom_min, mom_max, theta_min, theta_max
                    this_cut = ' && MCle_Energy>'+str(mom_min)+' && MCle_Energy<'+str(mom_max)+' && cos(MCNu_leptonTheta)>'+str(theta_min)+' && cos(MCNu_leptonTheta)<'+str(theta_max)
                    for t2_bin in range(len_theta):
                        h_temp = ROOT.TH2F("h_temp",'h_temp',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                        h_temp_r = ROOT.TH2F("h_temp_r",'h_temp_r',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                        globale.overlay_out.Draw('cos(TrackTheta):TrackMomMCS_mom'+'>>h_temp',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& numu_signal && track_end_uncontained)','',max_entries)
                        globale.overlay_out.Draw('cos(TrackTheta):TrackMomRange_mu'+'>>h_temp_r',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& numu_signal && !track_end_uncontained)','',max_entries)
                        h_temp.Add(h_temp_r)
                        for m2_bin in range(len_mom[t2_bin]):
                            #print bin_counter_true,bin_counter_reco
                            true_reco_matrix[para][uni][bin_counter_true][bin_counter_reco] = h_temp.GetBinContent(m2_bin+1,t2_bin+1)
                            #print bin_counter_true,bin_counter_reco, true_reco_matrix[bin_counter_true,bin_counter_reco]
                            bin_counter_reco+=1
                        true_reco_matrix[para][uni][bin_counter_true][bin_counter_reco] += h_temp.GetBinContent(len_mom[t2_bin]+1,t2_bin+1)
                        del h_temp
                        del h_temp_r
                        #bin_counter_reco+=1
                    bin_counter_true+=1
            # fill true overflow bin
            bin_counter_reco = 0
            this_cut = ' && MCle_Energy>2.5'
            for t2_bin in range(len_theta):
                h_temp = ROOT.TH2F("h_temp",'h_temp',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                h_temp_r = ROOT.TH2F("h_temp_r",'h_temp_r',len_mom[t2_bin],array('f',mom_bins[t2_bin]),len_theta,array('f',theta_bins))
                globale.overlay_out.Draw('cos(TrackTheta):TrackMomMCS_mom'+'>>h_temp',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& numu_signal && track_end_uncontained)','',max_entries)
                globale.overlay_out.Draw('cos(TrackTheta):TrackMomRange_mu'+'>>h_temp_r',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& numu_signal && !track_end_uncontained)','',max_entries)
                h_temp.Add(h_temp_r)
                for m2_bin in range(len_mom[t2_bin]):
                    true_reco_matrix[para][uni][sum_bins][bin_counter_reco] = h_temp.GetBinContent(m2_bin+1,t2_bin+1)
                    bin_counter_reco+=1
                true_reco_matrix[para][uni][sum_bins][bin_counter_reco] += h_temp.GetBinContent(len_mom[t2_bin]+1,t2_bin+1)
                del h_temp
                del h_temp_r


# In[ ]:


# save all genie true to reco matrix
h_true_reco = []
for para in range(num_para):
    h_true_reco.append([])
    for uni in range(num_universes):
        h_true_reco[para].append(ROOT.TH2F('h_true_reco['+str(para)+']['+str(uni)+']','Migration matrix',sum_bins+1,0,sum_bins+1,sum_bins+1,0,sum_bins+1))
        for i in range(sum_bins+1):
            for j in range(sum_bins+1):
                h_true_reco[para][uni].SetBinContent(i+1,j+1, true_reco_matrix[para][uni][i,j])
        h_true_reco[para][uni].Write('h_true_reco['+str(para)+']['+str(uni)+']')
if run:
    np.save(outputdir+'true_reco_matrix',true_reco_matrix)


# In[ ]:


# calculate cv smearing matrix
smearing_matrix_cv = np.zeros((sum_bins+1,sum_bins+1))
sum_reco_cv = true_reco_matrix_cv.sum(axis=1)
#print sum_reco
smearing_matrix_cv = true_reco_matrix_cv / (sum_reco_cv[:,None] + 1e-80)

np.save(outputdir+'smearing_matrix_cv',smearing_matrix_cv)


# In[ ]:


# save for the cv smearing matrix
h_smearing_cv = ROOT.TH2F("h_smearing_cv",'smearing matrix (cv)',sum_bins+1,0,sum_bins+1,sum_bins+1,0,sum_bins+1)
for i in range(sum_bins+1):
    for j in range(sum_bins+1):
        h_smearing_cv.SetBinContent(i+1,j+1, smearing_matrix_cv[i,j])

ROOT.gStyle.SetPaintTextFormat("0.1f");
ROOT.gStyle.SetPalette(55);
h_smearing_cv.SetContour(500);
h_smearing_cv.SetMaximum(1.0)
#c1.SetFrameFillColor(ROOT.TColor.GetColorPalette(0));
for i in range(1,sum_bins+1):
    h_smearing_cv.GetXaxis().SetBinLabel(i, str(i))
    h_smearing_cv.GetYaxis().SetBinLabel(i, str(i))
h_smearing_cv.GetXaxis().SetBinLabel(43, 'OF')
h_smearing_cv.GetYaxis().SetBinLabel(43, 'OF')

h_smearing_cv.SetXTitle("True bin number")
h_smearing_cv.SetYTitle("Reco bin number")
h_smearing_cv.GetYaxis().SetTitleSize(0.05)
h_smearing_cv.GetYaxis().SetTitleOffset(0.0)
h_smearing_cv.GetYaxis().SetLabelSize(0.02)
h_smearing_cv.GetXaxis().SetTitleSize(0.05)
h_smearing_cv.GetXaxis().SetLabelSize(0.02)
h_smearing_cv.GetXaxis().SetTitleOffset(1)

#h_smearing.Draw('colz')
h_smearing_cv.Write('h_smearing_cv')


# In[ ]:


smearing_matrix = np.zeros((num_para,num_universes,sum_bins+1,sum_bins+1))
sum_reco = []
for para in range(num_para):
    sum_reco.append([])
    for uni in range(num_universes):
        sum_reco[para].append(true_reco_matrix[para][uni].sum(axis=1))
        #print sum_reco
        smearing_matrix[para][uni] = true_reco_matrix[para][uni] / (sum_reco[para][uni][:,None] + 1e-80)

np.save(outputdir+'smearing_matrix',smearing_matrix)


# In[ ]:


h_smearing = []
for para in range(num_para):
    h_smearing.append([])
    for uni in range(num_universes):
        h_smearing[para].append(ROOT.TH2F('h_smearing['+str(para)+']['+str(uni)+']','smearing matrix',sum_bins+1,0,sum_bins+1,sum_bins+1,0,sum_bins+1))
        for i in range(sum_bins+1):
            for j in range(sum_bins+1):
                h_smearing[para][uni].SetBinContent(i+1,j+1, smearing_matrix[para][uni][i][j])

        ROOT.gStyle.SetPaintTextFormat("0.1f");
        ROOT.gStyle.SetPalette(55);
        h_smearing[para][uni].SetContour(500);
        h_smearing[para][uni].SetMaximum(1.0)
        #c1.SetFrameFillColor(ROOT.TColor.GetColorPalette(0));
        for i in range(1,sum_bins+1):
            h_smearing[para][uni].GetXaxis().SetBinLabel(i, str(i))
            h_smearing[para][uni].GetYaxis().SetBinLabel(i, str(i))
        h_smearing[para][uni].GetXaxis().SetBinLabel(43, 'OF')
        h_smearing[para][uni].GetYaxis().SetBinLabel(43, 'OF')

        h_smearing[para][uni].SetXTitle("True bin number")
        h_smearing[para][uni].SetYTitle("Reco bin number")
        h_smearing[para][uni].GetYaxis().SetTitleSize(0.05)
        h_smearing[para][uni].GetYaxis().SetTitleOffset(0.0)
        h_smearing[para][uni].GetYaxis().SetLabelSize(0.02)
        h_smearing[para][uni].GetXaxis().SetTitleSize(0.05)
        h_smearing[para][uni].GetXaxis().SetLabelSize(0.02)
        h_smearing[para][uni].GetXaxis().SetTitleOffset(1)

        #h_smearing.Draw('colz')
        h_smearing[para][uni].Write('h_smearing['+str(para)+']['+str(uni)+']')


# In[ ]:


# calculate cv smeared efficiency tilde
eff_tilde_cv = smearing_matrix_cv.dot(num_sel_cv)/(smearing_matrix_cv.dot(num_gen_cv)+1e-80)
h_1d_eff_tilde_cv = []
for t_bin in range(len_theta):
    h_1d_eff_tilde_cv.append( ROOT.TH1F('h_1d_eff_tilde_cv['+str(t_bin)+']',"reco efficiency (cv)",len_mom[t_bin],array('f',mom_bins[t_bin])) )
bin_counter = 0
for t_bin in range(len_theta):
    for m_bin in range(len_mom[t_bin]):
        h_1d_eff_tilde_cv[t_bin].SetBinContent(m_bin+1,eff_tilde_cv[bin_counter])
        h_1d_eff_tilde_cv[t_bin].SetBinError(m_bin+1,0)
        bin_counter += 1
    h_1d_eff_tilde_cv[t_bin].Write('h_1d_eff_tilde_cv['+str(t_bin)+']')
np.save(outputdir+'eff_tilde_cv',eff_tilde_cv)


# In[ ]:


# calculate smeared efficiency tilde
eff_tilde = []
h_1d_eff_tilde = []
for para in range(num_para):
    eff_tilde.append([])
    h_1d_eff_tilde.append([])
    for uni in range(num_universes):
        eff_tilde[para].append(smearing_matrix[para][uni].dot(num_sel[para][uni])/(smearing_matrix[para][uni].dot(num_gen[para][uni])+1e-80))
        h_1d_eff_tilde[para].append([])
        for t_bin in range(len_theta):
            h_1d_eff_tilde[para][uni].append( ROOT.TH1F('h_1d_eff_tilde['+str(para)+']['+str(uni)+']['+str(t_bin)+']',"reco efficiency",len_mom[t_bin],array('f',mom_bins[t_bin])) )
        bin_counter = 0
        for t_bin in range(len_theta):
            for m_bin in range(len_mom[t_bin]):
                #print eff_tilde[uni]
                #print '-------------------------'
                h_1d_eff_tilde[para][uni][t_bin].SetBinContent(m_bin+1,eff_tilde[para][uni][bin_counter])
                h_1d_eff_tilde[para][uni][t_bin].SetBinError(m_bin+1,0)
                bin_counter += 1
            h_1d_eff_tilde[para][uni][t_bin].Write('h_1d_eff_tilde['+str(para)+']['+str(uni)+']['+str(t_bin)+']')
np.save(outputdir+'eff_tilde',eff_tilde)


# In[ ]:



h_overlay_cv = []
h_temp_cv = []

for t_bin in range(len_theta):
    h_overlay_cv.append( ROOT.TH1F('h_overlay_cv['+str(t_bin)+']',"h_overlay (cv)",len_mom[t_bin],array('f',mom_bins[t_bin])) )
    h_temp_cv.append( ROOT.TH1F('h_temp_cv['+str(t_bin)+']','h_temp',len_mom[t_bin],array('f',mom_bins[t_bin])))
    
#bin_counter = 0
for t_bin in range(len_theta):
    this_cut = ' && cos(TrackTheta)>'+str(theta_bins[t_bin])+' && cos(TrackTheta)<'+str(theta_bins[t_bin+1])
    print this_cut

    globale.overlay_out.Draw('TrackMomMCS_mom'+'>>h_overlay_cv['+str(t_bin)+']',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& !numu_signal && track_end_uncontained)','',max_entries)
    globale.overlay_out.Draw('TrackMomRange_mu'+'>>h_temp_cv['+str(t_bin)+']',weight_name+'*TunedCentralValue_Genie*('+cut+this_cut+'&& !numu_signal && !track_end_uncontained)','',max_entries)
    h_overlay_cv[t_bin].Add(h_temp_cv[t_bin])
    h_overlay_cv[t_bin].Scale(scale[overlay])
    
    h_overlay_cv[t_bin].Write('h_overlay_cv['+str(t_bin)+']')


# In[ ]:



h_overlay = []
h_temp = []

bkg_all = np.zeros((num_para,num_universes,sum_bins+1))

for para in range(num_para):
    h_overlay.append([])
    h_temp.append([])
    for uni in range(num_universes):
        h_overlay[para].append([])
        h_temp[para].append([])
        for t_bin in range(len_theta):
            h_overlay[para][uni].append( ROOT.TH1F('h_overlay['+str(para)+']['+str(uni)+']['+str(t_bin)+']',"h_overlay",len_mom[t_bin],array('f',mom_bins[t_bin])) )
            h_temp[para][uni].append( ROOT.TH1F('h_temp['+str(para)+']['+str(uni)+']['+str(t_bin)+']','h_temp',len_mom[t_bin],array('f',mom_bins[t_bin])))

    #bin_counter = 0
for para in range(num_para):
    for uni in range(num_universes):
        bin_counter = 0
        for t_bin in range(len_theta):
            this_cut = ' && cos(TrackTheta)>'+str(theta_bins[t_bin])+' && cos(TrackTheta)<'+str(theta_bins[t_bin+1])
            print this_cut

            globale.overlay_out.Draw('TrackMomMCS_mom'+'>>h_overlay['+str(para)+']['+str(uni)+']['+str(t_bin)+']',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& !numu_signal && track_end_uncontained)','',max_entries)
            globale.overlay_out.Draw('TrackMomRange_mu'+'>>h_temp['+str(para)+']['+str(uni)+']['+str(t_bin)+']',weight_name+'*'+weight_list[para]+'['+str(uni)+']*('+cut+this_cut+'&& !numu_signal && !track_end_uncontained)','',max_entries)
            h_overlay[para][uni][t_bin].Add(h_temp[para][uni][t_bin])
            h_overlay[para][uni][t_bin].Scale(scale[overlay])
            
            for m_bin in range(len_mom[t_bin]):
                bkg_all[para][uni][bin_counter] = h_overlay[para][uni][t_bin].GetBinContent(m_bin)
                bin_counter += 1

            h_overlay[para][uni][t_bin].Write('h_overlay['+str(para)+']['+str(uni)+']['+str(t_bin)+']')

h_bkg_all = []
for para in range(num_para):
    h_bkg_all.append([])
    for uni in range(num_universes):
        h_bkg_all[para].append(ROOT.TH1F('h_bkg_all['+str(para)+']['+str(uni)+']','Background',sum_bins+1,0,sum_bins+1))
        for i in range(sum_bins+1):
            h_bkg_all[para][uni].SetBinContent(i+1, bkg_all[para][uni][i])
        h_bkg_all[para][uni].Write('h_bkg_all['+str(para)+']['+str(uni)+']')
        


# In[ ]:


# close the root file with all the histos
RootFile.Close()
#RootFile = ROOT.TFile(output_filedir+"Other_Genie_2D.root","update");


# In[ ]:




