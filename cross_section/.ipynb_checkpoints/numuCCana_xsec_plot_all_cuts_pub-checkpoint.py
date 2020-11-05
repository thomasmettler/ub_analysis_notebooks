#!/usr/bin/env python
# coding: utf-8

# In[4]:


from __future__ import division
import imp
import uproot
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
from array import array
import collections

inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/'
outputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/event_rates_all_pub_2/'
#input_filedir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/'
lib_function_dir = '/home/tmettler/Desktop/uBoone/do_plots/'

# helper functions
globale = imp.load_source('globale',lib_function_dir+'globale.py')
NuCC = imp.load_source('NuCC_function',lib_function_dir+'NuCC_function.py')
NuCC_w = imp.load_source('NuCCWeight_function',lib_function_dir+'NuCCWeight_function.py')

filename_overlay = 'NuCCana_overlay_V26_weight.root'

#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_CV.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_dEdx.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_LYAttenuation.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_LYdown.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_LYRayleigh.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_recomb2.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_SCE.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_WireAngleXZ.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_WireAngleYZ.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_WireModX.root'
#filename_overlay = 'detector_variation_reweight/NuCCana_ovleray_detsys_reweight_WireModYZ.root'

#outputdir = inputdir+'xsec_detsys/CV/' 
#outputdir = inputdir+'xsec_detsys/dEdx/' 
#outputdir = inputdir+'xsec_detsys/LYAttenuation/' 
#outputdir = inputdir+'xsec_detsys/LYdown/' 
#outputdir = inputdir+'xsec_detsys/LYRayleigh/' 
#outputdir = inputdir+'xsec_detsys/recomb2/' 
#outputdir = inputdir+'xsec_detsys/SCE/' 
#outputdir = inputdir+'xsec_detsys/WireAngleXZ/' 
#outputdir = inputdir+'xsec_detsys/WireAngleYZ/' 
#outputdir = inputdir+'xsec_detsys/WireModX/' 
#outputdir = inputdir+'xsec_detsys/WireModYZ/' 


# In[5]:


#!jupyter nbconvert --to script numuCCana_xsec_crosscheck.ipynb


# In[6]:


# initialte ROOT default canvas
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.05)
c1.SetBottomMargin(0.14)


# # CRTinBNB tagger, Overlay or MC

# In[7]:


# Load input files
#inputdir = '/home/thomasm/numuCC/V33/10kevents/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir)

#RootFile = ROOT.TFile(output_filedir+"xsec_cross_check.root","RECREATE");


#filename_overlay = 'NuCCana_overlay_V26_weight.root'

filename_data = 'NuCCana_data_V25.root'
filename_ext = 'NuCCana_ext_V25_G1.root'
filename_dirt = 'NuCCana_dirt_V26_weight.root'
    
tree_name = 'numuCCAna'


# In[8]:


#Open all the trees of the four files (data, ext, dirt, overlay)

data, ext, dirt, overlay = NuCC.openTrees(inputdir, filename_data, filename_ext, filename_dirt, filename_overlay, tree_name)
NuCC.printNumberOfEntries(data,ext,dirt,overlay)

pot_overlay = NuCC.getPOT(inputdir,filename_overlay,tree_name)
pot_dirt =  NuCC.getPOT(inputdir,filename_dirt,tree_name)
#V25 files
#pot_data =     8.649e+18  # best with tor875
#data_trigger = 2220362.0 #1854495.0 #4743794 # 1987072.0 # E1DCNT_wcut

pot_data =     7.644e+18  # best with tor875
data_trigger = 1838700.0 #1854495.0 #4743794 # 1987072.0 # E1DCNT_wcut

ext_trigger =  85768579.0  #2120135 #5685315 # EXT


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


# In[9]:


if 1:
    filename_overlay = filename_overlay+'out4.root'
    #filename_overlay = 'NuCCana_overlay_points_1kev.rootout2.root'
    filename_data = filename_data+'out4.root'
    filename_ext = filename_ext+'out4.root'
    filename_dirt = filename_dirt+'out4.root'
    
    #filename_overlay = 'NuCCana_overlay_Bin4.root'+'out4.root'
    #filename_data = 'NuCCana_data_Bin4.root'+'out4.root'
    #filename_ext = 'NuCCana_ext_Bin4.root'+'out4.root'
    #filename_dirt = 'NuCCana_dirt_Bin4.root'+'out4.root'

    tree_name = 't_out'

    data_out, ext_out, dirt_out, overlay_out = NuCC.openTreesOut(inputdir, filename_data, filename_ext, filename_dirt, filename_overlay, tree_name)
    NuCC.printNumberOfEntries(data_out,ext_out,dirt_out,overlay_out)

    sample_out = [data_out,ext_out,overlay_out,dirt_out]
    scale_out = {data_out:1.0,ext_out:1.0,overlay_out:1.0,dirt_out:1.0}
    name_out = {data_out:'data',ext_out:'ext',overlay_out:'overlay',dirt_out:'dirt'}

    scale_out[data_out], scale_out[ext_out], scale_out[dirt_out], scale_out[overlay_out] = NuCC.calculateScale(data_trigger, ext_trigger, pot_data, pot_dirt, pot_overlay)
    scale_out[dirt_out] = scale_out[dirt_out]
    scale_out[overlay_out] = scale_out[overlay_out]


# In[10]:


##### flux and number of tragets parameters###
#flux = 1.16859e11/1.592e20 # flux per POT per cm2
flux = 7.3789785277e-10
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


# In[11]:


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


# In[12]:


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


# In[13]:

numu_signal = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG==13 && MCTrackPurity>0.5' # numu CC signal definition
numu_nomuon = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG!=13 && MCTrackPurity>0.5' # numu CC signal definition
numu_nue = 'fidVol && MCNu_CCNC==0 && abs(MCNu_PDG)!=14 && MCTrackPurity>0.5' # e.g anti nu or nue
numu_antinu = 'fidVol && MCNu_CCNC==0 && MCNu_PDG==-14 && MCTrackPurity>0.5' # e.g anti nu or nue

numu_nc = 'fidVol && MCNu_CCNC==1 && MCTrackPurity>0.5' # nutral current
numu_ov = 'fidVol && !MCfidVol && MCTrackPurity>0.5' # out of fiducial

numu_cosmic = 'fidVol && MCTrackPurity<0.5' #low purity

#overlay_signals = ['numu_signal','numu_nue','numu_antinu','numu_nc','numu_ov','numu_cosmic']
overlay_signals = ['numu_cosmic','numu_ov','numu_nc','numu_antinu','numu_nue','numu_nomuon','numu_signal']


for x in sample_out:
    x.SetAlias('numu_signal',numu_signal)
    x.SetAlias('numu_nomuon',numu_nomuon)
    x.SetAlias('numu_nue',numu_nue)
    x.SetAlias('numu_antinu',numu_antinu)
    x.SetAlias('numu_nc',numu_nc)
    x.SetAlias('numu_ov',numu_ov)
    x.SetAlias('numu_cosmic',numu_cosmic)



# Load the global variables for access of functions
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)
#NuCC.printGlobal()
# In[12]:


#track_start_border_x = '(TrackStart_x_sce <(-1.55+1) || TrackStart_x_sce > (254.8-1))'
#track_end_border_x = '(TrackEnd_x_sce <(-1.55+1) || TrackEnd_x_sce > (254.8-1))'
#track_start_border_y = '(TrackStart_y_sce <(-115.53+1) || TrackStart_y_sce > (117.47-1))'
#track_end_border_y = '(TrackEnd_y_sce <(-115.53+1) || TrackEnd_y_sce > (117.47-1))'
#track_start_border_z = '(TrackStart_z_sce <(0.1+1) || TrackStart_z_sce > (1036.9-1))'
#track_end_border_z = '(TrackEnd_z_sce <(0.1+1) || TrackEnd_z_sce > (1039.9-1))'

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



# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
name = '_07'
side_right = 'right'
side_left = 'left'
folder = ''
#folder = 'All_cuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)
print 'All cuts'
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]' ,0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon'
name = '_01'
side_right = 'right'
side_left = 'left'
#folder = 'Precut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)
print 'preselection'
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8'
name = '_02'
side_right = 'right'
side_left = 'left'
#folder = 'Quality_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut'
name = '_03'
side_right = 'right'
side_left = 'left'
#folder = 'Crt_cuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8'
name = '_04'
side_right = 'right'
side_left = 'left'
#folder = 'TrackScore_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20'
name = '_05'
side_right = 'right'
side_left = 'left'
#folder = 'Tracklength_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78'
name = '_06'
side_right = 'right'
side_left = 'left'
#folder = 'TrackPID_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_MCC8_pub(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_MCC8_pub(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


# MCC8 plots
# Define signals as Marcos MCC8 analysis
'''
numu_signal_cont = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG==13 && MCTrackPurity>0.5 && !track_end_uncontained' # numu CC signal definition
numu_signal_uncont = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPDG==13 && MCTrackPurity>0.5 && track_end_uncontained' # numu CC signal definition

numu_nc_other = 'fidVol && MCfidVol && MCNu_CCNC==1 && MCTrackPDG!=211 && MCTrackPDG!=2212' # nutral current
numu_nc_pion = 'fidVol && MCfidVol && MCNu_CCNC==1 && MCTrackPDG==211' # nutral current
numu_nc_proton = 'fidVol && MCfidVol && MCNu_CCNC==1 && MCTrackPDG==2212' # nutral current

numu_cosmic_cont = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPurity<0.5 && !track_end_uncontained' #low purity
numu_cosmic_uncont = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==14 && MCTrackPurity<0.5 && track_end_uncontained' #low purity

numu_ov_cont = 'fidVol && !MCfidVol && !track_end_uncontained' # out of fiducial
numu_ov_uncont = 'fidVol && !MCfidVol && track_end_uncontained' # out of fiducial

numu_nue = 'fidVol && MCfidVol && MCNu_CCNC==0 && abs(MCNu_PDG)!=14' # e.g anti nu or nue
numu_antinu = 'fidVol && MCfidVol && MCNu_CCNC==0 && MCNu_PDG==-14' # e.g anti nu or nue

#tot_num_fidVol = num_fidVol[ext]+num_fidVol[dirt]+num_fidVol[overlay]
#overlay_signals = ['numu_signal_cont','numu_signal_uncont','numu_nue','numu_antinu','numu_nc_other','numu_nc_pion',\
#                   'numu_nc_proton', 'numu_ov_cont', 'numu_ov_uncont', 'numu_cosmic_cont', 'numu_cosmic_uncont']

overlay_signals = ['numu_cosmic_uncont','numu_cosmic_cont','numu_ov_uncont','numu_ov_cont','numu_nc_proton','numu_nc_pion',\
                   'numu_nc_other', 'numu_antinu', 'numu_nue', 'numu_signal_uncont', 'numu_signal_cont']


for x in sample_out:
    x.SetAlias('numu_signal_cont',numu_signal_cont)
    x.SetAlias('numu_signal_uncont',numu_signal_uncont)
    x.SetAlias('numu_nc_other',numu_nc_other)
    x.SetAlias('numu_nc_pion',numu_nc_pion)
    x.SetAlias('numu_nc_proton',numu_nc_proton)
    x.SetAlias('numu_cosmic_cont',numu_cosmic_cont)
    x.SetAlias('numu_cosmic_uncont',numu_cosmic_uncont)
    x.SetAlias('numu_ov_cont',numu_ov_cont)
    x.SetAlias('numu_ov_uncont',numu_ov_uncont)
    x.SetAlias('numu_nue',numu_nue)
    x.SetAlias('numu_antinu',numu_antinu)
    
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)
    


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
name = '_all_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_All_cuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78'
name = '_trackPID_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_TrackPID_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20'
name = '_tracklength_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_TrackLength_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8'
name = '_trackscore_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_TrackScore_cuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut'
name = '_crt_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_Crt_cuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8'
name = '_quality_cut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_Qualtity_cut/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:


cut = 'fidVol && muon'
name = '_precut'
side_right = 'right'
side_left = 'left'
folder = 'MCC8_Precuts/'
outputdir_png, outputdir_root,outputdir_pdf = NuCC.prepareOutput2(outputdir+folder)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)

NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vx_sce','EventWeight','#nu_{#mu}^{reco} vertex X position [cm]',-10,270,20,'NuVx_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vy_sce','EventWeight','#nu_{#mu}^{reco} vertex Y position [cm]',-120,120,20,'NuVy_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8_V(cut,'Nu_Vz_sce','EventWeight','#nu_{#mu}^{reco} vertex Z position [cm]',-50,1050,20,'NuVz_sce'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackMomMCS_mom','EventWeight','p_{#mu}^{reco} [GeV]',0,1.5,20,'TrackMom'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'cos(TrackTheta)','EventWeight','cos(#theta_{#mu}^{reco})',-1,1,20,'cosTheta'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPhi','EventWeight','#phi_{#mu}^{reco}',-3.15,3.15,20,'phi'+name,side_left)

NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackLength','EventWeight','Candidate Track Length [cm]',0,600,50,'TrackLength'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chiproton','EventWeight','Track PID proton',0,350,20,'PIDproton'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon','EventWeight','Track PID muon',0,60,20,'PIDmuon'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chipion','EventWeight','Track PID pion',0,60,20,'PIDpion'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut+'&& TrackPID_chimuon!=-999','TrackPID_chimuon/TrackPID_chiproton','EventWeight','Track PID muon/proton ration',0,0.25,20,'PIDmuprot_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackPID_chimuon/TrackPID_chipion','EventWeight','Track PID muon/pion ration',0.4,1.2,20,'PIDmupion_ration'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'nr_crthit_top','EventWeight','Number of CRT hits in top',0,4,4,'CRTTophit'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'crthit_vertex_zcut','EventWeight','Number of CRT hits upstream of vertex',0,4,4,'CRTVertexcut'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut,'TrackScore','EventWeight','Track Score',0.5,1,20,'TrackScore'+name,side_left)
NuCC_w.make_stacked_histo_weight_MCC8(cut,'NuScore','EventWeight','Topological score',0,1,20,'NuScore'+name,'right')
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NuTracks','EventWeight','Track Multiplicity',0,10,10,'NuTracks'+name,side_right)
NuCC_w.make_stacked_histo_weight_MCC8(cut + ' && 1','NumPfp','EventWeight','Particle Multiplicity',0,10,10,'NumPfp'+name,side_right)


# In[ ]:





# In[ ]:





# In[14]:


outputdir = outputdir+'Efficienies/'
outputdir_png, outputdir_root = NuCC.prepareOutput(outputdir)
NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root,outputdir_pdf)


# In[ ]:


nenner = 'numu_true'
cut = 'fidVol && muon'
zahler = cut+' && numu_signal'
NuCC_w.plot_eff_w(nenner,zahler,cut,'precut','precut')


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
zahler = cut+' && numu_signal'
name = 'all_cuts'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78'
zahler = cut+' && numu_signal'
name = 'trackPID_cut'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20'
zahler = cut+' && numu_signal'
name = 'trackLength_cut'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8'
zahler = cut+' && numu_signal'
name = 'trackScore_cut'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut'
zahler = cut+' && numu_signal'
name = 'crt_cuts'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8'
zahler = cut+' && numu_signal'
name = 'quality_cuts'
NuCC_w.plot_eff_w(nenner,zahler,cut,name,name)


# In[ ]:





# In[ ]:





# In[ ]:

'''



