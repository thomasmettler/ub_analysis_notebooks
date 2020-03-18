#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
helper = imp.load_source('helpfunction','/home/tmettler/Desktop/uBoone/do_plots/helpfunction.py')
globale = imp.load_source('globale','/home/tmettler/Desktop/uBoone/do_plots/globale.py')
#import /home/tmettler/Desktop/uBoone/do_plots/helpfunction as helper
import collections
NuCC = imp.load_source('NuCC_function','/home/tmettler/Desktop/uBoone/do_plots/NuCC_function.py')
NuCC_w = imp.load_source('NuCCWeight_function','/home/tmettler/Desktop/uBoone/do_plots/NuCCWeight_function.py')


# In[13]:


#get_ipython().system(u'jupyter nbconvert --to script numuCCana_muon_cut_V33.ipynb')


# In[3]:


# initialte ROOT default canvas
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# # CRTinBNB tagger, Overlay or MC

# In[4]:


# Load input files
#inputdir = '/home/thomasm/numuCC/V33/10kevents/'
inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/'
outputdir = inputdir+'muon_cut/' 
outputdir_png, outputdir_root = NuCC.prepareOutput(outputdir)

#filename_overlay = 'NuCCana_Overlay_V18_V4.root'
filename_overlay = 'NuCCana_overlay_V26_weight.root'
#filename_overlay = 'NuCCana_overlay_v25_weighted_all.root'
filename_data = 'NuCCana_data_V22.root'
filename_ext = 'NuCCana_ext_V25.root'
filename_dirt = 'NuCCana_dirt_V26_weight.root'
    
tree_name = 'numuCCAna'


# In[5]:


#Open all the trees of the four files (data, ext, dirt, overlay)

data, ext, dirt, overlay = NuCC.openTrees(inputdir, filename_data, filename_ext, filename_dirt, filename_overlay, tree_name)
NuCC.printNumberOfEntries(data,ext,dirt,overlay)

pot_overlay = NuCC.getPOT(inputdir,filename_overlay,tree_name)
pot_dirt =  NuCC.getPOT(inputdir,filename_dirt,tree_name)

#tree_pot = ROOT.TChain(tree_name + "/pottree","pottree")
#tree_pot.Add( inputdir + filename_overlay)
#tot_pot = 0.0
#for entry in tree_pot:
#    tot_pot = tot_pot+entry.pot
#print 'tot_pot', tot_pot

#V25 files
pot_data =     8.649e+18  # best with tor875
data_trigger = 1854495.0 #4743794 # 1987072.0 # E1DCNT_wcut
ext_trigger =  18997529.0  #2120135 #5685315 # EXT

#V17 files
#pot_data =     8.039e+18 # tor860_wcut
#data_trigger = 2326285.0 # E1DCNT
#ext_trigger =  25882118.0  #2120135 #5685315 # EXT
#pot_overlay = 95886250000000000.000000

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
    filename_overlay = filename_overlay+'out33.root'
    #filename_overlay = 'NuCCana_overlay_points_1kev.rootout2.root'
    filename_data = filename_data+'out33.root'
    filename_ext = filename_ext+'out33.root'
    filename_dirt = filename_dirt+'out33.root'

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
    #x.SetAlias('muon','(track_key == key_muon)')
    
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


# In[8]:


for x in sample_out:
    x.SetAlias('muon','(track_key == key_muon)')
    #x.SetAlias('muon','(key_muprotratio==track_key)')
    
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
    #num_fidVol[x] = x.GetEntries('fidVol && muon')*scale[x]
    
#tot_num_fidVol = num_fidVol[ext]+num_fidVol[dirt]+num_fidVol[overlay]
#overlay_signals = {'numu_signal','numu_nomu','numu_lowpur','numu_nc','numu_ov','numu_other'}


# In[9]:


# Load the global variables for access of functions

NuCC.loadGlobal(data,ext,dirt,overlay,data_out,ext_out,dirt_out,overlay_out,scale,scale_out,tot_num_fidVol,overlay_signals,sample,sample_out, name,name_out, outputdir_png, outputdir_root)
#NuCC.printGlobal()


# In[10]:


# initialte ROOT default canvas
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# In[11]:


track_start_border_x = '(TrackStart_x_sce <(-1.55+1) || TrackStart_x_sce > (254.8-1))'
track_end_border_x = '(TrackEnd_x_sce <(-1.55+1) || TrackEnd_x_sce > (254.8-1))'
track_start_border_y = '(TrackStart_y_sce <(-115.53+1) || TrackStart_y_sce > (117.47-1))'
track_end_border_y = '(TrackEnd_y_sce <(-115.53+1) || TrackEnd_y_sce > (117.47-1))'
track_start_border_z = '(TrackStart_z_sce <(0.1+1) || TrackStart_z_sce > (1036.9-1))'
track_end_border_z = '(TrackEnd_z_sce <(0.1+1) || TrackEnd_z_sce > (1039.9-1))'

track_end_uncontained = '(' + track_end_border_x + ' || ' + track_end_border_y + ' || ' + track_end_border_z+ ')'


data.SetAlias("track_end_uncontained",track_end_uncontained)
ext.SetAlias("track_end_uncontained",track_end_uncontained)
overlay.SetAlias("track_end_uncontained",track_end_uncontained)
dirt.SetAlias("track_end_uncontained",track_end_uncontained)
data_out.SetAlias("track_end_uncontained",track_end_uncontained)
ext_out.SetAlias("track_end_uncontained",track_end_uncontained)
overlay_out.SetAlias("track_end_uncontained",track_end_uncontained)
dirt_out.SetAlias("track_end_uncontained",track_end_uncontained)

data.SetAlias("crt_cut","(abs(crtt0_time-(69000-crt_trig_corr_med)/1000-4)<0.8 || crtt0_time==-1)")
ext.SetAlias("crt_cut","(abs(crtt0_time-(69000-crt_trig_corr_med)/1000-3.57+3.195-4)<0.8 || crtt0_time==-1)")
overlay.SetAlias("crt_cut","(abs(crtt0_time-4)<0.8 || crtt0_time==-1)")
dirt.SetAlias("crt_cut","(abs(crtt0_time-4)<0.8 || crtt0_time==-1)")
data_out.SetAlias("crt_cut","(abs(crtt0_time-(69000-crt_trig_corr_med)/1000-4)<0.8 || crtt0_time==-1)")
ext_out.SetAlias("crt_cut","(abs(crtt0_time-(69000-crt_trig_corr_med)/1000-3.57+3.195-4)<0.8 || crtt0_time==-1)")
overlay_out.SetAlias("crt_cut","(abs(crtt0_time-4)<0.8 || crtt0_time==-1)")
dirt_out.SetAlias("crt_cut","(abs(crtt0_time-4)<0.8 || crtt0_time==-1)")

crt_tom_cut = 'nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut'

data.SetAlias("crt_tom_cut",crt_tom_cut)
ext.SetAlias("crt_tom_cut",crt_tom_cut)
overlay.SetAlias("crt_tom_cut",crt_tom_cut)
dirt.SetAlias("crt_tom_cut",crt_tom_cut)
data_out.SetAlias("crt_tom_cut",crt_tom_cut)
ext_out.SetAlias("crt_tom_cut",crt_tom_cut)
overlay_out.SetAlias("crt_tom_cut",crt_tom_cut)
dirt_out.SetAlias("crt_tom_cut",crt_tom_cut)


# In[12]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
name = '_muon_cut'
start = time.time()
#NuCC_w.make_stacked_histo_weight(cut + ' && TrackTheta!=-9','cos(TrackTheta)','EventWeight','cos(theta)',-1,1,20,'costhetaw'+name)
#end = time.time()
#print(end - start)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vx_sce','EventWeight','Reco. Nu vertex X position [cm]',-10,270,20,'NuVx_sce'+name)
end = time.time()
print(end - start)
#NuCC.make_stacked_histo_plot_out(cut + ' && TrackTheta!=-9','cos(TrackTheta)','cos(theta)',-1,1,20,'costheta'+name)


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
name = '_muon_cut'
start = time.time()
NuCC_w.make_stacked_histo_weight(cut + ' && TrackTheta!=-9','TrackPhi','EventWeight','Track Phi',-3.141,3.141,20,'phi'+name)
end = time.time()
print(end - start)
NuCC_w.make_stacked_histo_weight(cut + ' && 1','TrackMomMCS_mom','EventWeight','TrackMomMCS_mom',0,1.5,20,'TrackMomMCS_mom'+name)
end = time.time()
print(end - start)
#NuCC_w.make_stacked_histo_weight(cut,'Nu_Vx_sce','EventWeight','Reco. Nu vertex X position [cm]',-10,270,20,'NuVx_sce'+name)
#end = time.time()
#print(end - start)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vy_sce','EventWeight','Reco. Nu vertex Y position [cm]',-120,120,20,'NuVy_sce'+name)
end = time.time()
print(end - start)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vz_sce','EventWeight','Reco. Nu vertex Z position [cm]',-10,1040,20,'NuVz_sce'+name)
end = time.time()
print(end - start)


# In[ ]:


NuCC_w.printEff_w('fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1')


# In[74]:


NuCC_w.printEff_out('fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1')


# In[49]:


title = 'Scale test'
name = 'CRTX'
xstart = -1.1
xend = 1
xbin = 1
cut = '1 && '
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
h_pos = ROOT.TH1F("h_pos",title,xbin,xstart,xend)
h_pos_ext = ROOT.TH1F("h_pos_ext",title,xbin,xstart,xend)
h_pos_MC = ROOT.TH1F("h_pos_MC",title,xbin,xstart,xend)
h_pos_dirt = ROOT.TH1F("h_pos_dirt",title,xbin,xstart,xend)
h_pos_all = ROOT.TH1F("h_pos_all",title,xbin,xstart,xend)
data_out.Draw('cos(TrackTheta)>>h_pos',cut+'1')
ext_out.Draw('cos(TrackTheta)>>h_pos_ext',cut+'1')
overlay_out.Draw('cos(TrackTheta)>>h_pos_MC',cut+'1')
dirt_out.Draw('cos(TrackTheta)>>h_pos_dirt',cut+'1')
h_pos.SetXTitle("CRT hit - MC vertex in x [cm]")
h_pos.SetYTitle("Number of entries")
legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
legend.AddEntry(h_pos,"all","l");
legend.AddEntry(h_pos_ext,"Cosmic","l");
legend.AddEntry(h_pos_MC,"MC","l");
legend.AddEntry(h_pos_dirt,"dirt","l");

h_pos_ext.Sumw2()
h_pos_ext.Scale(scale[ext])
h_pos_MC.Sumw2()
h_pos_MC.Scale(scale[overlay])
h_pos_dirt.Sumw2()
h_pos_dirt.Scale(scale[dirt])
h_pos.SetLineColor(1)
h_pos_ext.SetLineColor(3)
h_pos_MC.SetLineColor(2)
h_pos_dirt.SetLineColor(4)
h_pos_dirt.SetLineWidth(3)
h_pos.SetLineWidth(3)
h_pos_ext.SetLineWidth(3)
h_pos_MC.SetLineWidth(3)
h_pos.Draw()
h_pos_ext.Draw('same')
h_pos_MC.Draw('same')
h_pos_dirt.Draw('same')
h_pos_all.Add(h_pos_ext,1)
h_pos_all.Add(h_pos_dirt,1)
h_pos_all.Add(h_pos_MC,1)
h_pos_all.Draw('same')
legend.Draw('same');
c1.Draw()
c1.SaveAs(outputdir_png + "h_pos_"+name+".png")
c1.SaveAs(outputdir_root + "h_pos_"+name+".root")


# In[ ]:



t = ROOT.TH1F('t','t',20,0,10)
t.Fill(1,0.8)
t.Fill(1.4,0.9)
t.Fill(2.4,0.6)
t.Fill(2,0.9)
t.Fill(3,0.9)
t.Fill(8,1)
t.Scale(2)
print 'Entries: ', t.GetEntries()
print 'Entries: ', t.GetSumOfWeights()

t.Draw()
c1.Draw()


# In[76]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
NuCC.printEff_out(cut)


# In[ ]:





# In[153]:


cut = 'fidVol && muon'
NuCC_w.printEff_out(cut)


# In[13]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
NuCC_w.printEff_out(cut)


# In[138]:


NuCC.pdg_content33_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1')


# In[16]:


cut = 'fidVol && muon'
NuCC_w.printEff_out(cut)


# In[29]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
name = '_muon_cut'
NuCC_w.make_stacked_histo_weight(cut + ' && TrackTheta!=-9','cos(TrackTheta)','EventWeight','cos(theta)',-1,1,20,'costhetaw'+name)
NuCC_w.make_stacked_histo_weight(cut + ' && TrackTheta!=-9','TrackPhi','EventWeight','Track Phi',-3.141,3.141,20,'phi'+name)
NuCC_w.make_stacked_histo_weight(cut + ' && 1','TrackMomMCS_mom','EventWeight','TrackMomMCS_mom',0,1.5,20,'TrackMomMCS_mom'+name)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vx_sce','EventWeight','Reco. Nu vertex X position [cm]',-10,270,20,'NuVx_sce'+name)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vy_sce','EventWeight','Reco. Nu vertex Y position [cm]',-120,120,20,'NuVy_sce'+name)
NuCC_w.make_stacked_histo_weight(cut,'Nu_Vz_sce','EventWeight','Reco. Nu vertex Z position [cm]',-10,1040,20,'NuVz_sce'+name)

#NuCC.make_stacked_histo_plot_out(cut + ' && TrackTheta!=-9','cos(TrackTheta)','cos(theta)',-1,1,20,'costheta'+name)


# In[ ]:





# In[57]:





# In[110]:


nenner = 'numu_true'
zahler = 'fidVol && muon && numu_signal'
cut = 'fidVol && muon'
NuCC_w.plot_eff_out(nenner,zahler,cut,'initial','initial')


# In[ ]:


nenner = 'numu_true'
zahler = 'fidVol && muon && numu_signal'
cut = 'fidVol && muon'
NuCC_w.plot_eff_out(nenner,zahler,cut,'initial','initial')


# In[155]:


nenner = 'fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackPID_chiproton>78 && NuScore>0.1 && numu_true'
zahler = 'fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1 && numu_true'
cut = 'fidVol && muon && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1'
NuCC.plot_eff_outV2(nenner,zahler,cut,'tracklength','tracklength')


# In[ ]:


NuCC.pdg_content33_out('fidVol && muon && TrackLength>32 && nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut')


# In[ ]:


crt_tom_cut = 'nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut'

data.SetAlias("crt_tom_cut",crt_tom_cut)
ext.SetAlias("crt_tom_cut",crt_tom_cut)
overlay.SetAlias("crt_tom_cut",crt_tom_cut)
dirt.SetAlias("crt_tom_cut",crt_tom_cut)
data_out.SetAlias("crt_tom_cut",crt_tom_cut)
ext_out.SetAlias("crt_tom_cut",crt_tom_cut)
overlay_out.SetAlias("crt_tom_cut",crt_tom_cut)
dirt_out.SetAlias("crt_tom_cut",crt_tom_cut)


# In[ ]:


NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8')


# In[ ]:


#NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8 && TrackLength>20')


# In[ ]:


NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8 && TrackLength>32')
# 20cm = 150MeV MC momentum, 32 = 180 MeV, 25cm = 160MeV, 15cm = 130MeV, 10 = 115MeV, 8cm=100MeV


# In[ ]:


NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<24')


# In[ ]:


NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30')


# In[ ]:


NuCC.printEff_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1')


# In[ ]:


NuCC.pdg_content33_out('fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1')


# In[ ]:


cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
name = '_muon_cut'
NuCC.make_stacked_histo_plot_out(cut + ' && TrackTheta!=-9','cos(TrackTheta)','cos(theta)',-1,1,20,'costheta'+name)
NuCC.make_stacked_histo_plot_out(cut,'Nu_Vz_sce','Reco. Nu vertex Z position [cm]',-20,1020,50,'NuVz_sce'+name)
#NuCC.make_stacked_histo_plot_out(cut,'TrackTheta','Total NuScore of event',0,1,20,'NuScore'+name)


# In[ ]:


nenner = 'numu_true'
zahler = 'fidVol && muon && numu_signal'
cut = 'fidVol && muon'
NuCC.plot_eff_outV2(nenner,zahler,cut,'initial','initial')


# In[ ]:


nenner = 'numu_true'
zahler = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30 && numu_true'
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30'
NuCC.plot_eff_outV2(nenner,zahler,cut,'crt_muon_tot','crt_muon_tot')


# In[ ]:


nenner = 'fidVol && muon && numu_true && crt_tom_cut'
zahler = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30 && numu_true'
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30'
NuCC.plot_eff_outV2(nenner,zahler,cut,'crt_muon_only','crt_muon_only')


# In[ ]:


nenner = 'fidVol && muon && numu_true'
zahler = 'fidVol && muon && TrackLength>8 && crt_tom_cut && numu_true'
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut'
NuCC.plot_eff_outV2(nenner,zahler,cut,'crt_only','crt_only')


# In[ ]:


nenner = 'numu_true'
zahler = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1 && numu_true'
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chimuon<30 && NuScore>0.1'
NuCC.plot_eff_outV2(nenner,zahler,cut,'all','all')


# In[ ]:


nenner = 'numu_true && fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30'
zahler = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30 && NuScore>0.1 && numu_true'
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>32 && TrackPID_chimuon<30 && NuScore>0.1'
NuCC.plot_eff_outV2(nenner,zahler,cut,'nu_only','nu_only')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# initialte ROOT default canvas
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# In[ ]:


value = 'TrackLength'
xbins = 30
xstart = 0
xend = 30
step = (xend-xstart)/xbins
title = 'Track length [cm]'

h_eff =  ROOT.TH1F("h_eff",title,xbins,xstart,xend)
h_pur = ROOT.TH1F("h_pur",title,xbins,xstart,xend)
h_effpur = ROOT.TH1F("h_effpur",title,xbins,xstart,xend)
for i in range(xbins): 
    cut = 'muon && fidVol && TrackPID_chimuon<30 && NuScore>0.1 && crt_tom_cut && TrackScore>0.8 && ('+value+'>('+str(xstart)+'+'+str(step)+'*'+str(i)+'))'
    #print cut
    pur = overlay_out.GetEntries(cut+' && numu_signal')*scale[overlay]*100.0/(NuCC.getTotNum_out(cut)+0.00001)
    eff = overlay_out.GetEntries(cut+' && numu_signal')*100.0/overlay_out.GetEntries("numu_true")
    h_eff.SetBinContent(i+1,eff)
    h_pur.SetBinContent(i+1,pur)
    h_effpur.SetBinContent(i+1,eff*pur/100.0)
    #print eff ,pur, eff*pur,i

legend = ROOT.TLegend(0.15,0.75,0.48,0.9);
legend.AddEntry(h_eff,"efficiency","f");
legend.AddEntry(h_pur,"purity","f");
legend.AddEntry(h_effpur,"eff * pur","f");

h_eff.SetLineColor(1)
h_eff.SetMaximum(100)
h_eff.SetMinimum(0)
h_pur.SetLineColor(2)
h_effpur.SetLineColor(3)
h_pur.SetLineWidth(3)
h_eff.SetLineWidth(3)
h_effpur.SetLineWidth(3)
h_eff.Draw()
legend.Draw()
h_pur.Draw('same')
h_effpur.Draw('same')
h_eff.Draw('same')
c1.SetLogy(0)
c1.Draw()
c1.SaveAs(outputdir_png+ value +".png");
c1.SaveAs(outputdir_root+ value+ ".root");

print 'Variation of cut on: '+value
print 'Maximum value: ', h_effpur.GetMaximum()
print 'Bin: ', h_effpur.GetMaximumBin(),' -> Value: ', xstart+(h_effpur.GetMaximumBin()-1)*step


# In[ ]:


value = 'TrackPID_chimuon'
xbins = 30
xstart = 0
xend = 30
step = (xend-xstart)/xbins
title = 'Track PID chi square'

h_eff =  ROOT.TH1F("h_eff",title,xbins,xstart,xend)
h_pur = ROOT.TH1F("h_pur",title,xbins,xstart,xend)
h_effpur = ROOT.TH1F("h_effpur",title,xbins,xstart,xend)
for i in range(xbins): 
    cut = 'muon && fidVol && crt_tom_cut && TrackScore>0.8 &&         TrackLength>20 && NuScore>0.1 && ('+value+'<('+str(xstart)+'+'+str(step)+'*'+str(i)+'))'
    #print cut
    pur = overlay_out.GetEntries(cut+' && numu_signal')*scale[overlay]*100.0/(NuCC.getTotNum_out(cut)+0.00001)
    eff = overlay_out.GetEntries(cut+' && numu_signal')*100.0/overlay_out.GetEntries("numu_true")
    h_eff.SetBinContent(i+1,eff)
    h_pur.SetBinContent(i+1,pur)
    h_effpur.SetBinContent(i+1,eff*pur/100.0)
    #print eff ,pur, eff*pur,i

legend = ROOT.TLegend(0.15,0.75,0.48,0.9);
legend.AddEntry(h_eff,"efficiency","f");
legend.AddEntry(h_pur,"purity","f");
legend.AddEntry(h_effpur,"eff * pur","f");

h_eff.SetLineColor(1)
h_eff.SetMaximum(100)
h_eff.SetMinimum(0)
h_pur.SetLineColor(2)
h_effpur.SetLineColor(3)
h_pur.SetLineWidth(3)
h_eff.SetLineWidth(3)
h_effpur.SetLineWidth(3)
h_eff.Draw()
legend.Draw()
h_pur.Draw('same')
h_effpur.Draw('same')
h_eff.Draw('same')
c1.SetLogy(0)
c1.Draw()
c1.SaveAs(outputdir_png+ value +".png");
c1.SaveAs(outputdir_root+ value+ ".root");

print 'Variation of cut on: '+value
print 'Maximum value: ', h_effpur.GetMaximum()
print 'Bin: ', h_effpur.GetMaximumBin(),' -> Value: ', xstart+(h_effpur.GetMaximumBin()-1)*step


# In[ ]:


value = 'NuScore'
xbins = 20
xstart = 0
xend = 0.5
step = (xend-xstart)/xbins
title = 'Track length [cm]'

h_eff =  ROOT.TH1F("h_eff",title,xbins,xstart,xend)
h_pur = ROOT.TH1F("h_pur",title,xbins,xstart,xend)
h_effpur = ROOT.TH1F("h_effpur",title,xbins,xstart,xend)
for i in range(xbins): 
    cut = 'muon && fidVol && crt_tom_cut && TrackScore>0.8 && ('+value+'>('+str(xstart)+'+'+str(step)+'*'+str(i)+'))'
    #print cut
    pur = overlay_out.GetEntries(cut+' && numu_signal')*scale[overlay]*100.0/(NuCC.getTotNum_out(cut)+0.00001)
    eff = overlay_out.GetEntries(cut+' && numu_signal')*100.0/overlay_out.GetEntries("numu_true")
    h_eff.SetBinContent(i+1,eff)
    h_pur.SetBinContent(i+1,pur)
    h_effpur.SetBinContent(i+1,eff*pur/100.0)
    #print eff ,pur, eff*pur,i

legend = ROOT.TLegend(0.15,0.75,0.48,0.9);
legend.AddEntry(h_eff,"efficiency","f");
legend.AddEntry(h_pur,"purity","f");
legend.AddEntry(h_effpur,"eff * pur","f");

h_eff.SetLineColor(1)
h_eff.SetMaximum(100)
h_eff.SetMinimum(0)
h_pur.SetLineColor(2)
h_effpur.SetLineColor(3)
h_pur.SetLineWidth(3)
h_eff.SetLineWidth(3)
h_effpur.SetLineWidth(3)
h_eff.Draw()
legend.Draw()
h_pur.Draw('same')
h_effpur.Draw('same')
h_eff.Draw('same')
c1.SetLogy(0)
c1.Draw()
c1.SaveAs(outputdir_png+ value +".png");
c1.SaveAs(outputdir_root+ value+ ".root");

print 'Variation of cut on: '+value
print 'Maximum value: ', h_effpur.GetMaximum()
print 'Bin: ', h_effpur.GetMaximumBin(),' -> Value: ', xstart+(h_effpur.GetMaximumBin()-1)*step


# In[ ]:





# In[ ]:





# In[ ]:




