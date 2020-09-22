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
from array import array
import collections


# In[ ]:


#get_ipython().system(u'jupyter nbconvert --to script RooFit_model.ipynb')


# In[2]:


def draw_adding():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()

    pot_latex = ROOT.TLatex(.10, .92,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/25.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    return prelim, pot_latex

def draw_sim():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    
    return prelim


# In[3]:


from array import array
mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
binnum = len(mom_bins) - 1
h_xsec_mom = ROOT.TH1F('h_xsec_mom','h_xsec_mom',binnum,array('f',mom_bins))


# In[4]:


# calculate total flux integrated cross section:
cut = 'fidVol && muon && TrackLength>8 && crt_tom_cut && TrackScore>0.8                && TrackLength>20 && TrackPID_chiproton>78 && NuScore>0.1' #\
                #&& MCle_Energy>0 && MCle_Energy<2.5 && TrackMomMCS_mom>0 && TrackMomMCS_mom<2.5'


# In[5]:


# initialte ROOT default canvas
#ROOT.gStyle.SetOptStat(0)
#c1 = ROOT.TCanvas("c1","c1",1600,1200)
#c1.SetGrid(1)
#c1.SetLeftMargin(0.14)
#c1.SetRightMargin(0.05)
#c1.SetBottomMargin(0.14)


# In[6]:


#load data
#inputdir = '/home/tmettler/Desktop/ub_data/mcc9.1/v08_00_00_33/V08_00_00_35/fitting/tutorial/RooFit-tutorial/hists/'
inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/xsec_momentum_rooFit/'
f_mom = ROOT.TFile.Open(inputdir+"xsec_momentum_fit.root", 'read')

h_true = f_mom.Get('mom_truth')
h_temp = []
for i in range(binnum):
    h_temp.append(f_mom.Get('mom_reco_'+str(i)))
h_data = f_mom.Get('data_reco') #data_reco mc_reco
h_background = f_mom.Get('mom_bkg_reco')
h_signal = f_mom.Get('mom_truth_sig')

#h_true.Draw()
#c1.Draw()


# # start fitting using RooFit and HistFactory class

# In[7]:


#create the histfactory model
outputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/xsec_momentum_rooFit/' 

meas = ROOT.RooStats.HistFactory.Measurement("ProfiledUnfolding", "ProfiledUnfolding"); 
meas.SetOutputFilePrefix(outputdir+"workspaces/ProfiledUnfolding/tut");
try:
    os.stat(outputdir+"workspaces/ProfiledUnfolding/")
except:
    os.mkdir(outputdir+"workspaces/ProfiledUnfolding/")
meas.SetExportOnly(1);

# NOT SURE ABOUT THIS: This scales the histogram content, which already includes lumi, so set to 1
meas.SetLumi(1.00);
meas.SetLumiRelErr(0.0002);


# In[8]:


#CHECK WHATS THIS: create the SR
chan_ee = ROOT.RooStats.HistFactory.Channel("mom");
chan_ee.SetData(h_data);


# In[9]:


for i in range(binnum):
    mom_XS = ROOT.RooStats.HistFactory.NormFactor()
    mom_XS.SetName(ROOT.Form('Normalization_mom_'+str(i)));
    mom_XS.SetHigh(h_true.GetBinContent(i+1)*50000) # maximum value it can take
    mom_XS.SetLow(0) # minimum value it can take
    mom_XS.SetVal(h_true.GetBinContent(i+1)) #;//*rndm.Gaus(1,0.05)); // nominal value (randomize the initial value a bit)
    meas.AddPOI(ROOT.Form('Normalization_mom_'+str(i)))
    
    # add the signal samples
    sample = ROOT.RooStats.HistFactory.Sample(ROOT.Form('mom_acceptance_template_'+str(i)))
    sample.SetNormalizeByTheory(False)        
    #sample.SetNormalizeByTheory(True)
    sample.SetHisto(h_temp[i])
    sample.AddNormFactor(mom_XS)
    
    
    signal_shape = ROOT.RooStats.HistFactory.HistoSys()
    signal_shape.SetName("signal_shape")
    h_high = h_temp[i].Clone()
    h_high.Scale(1.1)
    
    h_low = h_temp[i].Clone()
    h_low.Scale(0.9)
    
    signal_shape.SetHistoHigh( h_high )
    signal_shape.SetHistoLow( h_low )
    sample.AddHistoSys( signal_shape )
    
    chan_ee.AddSample(sample)


# In[10]:


sample_nonFid = ROOT.RooStats.HistFactory.Sample("ZnonFid");
sample_nonFid.SetHisto(h_background);

bkg_shape = ROOT.RooStats.HistFactory.HistoSys()
bkg_shape.SetName('bkg_shape')
h_high = h_background.Clone()
h_high.Scale(1.05)

h_low = h_background.Clone()
h_low.Scale(0.95)
    
bkg_shape.SetHistoHigh( h_high )
bkg_shape.SetHistoLow( h_low )
sample_nonFid.AddHistoSys( bkg_shape )

chan_ee.AddSample(sample_nonFid);


# In[ ]:


# add the single region to the measurement
meas.AddChannel(chan_ee)

# make the workspace
ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)


# In[ ]:


print 'done'


# In[ ]:





# In[ ]:





# In[ ]:




