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


#!jupyter nbconvert --to script RooFit_fit.ipynb


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
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.05)
c1.SetBottomMargin(0.14)


# In[6]:


#load data
inputdir = '/home/tmettler/Desktop/ub_data/mcc9.1/v08_00_00_33/V08_00_00_35/fitting/tutorial/RooFit-tutorial/hists/'
f_mom = ROOT.TFile.Open(inputdir+"xsec_histos.root", 'read')

h_true = f_mom.Get('mom_truth_sig')
h_temp = []
for i in range(binnum):
    h_temp.append(f_mom.Get('mom_reco_'+str(i)))
h_data = f_mom.Get('data_reco')
h_background = f_mom.Get('mom_truth_sig')
h_signal = f_mom.Get('mom_bkg_reco')

#h_true.Draw()
#c1.Draw()


# # start fitting using RooFit and HistFactory class

# # Start fitting

# In[7]:



outputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted/xsec_momentum_rooFit/' 

f= ROOT.TFile(ROOT.Form(outputdir+"workspaces/ProfiledUnfolding/tut_combined_ProfiledUnfolding_model.root"));
w = f.Get("combined");
#if (!w):
#    print "ERROR::Workspace doesn't exist! Check file name"


# In[8]:


mc = w.obj("ModelConfig");
data = w.data("obsData");


# In[9]:


# Configure MINUIT
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2");
ROOT.Math.MinimizerOptions.SetDefaultStrategy(0);
ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(1);


# In[10]:


# Put the NPs and POIs into one set to send to the NLL
params = ROOT.RooArgSet(mc.GetNuisanceParameters(),mc.GetParametersOfInterest());


# In[11]:


# Build the NLL
nll =  mc.GetPdf().createNLL(data, ROOT.RooFit.Constrain(params),                ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),                ROOT.RooFit.Offset(1))


# In[12]:


def minimize(fcn, save=0, retry_mode=3,ret_status=0):
  #Grab the default minimizer options
  printLevel = ROOT.Math.MinimizerOptions.DefaultPrintLevel()
  msglevel = ROOT.RooMsgService.instance().globalKillBelow();
  if (printLevel < 0):
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL);
  
  strat = ROOT.Math.MinimizerOptions.DefaultStrategy();
  save_strat = strat;

  #Configure the minimizer
  minim = ROOT.RooMinimizer(fcn);
  minim.optimizeConst(2);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);
  #minim.setProfile(1);
  minim.setMinimizerType(str(ROOT.Math.MinimizerOptions.DefaultMinimizerType()));

  #Do the minimization here!
  status = minim.minimize(str(ROOT.Math.MinimizerOptions.DefaultMinimizerType()), str(ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()))
  
  #In case the fit failed, impose a semi-smart retry strategy

  #Here, increase the MINUIT strategy from 0 through 2 until the fit succeeds (or give up after 2)
  if (retry_mode == 0):
    # up the strategy
    if (status != 0 and status != 1 and strat < 2):
      strat+=1
      print "Fit failed with status ", status, ". Retrying with strategy ", strat
      minim.setStrategy(strat);
      status = minim.minimize(str(ROOT.Math.MinimizerOptions.DefaultMinimizerType()), str(ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()))

    
    # up the strategy
    if (status != 0 and status != 1 and strat < 2):
      strat+=1
      print "Fit failed with status ", status, ". Retrying with strategy ", strat
      minim.setStrategy(strat);
      status = minim.minimize(str(ROOT.Math.MinimizerOptions.DefaultMinimizerType()), str(ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()));
    

  #Here, just retry N times with the default strategy
  else:
    for i in range(retry_mode):
      if (status == 0 or status == 1):
        break
      print "Fit failed with status " , status , ". Retrying with strategy " , strat
      minim.setStrategy(strat);
      status = minim.minimize(ROOT.Math.MinimizerOptions.DefaultMinimizerType().str(), ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo().str())
  
  
  #Reset the global configuration to the previous one
  if (printLevel < 0):
    ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)
  ROOT.Math.MinimizerOptions.SetDefaultStrategy(save_strat)
  

  #Save the RooFitResult if asked to
  fitresult = 0
  if (save):
    fitresult = minim.save(ROOT.Form("fitresult_"+str(fcn.GetName())), ROOT.Form("fitresult_"+str(fcn.GetName())))
  

  #Save the fit status if the pointer is valid
  if (ret_status):
    ret_status = status;

  #Print the status of the fit
  if (status != 0 and status != 1):
    print  "Fit failed with status " , status
  #else cout << "Fit succeeded with status " << status << endl;      

  return fitresult


# In[ ]:





# In[ ]:


#Do the minimization and save the RooFitResult
#Definition in macros/minimize.C:
#minimize = ROOT.RooFitResult(fcn, save=0, retry_mode=3, ret_status=0)
result = minimize(nll, 1, 0);
result.SetName("result");
result.SaveAs(ROOT.Form(outputdir+"workspaces/ProfiledUnfolding/fit.root"))

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




