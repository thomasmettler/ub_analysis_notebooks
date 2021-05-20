#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import stuff you eventually need
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
from array import array


# In[2]:


# to convert the notebook into a python script...
#!jupyter nbconvert --to script FF_event_rate_tutorial.ipynb


# In[3]:


# define input (file) and output path (figures, files ..)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputdir = 'path to your input file' # please give here path to inputfile
outputdir = inputdir+'your_output/' 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# make output dir if not existing
try:
    os.stat(outputdir)
except:
    os.mkdir(outputdir)
# define root output file for histograms.
# Careful: RECREATE deletes existing files with same name
RootFile = ROOT.TFile(outputdir+"tutorial_histo.root","RECREATE");


# In[4]:


# for later if syst. are precalculated
# in the end change to calc_new = 1 and num_uni_override =1000
calc_new = 1 # in case the histograms are precalculated, the calculation takes long!
num_uni_override = 10 # change this to 1000 when ready with everything and calc_new to 0
path_sys = inputdir+'your precalculation/'


# In[5]:


# load TTree with toy data
input_file = inputdir+'toyData.root'

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# load tree as a TChain
data = 
data.Add(input_file)
print data.GetEntries()
#same for meta data
meta_data = 
meta_data.Add(input_file)
print meta_data.GetEntries()


# In[6]:


# Read meta data
for evt in meta_data:
    flux = evt.flux
    nTargets = evt.nTargets
    potData = evt.potData
    potMC = evt.potMC
    nUniverses = evt.nUniverses

# if num of universes wants to be limited...
nUniverses = num_uni_override

print 'Flux: ', flux
print 'Number of targets: ', nTargets
print 'POT data: ', potData
print 'POT MC: ', potMC
print 'Number universes: ', nUniverses


# In[7]:


# initialte ROOT default canvas for plotting
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.18)
c1.SetBottomMargin(0.14)


# In[8]:


# define some functions for 'preliminary' in plots etc.
def draw_adding():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);
    prelim.Draw()

    pot_latex = ROOT.TLatex(.10, .92,'Accumulated POT: '+str(potData)) 
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
    prelim.Draw()
    
    return prelim


# In[9]:


#function takes np array and returns the THXF
def arr2plot(this_arr):
    this_array = np.copy(this_arr)

    n_bins = this_array.shape
    dim = len(n_bins)
    if dim == 1:
        h_this = ROOT.TH1F('h_this',"Title",n_bins[0],0,n_bins[0])
        for i in range(n_bins[0]):
            h_this.SetBinContent(i+1,this_array[i])
    if dim == 2:
        h_this = ROOT.TH2F('h_this',"Title",n_bins[0],0,n_bins[0],n_bins[1],0,n_bins[1])
        for i in range(n_bins[0]):
            for j in range(n_bins[0]):
                h_this.SetBinContent(i+1,j+1,this_array[i][j])
    main_plot(h_this)        
    return h_this

#function takes THXF and the dimension and returns np array
def plot2arr(h_this_tmp, dim):
    h_this = h_this_tmp.Clone()
    if dim == 1:
        n_bins = h_this.GetNbinsX()
        arr = np.zeros((n_bins))
        for i in range(n_bins):
            arr[i] = h_this.GetBinContent(i+1)
    if dim == 2:
        n_bins = np.zeros((2))
        n_bins[0] = int(h_this.GetNbinsX())
        n_bins[1] = int(h_this.GetNbinsY())
        n_bins.astype(np.int64)
        arr = np.zeros((np.int64(n_bins[0]),np.int64(n_bins[1])))
        for i in range(np.int64(n_bins[0])):
            for j in range(np.int64(n_bins[1])):
                arr[i][j] = h_this.GetBinContent(i+1,j+1)
    return arr

# some formating / plot cosmetic
def main_plot(h_this):
    h_this.GetYaxis().SetTitleSize(0.05)
    h_this.GetYaxis().SetTitleOffset(0.0)
    h_this.GetYaxis().SetLabelSize(0.05)
    h_this.GetXaxis().SetTitleSize(0.05)
    h_this.GetXaxis().SetLabelSize(0.05)
    h_this.GetXaxis().SetTitleOffset(1)
    h_this.SetLineColor(ROOT.kBlack)
    h_this.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)
    return


# In[10]:


# define your irregular binning
# bins = [ -1.00, -0.50, 0.00, 0.28, 0.47, 0.63, 0.765, 0.865, 0.935, 1.00 ] # example of irregular binning
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# number of bins
binnum = len(bins) - 1
print binnum


# In[2]:


# plot data event rate and prediction event rate before cut

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# define the histograms h_data and h_mc
h_data = 
h_mc = 

# fill it with data.Draw()
data.Draw(%%%)
data.Draw(%%%)
# scale mc
h_mc.Scale(%%%)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# cosmetics
ROOT.gStyle.SetEndErrorSize(5)
h_data.SetXTitle("X reconstructed [unit]")
h_data.SetYTitle("# entries per bin")
main_plot(h_data)
h_data.SetMaximum(h_data.GetMaximum()*1.2)

h_mc.SetLineColor(ROOT.kRed)
h_mc.SetLineWidth(4)

# plot the data and MC histogram
h_data.Draw('E1')
h_mc.Draw('same hist')

# add the legend
legend = ROOT.TLegend(0.15,0.65,0.3,0.9)
legend.AddEntry(h_data,'data',"lep");
legend.AddEntry(h_mc,'MC',"l");
legend.Draw()

# add the premil test
prelim, pot_tex = draw_adding()
prelim.Draw()
pot_tex.Draw()

c1.Draw()
c1.SaveAs(outputdir + "h_eventrate_nocut.png") # png for slides
c1.SaveAs(outputdir + "h_eventrate_nocut.root") # root files for analysis/small changes
c1.SaveAs(outputdir + "h_eventrate_nocut.pdf") # pdf for documents

h_data.Write('h_data_nocut')
h_mc.Write('h_mc_nocut')


# In[12]:


# plot data event rate and prediction event rate after cut
# same as before but with isSelected

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_data 
h_mc = 
data.Draw(%%%)
data.Draw(%%%)
h_mc.Scale(%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cosmetics
ROOT.gStyle.SetEndErrorSize(5)
h_data.SetXTitle("X reconstructed [unit]")
h_data.SetYTitle("# entries per bin")
main_plot(h_data)
h_data.SetMaximum(h_data.GetMaximum()*1.2)

h_mc.SetLineColor(ROOT.kRed)
h_mc.SetLineWidth(4)

# plot the data and MC histogram
h_data.Draw('E1')
h_mc.Draw('same hist')

# add the legend
legend = ROOT.TLegend(0.15,0.65,0.3,0.9)
legend.AddEntry(h_data,'data',"lep");
legend.AddEntry(h_mc,'MC',"l");
legend.Draw()

# add the premil test
prelim, pot_tex = draw_adding()
prelim.Draw()
pot_tex.Draw()

c1.Draw()
c1.SaveAs(outputdir + "h_eventrate_selected.png") # png for slides
c1.SaveAs(outputdir + "h_eventrate_selected.root") # root files for analysis/small changes
c1.SaveAs(outputdir + "h_eventrate_selected.pdf") # pdf for documents

h_data.Write('h_data')
h_mc.Write('h_mc')


# In[13]:


# get efficiency
weight_name = '1' # ub_tune_weight or whatever...

#define nom and denom 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nenner_cut = 
zahler_cut = 
#define the histograms
h_sel = 
h_gen = 
#fill with data.Draw
data.Draw(%%%)
data.Draw(%%%)

# Calculate TEfficiency
eff =  ROOT.TEfficiency(%%%) # for correct efficiency error bars (asymmetric)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
eff.SetConfidenceLevel(0.68)
eff.SetTitle('Effieciency')
eff.Draw("AP")
ROOT.gPad.Update()
graph = eff.GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1)
graph.SetLineWidth(2)
graph.GetXaxis().SetTitle("Truth X")
graph.GetYaxis().SetTitle("Signal efficiency")
graph.GetYaxis().SetTitleSize(0.05)
graph.GetYaxis().SetTitleOffset(0.0)
graph.GetYaxis().SetLabelSize(0.05)
graph.GetXaxis().SetTitleSize(0.05)
graph.GetXaxis().SetLabelSize(0.05)
graph.GetXaxis().SetTitleOffset(1)
graph.SetLineColor(ROOT.kBlack)
graph.SetLineWidth(4)
ROOT.gStyle.SetEndErrorSize(5)
graph.Draw("AP")
prelim = draw_sim()
prelim.Draw()
c1.Draw()
c1.SaveAs(outputdir+ "h_efficiency.png")
#c1.SaveAs(outputdir + "h_efficiency.root")
#c1.SaveAs(outputdir + "h_efficiency.pdf")

eff.Write("h_eff_true") # TEfficiency plot
h_eff = h_sel.Clone() # makes a copy
h_eff.Divide(h_gen) # efficiency as a simple TH1F
h_eff.Write('h_eff')

# store scaled selected and generated histograms
h_sel.Scale(potData/potMC)
h_gen.Scale(potData/potMC)
h_sel.Write('h_true_sel') # Selected true events (scaled)
h_gen.Write('h_true_gen') # gen true event (scaled)

# make a numpy array with the efficiency
eff_arr = plot2arr(h_eff,1)


# In[14]:


# plot migration matrix with arbitary small binning

#define 2d histogram
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_migration = 
# fill it
data.Draw(%%%)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cosmetics
h_migration.SetXTitle("reco X")
h_migration.SetYTitle("true X")
h_migration.GetYaxis().SetTitleSize(0.05)
h_migration.GetYaxis().SetTitleOffset(0.0)
h_migration.GetYaxis().SetLabelSize(0.05)
h_migration.GetXaxis().SetTitleSize(0.05)
h_migration.GetXaxis().SetLabelSize(0.05)
h_migration.GetXaxis().SetTitleOffset(1)
h_migration.SetLineColor(ROOT.kBlack)
h_migration.SetLineWidth(4)
ROOT.gStyle.SetEndErrorSize(5)
c1.SetRightMargin(0.1)
h_migration.Draw("colz")
prelim = draw_sim()
prelim.Draw()
#h_migration.Draw("same text")
c1.Draw()
#c1.SaveAs(outputdir + "h2_RecoVsTrue_fine.root")
c1.SaveAs(outputdir + "h2_RecoVsTrue_fine.png")
#c1.SaveAs(outputdir + "h2_RecoVsTrue_fine.pdf")

#store it
h_migration.Write("h2_RecoVsTrue_fine")


# In[15]:


# plot migration matrix now with specified binning

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# define it
h_migration = 
# fill it
data.Draw(%%%)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#cosmetics
h_migration.SetXTitle("true X")
h_migration.SetYTitle("reco X")
h_migration.GetYaxis().SetTitleSize(0.05)
h_migration.GetYaxis().SetTitleOffset(0.0)
h_migration.GetYaxis().SetLabelSize(0.05)
h_migration.GetXaxis().SetTitleSize(0.05)
h_migration.GetXaxis().SetLabelSize(0.05)
h_migration.GetXaxis().SetTitleOffset(1)

ROOT.gStyle.SetPaintTextFormat('0.0f')
c1.SetRightMargin(0.1)

h_migration.Draw("colz")
h_migration.Draw("text same")
prelim = draw_sim()
prelim.Draw()
c1.Draw()

#c1.SaveAs(outputdir + "h2_RecoVsTrue.root")
c1.SaveAs(outputdir + "h2_RecoVsTrue.png")
#c1.SaveAs(outputdir + "h2_RecoVsTrue.pdf")

#save it
h_migration.Write("h2_RecoVsTrue")


# In[16]:


# calculate smearing matrix from migration matrix
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# make np array
migration = plot2arr(%%%,2)
# initialize smearing matrix
smearing_matrix = np.zeros((%%%,%%%))
# prepare for line normalization
sum_reco = migration.sum(axis=1)
# normalize generated events to 1 (careful with division by 0)
smearing_matrix = get_ipython().magic(u'%%')
# fill smearing matirx and divide by efficiency
for i in range(binnum):
    for j in range(binnum):
        smearing_matrix[i][j] = %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[17]:


# plot the smearing matrix

# get the histogram from the array
h_smearing_matrix = arr2plot(smearing_matrix)
#cosmetics
h_smearing_matrix.SetXTitle("True bin i")
h_smearing_matrix.SetYTitle("Reco bin j")
h_smearing_matrix.Draw('colz')
h_smearing_matrix.Draw('same text')
ROOT.gStyle.SetPaintTextFormat('0.3f')
prelim = draw_sim()
prelim.Draw()
c1.Draw()
#c1.SaveAs(outputdir + "h_smearing_matrix.root")
c1.SaveAs(outputdir + "h_smearing_matrix.png")
#c1.SaveAs(outputdir + "h_smearing_matrix.pdf")

#store it
h_smearing_matrix.Write("h_smearing_matrix")


# In[18]:


# check your result:
# get true generated event rate, smear it -> should be the same as reco event rate

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# get true as an array from histogram
true_evrate = plot2arr(%%%,1)
# smear it
smear_true = 
# initialize histogram and fill it from the smeared array
h_smear_true = get_ipython().magic(u'%%')
for i in range(binnum):
    h_smear_true.SetBinContent(i+1,smear_true[i])
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
# get reconstructed signal as cross check
# define histo
h_sel_reco = ROOT.TH1F("h_sel_reco",'h_sel_reco',binnum,array('f',bins))
# fill it
data.Draw('xReco>>h_sel_reco','isSignal && isSelected && !isData')
#scale it
h_sel_reco.Scale(potData/potMC)
#make array
sel_reco = plot2arr(h_sel_reco,1)

# if you want to look at the both plots..
#h_sel_reco.SetLineColor(ROOT.kBlack)
#h_sel_reco.Draw('E1')
#h_smear_true.Draw('same')
#c1.Draw()


# In[19]:


# get background prediction

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# define histogram
h_bkg = get_ipython().magic(u'%%')
# fill it
data.Draw(%%%)
# scale it
%%%
# generate an array
%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[20]:


# summary:
# we have data event, background (in reco) and the smearing matrix
# also we have the true predicted event rate
# for all we have it as ROOT histogram and np array

# data: h_data; data_arr = plot2arr(h_data)
# bkg: h_bkg; bkg
# true: h_gen; true_evrate
# smearing matrix: h_smearing_matrix; smearing_matrix

# -> so any prediction (true event rate) can be smeared + bkg and then compared to data
# still missing: uncertainties
# what is needed: varied background predictions + varied smearing matrixes
# prepeare both and store them, they will be used for each new prediction/model comparison


# In[21]:


# get file with histograms
if(calc_new==0):
    f_sys = ROOT.TFile.Open(path_sys+"tutorial_histo.root", 'read')


# In[22]:


# vary backgrounds:

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define array with number of universes and bin numbers
bkg_err = get_ipython().magic(u'')
# initialize list of histograms
h_bkg_err = [] # all histograms (for eventual plotting)
# if the syst. are not precalculated
if(calc_new==1): # this takes long
    # loop over universes
    for uni in range(nUniverses):
        # initialize histogram
        h_bkg_tmp = 
        # fill it
        data.Draw(%%%)
        # scale it
        get_ipython().magic(u'%% %')
        # fill array
        bkg_err[uni] = %%%
        # add histogram to list
        h_bkg_err.append(%%%)
        # store it for later use...
        h_bkg_tmp.Write('h_bkg_err['+str(uni)+']')
        del h_bkg_tmp
else: # if precalculated...
    for uni in range(nUniverses):
        h_bkg_err.append(f_sys.Get('h_bkg_err['+str(uni)+']'))
        bkg_err[uni] = plot2arr(h_bkg_err[uni],1)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[23]:


# vary efficiency:
# same as before but now with the efficiency
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eff_err = get_ipython().magic(u'%%')
h_eff_err = [] # all histograms (for eventual plotting)
if(calc_new==1): # this takes long
    for uni in range(nUniverses):
        # calculate each efficiency
        nenner_cut = get_ipython().magic(u'%%')
        zahler_cut = get_ipython().magic(u'%%')
        h_sel_tmp = get_ipython().magic(u'%%')
        h_gen_tmp = get_ipython().magic(u'%%')
        data.Draw(%%%)
        data.Draw(%%%)
        h_sel_tmp.Divide(h_gen_tmp)
        h_eff_err.append(h_sel_tmp)
        eff_err[uni] = %%%
        h_sel_tmp.Write('h_eff_err['+str(uni)+']')
        del h_sel_tmp
else:
    for uni in range(nUniverses):
        h_eff_err.append(f_sys.Get('h_eff_err['+str(uni)+']'))
        eff_err[uni] = plot2arr(h_eff_err[uni],1)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[24]:


# vary smearing:
# same as before but now with the smearing matrix (needs varied efficiecny)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smear_err = get_ipython().magic(u'%%')
h_smear_err = [] # all histograms (for eventual plotting)

if(calc_new==1): # takes long
    for uni in range(nUniverses):
        h_smear_tmp = get_ipython().magic(u'%%')
        data.Draw(%%%)

        migration = plot2arr(h_smear_tmp,2)
        smearing_matrix = get_ipython().magic(u'%%%')
        sum_reco = migration.sum(axis=1)
        smearing_matrix = get_ipython().magic(u'%%%')
        for i in range(binnum):
            for j in range(binnum):
                smearing_matrix[i][j] = %%%

        h_smear_err.append(arr2plot(smearing_matrix))
        h_smear_err[uni].Write('h_smear_err['+str(uni)+']')
else:
    for uni in range(nUniverses):
        h_smear_err.append(f_sys.Get('h_smear_err['+str(uni)+']'))
        smear_err[uni] = plot2arr(h_smear_err[uni],2)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[25]:


# calculate syst. covar from varied bkg and smearing histograms

# function takes true event rate and calculates the variation from the universes
# needs filled list of h_smear_err and h_bkg_err from the calculations before!
def return_covar(h_true):
    # defince CV values
    h_true_cv = h_true.Clone()
    h_bkg_cv = h_bkg.Clone()   # idealy from a file 
    h_smear_cv = h_smearing_matrix.Clone() # idealy from a file 
    
    num_bins = h_smearing_matrix.GetNbinsX() # number of bins in the histograms
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # generate numpy arrays from histograms
    smear_mat = get_ipython().magic(u'%%')
    true_vec = get_ipython().magic(u'%%')
    bkg_vec = get_ipython().magic(u'%%')

    #calculate the semared true
    reco_vec = get_ipython().magic(u'%%')
    
    # for the storage of smeared true and background varied distributions as arrays
    reco_err = get_ipython().magic(u'%%')
    bkg_err = get_ipython().magic(u'%%')
    
    # also for the residuals
    residual_vec = get_ipython().magic(u'%%')
    
    # loop over the variations
    for uni in range(nUniverses):
        # get varied smearing and background distributions 
        # get varied histogram from list
        h_smear_this = get_ipython().magic(u'%%')
        # make array
        smear_this = get_ipython().magic(u'%%')
        # smear true to reco
        reco_err[uni] = %%%
        # get background histogram
        h_bkg_this = get_ipython().magic(u'%%')
        # make array
        bkg_err[uni] =%%%
        
        #calculated residual
        residual_vec[uni] = %%%

    # calculate covariance
    # initialize it
    covar = get_ipython().magic(u'%%')
    #loop over universes
    for uni in range(nUniverses):
        # loop over all bins in 2d
        for j in range(num_bins):
            for k in range(num_bins):
                #calculate the covariance matrix
                covar[j][k] += %%%
    covar = covar/nUniverses
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #return the covariance matrix
    return covar


# In[26]:


# get sys. covar with your model prediction as an event rate
covar = return_covar(h_gen)


# In[27]:


# calculate statistical error on prediction

# takes the prediction histograms as an input
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def return_stat_covar(h_true_func,h_bkg_func):
    # get number of bins
    n_bins = get_ipython().magic(u'%%')
    # initialize covar matrix
    covar_stat = get_ipython().magic(u'%%')
    
    # get array from the input histograms
    true_vec = get_ipython().magic(u'%%')
    bkg_vec = get_ipython().magic(u'%%')
    #get the semaring matix (!has to be defined before in the code or loaded from a file)
    h_smear = h_smearing_matrix.Clone() # idealy from a file 
    # as an array
    smearing = get_ipython().magic(u'%%')
    # smear true to reco
    reco_vec = get_ipython().magic(u'%%')
    
    #fill the diagonal with stat err, rest is zero
    for i in range(n_bins):
        covar_stat[i][i] = %%%
    
    #return the covariance marix
    return covar_stat
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In[28]:


# Get stat. covar
stat_covar = return_stat_covar(h_gen,h_bkg)


# In[29]:


# plot event rate comparison and calculate chi2
def eventrate_comparison(h_data_func,h_bkg_func,h_true_func,this_covar, filename,model_name):
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.1)
    c1.cd()

    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)

    main_plot(h_data_func)
    
    true_vec = plot2arr(h_true_func,1)
    h_smear = h_smearing_matrix.Clone() # idealy from a file 
    smearing = plot2arr(h_smear,2)
    reco_vec = true_vec.dot(smearing)    
    h_reco_func = h_true_func.Clone() # in order to keep binning
    for i in range(h_reco_func.GetNbinsX()):
        h_reco_func.SetBinContent(i+1,reco_vec[i])
    
    hs = ROOT.THStack("hs","");
    h_bkg_func.SetFillColor(ROOT.kGray)
    h_bkg_func.SetLineColor(ROOT.kGray)
    h_reco_func.SetFillColor(ROOT.kRed)
    h_reco_func.SetLineColor(ROOT.kRed)
    hs.Add(h_bkg_func)
    hs.Add(h_reco_func)

    h_data_func.SetXTitle("Reco bin number")
    h_data_func.SetYTitle("Number of events")

    h_tot = h_reco_func.Clone()
    h_tot.Add(h_bkg_func)
    
    num_bins = h_true_func.GetNbinsX() # get eventually better than hardcode
    #pred_vec = FF.plot2arr(h_tot,1)
    
    h_tot.SetFillColor(ROOT.kBlack)
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetFillStyle(3004)
    
    for i in range(num_bins):
        h_tot.SetBinError(i+1,math.sqrt(this_covar[i][i]))
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # calculate CHi2
    # data distribution as array
    dat_vec = get_ipython().magic(u'%%')
    # prediction as array
    tot_vec = get_ipython().magic(u'%%')
    # invert the covariance matrix
    inv_covar = get_ipython().magic(u'%%')
    # take the difference between the data and the prediction
    diff_vec = get_ipython().magic(u'%%')
    chi2 = 0.0
    # calculate the chi2 value
    for i in range(num_bins):
        for j in range(num_bins):
            chi2 += %%%

    print 'Chi2 = ',chi2
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h_data_func.SetMaximum(-1111)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    h_data_func.Draw('E1')
    hs.Draw('same hist')
    h_tot.Draw('same E2')
    
    legend = ROOT.TLegend(0.15,0.45,0.5,0.9)
    legend.AddEntry(h_data_func,'Data + stat.',"lep");
    legend.AddEntry(h_tot,model_name+' + bkg: #chi^{2}:'+'{:04.1f})'.format(chi2),"f");
    legend.AddEntry(h_reco_func,'Signal',"f")
    legend.AddEntry(h_bkg_func,'background',"f")
    legend.Draw()
    
    h_data_func.Draw('E1 same')

    prelim, pot_tex = draw_adding()
    prelim.Draw()
    pot_tex.Draw()

    h_data_func.GetYaxis().SetTitleOffset(1);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);

    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    h_tot_noerr = h_tot.Clone()
    for i in range(43):
        h_tot_noerr.SetBinError(i,0)
    h_xsec_data_ratio = h_data_func.Clone()
    h_xsec_data_ratio.Divide(h_tot_noerr)
    h_mc_ratio = h_tot.Clone()
    h_mc_ratio.Divide(h_tot_noerr)

    h_xsec_data_ratio.SetYTitle('Data/(Ext+MC)')
    h_xsec_data_ratio.GetYaxis().SetTitleSize(0.1)
    h_xsec_data_ratio.GetYaxis().SetTitleOffset(0.3)
    h_xsec_data_ratio.GetYaxis().SetLabelSize(0.07)
    h_xsec_data_ratio.GetXaxis().SetTitleSize(0.15)
    h_xsec_data_ratio.GetXaxis().SetLabelSize(0.15)
    h_xsec_data_ratio.GetXaxis().SetTitleOffset(1)
    h_xsec_data_ratio.SetMaximum(1.5)
    h_xsec_data_ratio.SetMinimum(0.5)
    h_xsec_data_ratio.Draw('E1')
    h_mc_ratio.Draw('E2 same')

    c1.Draw()
    c1.SaveAs(outputdir + filename+".png")
    #c1.SaveAs(outputdir + filename+".root")
    #c1.SaveAs(outputdir + filename+".pdf")

    return


# In[30]:


# now try new model
# this means: error recalculation

# define new model:
h_true_m2 = h_gen.Clone()
true_m2 = np.zeros((h_true_m2.GetNbinsX()))
for i in range( h_true_m2.GetNbinsX()):
    h_true_m2.SetBinContent(i+1,h_gen.GetBinContent(i+1)*0.9)
    true_m2[i] = h_true_m2.GetBinContent(i+1)


# In[31]:


# calculation of the uncertainties using model 2
covar_m2 = return_covar(h_true_m2)
stat_covar_m2 = return_stat_covar(h_true_m2,h_bkg)


# In[32]:


#plot the event rates
eventrate_comparison(h_data,h_bkg,h_gen,covar+stat_covar, 'event_rate','default')
eventrate_comparison(h_data,h_bkg,h_true_m2,covar_m2+stat_covar_m2, 'event_rate_m2','model 2')


# In[ ]:


RootFile.Close()


# In[ ]:





# In[ ]:





# In[ ]:




