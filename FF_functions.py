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
from array import array
import collections

from array import array
#mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
#mom_bins = [ 0.00, 0.225, 0.28, 0.33, 0.39, 0.52, 0.78, 1.21, 2.5]

## Get bining for binwidths...#########################################################
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
#########################################################################################

## define flux and Number of targets ####################################################
NumberOfFiles = 5000
POT_file = 5e8
activeVol_Area = 128.175*2*116.5*2
pot_data = 2.144e+20

path = '/home/tmettler/Desktop/ub_data/mcc9.1/v08_00_00_33/V08_00_00_35/Flux/zarko/'
f_flux_int = ROOT.TFile.Open(path+"MCC9_FluxHist_volTPCActive.root", 'read')
h_flux_cv = f_flux_int.Get("hEnumu_cv")
flux_cv = h_flux_cv.Integral(-1,201)/(NumberOfFiles*POT_file*activeVol_Area)*pot_data

roh_data = 1.3836 #data denisity g/cm3
roh_MC = 1.3954 # MC denisity g/cm3
mol = 39.95 # g for argon
N_A = 6.022140857e23 # molec/mol avogadro zahl
N_nucleons = 40.0
V_fid = ((254.8-10)-(-1.55+10))*((117.47-10)-(-115.53+10))*((1036.9-50)-(-0.1+10))

N_tot = roh_data*N_A*N_nucleons*V_fid/mol

print 'N_tot = ',N_tot,'  flux= ',flux_cv, 'for pot: ', pot_data

#########################################################################################


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
    
def xsec2rate(h_true):
    #xsec to event rate:
    # need pot/flux, xsec
    #h_rate = h_true.Clone()
    h_rate = ROOT.TH1F('h_rate',"h_rate",43,0,43)
    bin_counter = 0
    for t_bin in range(len_theta): #len_theta
        #h_temp = f_cv.Get('h_data['+str(t_bin)+']')
        for m_bin in range(len_mom[t_bin]):
            bin_counter+=1
            bin_width =  mom_bins[t_bin][m_bin+1] - mom_bins[t_bin][m_bin]
            bin_width_theta = theta_bins[t_bin+1] - theta_bins[t_bin]
            h_rate.SetBinContent(bin_counter,h_true.GetBinContent(bin_counter)*(bin_width*bin_width_theta))
            h_rate.SetBinError(bin_counter,h_true.GetBinError(bin_counter)*(bin_width*bin_width_theta))
    h_rate.SetBinContent(43,1e-80)
    h_rate.SetBinError(43,1e-80)
    h_rate.Scale(N_tot*flux_cv)
    return h_rate

#function takes np array and returns the THXF
def arr2plot(this_arr):
    this_array = np.copy(this_arr)
    #tmp = ROOT.gROOT.FindObject("h_this")
    #del tmp
    n_bins = this_array.shape
    dim = len(n_bins)
    if dim == 1:
        h_this1 = ROOT.TH1F('h_this1',"Title",n_bins[0],0,n_bins[0])
        for i in range(n_bins[0]):
            h_this1.SetBinContent(i+1,this_array[i])
    if dim == 2:
        h_this1 = ROOT.TH2F('h_this1',"Title",n_bins[0],0,n_bins[0],n_bins[1],0,n_bins[1])
        for i in range(n_bins[0]):
            for j in range(n_bins[0]):
                h_this1.SetBinContent(i+1,j+1,this_array[i][j])
    h_this1.SetDirectory(0)
    main_plot(h_this1)  
    #del h_this      
    return h_this1

#function takes THXF and the dimension and returns np array
def plot2arr(h_this_tmp, dim):
    #h_this = h_this_tmp.Clone()
    if dim == 1:
        n_bins = h_this_tmp.GetNbinsX()
        arr = np.zeros((n_bins))
        for i in range(n_bins):
            arr[i] = h_this_tmp.GetBinContent(i+1)
    if dim == 2:
        n_bins = np.zeros((2))
        n_bins[0] = int(h_this_tmp.GetNbinsX())
        n_bins[1] = int(h_this_tmp.GetNbinsY())
        n_bins.astype(np.int64)
        #print n_bins
        arr = np.zeros((np.int64(n_bins[0]),np.int64(n_bins[1])))
        for i in range(np.int64(n_bins[0])):
            for j in range(np.int64(n_bins[1])):
                arr[i][j] = h_this_tmp.GetBinContent(i+1,j+1)
    #del h_this
    return arr

def histBkg(hist):
    h_bin_bkg = hist.Clone()
    bin_counter=0
    for t_bin in range(len_theta): #len_theta
            for m_bin in range(len_mom[t_bin]):
                bin_counter+=1
                h_bin_bkg.SetBinContent(bin_counter,0)
                if t_bin%2==0:
                    h_bin_bkg.SetBinContent(bin_counter,10000)
    h_bin_bkg.SetLineWidth(0)
    h_bin_bkg.SetFillColorAlpha(ROOT.kGray,0.2)
    return h_bin_bkg
    
def plot_err_array(this_arr):
    this_array = np.copy(this_arr)
    n_bins = this_array.shape
    dim = len(n_bins)
    if dim == 1:
        h_this_err = ROOT.TH1F('h_this_err',"Title",n_bins[0]-1,0,n_bins[0]-1)
        for i in range(n_bins[0]-1):
            h_this_err.SetBinContent(i+1,1)
            h_this_err.SetBinError(i+1,this_array[i])
    if dim == 2:
        h_this_err = ROOT.TH1F('h_this_err',"Title",n_bins[0]-1,0,n_bins[0]-1)
        for i in range(n_bins[0]-1):
            h_this_err.SetBinContent(i+1,1)
            h_this_err.SetBinError(i+1,math.sqrt(this_array[i][i]))
                                 
    main_plot(h_this_err) 
    h_this_err.SetXTitle("Reco bin number")
    h_this_err.SetYTitle("frac error")
    h_this_err.SetDirectory(0)
    return h_this_err
  
  
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

def cov2corr(this_covar):
    n_bins = this_covar.shape
    this_corr = np.copy(this_covar)
    for i in range(n_bins[0]):
        for j in range(n_bins[0]):
            this_corr[i][j] = this_covar[i][j]/(math.sqrt(this_covar[i][i]*this_covar[j][j]))

    return this_corr
        

def frac2cov(frac, pred):
  num_bins = frac.shape
  covar = np.zeros((num_bins[0],num_bins[0]))
  for i in range(num_bins[0]):
    for j in range(num_bins[0]):
      covar[i][j] = frac[i][j]*pred[i]*pred[j]
  
  
  return covar













