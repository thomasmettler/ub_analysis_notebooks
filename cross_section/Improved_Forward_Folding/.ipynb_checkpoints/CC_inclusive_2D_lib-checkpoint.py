
#import stuff you eventually need
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
from array import array



inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted_improved/My_Measurement/' # please give here path to inputfile
outputdir = inputdir+'plots/' 
# change this to a single file:
f_in = ROOT.TFile.Open(inputdir+'FF_detsys.root', 'read')
f_flux = ROOT.TFile.Open(inputdir+'FF_flux.root', 'read')
f_flux_int = ROOT.TFile.Open(inputdir+'FF_neutrino_flux.root', 'read')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.1)
c1.SetBottomMargin(0.14)

#define binning
mom_bins = {}
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

def main_plot(h_this):
    h_this.GetYaxis().SetTitleSize(0.05)
    h_this.GetYaxis().SetTitleOffset(0.0)
    h_this.GetYaxis().SetLabelSize(0.05)
    h_this.GetXaxis().SetTitleSize(0.05)
    h_this.GetXaxis().SetLabelSize(0.05)
    h_this.GetXaxis().SetTitleOffset(1)
    h_this.SetLineColor(ROOT.kBlack)
    h_this.SetLineWidth(2)
    ROOT.gStyle.SetEndErrorSize(5)
    return

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
    h_this.SetDirectory(0)        
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
        #print n_bins
        arr = np.zeros((np.int64(n_bins[0]),np.int64(n_bins[1])))
        for i in range(np.int64(n_bins[0])):
            for j in range(np.int64(n_bins[1])):
                arr[i][j] = h_this.GetBinContent(i+1,j+1)
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
        h_this = ROOT.TH1F('h_this',"Title",n_bins[0]-1,0,n_bins[0]-1)
        for i in range(n_bins[0]-1):
            h_this.SetBinContent(i+1,1)
            h_this.SetBinError(i+1,this_array[i])
    if dim == 2:
        h_this = ROOT.TH1F('h_this',"Title",n_bins[0]-1,0,n_bins[0]-1)
        for i in range(n_bins[0]-1):
            h_this.SetBinContent(i+1,1)
            h_this.SetBinError(i+1,math.sqrt(this_array[i][i]))
                                 
    main_plot(h_this) 
    h_this.SetXTitle("Reco bin number")
    h_this.SetYTitle("frac error")
    h_this.SetDirectory(0)        
    return h_this
  
  
def draw_adding(this_pot):
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()

    pot_latex = ROOT.TLatex(.10, .92,'Accumulated POT: '+str(this_pot)) 
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


# In[6]:


def smear_plot(h_true, h_smearing):
    if h_true.GetNbinsX() != 43:
        h_true_this = ROOT.TH1F('h_true_this',"h_true_this",43,0,43)
        for i in range(43):
            h_true_this.SetBinContent(i+1,h_true.GetBinContent(i+1))
    else:
        h_true_this = h_true.Clone()
    h_smear = h_smearing.Clone()
    smear_mat = plot2arr(h_smear,2)
    eff_vec = plot2arr(h_eff,1)
    true_vec = plot2arr(h_true_this,1)
    #smear_mat = smear_mat.T*eff_vec
    #reco_vec = smear_mat.dot(true_vec)
    reco_vec = true_vec.dot(smear_mat)
    h_reco = arr2plot(reco_vec)
    h_true_this.SetDirectory(0)        
    return h_reco

def smear_plot(h_true):
    h_smearing = f_in.Get('h_smear_cv')
    if h_true.GetNbinsX() != 43:
        h_true_this = ROOT.TH1F('h_true_this',"h_true_this",43,0,43)
        for i in range(43):
            h_true_this.SetBinContent(i+1,h_true.GetBinContent(i+1))
    else:
        h_true_this = h_true.Clone()
    h_smear = h_smearing.Clone()
    smear_mat = plot2arr(h_smear,2)
    #eff_vec = plot2arr(h_eff,1)
    true_vec = plot2arr(h_true_this,1)
    #smear_mat = smear_mat.T*eff_vec
    #reco_vec = smear_mat.dot(true_vec)
    reco_vec = true_vec.dot(smear_mat)
    h_reco = arr2plot(reco_vec)
    h_true_this.SetDirectory(0)        
    return h_reco

# In[7]:


def return_detsys_covar(h_true,this_pot):
    print 'Doing detector systematics'
    para = ['dedx','LYatt','LYdown','LYray','recomb2','sce','waxz','wayz','wmx','wmyz']
    
    h_true_cv = h_true.Clone()
    h_bkg_cv = f_in.Get('h_bkg_detcv')
    h_bkg_cv.Scale(this_pot/9.457e+18)
    h_smear_cv = f_in.Get('h_smear_detcv')
    
    h_ext = f_in.Get('h_ext_cv')
    h_dirt = f_in.Get('h_dirt_cv')
    h_ext.Scale(this_pot/2.144e+20)
    h_dirt.Scale(this_pot/2.144e+20)
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt,1)
    
    num_bins = h_smear_cv.GetNbinsX()
    smear_mat = plot2arr(h_smear_cv,2)
    true_vec = plot2arr(h_true_cv,1)
    bkg_vec = plot2arr(h_bkg_cv,1)

    reco_vec = true_vec.dot(smear_mat)

    bkg_vec_det = np.zeros((len(para),num_bins))
    reco_vec_det = np.zeros((len(para),num_bins))
    smear_mat_det = np.zeros((len(para),num_bins,num_bins))

    res_vec_det = np.zeros((len(para),num_bins)) # residual
    res_vec = np.zeros((num_bins)) # quadratic sum

    h_bkg_det = []
    h_smear_det = []
    #print 'go into loop'
    for i,x in enumerate(para):
        h_bkg_det.append(f_in.Get('h_bkg_'+x))
        h_bkg_det[i].Scale(this_pot/9.457e18)
        h_smear_det.append(f_in.Get('h_smear_'+x))
        
        bkg_vec_det[i] = plot2arr(h_bkg_det[i],1)
        smear_mat_det[i] = plot2arr(h_smear_det[i],2)
        
        reco_vec_det[i]= true_vec.dot(smear_mat_det[i])
        #print reco_vec[1],bkg_vec[1], reco_vec_det[i][1] , bkg_vec_det[i][1]
        res_vec_det[i] = ((reco_vec+bkg_vec) - (reco_vec_det[i] + bkg_vec_det[i]))/(reco_vec+bkg_vec+ext_vec+dirt_vec)

        for j in range(num_bins):
            res_vec[j] += res_vec_det[i][j]*res_vec_det[i][j]

    frac_covar_det = np.zeros((num_bins,num_bins))
    for i,x in enumerate(para):
        for j in range(num_bins):
            for k in range(num_bins):
                frac_covar_det[j][k] += res_vec_det[i][j] * res_vec_det[i][k]

    return frac_covar_det

def return_other_covar(h_true,this_pot):
    print 'Doing Genie single variation uncertainties'
    
    weight_list = [ 'AxFFCCQEshape_UBGenie' , 'DecayAngMEC_UBGenie', 'NormCCCOH_UBGenie', 'NormNCCOH_UBGenie',                    'RPA_CCQE_UBGenie','ThetaDelta2NRad_UBGenie','Theta_Delta2Npi_UBGenie',                    'VecFFCCQEshape_UBGenie','XSecShape_CCMEC_UBGenie']
    w_list = [0,1,2,3,4,5,6,7,8]
    
    h_true_cv = h_true.Clone()
    num_bins = h_true_cv.GetNbinsX()
    frac_covar = np.zeros((num_bins,num_bins))
    
    h_ext = f_in.Get('h_ext_cv')
    h_dirt = f_in.Get('h_dirt_cv')
    h_ext.Scale(this_pot/2.144e+20)
    h_dirt.Scale(this_pot/2.144e+20)
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt,1)
    true_vec = plot2arr(h_true_cv,1)
        
    for w in w_list:
        #take first uiverse as cv
        h_bkg_cv = f_in.Get('h_bkg_other['+str(w)+'][0]')
        h_bkg_cv.Scale(this_pot/7.644e+18)
        h_smear_cv = f_in.Get('h_smear_other['+str(w)+'][0]')

        h_bkg_1 = f_in.Get('h_bkg_other['+str(w)+'][1]')
        h_bkg_1.Scale(this_pot/7.644e+18)
        h_smear_1 = f_in.Get('h_smear_other['+str(w)+'][1]')
        
        bkg_vec = plot2arr(h_bkg_cv,1)
        smear_mat = plot2arr(h_smear_cv,2)
        bkg_vec1 = plot2arr(h_bkg_1,1)
        smear_mat1 = plot2arr(h_smear_1,2)

        reco_vec = true_vec.dot(smear_mat)
        reco_vec1 = true_vec.dot(smear_mat1)

        res_vec_other = ((reco_vec+bkg_vec) - (reco_vec1 + bkg_vec1))/(reco_vec+bkg_vec+ext_vec+dirt_vec)
        frac_covar_other = np.zeros((num_bins,num_bins))
        for j in range(num_bins):
            for k in range(num_bins):
                frac_covar_other[j][k] += res_vec_other[j] * res_vec_other[k]

        frac_covar += frac_covar_other
    return frac_covar

def return_all_covar(h_true,this_pot):
    print 'Doing Genie multiverse systematics'
    h_true_cv = h_true.Clone()
    h_bkg_cv = f_in.Get('h_bkg_all_cv')
    
    h_bkg_cv.Scale(this_pot/7.644e+18)
    h_smear_cv = f_in.Get('h_smearing_all_cv')
    
    h_ext = f_in.Get('h_ext_cv')
    h_dirt = f_in.Get('h_dirt_cv')
    h_ext.Scale(this_pot/2.144e+20)
    h_dirt.Scale(this_pot/2.144e+20)
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt,1)
    true_vec = plot2arr(h_true_cv,1)

    num_bins = h_smear_cv.GetNbinsX()
    smear_mat = plot2arr(h_smear_cv,2)
    bkg_vec = plot2arr(h_bkg_cv,1)

    reco_vec = true_vec.dot(smear_mat)
    
    num_universe = 600

    bkg_vec_all = np.zeros((num_universe,num_bins))
    reco_vec_all = np.zeros((num_universe,num_bins))
    smear_mat_all = np.zeros((num_universe,num_bins,num_bins))
    res_vec_all = np.zeros((num_universe,num_bins))

    h_bkg_all = []
    h_smear_all = []
    for uni in range(num_universe):
        #do something
        h_bkg_all.append(f_in.Get('h_bkg_all['+str(uni)+'];1'))
        h_bkg_all[uni].Scale(this_pot/7.644e+18)
        h_smear_all.append(f_in.Get('h_smearing_all['+str(uni)+']'))
 
        bkg_vec_all[uni] = plot2arr(h_bkg_all[uni],1)
        smear_mat_all[uni] = plot2arr(h_smear_all[uni],2)
        
        reco_vec_all[uni]= true_vec.dot(smear_mat_all[uni])
        #print uni
        #bin_ = 36
        #print reco_vec[bin_],bkg_vec[bin_],reco_vec_all[uni][bin_],bkg_vec_all[uni][bin_]
        res_vec_all[uni] = ((reco_vec+bkg_vec) - (reco_vec_all[uni] + bkg_vec_all[uni]))/(reco_vec+bkg_vec+ext_vec+dirt_vec)

    frac_covar_all = np.zeros((num_bins,num_bins))
    for uni in range(num_universe):
        for j in range(num_bins):
            for k in range(num_bins):
                frac_covar_all[j][k] += res_vec_all[uni][j] * res_vec_all[uni][k]
    frac_covar_all = frac_covar_all/num_universe
    
    del h_bkg_all
    del h_smear_all
    
    return frac_covar_all,bkg_vec_all


def return_flux_covar(h_true,this_pot):
    print 'Doing flux systematics'
    # give h_true_cv and POT since bkg is for specific POT
    weight_list = [ 'expskin_FluxUnisim', 'horncurrent_FluxUnisim', 'kminus_PrimaryHadronNormalization', 'kplus_PrimaryHadronFeynmanScaling',        'kzero_PrimaryHadronSanfordWang', 'nucleoninexsec_FluxUnisim', 'nucleonqexsec_FluxUnisim', 'nucleontotxsec_FluxUnisim',        'piminus_PrimaryHadronSWCentralSplineVariation', 'pioninexsec_FluxUnisim', 'pionqexsec_FluxUnisim', 'piontotxsec_FluxUnisim',        'piplus_PrimaryHadronSWCentralSplineVariation' ]

    h_true_cv = h_true.Clone()
    h_bkg_cv = f_flux.Get('h_bkg_flux_cv')
    h_smear_cv = f_flux.Get('h_smear_flux_cv')
    
    NumberOfFiles = 5000
    POT_file = 5e8
    activeVol_Area = 128.175*2*116.5*2
    h_flux_cv = f_flux_int.Get("h_flux_cv")
    flux_cv = h_flux_cv.Integral(-1,201)/(NumberOfFiles*POT_file*activeVol_Area)*this_pot

    h_bkg_cv.Scale(this_pot/7.644e+18)
    
    h_ext = f_in.Get('h_ext_cv')
    h_dirt = f_in.Get('h_dirt_cv')
    h_ext.Scale(this_pot/2.144e+20)
    h_dirt.Scale(this_pot/2.144e+20)
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt,1)
    true_vec = plot2arr(h_true_cv,1)

    num_bins = h_smear_cv.GetNbinsX()
    smear_mat = plot2arr(h_smear_cv,2)
    bkg_vec = plot2arr(h_bkg_cv,1)

    reco_vec = true_vec.dot(smear_mat)

    num_bins = h_smear_cv.GetNbinsX()

    num_universe = 1000

    bkg_vec_flux = np.zeros((len(weight_list),num_universe,num_bins))
    reco_vec_flux = np.zeros((len(weight_list),num_universe,num_bins))
    smear_mat_flux = np.zeros((len(weight_list),num_universe,num_bins,num_bins))

    res_vec_flux = np.zeros((len(weight_list),num_universe,num_bins))
    #h_bkg_flux = []
    #h_smear_flux = []
    for i,x in enumerate(weight_list):
        #h_bkg_flux.append([])
        #h_smear_flux.append([])
        for uni in range(num_universe):
            if uni%100 ==0:
                print 'Variation: ',i,' At universe: ',uni
            #do something
            #h_bkg_flux[i].append(f_flux.Get('h_bkg_flux['+str(i)+']['+str(uni)+']'))
            #h_bkg_flux[i][uni].Scale(this_pot/7.644e+18)
            h_bkg_tmp = f_flux.Get('h_bkg_flux['+str(i)+']['+str(uni)+']')
            h_bkg_tmp.Scale(this_pot/7.644e+18)
            #h_smear_flux[i].append(f_flux.Get('h_smear_flux['+str(i)+']['+str(uni)+']'))
            h_smear_tmp = f_flux.Get('h_smear_flux['+str(i)+']['+str(uni)+']')
            
            h_flux_var = f_flux_int.Get('h_flux_flux['+str(i)+']['+str(uni)+']')
            flux_var = h_flux_var.Integral(-1,201)/(NumberOfFiles*POT_file*activeVol_Area)*this_pot
            #print flux_cv, flux_var, (flux_var-flux_cv)/flux_var 
            bkg_vec_flux[i][uni] = plot2arr(h_bkg_tmp,1)
            smear_mat_flux[i][uni] = plot2arr(h_smear_tmp,2)
            
            reco_vec_flux[i][uni]= true_vec.dot(smear_mat_flux[i][uni])
            res_vec_flux[i][uni] = ((reco_vec+bkg_vec) - (reco_vec_flux[i][uni] + bkg_vec_flux[i][uni])*flux_var/flux_cv)/((reco_vec+bkg_vec+ext_vec+dirt_vec))
            
    frac_covar_flux_vec = np.zeros((len(weight_list),num_bins,num_bins))
    for i,x in enumerate(weight_list):
        for uni in range(num_universe):
            for j in range(num_bins):
                for k in range(num_bins):
                    frac_covar_flux_vec[i][j][k] += res_vec_flux[i][uni][j] * res_vec_flux[i][uni][k]
        frac_covar_flux_vec[i] = frac_covar_flux_vec[i]/num_universe

    frac_covar_flux = np.zeros((num_bins,num_bins))
    for i,x in enumerate(weight_list):
        frac_covar_flux += frac_covar_flux_vec[i]
    return frac_covar_flux


def return_crt_covar(h_true,this_pot):
    print 'Doing crt systematics'
    #print 'go into function'

    h_true_cv = h_true.Clone()
    h_bkg_cv = f_in.Get('h_bkg_cv')
    h_bkg_crt = f_in.Get('h_bkg_cv_nocrt')
    h_dirt_cv = f_in.Get('h_dirt_cv')
    h_dirt_crt = f_in.Get('h_dirt_cv_nocrt')
    h_ext = f_in.Get('h_ext_cv')

    h_bkg_cv.Scale(this_pot/2.144e+20)
    h_bkg_crt.Scale(this_pot/2.144e+20)
    h_dirt_cv.Scale(this_pot/2.144e+20)
    h_dirt_crt.Scale(this_pot/2.144e+20)
    h_ext.Scale(this_pot/2.144e+20)
    
    h_smear_cv = f_in.Get('h_smear_cv')
    h_smear_crt = f_in.Get('h_smear_cv_nocrt')
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt_cv,1)
    dirt_vec_crt = plot2arr(h_dirt_crt,1)
    bkg_vec = plot2arr(h_bkg_cv,1)
    bkg_vec_crt = plot2arr(h_bkg_crt,1)
    smear_mat = plot2arr(h_smear_cv,2)
    smear_mat_crt = plot2arr(h_smear_crt,2)

    true_vec = plot2arr(h_true_cv,1)

    num_bins = h_smear_cv.GetNbinsX()
    reco_vec = true_vec.dot(smear_mat)
    reco_vec_crt = true_vec.dot(smear_mat_crt)

    delta = bkg_vec+dirt_vec+reco_vec - (bkg_vec_crt+dirt_vec_crt+reco_vec_crt)
    
    delta = delta*0.1
    tot = bkg_vec+dirt_vec+reco_vec+ext_vec
    
    frac_covar_det = np.zeros((num_bins,num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            frac_covar_det[i][j] = (delta[i]*delta[j])/(tot[i]*tot[j])
            
    return frac_covar_det

def return_dirt_covar(h_true,this_pot):
    print 'Doing dirt systematics'
    n_bins = 43#h_true.GetNbinsX()
    frac_covar_dirt = np.zeros((n_bins,n_bins))
    
    h_true_cv = h_true.Clone()
    h_bkg_cv = f_in.Get('h_bkg_cv')
    h_dirt_cv = f_in.Get('h_dirt_cv')
    h_ext = f_in.Get('h_ext_cv')

    h_bkg_cv.Scale(this_pot/2.144e+20)
    h_dirt_cv.Scale(this_pot/2.144e+20)
    h_ext.Scale(this_pot/2.144e+20)
    
    h_smear_cv = f_in.Get('h_smear_cv')
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt_cv,1)
    bkg_vec = plot2arr(h_bkg_cv,1)
    smear_mat = plot2arr(h_smear_cv,2)  
    true_vec = plot2arr(h_true_cv,1)
    
    reco_vec = true_vec.dot(smear_mat)
    
    all_vec = ext_vec+dirt_vec+bkg_vec+reco_vec
    
    if n_bins == 43:
        n_bins-=1
    for i in range(n_bins):
        for j in range(n_bins):
            #print all_vec[i], dirt_vec[i]
            frac_covar_dirt[i][j] = (dirt_vec[i])*(dirt_vec[j])/(all_vec[i]*all_vec[j])
    
    
    return frac_covar_dirt


def predict(h_true,this_pot):
    h_true_cv = h_true.Clone()
    h_bkg_cv = f_in.Get('h_bkg_cv')
    h_dirt_cv = f_in.Get('h_dirt_cv')
    h_ext = f_in.Get('h_ext_cv')

    h_bkg_cv.Scale(this_pot/2.144e+20)
    h_dirt_cv.Scale(this_pot/2.144e+20)
    h_ext.Scale(this_pot/2.144e+20)
    
    h_smear_cv = f_in.Get('h_smear_cv')
    
    ext_vec = plot2arr(h_ext,1)
    dirt_vec = plot2arr(h_dirt_cv,1)
    bkg_vec = plot2arr(h_bkg_cv,1)
    smear_mat = plot2arr(h_smear_cv,2)  
    true_vec = plot2arr(h_true_cv,1)
    
    reco_vec = true_vec.dot(smear_mat)
    
    pred_vec = reco_vec + dirt_vec+ bkg_vec + ext_vec
    h_pred = arr2plot(pred_vec)
    return pred_vec, ext_vec, dirt_vec, bkg_vec, reco_vec
    


# In[14]:


def return_stat_covar(h_true,this_pot):
    print 'Doing stat. systematics'
    all_vec = predict(h_true,this_pot)
    pred_vec = all_vec[0]
    n_bins = h_true.GetNbinsX()
    frac_covar_stat = np.zeros((n_bins,n_bins))
    if n_bins == 43:
        n_bins-=1
    for i in range(n_bins):
        #print err ,val
        frac_covar_stat[i][i] = abs(pred_vec[i]) / (pred_vec[i]*pred_vec[i])
    
    
    return frac_covar_stat

def return_pot_covar(h_true,this_pot):
    print 'Doing pot systematics'
    pred_vec, ext_vec, dirt_vec, bkg_vec, reco_vec = predict(h_true,this_pot)
    n_bins = h_true.GetNbinsX()
    pred_vec_2 = ext_vec + 1.02*(dirt_vec + bkg_vec + reco_vec)
    
    frac_covar_pot = np.zeros((n_bins,n_bins))
    if n_bins == 43:
        n_bins-=1
    for i in range(n_bins):
        for j in range(n_bins):
        #print err ,val
            frac_covar_pot[i][j] = (pred_vec[i] - pred_vec_2[i])*(pred_vec[j] - pred_vec_2[j])/(pred_vec[i] * pred_vec[j])
    
    
    return frac_covar_pot


def eventrate_comparison(h_data_func,h_true_func,this_pot,this_frac_tot, filename,model_name):
    print 'Making event rate comparison'
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
    
    pred_vec, ext_vec, dirt_vec, bkg_vec, reco_vec = predict(h_true_func,this_pot)
    h_tot = arr2plot(pred_vec[0:-1])
    h_ext_func = arr2plot(ext_vec)
    h_dirt_func = arr2plot(dirt_vec)
    h_bkg_func = arr2plot(bkg_vec)
    h_reco_func = arr2plot(reco_vec)
    
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_bkg_func.SetFillColor(ROOT.kGray)
    h_bkg_func.SetLineColor(ROOT.kGray)
    h_reco_func.SetFillColor(ROOT.kRed)
    h_reco_func.SetLineColor(ROOT.kRed)
    
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    hs.Add(h_bkg_func)
    hs.Add(h_reco_func)

    h_data_func.SetXTitle("Reco bin number")
    h_data_func.SetYTitle("Number of events")

    num_bins = h_true_func.GetNbinsX() # get eventually better than hardcode    
    
    h_tot.SetFillColor(ROOT.kBlack)
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetFillStyle(3004)
    
    # calculate CHi2
    dat_vec = plot2arr(h_data_func,1)
    #dat_vec = dat_vec[0:-1]
    tot_vec = pred_vec[0:-1]
    covar = np.zeros((num_bins-1,num_bins-1))
    
    covar_data = np.zeros((num_bins-1,num_bins-1))
    for i in range(num_bins-1):
        for j in range(num_bins-1):
            covar[i][j] = this_frac_tot[i][j]*tot_vec[i]*tot_vec[j]
            
    for i in range(num_bins):
        h_tot.SetBinError(i+1,math.sqrt(this_frac_tot[i][i])*h_tot.GetBinContent(i+1))
        
    inv_covar = np.linalg.inv(covar)
    diff_vec = (dat_vec - tot_vec)

    chi2 = 0.0
    for i in range(42):
        for j in range(42):
            chi2 += diff_vec[i]*inv_covar[i][j]*diff_vec[j]
        
    print 'Chi2 = ',chi2
    h_data_func.SetMaximum(-1111)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.2)
    h_data_func.Draw('E1')
    h_histo = histBkg(h_true_func)
    #h_histo.SetFillColorAlpha(ROOT.kRed,0.2)
    h_histo.Draw('hist same')
    hs.Draw('same hist')
    h_tot.Draw('same E2')
    
    legend = ROOT.TLegend(0.13,0.45,0.47,0.9)
    legend.AddEntry(h_data_func,'Data + stat.',"lep");
    legend.AddEntry(h_tot,model_name+' + bkg: #chi^{2}:'+'{:0.1f}'.format(chi2),"f");
    legend.AddEntry(h_reco_func,'#nu_{#mu} CC',"f")
    legend.AddEntry(h_bkg_func,'Background, in-beam',"f")
    legend.AddEntry(h_dirt_func,'Dirt',"f")
    legend.AddEntry(h_ext_func,'Cosmic, off-beam',"f")
    #legend.AddEntry(h_tot,'Syst. error, #chi^{2}:'+'{:04.1f})'.format(chi2),"f");
    legend.Draw()
    
    h_data_func.Draw('E1 same')

    prelim, pot_tex = draw_adding(this_pot)
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
    h_xsec_data_ratio.SetMaximum(-1111)
    h_xsec_data_ratio.Draw('E1')
    h_histo.Draw('hist same')
    h_mc_ratio.Draw('E2 same')

    c1.Draw()
    c1.SaveAs(outputdir + filename+".png")
    c1.SaveAs(outputdir + filename+".root")
    c1.SaveAs(outputdir + filename+".pdf")

    return


def my_chi2(h_data, h_true_func,this_pot, frac):
    data_vec = plot2arr(h_data,1)
    all_vec = predict(h_true_func,this_pot)
    pred_vec = all_vec[0]
    covar = np.zeros((42,42))
    for i in range(42):
        for j in range(42):
            covar[i][j] = frac[i][j]*pred_vec[i]*pred_vec[j]
            
    inv_covar = np.linalg.inv(covar)
    chi_vec = 0
    for i in range(42):
        for j in range(42):
            chi2_tmp = (data_vec[i]-pred_vec[i])*inv_covar[i][j]*(data_vec[j]-pred_vec[j])
            chi_vec +=chi2_tmp
    return chi_vec

def chi2_minus1(h_data, h_true_func,this_pot, frac):
    data_vec = plot2arr(h_data,1)
    all_vec = predict(h_true_func,this_pot)
    pred_vec = all_vec[0]    
    covar = np.zeros((42,42))
    for i in range(42):
        for j in range(42):
            covar[i][j] = frac[i][j]*pred_vec[i]*pred_vec[j]
            
    inv_covar = np.linalg.inv(covar)
    chi_vec = np.zeros((42))
    for z in range(42):
        for i in range(42):
            for j in range(42):
                if i != z and j!=z:
                    chi2_tmp = (data_vec[i]-pred_vec[i])*inv_covar[i][j]*(data_vec[j]-pred_vec[j])
                    chi_vec[z] +=chi2_tmp
    return chi_vec


