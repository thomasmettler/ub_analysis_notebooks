#import stuff you eventually need
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
from array import array

import imp
CC = imp.load_source('CC_inclusive_2D_lib','CC_inclusive_2D_lib.py')
FF = imp.load_source('FF_functions','/home/tmettler/Desktop/uBoone/do_plots/'+'FF_functions.py')


inputdir = '/home/tmettler/Desktop/v08_00_00_33/V08_00_00_35/weighted_improved/My_Measurement/' # please give here path to inputfile
outputdir = inputdir+'plots/' 
# make output dir if not existing
try:
    os.stat(outputdir)
except:
    os.mkdir(outputdir)
RootFile = ROOT.TFile(outputdir+"result_histo.root","RECREATE");


f = ROOT.TFile.Open(inputdir+'FF_detsys.root', 'read')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c1 = ROOT.TCanvas("c1","c1",1600,1200)
c1.SetGrid(1)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.1)
c1.SetBottomMargin(0.14)

h_data = f.Get('h_data')
h_guB = f.Get('h_gen_cv')


f_gen = ROOT.TFile.Open(inputdir+"FF_generators.root", 'read')

h_xsec_g3 = f_gen.Get('h_rate_g3')
h_xsec_g2 = f_gen.Get('h_rate_g2')
h_xsec_gibuu = f_gen.Get('h_rate_gibuu')
h_xsec_nuwro = f_gen.Get('h_rate_nuwro')
h_xsec_neut = f_gen.Get('h_rate_neut')


h_true = h_xsec_nuwro.Clone()
model = 'NuWro'



frac_det = CC.return_detsys_covar(h_true,2.144e+20)
frac_other = CC.return_other_covar(h_true,2.144e+20)
frac_det_all,bkg_all = CC.return_all_covar(h_true,2.144e+20)
frac_det_flux = CC.return_flux_covar(h_true,2.144e+20)
frac_crt = CC.return_crt_covar(h_true,2.144e+20)
frac_dirt = CC.return_dirt_covar(h_true,2.144e+20)
frac_pot = CC.return_pot_covar(h_true,2.144e+20)
frac_stat = CC.return_stat_covar(h_true,2.144e+20)

np.save(outputdir+'frac_det_'+model,frac_det)
np.save(outputdir+'frac_other_'+model,frac_other)
np.save(outputdir+'frac_det_all_'+model,frac_det_all)
#np.save(outputdir+'bkg_all_'+model,bkg_all)
np.save(outputdir+'frac_det_flux_'+model,frac_det_flux)
np.save(outputdir+'frac_crt_'+model,frac_crt)
np.save(outputdir+'frac_dirt_'+model,frac_dirt)
np.save(outputdir+'frac_stat_'+model,frac_stat)
np.save(outputdir+'frac_pot_'+model,frac_pot)

frac_tot = frac_det+frac_other+frac_det_all+frac_det_flux+frac_crt+frac_dirt+frac_stat+frac_pot

np.save(outputdir+'frac_tot_'+model,frac_tot)

chi2 = CC.my_chi2(h_data, h_true,2.144e+20,frac_tot)
print chi2

h_frac = CC.plot_err_array(frac_tot)
h_frac.Draw()
#h_histo = CC.histBkg(h_frac)
#h_histo.Draw('hist same')
c1.Draw()
c1.SaveAs(outputdir + 'Frac_Tot'+model+".png")

CC.eventrate_comparison(h_data,h_true,2.144e+20,frac_tot, 'Event_rate_'+model,model)



