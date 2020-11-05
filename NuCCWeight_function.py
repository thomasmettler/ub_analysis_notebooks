import uproot
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
import imp
from array import array
globale = imp.load_source('globale','/home/tmettler/Desktop/uBoone/do_plots/globale.py')


#pot_data = 8.649e+18
pot_data = 7.644e+18

def draw_adding():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()

    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
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
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    
    return prelim

def draw_adding_ratio():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()

    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data_draw)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    return prelim, pot_latex

def draw_sim_ratio():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    
    return prelim

# @input: ttree with events, variable to plot, the sistogram, the selection cuts...
def fill_histo(tree,variable,histo,cut):
    sel = tree.CopyTree(cut)
    #if variable == 'cos(TrackTheta)':
    #    var = math.cos(entry.TrackTheta)
    for entry in sel:
        for j in range(entry.TrackTheta.size()):
            if entry.track_key.at(j) == entry.key_muon:
                if variable == 'cos(TrackTheta)':
                    histo.Fill(math.cos(entry.TrackTheta.at(j)),entry.EventWeight*entry.TunedCentralValue_Genie) 
                    #print 'Event weight: ',entry.EventWeight
                if variable == 'TrackPhi':
                    histo.Fill(entry.TrackPhi.at(j),entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'TrackMomMCS_mom':
                    histo.Fill(entry.TrackMomMCS_mom.at(j),entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'TrackLength':
                    histo.Fill(entry.TrackLength.at(j),entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'Nu_Vx_sce':
                    histo.Fill(entry.Nu_Vx_sce,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'Nu_Vy_sce':
                    histo.Fill(entry.Nu_Vy_sce,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'Nu_Vz_sce':
                    histo.Fill(entry.Nu_Vz_sce,entry.EventWeight*entry.TunedCentralValue_Genie)
                # and now the efficiency variables:
                if variable == 'MCNu_Energy':
                    histo.Fill(entry.MCNu_Energy,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'cos(MCNu_leptonTheta)':
                    histo.Fill(math.cos(entry.MCNu_leptonTheta),entry.EventWeight*entry.TunedCentralValue_Genie)
                #if variable == 'TrackPhi':
                #    histo.Fill(entry.TrackPhi,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'MCNu_Vx':
                    histo.Fill(entry.MCNu_Vx,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'MCNu_Vy':
                    histo.Fill(entry.MCNu_Vy,entry.EventWeight*entry.TunedCentralValue_Genie)
                if variable == 'MCNu_Vz':
                    histo.Fill(entry.MCNu_Vz,entry.EventWeight*entry.TunedCentralValue_Genie)
                #if variable == 'Nu_Vz_sce':
                #    histo.Fill(entry.Nu_Vz_sce,entry.EventWeight*entry.TunedCentralValue_Genie)
    sel.Delete()
    return histo

def make_stacked_histo_weight(cut,variable,weight,title,xstart,xend,xbins,file_name):
    #initialize the 1d histograms
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',cut,'')
    
    cut = cut +' && '
    for x in globale.overlay_signals:
        #histo = x
        #globale.overlay_out.Draw(variable+'>>'+histo,cut+x,'')
        h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    legend = ROOT.TLegend(0.75,0.65,0.9,0.9); #LEGEND RIGHT
    #legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
    legend.AddEntry(h_data_func,"data","lep");
    legend.AddEntry(h_ext_func,"ext","f");
    legend.AddEntry(h_dirt_func,"dirt","f");
    for x in globale.overlay_signals:
        legend.AddEntry(h_overlay_func[x],x,"f");
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    for i,x in enumerate(globale.overlay_signals):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        mc_events = mc_events+h_overlay_func[x].GetSumOfWeights()#*globale.scale[globale.overlay]
        hs.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    mc_events = mc_events + h_dirt_func.GetEntries()*globale.scale[globale.dirt]
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.Draw('E')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E same')

    # second pad
    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
        
    return normalization

def make_stacked_histo_weightV2(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    #mom_bins = [ -1.00, -0.50, 0.00, 0.28, 0.47, 0.63, 0.765, 0.865, 0.935, 1.00 ]
    mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
    binnum = len(mom_bins) - 1
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    '''h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))'''
    overlay_signals_name = ['NC','OUTFV','#bar{#nu_{#mu}},#nu_{e}, #bar{#nu_{e}} CC','#nu_{#mu} wrong matched','#nu_{#mu} (bad matched)','#nu_{#mu} CC']
    #['#nu_{#mu} CC (stopping #mu)','#nu_{#mu} CC (other)','#nu_{e}, #bar{#nu_{e}}','#bar{#nu_{#mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping #mu)', 'OUTFV (other)', 'Cosmic (stopping #mu)', 'Cosmic (other)']
    #overlay_signals_name = overlay_signals_name.reverse()
    overlay_signals_func = ['numu_nc','numu_ov','numu_other','numu_nomu','numu_lowpur','numu_signal']
    
    
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in overlay_signals_func:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(overlay_signals_func):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        #h_overlay_func[x].SetFillColor((2*i+31))
        #h_overlay_func[x].SetLineColor((2*i+31))
        h_overlay_func[x].SetLineColor(1)
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.1,0.65,0.6,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.4,0.65,0.9,0.9); #LEGEND RIGHT
    legend.SetNColumns(2)
    data_name = 'Data: {0:0.1f}'.format(data_events)
    ext_name = 'Cosmic: {0:0.1f}'.format(ext_events)
    dirt_name = 'Dirt: {0:0.1f}'.format(dirt_events)
    #['numu_nomu','numu_lowpur','numu_nc','numu_ov','numu_other','numu_signal']
    h_overlay_func['numu_nomu'].SetFillColor(ROOT.kGray+3)
    h_overlay_func['numu_nomu'].SetLineColor(ROOT.kGray+3)
    h_overlay_func['numu_lowpur'].SetFillColor(ROOT.kGray+2)
    h_overlay_func['numu_lowpur'].SetLineColor(ROOT.kGray+2)
    h_overlay_func['numu_nc'].SetFillColor(ROOT.kOrange)
    h_overlay_func['numu_nc'].SetLineColor(ROOT.kOrange)
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_ov'].SetLineColor(ROOT.kBlue)
    h_overlay_func['numu_other'].SetFillColor(ROOT.kGreen)
    h_overlay_func['numu_other'].SetLineColor(ROOT.kGreen)
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kGray)
        
    legend.AddEntry(h_data_func,data_name,"lep");
    legend.AddEntry(h_ext_func,ext_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f");
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.1)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)

    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_onoff(cut,variable,weight,title,xstart,xend,xbins,file_name,side,bin0,bin1):

    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    '''h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))'''
    overlay_signals_name = ['NC','OUTFV','#bar{#nu_{#mu}},#nu_{e}, #bar{#nu_{e}} CC','#nu_{#mu} wrong matched','#nu_{#mu} (bad matched)','#nu_{#mu} CC']
    #['#nu_{#mu} CC (stopping #mu)','#nu_{#mu} CC (other)','#nu_{e}, #bar{#nu_{e}}','#bar{#nu_{#mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping #mu)', 'OUTFV (other)', 'Cosmic (stopping #mu)', 'Cosmic (other)']
    #overlay_signals_name = overlay_signals_name.reverse()
    overlay_signals_func = ['numu_nc','numu_ov','numu_other','numu_nomu','numu_lowpur','numu_signal']
    
    
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in overlay_signals_func:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(overlay_signals_func):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        #h_overlay_func[x].SetFillColor((2*i+31))
        #h_overlay_func[x].SetLineColor((2*i+31))
        h_overlay_func[x].SetLineColor(1)
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.1,0.65,0.6,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.4,0.65,0.9,0.9); #LEGEND RIGHT
    legend.SetNColumns(2)
    data_name = 'Data: {0:0.1f}'.format(data_events)
    ext_name = 'Cosmic: {0:0.1f}'.format(ext_events)
    dirt_name = 'Dirt: {0:0.1f}'.format(dirt_events)
    #['numu_nomu','numu_lowpur','numu_nc','numu_ov','numu_other','numu_signal']
    h_overlay_func['numu_nomu'].SetFillColor(ROOT.kGray+3)
    h_overlay_func['numu_nomu'].SetLineColor(ROOT.kGray+3)
    h_overlay_func['numu_lowpur'].SetFillColor(ROOT.kGray+2)
    h_overlay_func['numu_lowpur'].SetLineColor(ROOT.kGray+2)
    h_overlay_func['numu_nc'].SetFillColor(ROOT.kOrange)
    h_overlay_func['numu_nc'].SetLineColor(ROOT.kOrange)
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_ov'].SetLineColor(ROOT.kBlue)
    h_overlay_func['numu_other'].SetFillColor(ROOT.kGreen)
    h_overlay_func['numu_other'].SetLineColor(ROOT.kGreen)
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kGray)
        
    legend.AddEntry(h_data_func,data_name,"lep");
    legend.AddEntry(h_ext_func,ext_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f");
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.1)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)

    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(-999)
    h_div_func.GetXaxis().SetLabelOffset(0.02)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetBinLabel(1, bin0)
    h_div_func.GetXaxis().SetBinLabel(2, bin1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization


def make_stacked_histo_weight_MCC8(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    mom_bins = [ -1.00, -0.50, 0.00, 0.28, 0.47, 0.63, 0.765, 0.865, 0.935, 1.00 ]
    #mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
    binnum = len(mom_bins) - 1
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    #h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    #h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    #h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    #h_ext_func.SetFillColor(2)
    #h_ext_func.SetLineColor(1)
    #h_dirt_func.SetFillColor(42)
    #h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    
    h_data_func.SetMarkerStyle(ROOT.kFullCircle);
    h_data_func.SetMarkerSize(0.9);

    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(globale.overlay_signals):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # use marcos colors:
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    
    h_overlay_func['numu_signal_cont'].SetFillColor(ROOT.kRed)
    h_overlay_func['numu_signal_cont'].SetLineColor(ROOT.kRed)
    
    h_overlay_func['numu_signal_uncont'].SetFillColor(ROOT.kRed+2)
    h_overlay_func['numu_signal_uncont'].SetLineColor(ROOT.kRed+2)
    
    h_overlay_func['numu_nue'].SetFillColor(ROOT.kMagenta+1)
    h_overlay_func['numu_nue'].SetLineColor(ROOT.kMagenta+1)
    
    h_overlay_func['numu_antinu'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['numu_antinu'].SetLineColor(ROOT.kOrange-3)
    
    h_overlay_func['numu_nc_other'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_nc_other'].SetLineColor(ROOT.kGray)
    
    h_overlay_func['numu_nc_pion'].SetFillColor(ROOT.kGray+1)
    h_overlay_func['numu_nc_pion'].SetLineColor(ROOT.kGray+1)
    
    h_overlay_func['numu_nc_proton'].SetFillColor(ROOT.kGray+2)
    h_overlay_func['numu_nc_proton'].SetLineColor(ROOT.kGray+2)
    
    h_overlay_func['numu_cosmic_cont'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_cosmic_cont'].SetLineColor(ROOT.kBlue)
    
    h_overlay_func['numu_cosmic_uncont'].SetFillColor(ROOT.kBlue+2)
    h_overlay_func['numu_cosmic_uncont'].SetLineColor(ROOT.kBlue+2)
    
    h_overlay_func['numu_ov_cont'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['numu_ov_cont'].SetLineColor(ROOT.kGreen+2)
    
    h_overlay_func['numu_ov_uncont'].SetFillColor(ROOT.kGreen+3)
    h_overlay_func['numu_ov_uncont'].SetLineColor(ROOT.kGreen+3)
    
    
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.16,0.30,0.42,0.9) # LEGEND LEFT
        #legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        #legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    #legend.SetNColumns(1)
    #data_name = 'Data: {0:0.1f}'.format(data_events)
    data_name = 'Data'
    ext_name = 'Data (Beam off): {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    dirt_name = 'Dirt: {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))
    #overlay_signals_name = ['\nu_{\mu} CC (stopping \mu)','\nu_{\mu} CC (other)','\nu_e, \overline{\nu_e}','\overline{\nu_{\mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping \mu)', 'OUTFV (other)', 'Cosmic (stopping \mu)', 'Cosmic (other)']
    overlay_signals_name = ['#nu_{#mu} CC (stopping #mu)','#nu_{#mu} CC (other)','#nu_{e}, #bar{#nu_{e}}','#bar{#nu_{#mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping #mu)', 'OUTFV (other)', 'Cosmic (stopping #mu)', 'Cosmic (other)']
    
    #overlay_signals_name = ['Cosmic (other)','Cosmic (stopping \mu)','OUTFV (other)','OUTFV (stopping \mu)', 'NC (proton)', 'NC (pion)', 'NC (other)', '\overline{\nu_{\mu}} CC', '\nu_e, \overline{\nu_e}', '\nu_{\mu} CC (other)', '\nu_{\mu} CC (stopping \mu)']
    
    
    #overlay_signals_name = overlay_signals_name.reverse()
    overlay_signals_func = ['numu_signal_cont','numu_signal_uncont','numu_nue','numu_antinu','numu_nc_other','numu_nc_pion',\
                   'numu_nc_proton', 'numu_ov_cont', 'numu_ov_uncont', 'numu_cosmic_cont', 'numu_cosmic_uncont']

    
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f")
    legend.AddEntry(h_ext_func,ext_name,"f")
    legend.AddEntry(h_data_func,data_name,"lep")
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)
    
    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    #if side != 'none':
        #legend.Draw('same')

        
        

    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_png+ file_name + ".pdf")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_weight_MCC8_V(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    
    overlay_signals_name = ['Background (MC)','Signal']
    
    overlay_signals_func = ['!numu_signal','numu_signal']
    
    
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in overlay_signals_func:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_data_func.SetLineWidth(1)
    
    h_data_func.SetMarkerStyle(ROOT.kFullCircle);
    h_data_func.SetMarkerSize(0.9);

    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(overlay_signals_func):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # use marcos colors:
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kRed-3)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kRed-3)
    
    h_overlay_func['!numu_signal'].SetFillColor(ROOT.kBlue-3)
    h_overlay_func['!numu_signal'].SetLineColor(ROOT.kBlue-3)
    
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.16,0.65,0.42,0.9) # LEGEND LEFT
        #legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        #legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
        legend = ROOT.TLegend(0.6,0.65,0.9,0.9)
    #legend.SetNColumns(1)
    #data_name = 'Data: {0:0.1f}'.format(data_events)
    data_name = 'Data'
    ext_name = 'Background (Beam off): {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    dirt_name = 'Background (Dirt): {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))
    
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f")
    legend.AddEntry(h_ext_func,ext_name,"f")
    legend.AddEntry(h_data_func,data_name,"lep")
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    max = h_data_func.GetMaximum()
    h_data_func.SetMaximum(max*1.3)
    
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)
    
    
    
    h_data_func.Draw('E1')
    #legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    if side != 'none':
        legend.Draw('same')
        
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();

    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_png+ file_name + ".pdf")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_particle(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    #mom_bins = [ -1.00, -0.50, 0.00, 0.28, 0.47, 0.63, 0.765, 0.865, 0.935, 1.00 ]
    mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
    binnum = len(mom_bins) - 1
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    '''h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))'''
    overlay_signals_name = ['OUTFV','low purity','other particle','proton','#pi^-','#pi^+','#mu']

    PDG_muon = 'fidVol && MCfidVol && MCTrackPDG==13 && MCTrackPurity>0.5' # numu CC signal definition
    PDG_pion = 'fidVol && MCfidVol && MCTrackPDG==211 && MCTrackPurity>0.5' # numu CC signal definition
    PDG_pionMinus = 'fidVol && MCfidVol && MCTrackPDG==-211 && MCTrackPurity>0.5' # numu CC signal definition
    PDG_proton = 'fidVol && MCfidVol && MCTrackPDG==2212 && MCTrackPurity>0.5' # numu CC signal definition
    PDG_other = 'fidVol && MCfidVol && MCTrackPDG!=13 && MCTrackPDG!=211 && MCTrackPDG!=-211 && MCTrackPDG!=2212 && MCTrackPurity>0.5' # numu CC signal definition
    numu_lowpur = 'fidVol && MCfidVol && MCTrackPurity<0.5' #low purity
    numu_ov = 'fidVol && !MCfidVol' # out of fiducial
    
    globale.overlay_out.SetAlias('PDG_muon',PDG_muon)
    globale.overlay_out.SetAlias('PDG_pion',PDG_pion)
    globale.overlay_out.SetAlias('PDG_pionMinus',PDG_pionMinus)
    globale.overlay_out.SetAlias('PDG_proton',PDG_proton)
    globale.overlay_out.SetAlias('PDG_other',PDG_other)
    globale.overlay_out.SetAlias('numu_lowpur',numu_lowpur)
    globale.overlay_out.SetAlias('numu_ov',numu_ov)
    
    overlay_signals_func = ['numu_ov','numu_lowpur','PDG_other','PDG_proton','PDG_pionMinus','PDG_pion','PDG_muon']
    
    
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in overlay_signals_func:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_ext_func.SetFillStyle(3004)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_dirt_func.SetFillStyle(3004)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(overlay_signals_func):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetLineColor(1)
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.1,0.65,0.6,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.4,0.65,0.9,0.9); #LEGEND RIGHT
    legend.SetNColumns(2)
    data_name = 'Data: {0:0.1f}'.format(data_events)
    ext_name = 'Cosmic: {0:0.1f}'.format(ext_events)
    dirt_name = 'Dirt: {0:0.1f}'.format(dirt_events)
    #['numu_ov','numu_lowpur','PDG_other','PDG_proton','PDG_pionMinus','PDG_pion','PDG_muon']
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kBlue+1)
    h_overlay_func['numu_lowpur'].SetFillColor(ROOT.kCyan+1)
    h_overlay_func['PDG_pionMinus'].SetFillColor(ROOT.kPink+10)
    h_overlay_func['PDG_proton'].SetFillColor(ROOT.kRed+1)
    h_overlay_func['PDG_pion'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['PDG_muon'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['PDG_other'].SetFillColor(ROOT.kGray+1)
    
    '''h_overlay_func['numu_ov'].SetLineColor(ROOT.kBlue+1)
    h_overlay_func['numu_lowpur'].SetLineColor(ROOT.kCyan+1)
    h_overlay_func['PDG_pionMinus'].SetLineColor(ROOT.kPink+10)
    h_overlay_func['PDG_proton'].SetLineColor(ROOT.kRed+1)
    h_overlay_func['PDG_pion'].SetLineColor(ROOT.kOrange-3)
    h_overlay_func['PDG_muon'].SetLineColor(ROOT.kGreen+2)
    h_overlay_func['PDG_other'].SetLineColor(ROOT.kGray+1)'''
    
    h_overlay_func['numu_ov'].SetLineColor(1)
    h_overlay_func['numu_lowpur'].SetLineColor(1)
    h_overlay_func['PDG_pionMinus'].SetLineColor(1)
    h_overlay_func['PDG_proton'].SetLineColor(1)
    h_overlay_func['PDG_pion'].SetLineColor(1)
    h_overlay_func['PDG_muon'].SetLineColor(1)
    h_overlay_func['PDG_other'].SetLineColor(1)
    
    legend.AddEntry(h_data_func,data_name,"lep");
    legend.AddEntry(h_ext_func,ext_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f");
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.1)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)

    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    #for x in globale.overlay_signals:
    #    h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization
def make_stacked_histo_MCC8_pub(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):

    #mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
    #binnum = len(mom_bins) - 1
    #initialize the 1d histograms
    #h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    #h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    #h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))
    
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_data_func.SetLineWidth(1)
    
    h_data_func.SetMarkerStyle(ROOT.kFullCircle);
    h_data_func.SetMarkerSize(0.9);

    #scale the histograms
    h_data_func.Sumw2()
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(globale.overlay_signals):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # use marcos colors:
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kRed)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kRed)
    
    h_overlay_func['numu_nomuon'].SetFillColor(ROOT.kRed+2)
    h_overlay_func['numu_nomuon'].SetLineColor(ROOT.kRed+2)
    
    h_overlay_func['numu_nue'].SetFillColor(ROOT.kMagenta+1)
    h_overlay_func['numu_nue'].SetLineColor(ROOT.kMagenta+1)
    
    h_overlay_func['numu_antinu'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['numu_antinu'].SetLineColor(ROOT.kOrange-3)
    
    h_overlay_func['numu_nc'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_nc'].SetLineColor(ROOT.kGray)
    
    h_overlay_func['numu_cosmic'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_cosmic'].SetLineColor(ROOT.kBlue)
    
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['numu_ov'].SetLineColor(ROOT.kGreen+2)

    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.15,0.30,0.5,0.9) # LEGEND LEFT
        #legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        #legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
        legend = ROOT.TLegend(0.55,0.30,0.9,0.9)
    if variable == 'TrackPID_chiproton' or variable == 'NuScore' or variable == 'TrackPhi' or variable == 'Nu_Vx_sce' or variable == 'Nu_Vy_sce' or variable == 'Nu_Vz_sce':
        legend = ROOT.TLegend(0.15,0.68,0.85,0.9) # LEGEND LEFT
        legend.SetNColumns(2)
    #data_name = 'Data: {0:0.1f}'.format(data_events)
    data_name = 'Data (Beam-on, stat. only)'#: {0:0.1f}'.format(data_events)
    #ext_name = 'Data (Beam off):'+' {0:0.1f}'.format(ext_events)+', {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    #dirt_name = 'Dirt:'+' {0:0.1f}'.format(dirt_events)+', {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))
    ext_name = 'Data (Beam-off): {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))# +' {0:0.1f}'.format(ext_events)
    dirt_name = 'Dirt: {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))#+' {0:0.1f}'.format(dirt_events)

    overlay_signals_name = ['#nu_{#mu} CC (signal)','#nu_{#mu} CC (not #mu)','#nu_{e}, #bar{#nu_{e}} CC','#bar{#nu_{#mu}} CC','NC', 'OUTFV', 'Cosmic']
    overlay_signals_func = ['numu_signal','numu_nomuon','numu_nue','numu_antinu','numu_nc','numu_ov','numu_cosmic']

    for i,x in enumerate(overlay_signals_func):
        #ov_name = '%20s'%overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        #ov_name = overlay_signals_name[i]+':\t {0:0.1f}'.format(mc_event_list[x])+',\t {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))# + ': {0:0.1f}'.format(mc_event_list[x])
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f")
    legend.AddEntry(h_ext_func,ext_name,"f")
    legend.AddEntry(h_data_func,data_name,"lep")
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    if variable == 'TrackPID_chiproton' or variable == 'TrackPhi' or variable == 'Nu_Vx_sce' or variable == 'Nu_Vy_sce' or variable == 'Nu_Vz_sce':
        h_data_func.SetMaximum(h_data_func.GetMaximum()*1.5)
    else:
        h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(2)
    ROOT.gStyle.SetEndErrorSize(9)
    
    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')

    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().CenterTitle()
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization
    
def make_stacked_histo_MCC8_pub_impMom(cut,weight,title,xstart,xend,xbins,file_name,side='right'):

    #mom_bins = [ 0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50 ]
    #binnum = len(mom_bins) - 1
    #initialize the 1d histograms
    #h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    #h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    #h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))
    
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {} # make an array of histograms for the different interactions
    h_data_func_r = ROOT.TH1F("h_data_func_r",title,xbins,xstart,xend)
    h_ext_func_r = ROOT.TH1F("h_ext_func_r",title,xbins,xstart,xend)
    h_dirt_func_r = ROOT.TH1F("h_dirt_func_r",title,xbins,xstart,xend)
    h_overlay_func_r = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        h_overlay_func_r[x] = ROOT.TH1F(x+'_r',title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
        
    globale.data_out.Draw('TrackMomMCS_mom'+'>>h_data_func',cut+'&& track_end_uncontained','')
    globale.ext_out.Draw('TrackMomMCS_mom'+'>>h_ext_func',cut+'&& track_end_uncontained','')
    globale.dirt_out.Draw('TrackMomMCS_mom'+'>>h_dirt_func',weight_name+'*('+cut+'&& track_end_uncontained'+')','')
    
    globale.data_out.Draw('TrackMomRange_mu'+'>>h_data_func_r',cut+'&& !track_end_uncontained','')
    globale.ext_out.Draw('TrackMomRange_mu'+'>>h_ext_func_r',cut+'&& !track_end_uncontained','')
    globale.dirt_out.Draw('TrackMomRange_mu'+'>>h_dirt_func_r',weight_name+'*('+cut+'&& !track_end_uncontained'+')','')

    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw('TrackMomMCS_mom'+'>>'+histo,weight_name+'*('+cut+x+'&& track_end_uncontained'+')','')
        globale.overlay_out.Draw('TrackMomRange_mu'+'>>'+histo+'_r',weight_name+'*('+cut+x+'&& !track_end_uncontained'+')','')
    # prepare the stacked histogram
    h_data_func.Add(h_data_func_r)
    h_ext_func.Add(h_ext_func_r)
    h_dirt_func.Add(h_dirt_func_r)
    
    for x in globale.overlay_signals:
        h_overlay_func[x].Add(h_overlay_func_r[x])
    
    
    hs = ROOT.THStack("hs","");
    h_data_func.SetLineWidth(1)
    
    h_data_func.SetMarkerStyle(ROOT.kFullCircle);
    h_data_func.SetMarkerSize(0.9);

    #scale the histograms
    h_data_func.Sumw2()
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(globale.overlay_signals):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # use marcos colors:
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kRed)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kRed)
    
    h_overlay_func['numu_nomuon'].SetFillColor(ROOT.kRed+2)
    h_overlay_func['numu_nomuon'].SetLineColor(ROOT.kRed+2)
    
    h_overlay_func['numu_nue'].SetFillColor(ROOT.kMagenta+1)
    h_overlay_func['numu_nue'].SetLineColor(ROOT.kMagenta+1)
    
    h_overlay_func['numu_antinu'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['numu_antinu'].SetLineColor(ROOT.kOrange-3)
    
    h_overlay_func['numu_nc'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_nc'].SetLineColor(ROOT.kGray)
    
    h_overlay_func['numu_cosmic'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_cosmic'].SetLineColor(ROOT.kBlue)
    
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['numu_ov'].SetLineColor(ROOT.kGreen+2)

    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.15,0.30,0.5,0.9) # LEGEND LEFT
        #legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        #legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
        legend = ROOT.TLegend(0.55,0.30,0.9,0.9)
    #data_name = 'Data: {0:0.1f}'.format(data_events)
    data_name = 'Data (Beam-on, stat. only)'#: {0:0.1f}'.format(data_events)
    #ext_name = 'Data (Beam off):'+' {0:0.1f}'.format(ext_events)+', {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    #dirt_name = 'Dirt:'+' {0:0.1f}'.format(dirt_events)+', {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))
    ext_name = 'Data (Beam-off): {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    dirt_name = 'Dirt: {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))

    overlay_signals_name = ['#nu_{#mu} CC (signal)','#nu_{#mu} CC (not #mu)','#nu_{e}, #bar{#nu_{e}} CC','#bar{#nu_{#mu}} CC','NC', 'OUTFV', 'Cosmic']
    overlay_signals_func = ['numu_signal','numu_nomuon','numu_nue','numu_antinu','numu_nc','numu_ov','numu_cosmic']

    for i,x in enumerate(overlay_signals_func):
        #ov_name = '%20s'%overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        #ov_name = overlay_signals_name[i]+':\t {0:0.1f}'.format(mc_event_list[x])+',\t {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f")
    legend.AddEntry(h_ext_func,ext_name,"f")
    legend.AddEntry(h_data_func,data_name,"lep")
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)

    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(2)
    ROOT.gStyle.SetEndErrorSize(9)
    
    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')

    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().CenterTitle()
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_onoff_pub(cut,variable,weight,title,xstart,xend,xbins,file_name,side,bin0,bin1):

    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    '''h_data_func = ROOT.TH1F("h_data_func",title,binnum,array('f',mom_bins))
    h_ext_func = ROOT.TH1F("h_ext_func",title,binnum,array('f',mom_bins))
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,binnum,array('f',mom_bins))'''
    #overlay_signals_name = ['NC','OUTFV','#bar{#nu_{#mu}},#nu_{e}, #bar{#nu_{e}} CC','#nu_{#mu} wrong matched','#nu_{#mu} (bad matched)','#nu_{#mu} CC']
    #['#nu_{#mu} CC (stopping #mu)','#nu_{#mu} CC (other)','#nu_{e}, #bar{#nu_{e}}','#bar{#nu_{#mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping #mu)', 'OUTFV (other)', 'Cosmic (stopping #mu)', 'Cosmic (other)']
    #overlay_signals_name = overlay_signals_name.reverse()
    #overlay_signals_func = ['numu_nc','numu_ov','numu_other','numu_nomu','numu_lowpur','numu_signal']
    overlay_signals_name = ['#nu_{#mu} CC (signal)','#nu_{#mu} CC (not #mu)','#nu_{e}, #bar{#nu_{e}} CC','#bar{#nu_{#mu}} CC','NC', 'OUTFV', 'Cosmic']
    overlay_signals_func = ['numu_cosmic','numu_ov','numu_nc','numu_antinu','numu_nue','numu_nomuon','numu_signal']
    
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
        #h_overlay_func[x] = ROOT.TH1F(x,title,binnum,array('f',mom_bins))
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in overlay_signals_func:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(overlay_signals_func):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        #h_overlay_func[x].SetFillColor((2*i+31))
        #h_overlay_func[x].SetLineColor((2*i+31))
        h_overlay_func[x].SetLineColor(1)
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.10,0.30,0.45,0.9) # LEGEND LEFT
        #legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        #legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    if variable == 'TrackPhi':
        legend = ROOT.TLegend(0.15,0.68,0.85,0.9) # LEGEND LEFT
        legend.SetNColumns(2)
    #data_name = 'Data: {0:0.1f}'.format(data_events)
    data_name = 'Data (Beam-on, stat. only)'#.format(data_events)
    #ext_name = 'Data (Beam off):'+' {0:0.1f}'.format(ext_events)+', {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    #dirt_name = 'Dirt:'+' {0:0.1f}'.format(dirt_events)+', {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))
    ext_name = 'Data (Beam off): {0:0.1f}%'.format(ext_events*100.0/(mc_events+ext_events))
    dirt_name = 'Dirt: {0:0.1f}%'.format(dirt_events*100.0/(mc_events+ext_events))

    overlay_signals_func = ['numu_signal','numu_nomuon','numu_nue','numu_antinu','numu_nc','numu_ov','numu_cosmic']

    for i,x in enumerate(overlay_signals_func):
        #ov_name = '%20s'%overlay_signals_name[i]+': {0:0.1f}'.format(mc_event_list[x])+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        #ov_name = overlay_signals_name[i]+':\t {0:0.1f}'.format(mc_event_list[x])+',\t {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100.0/(mc_events+ext_events))
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f")
    legend.AddEntry(h_ext_func,ext_name,"f")
    legend.AddEntry(h_data_func,data_name,"lep")
    
    h_ext_func.SetFillColor(ROOT.kBlue+2)
    h_ext_func.SetLineColor(ROOT.kBlue+2)
    h_ext_func.SetFillStyle(3004)
    
    h_dirt_func.SetLineColor(ROOT.kOrange+2);
    h_dirt_func.SetFillColor(ROOT.kOrange+2);
    
    h_overlay_func['numu_signal'].SetFillColor(ROOT.kRed)
    h_overlay_func['numu_signal'].SetLineColor(ROOT.kRed)
    
    h_overlay_func['numu_nomuon'].SetFillColor(ROOT.kRed+2)
    h_overlay_func['numu_nomuon'].SetLineColor(ROOT.kRed+2)
    
    h_overlay_func['numu_nue'].SetFillColor(ROOT.kMagenta+1)
    h_overlay_func['numu_nue'].SetLineColor(ROOT.kMagenta+1)
    
    h_overlay_func['numu_antinu'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['numu_antinu'].SetLineColor(ROOT.kOrange-3)
    
    h_overlay_func['numu_nc'].SetFillColor(ROOT.kGray)
    h_overlay_func['numu_nc'].SetLineColor(ROOT.kGray)
    
    h_overlay_func['numu_cosmic'].SetFillColor(ROOT.kBlue)
    h_overlay_func['numu_cosmic'].SetLineColor(ROOT.kBlue)
    
    h_overlay_func['numu_ov'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['numu_ov'].SetLineColor(ROOT.kGreen+2)
        
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.1)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.35,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.SetMaximum(h_data_func.GetMaximum()*1.3)
    
    h_data_func.GetYaxis().SetLabelSize(0.06);
    h_data_func.GetYaxis().SetTitleSize(0.06);
    h_data_func.GetYaxis().SetTitleOffset(0.8);
    h_data_func.GetXaxis().SetLabelOffset(999);
    h_data_func.GetXaxis().SetTitleOffset(999);
    h_data_func.GetXaxis().SetLabelSize(0);
    h_data_func.GetXaxis().SetTitleSize(0);
    pad1.SetBottomMargin(0.03);
    pad1.SetTopMargin(0.1)
    h_data_func.SetLineColor(ROOT.kBlack)
    h_data_func.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)

    h_data_func.Draw('E1')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    
    
    # second pad
    c1.cd()
    h_data_func.GetXaxis().SetLabelOffset(0);
    h_data_func.GetXaxis().SetTitleOffset(0);
    pad2 = ROOT.TPad('pad2','pad2',0,0,1,0.35)
    pad2.SetGrid(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.GetYaxis().SetTitleSize(0.1)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.07)
    h_div_func.GetXaxis().SetTitleSize(0.15)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetTitleOffset(-999)
    h_div_func.GetXaxis().SetLabelOffset(0.02)
    h_div_func.GetXaxis().SetLabelSize(0.15)
    h_div_func.GetXaxis().SetBinLabel(1, bin0)
    h_div_func.GetXaxis().SetBinLabel(2, bin1)
    
    h_div_func.Draw('E1')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_weight_pad2(h_div_func):
    #h_div_func.SetMinimum(0.5)
    #h_div_func.SetMaximum(1.5)
    h_div_func.SetMaximum(-1111)
    h_div_func.SetMinimum(-1111)
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    return h_div_func

def make_stacked_histo_weight_pad2_fine(h_div_func):
    h_div_func.SetMinimum(0.8)
    h_div_func.SetMaximum(1.2)
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    return h_div_func

def draw_prelim():
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/15.);
    prelim.SetTextAlign(32);
    #prelim.SetTextSize(0.04631579);
    prelim.Draw()
    return

def draw_pot():
    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/15.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    pot_latex.Draw();
    return

def eff_channel_histo(cut,variable,xstart,xend,xbins,axis_name,file_name,side='right'):
    #overlay_signals_func = ['other','DIS','COH','RES','MEC','QE']
    overlay_signals_func = ['QE','MEC','RES','DIS']
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    globale.overlay_out.SetAlias('DIS','MCNu_Interaction==2')
    globale.overlay_out.SetAlias('COH','(MCNu_Interaction==3 || MCNu_Interaction==4)')
    globale.overlay_out.SetAlias('RES','MCNu_Interaction==1')
    globale.overlay_out.SetAlias('MEC','MCNu_Interaction==10')
    globale.overlay_out.SetAlias('QE','MCNu_Interaction==0')
    globale.overlay_out.SetAlias('other','!(DIS || COH || RES || MEC || QE)')
    
    h_overlay_func = {} # make an array of histograms for the 
    h_overlay_func_gen = {} # make an array of histograms for the 
    #tot_events = 0.0
    #channel_event = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,'title',xbins,xstart,xend)
        h_overlay_func_gen[x] = ROOT.TH1F(x+'gen','title',xbins,xstart,xend)
    
    #hs = ROOT.THStack("hs",variable+';'+axis_name+'; Entries per bin');
    cut = cut +'&& numu_signal && '
    cut_true = 'numu_true && '
    if variable == 'MCTrackPhi':
        cut_true = cut_true+'MCTrackPhi!=-1 && MCTrackPhi!=0 &&'
    if variable == 'MCTrackMomentum':
        cut_true = cut_true+'MCTrackMomentum!=-1 && MCTrackMomentum!=0 && '
        
    for i,x in enumerate(overlay_signals_func):
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        globale.overlay_out.Draw(variable+'>>'+histo+'gen',weight_name+'*('+cut_true+x+')','')
        h_overlay_func[x].Divide(h_overlay_func_gen[x])
        h_overlay_func[x].GetYaxis().SetTitleSize(0.05)
        h_overlay_func[x].GetYaxis().SetTitleOffset(0.0)
        h_overlay_func[x].GetYaxis().SetLabelSize(0.05)
        h_overlay_func[x].GetXaxis().SetTitleSize(0.05)
        h_overlay_func[x].GetXaxis().SetLabelSize(0.05)
        h_overlay_func[x].GetXaxis().SetTitleOffset(1)
        h_overlay_func[x].SetYTitle("Efficiency")
        h_overlay_func[x].SetXTitle(axis_name)
        h_overlay_func[x].SetLineWidth(4)
        #hs.Add(h_overlay_func[x])
    
    #h_overlay_func['other'].SetFillColor(ROOT.kGray+1)
    h_overlay_func['DIS'].SetFillColor(ROOT.kBlue+1)
    #h_overlay_func['COH'].SetFillColor(ROOT.kPink+10)
    h_overlay_func['RES'].SetFillColor(ROOT.kRed+1)
    h_overlay_func['MEC'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['QE'].SetFillColor(ROOT.kGreen+2)
    #h_overlay_func['other'].SetLineColor(ROOT.kGray+1)
    h_overlay_func['DIS'].SetLineColor(ROOT.kBlue+1)
    #h_overlay_func['COH'].SetLineColor(ROOT.kPink+10)
    h_overlay_func['RES'].SetLineColor(ROOT.kRed+1)
    h_overlay_func['MEC'].SetLineColor(ROOT.kOrange-3)
    h_overlay_func['QE'].SetLineColor(ROOT.kGreen+2)
    
    if side == 'left':
        legend = ROOT.TLegend(0.15,0.65,0.4,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.7,0.15,0.95,0.4); #LEGEND RIGHT
        
    for i,x in enumerate(overlay_signals_func):
        #ov_name = overlay_signals_func[i]+': {0:0.2f}%'.format(channel_event[i]*100.0/(tot_events+1e-10))
        legend.AddEntry(h_overlay_func[x],overlay_signals_func[i],"f");
        
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);

    #hs.SetYTitle("Entries per bin")
    h_overlay_func['DIS'].SetMinimum(0)
    h_overlay_func['DIS'].SetMaximum(1)
    
    #hs.SetLineColor(ROOT.kBlack)
    #hs.SetLineWidth(4)
    ROOT.gStyle.SetEndErrorSize(5)

    #h_overlay_func['other'].Draw('E1')
    h_overlay_func['DIS'].Draw('E1')
    #h_overlay_func['COH'].Draw('E1 same')
    h_overlay_func['RES'].Draw('E1 same')
    h_overlay_func['MEC'].Draw('E1 same')
    h_overlay_func['QE'].Draw('E1 same')
    legend.Draw()
    
    prelim.Draw()
    
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_channels_eff_"+file_name+".png")
    c1.SaveAs(globale.outputdir_root + "h_channels_eff_"+file_name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_channels_eff_"+file_name + ".pdf")
    
    return 

def make_channel_histo(cut,variable,xstart,xend,xbins,axis_name,file_name,side='right'):
    #overlay_signals_func = ['other','DIS','COH','RES','MEC','QE']
    overlay_signals_func = ['QE','MEC','RES','COH','DIS','other']
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    globale.overlay_out.SetAlias('DIS','MCNu_Interaction==2')
    globale.overlay_out.SetAlias('COH','(MCNu_Interaction==3 || MCNu_Interaction==4)')
    globale.overlay_out.SetAlias('RES','MCNu_Interaction==1')
    globale.overlay_out.SetAlias('MEC','MCNu_Interaction==10')
    globale.overlay_out.SetAlias('QE','MCNu_Interaction==0')
    globale.overlay_out.SetAlias('other','!(DIS || COH || RES || MEC || QE)')
    
    h_overlay_func = {} # make an array of histograms for the 
    tot_events = 0.0
    channel_event = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for x in overlay_signals_func:
        h_overlay_func[x] = ROOT.TH1F(x,'title',xbins,xstart,xend)
    
    hs = ROOT.THStack("hs",variable+';'+axis_name+'; Entries per bin');
    cut = cut +' && '
    for i,x in enumerate(overlay_signals_func):
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        channel_event[i] = h_overlay_func[x].GetSumOfWeights()
        tot_events += channel_event[i]

        h_overlay_func[x].GetYaxis().SetTitleSize(0.05)
        h_overlay_func[x].GetYaxis().SetTitleOffset(0.0)
        h_overlay_func[x].GetYaxis().SetLabelSize(0.05)
        h_overlay_func[x].GetXaxis().SetTitleSize(0.05)
        h_overlay_func[x].GetXaxis().SetLabelSize(0.05)
        h_overlay_func[x].GetXaxis().SetTitleOffset(1)
        h_overlay_func[x].SetYTitle("Entries per bin")
        h_overlay_func[x].SetXTitle(axis_name)
        hs.Add(h_overlay_func[x])
    
    h_overlay_func['other'].SetFillColor(ROOT.kGray+1)
    h_overlay_func['DIS'].SetFillColor(ROOT.kBlue+1)
    h_overlay_func['COH'].SetFillColor(ROOT.kPink+10)
    h_overlay_func['RES'].SetFillColor(ROOT.kRed+1)
    h_overlay_func['MEC'].SetFillColor(ROOT.kOrange-3)
    h_overlay_func['QE'].SetFillColor(ROOT.kGreen+2)
    h_overlay_func['other'].SetLineColor(ROOT.kGray+1)
    h_overlay_func['DIS'].SetLineColor(ROOT.kBlue+1)
    h_overlay_func['COH'].SetLineColor(ROOT.kPink+10)
    h_overlay_func['RES'].SetLineColor(ROOT.kRed+1)
    h_overlay_func['MEC'].SetLineColor(ROOT.kOrange-3)
    h_overlay_func['QE'].SetLineColor(ROOT.kGreen+2)
    
    if side == 'left':
        legend = ROOT.TLegend(0.15,0.65,0.4,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.7,0.65,0.95,0.9); #LEGEND RIGHT
        
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_func[i]+': {0:0.2f}%'.format(channel_event[i]*100.0/(tot_events+1e-10))
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
        
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);

    pot_latex = ROOT.TLatex(.10, .91,'Accumulated POT: '+str(pot_data)) 
    pot_latex.SetTextFont(62);
    pot_latex.SetTextColor(ROOT.kGray+2);
    pot_latex.SetNDC();
    pot_latex.SetTextSize(1/25.);
    pot_latex.SetTextAlign(10) #;//left adjusted
    
    #hs.SetYTitle("Entries per bin")
    hs.SetMinimum(0)
    hs.SetMaximum(hs.GetMaximum()*1.1)
    #hs.SetLabelSize(0.07);
    hs.SetHistogram( ROOT.TH1F("hstot","",xbins,xstart,xend))
    hs.GetHistogram().GetXaxis().SetTitle(axis_name)
    hs.GetHistogram().GetYaxis().SetTitle("Entries per bin")
    hs.GetHistogram().GetYaxis().SetTitleSize(0.05)
    hs.GetHistogram().GetYaxis().SetTitleOffset(0.0)
    hs.GetHistogram().GetYaxis().SetLabelSize(0.05)
    hs.GetHistogram().GetXaxis().SetTitleSize(0.05)
    hs.GetHistogram().GetXaxis().SetLabelSize(0.05)
    hs.GetHistogram().GetXaxis().SetTitleOffset(1)
    #hs.SetLineColor(ROOT.kBlack)
    #hs.SetLineWidth(4)
    #ROOT.gStyle.SetEndErrorSize(5)

    hs.Draw('histo')
    legend.Draw()
    hs.Draw('histo same')
    prelim.Draw()
    pot_latex.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_channels_"+file_name+".png")
    c1.SaveAs(globale.outputdir_root + "h_channels_"+file_name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_channels_"+file_name + ".pdf")
    
    return 

def printEff_out(cut):
    print "Efficiency for cut: ", cut
    # calculate the number of events including weights
    ext_pass =      globale.ext_out.GetEntries(cut+' && fidVol')
    dirt_pass =     globale.dirt_out.GetEntries(cut+' && fidVol')
    overlay_pass =  globale.overlay_out.GetEntries(cut+' && fidVol')
    overlay_sig =   globale.overlay_out.GetEntries(cut+' && numu_signal')
    ext_all =       globale.ext_out.GetEntries('muon && fidVol')
    dirt_all =      globale.dirt_out.GetEntries('muon && fidVol')
    overlay_all =   globale.overlay_out.GetEntries('muon && fidVol')
    overlay_true =  globale.overlay_out.GetEntries('numu_true')
    all_pass =      getTotNum_out(cut) # takes weights into account
    
    # print passing and rejection per sample:
    print '\nKeep of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass * 100.0 / ext_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - ext_pass * 100.0 / ext_all )
    print 'Keep of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass * 100.0 / dirt_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - dirt_pass * 100.0 / dirt_all )
    print 'Keep of %7s'%'overlay'+':\t{0:0.2f}%'.format( overlay_pass * 100.0 / overlay_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - overlay_pass * 100.0 / overlay_all )
    
    #calculate efficiency and purity:
    efficiency_func = overlay_sig*100.0/overlay_true
    purity_func = overlay_sig*globale.scale[globale.overlay]*100.0/all_pass
    print '\nEfficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}\n'.format(purity_func*efficiency_func/100)
    
    #print amount of each sample
    print 'Amount of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass*globale.scale_out[globale.ext_out]*100.0/all_pass)
    print 'Amount of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass*globale.scale_out[globale.dirt_out]*100.0/all_pass)
    print 'Amount of %7s'%'overlay'+':\t{0:0.2f}%\n'.format( overlay_pass*globale.scale_out[globale.overlay_out]*100.0/all_pass)
    
    #print the amount per signal definition:
    for x in globale.overlay_signals:
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay_out.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0/all_pass)
    
    return efficiency_func, purity_func

def printEff_w(cut):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    print "Efficiency for cut: ", cut
    print 'Used weight: ', weight_name
    h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',10000,0,1000)
    
    h_test_func = ROOT.TH1F("h_test_func",'h_test_func',5,-1,2)
    
    #calculating the mean weights
    globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && fidVol','0')
    overlay_pass_weight = h_weight_func.GetMean()  
    print overlay_pass_weight
    #print overlay_pass_weight/
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+' && fidVol)','0')
    print h_test_func.GetSumOfWeights()
    print h_test_func.GetSumOfWeights()/h_test_func.GetEntries()
    
    
    globale.overlay_out.Draw(weight_name+'>>h_weight_func','muon &&fidVol','0')
    overlay_reco = h_weight_func.GetMean()  
    globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && numu_signal','0') # weights for signal definition
    overlay_signal = h_weight_func.GetMean()
    globale.overlay_out.Draw(weight_name+'>>h_weight_func','numu_true','0') # weights for true signal
    overlay_true = h_weight_func.GetMean()
    
    globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && fidVol','0')
    dirt_pass_weight = h_weight_func.GetMean()  
    globale.overlay_out.Draw(weight_name+'>>h_weight_func','muon && fidVol','0')
    dirt_all_weight = h_weight_func.GetMean()  
    
    # calculate the number of events including weights
    ext_pass =      globale.ext_out.GetEntries(cut+' && fidVol')
    dirt_pass =     globale.dirt_out.GetEntries(cut+' && fidVol')*dirt_pass_weight
    overlay_pass =  globale.overlay_out.GetEntries(cut+' && fidVol')*overlay_pass_weight
    overlay_sig =   globale.overlay_out.GetEntries(cut+' && numu_signal')*overlay_signal
    ext_all =       globale.ext_out.GetEntries('muon && fidVol')
    dirt_all =      globale.dirt_out.GetEntries('muon && fidVol')*dirt_all_weight
    overlay_all =   globale.overlay_out.GetEntries('muon && fidVol')*overlay_reco
    overlay_true =  globale.overlay_out.GetEntries('numu_true')*overlay_true
    all_pass =      getTotNum_out_w(cut) # takes weights into account
    
    # print passing and rejection per sample:
    print '\nKeep of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass * 100.0 / ext_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - ext_pass * 100.0 / ext_all )
    print 'Keep of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass * 100.0 / dirt_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - dirt_pass * 100.0 / dirt_all )
    print 'Keep of %7s'%'overlay'+':\t{0:0.2f}%'.format( overlay_pass * 100.0 / overlay_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - overlay_pass * 100.0 / overlay_all )
    
    #calculate efficiency and purity:
    efficiency_func = overlay_sig*100.0/overlay_true
    purity_func = overlay_sig*globale.scale[globale.overlay]*100.0/all_pass
    print '\nEfficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}\n'.format(purity_func*efficiency_func/100)
    
    #print amount of each sample
    print 'Amount of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass*globale.scale_out[globale.ext_out]*100.0/all_pass )
    print 'Amount of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass*globale.scale_out[globale.dirt_out]*100.0/all_pass )
    print 'Amount of %7s'%'overlay'+':\t{0:0.2f}%\n'.format( overlay_pass*globale.scale_out[globale.overlay_out]*100.0/all_pass )
    
    #print the amount per signal definition:
    h_weight_func2 = ROOT.TH1F("h_weight_func2",'h_weight_func2',10000,0,1000)
    for x in globale.overlay_signals:
        globale.overlay_out.Draw(weight_name+'>>h_weight_func2',cut+' && '+x,'0')
        this_weight2 = h_weight_func2.GetMean()
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay_out.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0*this_weight2/all_pass)
    
    return efficiency_func, purity_func
  
def printEff_w2(cut):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    print "Efficiency for cut: ", cut
    print 'Used weight: ', weight_name
    #h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',10000,0,1000)
    
    h_test_func = ROOT.TH1F("h_test_func",'h_test_func',5,-1,2)
    
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func','muon &&fidVol','0')
    #overlay_reco = h_weight_func.GetMean()  
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && numu_signal','0') # weights for signal definition
    #overlay_signal = h_weight_func.GetMean()
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func','numu_true','0') # weights for true signal
    #overlay_true = h_weight_func.GetMean()
    
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && fidVol','0')
    #dirt_pass_weight = h_weight_func.GetMean()  
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func','muon && fidVol','0')
    #dirt_all_weight = h_weight_func.GetMean() 
    
    
    
    # calculate the number of events including weights
    ext_pass =      globale.ext_out.GetEntries(cut+' && fidVol')
    globale.dirt_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+' && fidVol)','0')
    dirt_pass =     h_test_func.GetSumOfWeights()
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+' && fidVol)','0')
    overlay_pass =  h_test_func.GetSumOfWeights()
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+' && numu_signal)','0')
    overlay_sig =   h_test_func.GetSumOfWeights()
    ext_all =       globale.ext_out.GetEntries('muon && fidVol')
    globale.dirt_out.Draw('fidVol'+'>>h_test_func',weight_name+'*(muon && fidVol)','0')
    dirt_all =      h_test_func.GetSumOfWeights()
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*(muon && fidVol)','0')
    overlay_all =   h_test_func.GetSumOfWeights()
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*(numu_true)','0')
    overlay_true =  h_test_func.GetSumOfWeights()
    all_pass =      getTotNum_out_w2(cut) # takes weights into account
    
    # print passing and rejection per sample:
    print '\nKeep of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass * 100.0 / ext_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - ext_pass * 100.0 / ext_all )
    print 'Keep of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass * 100.0 / dirt_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - dirt_pass * 100.0 / dirt_all )
    print 'Keep of %7s'%'overlay'+':\t{0:0.2f}%'.format( overlay_pass * 100.0 / overlay_all ) + '\t reject:\t{0:0.2f}%'.format( 100.0 - overlay_pass * 100.0 / overlay_all )
    
    #calculate efficiency and purity:
    efficiency_func = overlay_sig*100.0/overlay_true
    purity_func = overlay_sig*globale.scale[globale.overlay]*100.0/all_pass
    print '\nEfficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}\n'.format(purity_func*efficiency_func/100)
    
    #print amount of each sample
    print 'Amount of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass*globale.scale_out[globale.ext_out]*100.0/all_pass )
    print 'Amount of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass*globale.scale_out[globale.dirt_out]*100.0/all_pass )
    print 'Amount of %7s'%'overlay'+':\t{0:0.2f}%\n'.format( overlay_pass*globale.scale_out[globale.overlay_out]*100.0/all_pass )
    
    #print the amount per signal definition:
    #h_weight_func2 = ROOT.TH1F("h_weight_func2",'h_weight_func2',10000,0,1000)
    test = 0.0
    h_test_func2 = ROOT.TH1F("h_test_func2",'h_test_func2',5,-1,2)
    for x in globale.overlay_signals:
        #globale.overlay_out.Draw(weight_name+'>>h_weight_func2',cut+' && '+x,'0')
        #this_weight2 = h_weight_func2.GetMean()
        globale.overlay_out.Draw('fidVol'+'>>h_test_func2',weight_name+'*('+cut+' && '+x+')','0')
        this_value = h_test_func2.GetSumOfWeights()
        #print this_value, all_pass, globale.scale[globale.overlay]
        test = test+this_value*globale.scale[globale.overlay]*100.0/all_pass
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(this_value*globale.scale[globale.overlay]*100.0/all_pass)
    print test
      
    #globale.overlay_out.Draw('fidVol'+'>>h_test_func2',weight_name+'*('+cut+' && !numu_cosmic && !numu_ov && !numu_nc && !numu_antinu && !numu_nue && !numu_nomuon && !numu_signal)','0')
    #this_value = h_test_func2.GetSumOfWeights()
    #print this_value
    #test = test+this_value*globale.scale[globale.overlay]*100.0/all_pass
    return efficiency_func, purity_func

# calculates the total number of events passing the cut
# It takes the EventWeight*TunedCentralValue_Genie weights into account (for overlay)
def getTotNum_out_w(cut):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    num_fidVol = 0.0
    h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',1000,0,100)
    num_fidVol = num_fidVol+globale.ext_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.ext_out]               # number of events in ext added
    #num_fidVol = num_fidVol+globale.dirt_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.dirt_out]             # number of events in dirt added
    globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut,'0')                               # calculating the mean weights for overlay
    weight = h_weight_func.GetMean()
    num_fidVol = num_fidVol+globale.overlay_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.overlay_out]*weight # weighted overlay events added
    globale.dirt_out.Draw(weight_name+'>>h_weight_func',cut,'0')                               # calculating the mean weights for overlay
    weight = h_weight_func.GetMean()
    num_fidVol = num_fidVol+globale.dirt_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.dirt_out]*weight # weighted dirt events added
    return num_fidVol
  
def getTotNum_out_w2(cut):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    num_fidVol = 0.0
    #h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',1000,0,100)
    h_test_func = ROOT.TH1F("h_test_func",'h_test_func',5,-1,2)

    num_fidVol = num_fidVol+globale.ext_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.ext_out]               # number of events in ext added
    #num_fidVol = num_fidVol+globale.dirt_out.GetEntries(cut + ' && fidVol')*globale.scale_out[globale.dirt_out]             # number of events in dirt added
    #globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut,'0')                               # calculating the mean weights for overlay
    #weight = h_weight_func.GetMean()
    globale.overlay_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+')','0')
    this_value = h_test_func.GetSumOfWeights()
    
    num_fidVol = num_fidVol+this_value*globale.scale_out[globale.overlay_out] #weighted overlay events added
    #globale.dirt_out.Draw(weight_name+'>>h_weight_func',cut,'0')                               # calculating the mean weights for overlay
    #weight = h_weight_func.GetMean()
    
    globale.dirt_out.Draw('fidVol'+'>>h_test_func',weight_name+'*('+cut+')','0')
    this_value = h_test_func.GetSumOfWeights()
    
    num_fidVol = num_fidVol+this_value*globale.scale_out[globale.dirt_out] # weighted dirt events added
    return num_fidVol

# calculates the total number of events passing the cut, ignoring any weights
def getTotNum_out(cut):
    num_fidVol = 0.0
    for x in globale.sample_out:
        if x != globale.data_out:
            num_fidVol = num_fidVol+x.GetEntries(cut + ' && fidVol')*globale.scale_out[x]
    return num_fidVol


def plot_eff_w(nenner_cut, zahler_cut, cut, name, title):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    # define weight histogram
    h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',10000,0,1000)
    #get weights of events with cut
    globale.overlay_out.Draw('EventWeight*TunedCentralValue_Genie'+'>>h_weight_func',cut+' && numu_signal','0')
    weight_cut = h_weight_func.GetMean()
    #get weights for all true events
    globale.overlay_out.Draw('EventWeight*TunedCentralValue_Genie'+'>>h_weight_func','numu_true','0')
    weight = h_weight_func.GetMean()
    #calculate efficiency
    efficiency = globale.overlay_out.GetEntries(cut+' && numu_signal')*100.0*weight_cut /(globale.overlay_out.GetEntries('numu_true')*weight+0.001)
    #purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100.0/getTotNum_out(cut)
    purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100*weight_cut /(getTotNum_out_w(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/25.);
    prelim.SetTextAlign(32);

    xstart = 0
    xend = 2
    xbin = 50
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")
    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino energy [GeV]")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_energy"+name + ".pdf")
    
    xstart = 0.00001
    xend = 2
    xbin = 50
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw('MCTrackMomentum>>h_init_eff_energy',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('MCTrackMomentum>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")
    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth muon momentum [GeV]")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_momentum"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_momentum"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_momentum"+name + ".pdf")

    xstart = -1
    xend = 1
    xbin = 50
    h_init_eff_theta = ROOT.TH1F("h_init_eff_theta",title,xbin,xstart,xend)
    h_init_eff_theta_1 = ROOT.TH1F("h_init_eff_theta_1",title,xbin,xstart,xend)
    #h_init_eff_theta = fill_histo(globale.overlay_out,'cos(MCNu_leptonTheta)',h_init_eff_theta,zahler_cut)
    #h_init_eff_theta_1 = fill_histo(globale.overlay_out,'cos(MCNu_leptonTheta)',h_init_eff_theta_1,nenner_cut)
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_theta,h_init_eff_theta_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth lepton cos(theta)")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_theta"+name+".pdf")

    xstart = -3.14159
    xend = 3.14159
    xbin = 50
    h_init_eff_phi = ROOT.TH1F("h_init_eff_phi",title,xbin,xstart,xend)
    h_init_eff_phi_1 = ROOT.TH1F("h_init_eff_phi_1",title,xbin,xstart,xend)
    #h_init_eff_phi = fill_histo(globale.overlay_out,'TrackPhi',h_init_eff_phi,zahler_cut)
    #h_init_eff_phi_1 = fill_histo(globale.overlay_out,'TrackPhi',h_init_eff_phi_1,nenner_cut)
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi',weight_name+'*('+zahler_cut+' && muon'+')')
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi_1',weight_name+'*('+nenner_cut+' && muon'+')')
    eff =  ROOT.TEfficiency(h_init_eff_phi,h_init_eff_phi_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Reco lepton phi [pi]")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_phi"+name+".pdf")

    xstart = -50
    xend = 300
    xbin = 50
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    #h_init_eff_V = fill_histo(globale.overlay_out,'MCNu_Vx',h_init_eff_V,zahler_cut)
    #h_init_eff_V_1 = fill_histo(globale.overlay_out,'MCNu_Vx',h_init_eff_V_1,nenner_cut)
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino vertex X [cm]")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_Vx"+name+".pdf")

    xstart = -150
    xend = 150
    xbin = 50
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    #h_init_eff_V = fill_histo(globale.overlay_out,'MCNu_Vy',h_init_eff_V,zahler_cut)
    #h_init_eff_V_1 = fill_histo(globale.overlay_out,'MCNu_Vy',h_init_eff_V_1,nenner_cut)
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino vertex Y [cm]")
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_Vy"+name+".pdf")

    xstart = -50
    xend = 1000
    xbin = 50
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    #h_init_eff_V = fill_histo(globale.overlay_out,'MCNu_Vz',h_init_eff_V,zahler_cut)
    #h_init_eff_V_1 = fill_histo(globale.overlay_out,'MCNu_Vz',h_init_eff_V_1,nenner_cut)
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V',weight_name+'*('+zahler_cut+')')
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V_1',weight_name+'*('+nenner_cut+')')
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle('Truth neutrino vertex Z [cm]')
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
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_Vz"+name+".pdf")
    return

def plot_eff_all_w(variable,xstart,xend,xbin,xaxis_name, name,side):
    title = 'efficiency'
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);

    nenner_cut = 'numu_true'
    precut = 'fidVol && muon && numu_signal'
    if(variable=='MCTrackPhi'):
        nenner_cut = nenner_cut+' && MCTrackPhi!=-1'
        precut = precut+' && MCTrackPhi!=-1'
    #quality_cut = ' && TrackLength>8'
    crt_cut = ' && crt_tom_cut'
    trackscore_cut = ' && TrackScore>0.8'
    tracklength_cut = ' && TrackLength>20'
    trackPID_cut = ' && TrackPID_chiproton>78'
    nuscore_cut = ' && NuScore>0.1'
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    '''eff_pre =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff_pre.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff_pre.SetConfidenceLevel(0.68)'''
    eff_pre = h_init_eff_energy.Clone()
    eff_pre.Divide(h_init_eff_energy_1)
    
    #precut = precut + quality_cut
    #globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    #globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    #eff_qc = h_init_eff_energy.Clone()
    #eff_qc.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt = h_init_eff_energy.Clone()
    eff_crt.Divide(h_init_eff_energy_1)
    
    precut = precut + trackscore_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trkscore = h_init_eff_energy.Clone()
    eff_trkscore.Divide(h_init_eff_energy_1)
    
    precut = precut + tracklength_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trklength = h_init_eff_energy.Clone()
    eff_trklength.Divide(h_init_eff_energy_1)
    
    precut = precut + trackPID_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trkPID = h_init_eff_energy.Clone()
    eff_trkPID.Divide(h_init_eff_energy_1)
    
    precut = precut + nuscore_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_nuscore = h_init_eff_energy.Clone()
    eff_nuscore.Divide(h_init_eff_energy_1)
    
    if side == 'left':
        legend = ROOT.TLegend(0.3,0.15,0.85,0.3) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    legend.SetNColumns(3)
    legend.AddEntry(eff_pre,'SliceID',"f")
    #legend.AddEntry(eff_qc,'Quality cut',"f")
    legend.AddEntry(eff_crt,'CRT cuts',"f")
    legend.AddEntry(eff_trkscore,'Track score',"f")
    legend.AddEntry(eff_trklength,'Track length',"f")
    legend.AddEntry(eff_trkPID,'Track PID',"f")
    legend.AddEntry(eff_nuscore,'Topological score',"f")
    
    eff_pre.SetMinimum(0)
    eff_pre.SetMaximum(1)
    eff_pre.SetLineWidth(1)
    eff_pre.GetXaxis().SetTitle(xaxis_name)
    eff_pre.GetYaxis().SetTitle("Signal efficiency")
    eff_pre.GetYaxis().SetTitleSize(0.05)
    eff_pre.GetYaxis().SetTitleOffset(0.0)
    eff_pre.GetYaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleSize(0.05)
    eff_pre.GetXaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleOffset(1)
    eff_pre.SetLineColor(ROOT.kOrange)
    eff_pre.SetLineWidth(3)
    ROOT.gStyle.SetEndErrorSize(3)
    eff_pre.Draw("E1")
    legend.Draw()
    eff_pre.SetLineColor(ROOT.kRed)
    eff_pre.Draw("E1 same")
    #eff_qc.SetLineColor(ROOT.kRed)
    #eff_qc.SetLineWidth(3)
    #eff_qc.Draw("E1 same")
    eff_crt.SetLineColor(ROOT.kBlue)
    eff_crt.SetLineWidth(3)
    eff_crt.Draw("E1 same")
    eff_trkscore.SetLineColor(ROOT.kGreen-2)
    eff_trkscore.SetLineWidth(3)
    eff_trkscore.Draw("E1 same")
    eff_trklength.SetLineColor(ROOT.kYellow)
    eff_trklength.SetLineWidth(3)
    eff_trklength.Draw("E1 same")
    eff_trkPID.SetLineColor(ROOT.kPink-9)
    eff_trkPID.SetLineWidth(3)
    eff_trkPID.Draw("E1 same")
    eff_nuscore.SetLineColor(ROOT.kBlack)
    eff_nuscore.SetLineWidth(3)
    eff_nuscore.Draw("E1 same")
    
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_"+name+".pdf")
    return

def plot_eff_all_rel_w(variable,xstart,xend,xbin,xaxis_name, name,side):
    title = 'efficiency'
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);

    nenner_cut = 'fidVol && muon && numu_signal'
    precut = 'fidVol && muon && numu_signal'
    if(variable=='MCTrackPhi'):
        nenner_cut = nenner_cut+' && MCTrackPhi!=-1'
        precut = precut+' && MCTrackPhi!=-1'
    #quality_cut = ' && TrackLength>8'
    crt_cut = ' && crt_tom_cut'
    trackscore_cut = ' && TrackScore>0.8'
    tracklength_cut = ' && TrackLength>20'
    trackPID_cut = ' && TrackPID_chiproton>78'
    nuscore_cut = ' && NuScore>0.1'
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    '''eff_pre =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff_pre.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff_pre.SetConfidenceLevel(0.68)'''
    eff_pre = h_init_eff_energy.Clone()
    eff_pre.Divide(h_init_eff_energy_1)
    
    #precut = precut + quality_cut
    #globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    #globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    #eff_qc = h_init_eff_energy.Clone()
    #eff_qc.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt = h_init_eff_energy.Clone()
    eff_crt.Divide(h_init_eff_energy_1)
    
    precut = precut + trackscore_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trkscore = h_init_eff_energy.Clone()
    eff_trkscore.Divide(h_init_eff_energy_1)
    
    precut = precut + tracklength_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trklength = h_init_eff_energy.Clone()
    eff_trklength.Divide(h_init_eff_energy_1)
    
    precut = precut + trackPID_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_trkPID = h_init_eff_energy.Clone()
    eff_trkPID.Divide(h_init_eff_energy_1)
    
    precut = precut + nuscore_cut
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_nuscore = h_init_eff_energy.Clone()
    eff_nuscore.Divide(h_init_eff_energy_1)
    
    if side == 'left':
        legend = ROOT.TLegend(0.3,0.15,0.85,0.3) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    legend.SetNColumns(3)
    #legend.AddEntry(eff_pre,'SliceID',"f")
    #legend.AddEntry(eff_qc,'Quality cut',"f")
    legend.AddEntry(eff_crt,'CRT cuts',"f")
    legend.AddEntry(eff_trkscore,'Track score',"f")
    legend.AddEntry(eff_trklength,'Track length',"f")
    legend.AddEntry(eff_trkPID,'Track PID',"f")
    legend.AddEntry(eff_nuscore,'Topological score',"f")
    
    if (eff_nuscore.GetMinimum() > 0.1):
      eff_pre.SetMinimum(eff_nuscore.GetMinimum()*0.7)
    else:
      eff_pre.SetMinimum(0.3)
    eff_pre.SetMaximum(1)
    eff_pre.SetLineWidth(1)
    eff_pre.GetXaxis().SetTitle(xaxis_name)
    eff_pre.GetYaxis().SetTitle("Signal efficiency")
    eff_pre.GetYaxis().SetTitleSize(0.05)
    eff_pre.GetYaxis().SetTitleOffset(0.0)
    eff_pre.GetYaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleSize(0.05)
    eff_pre.GetXaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleOffset(1)
    eff_pre.SetLineColor(ROOT.kOrange)
    eff_pre.SetLineWidth(3)
    ROOT.gStyle.SetEndErrorSize(3)
    eff_pre.Draw("E1")
    legend.Draw()
    eff_pre.SetLineColorAlpha(ROOT.kRed,0.0)
    eff_pre.Draw("E1 same")
    #eff_qc.SetLineColor(ROOT.kRed)
    #eff_qc.SetLineWidth(3)
    #eff_qc.Draw("E1 same")
    eff_crt.SetLineColor(ROOT.kBlue)
    eff_crt.SetLineWidth(3)
    eff_crt.Draw("E1 same")
    eff_trkscore.SetLineColor(ROOT.kGreen-2)
    eff_trkscore.SetLineWidth(3)
    eff_trkscore.Draw("E1 same")
    eff_trklength.SetLineColor(ROOT.kYellow)
    eff_trklength.SetLineWidth(3)
    eff_trklength.Draw("E1 same")
    eff_trkPID.SetLineColor(ROOT.kPink-9)
    eff_trkPID.SetLineWidth(3)
    eff_trkPID.Draw("E1 same")
    eff_nuscore.SetLineColor(ROOT.kBlack)
    eff_nuscore.SetLineWidth(3)
    eff_nuscore.Draw("E1 same")
    
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_rel_"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_rel_"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_eff_rel_"+name+".pdf")
    return
  
def plot_eff_crt_w(variable,xstart,xend,xbin,xaxis_name, name,side):
    title = 'efficiency'
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);

    nenner_cut = 'numu_true'
    precut = 'fidVol && muon && numu_signal'
    if(variable=='MCTrackPhi'):
        nenner_cut = nenner_cut+' && MCTrackPhi!=-1'
        precut = precut+' && MCTrackPhi!=-1'
    #quality_cut = ' && TrackLength>8'
    crt_cut_1 = ' && nr_crthit_top==0'
    crt_cut_2 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0'
    crt_cut_3 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0)'
    crt_cut_4 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut'
    
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')

    eff_pre = h_init_eff_energy.Clone()
    eff_pre.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_1
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_1 = h_init_eff_energy.Clone()
    eff_crt_1.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_2
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_2 = h_init_eff_energy.Clone()
    eff_crt_2.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_3
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_3 = h_init_eff_energy.Clone()
    eff_crt_3.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_4
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_4 = h_init_eff_energy.Clone()
    eff_crt_4.Divide(h_init_eff_energy_1)
    
    if side == 'left':
        legend = ROOT.TLegend(0.3,0.15,0.85,0.3) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    legend.SetNColumns(3)
    legend.AddEntry(eff_pre,'SliceID',"f")
    #legend.AddEntry(eff_qc,'Quality cut',"f")
    legend.AddEntry(eff_crt_1,'Top veto',"f")
    legend.AddEntry(eff_crt_2,'CRT-TPC Z pos',"f")
    legend.AddEntry(eff_crt_3,'CRT veto (cont.)',"f")
    legend.AddEntry(eff_crt_4,'CRT asso.',"f")

    
    eff_pre.SetMinimum(0)
    eff_pre.SetMaximum(1)
    eff_pre.SetLineWidth(1)
    eff_pre.GetXaxis().SetTitle(xaxis_name)
    eff_pre.GetYaxis().SetTitle("Signal efficiency")
    eff_pre.GetYaxis().SetTitleSize(0.05)
    eff_pre.GetYaxis().SetTitleOffset(0.0)
    eff_pre.GetYaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleSize(0.05)
    eff_pre.GetXaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleOffset(1)
    eff_pre.SetLineColor(ROOT.kOrange)
    eff_pre.SetLineWidth(3)
    ROOT.gStyle.SetEndErrorSize(3)
    eff_pre.Draw("E1")
    legend.Draw()
    eff_pre.SetLineColor(ROOT.kRed)
    eff_pre.Draw("E1 same")
    #eff_qc.SetLineColor(ROOT.kRed)
    #eff_qc.SetLineWidth(3)
    #eff_qc.Draw("E1 same")
    eff_crt_1.SetLineColor(ROOT.kBlue)
    eff_crt_1.SetLineWidth(3)
    eff_crt_1.Draw("E1 same")
    eff_crt_2.SetLineColor(ROOT.kGreen-2)
    eff_crt_2.SetLineWidth(3)
    eff_crt_2.Draw("E1 same")
    eff_crt_3.SetLineColor(ROOT.kYellow-1)
    eff_crt_3.SetLineWidth(3)
    eff_crt_3.Draw("E1 same")
    eff_crt_4.SetLineColor(ROOT.kPink-9)
    eff_crt_4.SetLineWidth(3)
    eff_crt_4.Draw("E1 same")

    
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_crt_eff_"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_crt_eff_"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_crt_eff_"+name+".pdf")
    return

def plot_eff_crt_rel_w(variable,xstart,xend,xbin,xaxis_name, name,side):
    title = 'efficiency'
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.14)
    
    prelim = ROOT.TLatex(0.9,0.93, "MicroBooNE Simulation Preliminary");
    prelim.SetTextFont(62);
    prelim.SetTextColor(ROOT.kGray+2);
    prelim.SetNDC();
    prelim.SetTextSize(1/20.);
    prelim.SetTextAlign(32);

    nenner_cut = 'fidVol && muon && numu_signal'
    precut = 'fidVol && muon && numu_signal'
    if(variable=='MCTrackPhi'):
        nenner_cut = nenner_cut+' && MCTrackPhi!=-1'
        precut = precut+' && MCTrackPhi!=-1'
    #quality_cut = ' && TrackLength>8'
    crt_cut_1 = ' && nr_crthit_top==0'
    crt_cut_2 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0'
    crt_cut_3 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0)'
    crt_cut_4 = ' && nr_crthit_top==0 && crthit_vertex_zcut==0 && (track_end_uncontained==1 || nr_crthit_beam_tres==0) && crt_cut'
    
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    #h_init_eff_energy = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy,zahler_cut)
    #h_init_eff_energy_1 = fill_histo(globale.overlay_out,'MCNu_Energy',h_init_eff_energy_1,nenner_cut)
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')

    eff_pre = h_init_eff_energy.Clone()
    eff_pre.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_1
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_1 = h_init_eff_energy.Clone()
    eff_crt_1.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_2
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_2 = h_init_eff_energy.Clone()
    eff_crt_2.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_3
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_3 = h_init_eff_energy.Clone()
    eff_crt_3.Divide(h_init_eff_energy_1)
    
    precut = precut + crt_cut_4
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy',weight_name+'*('+precut+')')
    globale.overlay_out.Draw(variable+'>>h_init_eff_energy_1',weight_name+'*('+nenner_cut+')')
    eff_crt_4 = h_init_eff_energy.Clone()
    eff_crt_4.Divide(h_init_eff_energy_1)
    
    if side == 'left':
        legend = ROOT.TLegend(0.3,0.15,0.85,0.3) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.6,0.30,0.9,0.9)
    legend.SetNColumns(2)
    #legend.AddEntry(eff_pre,'SliceID',"f")
    #legend.AddEntry(eff_qc,'Quality cut',"f")
    legend.AddEntry(eff_crt_1,'Top veto',"f")
    legend.AddEntry(eff_crt_2,'CRT-TPC Z pos',"f")
    legend.AddEntry(eff_crt_3,'CRT veto (cont.)',"f")
    legend.AddEntry(eff_crt_4,'CRT asso.',"f")

    if (eff_crt_4.GetMinimum() > 0.1):
      eff_pre.SetMinimum(eff_crt_4.GetMinimum()*0.9)
    else:
      eff_pre.SetMinimum(0.8)
    eff_pre.SetMaximum(1)
    eff_pre.SetLineWidth(1)
    eff_pre.GetXaxis().SetTitle(xaxis_name)
    eff_pre.GetYaxis().SetTitle("Signal efficiency")
    eff_pre.GetYaxis().SetTitleSize(0.05)
    eff_pre.GetYaxis().SetTitleOffset(0.0)
    eff_pre.GetYaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleSize(0.05)
    eff_pre.GetXaxis().SetLabelSize(0.05)
    eff_pre.GetXaxis().SetTitleOffset(1)
    eff_pre.SetLineColor(ROOT.kOrange)
    eff_pre.SetLineWidth(3)
    ROOT.gStyle.SetEndErrorSize(0)
    eff_pre.Draw("l")
    legend.Draw()
    eff_pre.SetLineColorAlpha(ROOT.kRed,0.0)
    eff_pre.Draw("l same")
    #eff_qc.SetLineColor(ROOT.kRed)
    #eff_qc.SetLineWidth(3)
    #eff_qc.Draw("E1 same")
    ROOT.gStyle.SetEndErrorSize(3)
    eff_crt_1.SetLineColor(ROOT.kBlue)
    eff_crt_1.SetLineWidth(3)
    eff_crt_1.Draw("E1 same")
    eff_crt_2.SetLineColor(ROOT.kGreen-2)
    eff_crt_2.SetLineWidth(3)
    eff_crt_2.Draw("E1 same")
    eff_crt_3.SetLineColor(ROOT.kYellow-1)
    eff_crt_3.SetLineWidth(3)
    eff_crt_3.Draw("E1 same")
    eff_crt_4.SetLineColor(ROOT.kPink-9)
    eff_crt_4.SetLineWidth(3)
    eff_crt_4.Draw("E1 same")

    
    prelim.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_crt_eff_rel_"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_crt_eff_rel_"+name+".root")
    c1.SaveAs(globale.outputdir_pdf+ "h_crt_eff_rel_"+name+".pdf")
    return

def pdg_content33_w(cut):
    xstart = -3000
    xend = 3000
    xbin = 6000
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay_out.Draw('MCTrackPDG>>h_1d',weight_name+'*('+cut+')','')
    h_1d.SetXTitle("Track PDG code")
    h_1d.SetYTitle("number of entries")
    tot = h_1d.GetSumOfWeights()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100.0/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100.0/tot),',\tnumber: ',particle[x+1][0]
    return

#calculate the bin content for a histogram in variable for all universes
def calc_genie_systematic(cut,variable,xstart,xend,xbin,num_universes = 100):
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_1d = {}
    for uni in range(num_universes):
        #print 'At universe: ',uni
        h_1d[uni] = ROOT.TH1F("h_1d["+str(uni)+']',"Track cos(theta)",xbin, xstart, xend)
        globale.overlay_out.Draw(variable+">>h_1d["+str(uni)+']',weight_name+'*All_Genie['+str(uni)+']*('+cut+')')
    bin_entry = np.zeros((xbin,num_universes))
    for i in range(xbin):
        for uni in range(num_universes):
            bin_entry[i][uni] = h_1d[uni].GetBinContent(i+1)
    #print bin_entry
    mean = np.mean(bin_entry,axis=1)
    std = np.std(bin_entry,axis=1)
    print mean.size
    print std.size  
         
    return mean, std, bin_entry

def make_stacked_histo_sys(cut,variable,weight,title,xstart,xend,xbins,file_name,side,this_std):
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_all_func = ROOT.TH1F("h_all_func",title,xbins,xstart,xend)
    h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    #variable = variable
    # Draw/Fill the histograms for all kind of interactions
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    #mofify this...
    globale.dirt_out.Draw(variable+'>>h_dirt_func',weight_name+'*('+cut+')','')

    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,weight_name+'*('+cut+x+')','')
        #h_overlay_func[x] = fill_histo(globale.overlay_out,variable,h_overlay_func[x],cut+x)
    # prepare the stacked histogram
    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(1)
    #scale the histograms
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    #fill the stacked histogram
    hs.Add(h_ext_func)
    h_all_func.Add(h_ext_func)
    hs.Add(h_dirt_func)
    h_all_func.Add(h_dirt_func)
    mc_events = 0
    mc_event_list = {}
    for i,x in enumerate(globale.overlay_signals):
        #h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        mc_event_list[x] = h_overlay_func[x].GetSumOfWeights()
        mc_events = mc_events+mc_event_list[x]
        hs.Add(h_overlay_func[x])
        h_all_func.Add(h_overlay_func[x])
    # calculate the data - MC ratio
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    dirt_events = h_dirt_func.GetSumOfWeights()
    mc_events = mc_events + dirt_events
    normalization = (data_events)/(mc_events+ext_events)
    print 'Normalization (data)/(mc +ext) = ', normalization
    if side == 'left':
        legend = ROOT.TLegend(0.15,0.65,0.45,0.9) # LEGEND LEFT
    else:
        legend = ROOT.TLegend(0.6,0.65,0.9,0.9); #LEGEND RIGHT
    legend.SetNColumns(2)
    data_name = 'data: {0:0.1f}'.format(data_events)
    ext_name = 'ext: {0:0.1f}'.format(ext_events)
    dirt_name = 'dirt: {0:0.1f}'.format(dirt_events)
    legend.AddEntry(h_data_func,data_name,"lep");
    legend.AddEntry(h_ext_func,ext_name,"f");
    legend.AddEntry(h_dirt_func,dirt_name,"f");
    for x in globale.overlay_signals:
        ov_name = x+': {0:0.1f}'.format(mc_event_list[x])
        legend.AddEntry(h_overlay_func[x],ov_name,"f");
    #prepare the canvas with thw pads
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    # first pad
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.Draw('E')
    #legend.Draw();
    hs.Draw('same hist')
    #h_all_func = hs.GetHistogram()
    for i in range(xbins):
        #test.SetBinContent(i+1,mean[i])
        h_all_func.SetBinError(i+1,this_std[i]*globale.scale[globale.overlay])
    h_all_func.SetFillStyle(3005)
    h_all_func.SetFillColor(1);
    h_all_func.Draw('E2 same')
    
    h_data_func.Draw('E same')
    legend.Draw('same')

    # second pad
    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    # Draw data - MC difference
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    h_all_func2 = h_all_func.Clone()
    h_all_func2.Divide(h_all_func)
    
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.Draw()
    h_all_func2.SetFillStyle(3005)
    h_all_func2.SetFillColor(1);
    h_all_func2.Draw('E2 same')
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    c1.SaveAs(globale.outputdir_pdf+ file_name + ".pdf")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization
