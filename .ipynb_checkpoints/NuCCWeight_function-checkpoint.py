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
globale = imp.load_source('globale','/home/tmettler/Desktop/uBoone/do_plots/globale.py')


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
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
        
    return normalization

def make_stacked_histo_weightV2(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
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
    hs.Add(h_dirt_func)
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
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization

def make_stacked_histo_weight_MCC8(cut,variable,weight,title,xstart,xend,xbins,file_name,side='right'):
    #initialize the 1d histograms
    weight_name = 'EventWeight*TunedCentralValue_Genie'
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
    ext_name = 'Data (Beam off): {0:0.1f}%'.format(ext_events*100/(mc_events+ext_events))
    dirt_name = 'Dirt: {0:0.1f}%'.format(dirt_events*100/(mc_events+ext_events))
    #overlay_signals_name = ['\nu_{\mu} CC (stopping \mu)','\nu_{\mu} CC (other)','\nu_e, \overline{\nu_e}','\overline{\nu_{\mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping \mu)', 'OUTFV (other)', 'Cosmic (stopping \mu)', 'Cosmic (other)']
    overlay_signals_name = ['#nu_{#mu} CC (stopping #mu)','#nu_{#mu} CC (other)','#nu_{e}, #bar{#nu_{e}}','#bar{#nu_{#mu}} CC','NC (other)','NC (pion)', 'NC (proton)', 'OUTFV (stopping #mu)', 'OUTFV (other)', 'Cosmic (stopping #mu)', 'Cosmic (other)']
    
    #overlay_signals_name = ['Cosmic (other)','Cosmic (stopping \mu)','OUTFV (other)','OUTFV (stopping \mu)', 'NC (proton)', 'NC (pion)', 'NC (other)', '\overline{\nu_{\mu}} CC', '\nu_e, \overline{\nu_e}', '\nu_{\mu} CC (other)', '\nu_{\mu} CC (stopping \mu)']
    
    
    #overlay_signals_name = overlay_signals_name.reverse()
    overlay_signals_func = ['numu_signal_cont','numu_signal_uncont','numu_nue','numu_antinu','numu_nc_other','numu_nc_pion',\
                   'numu_nc_proton', 'numu_ov_cont', 'numu_ov_uncont', 'numu_cosmic_cont', 'numu_cosmic_uncont']

    
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100/(mc_events+ext_events))
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
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    h_data_func.Draw('E1')
    #legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    if side != 'none':
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
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.SetMarkerStyle(ROOT.kFullCircle);
    h_div_func.SetMarkerSize(0.9);
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
    ext_name = 'Background (Beam off): {0:0.1f}%'.format(ext_events*100/(mc_events+ext_events))
    dirt_name = 'Background (Dirt): {0:0.1f}%'.format(dirt_events*100/(mc_events+ext_events))
    
    for i,x in enumerate(overlay_signals_func):
        ov_name = overlay_signals_name[i]+': {0:0.1f}%'.format(mc_event_list[x]*100/(mc_events+ext_events))
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
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    # draw fisrt histogram with data points and stacked ext and MC
    h_data_func.SetYTitle("Entries per bin")
    h_data_func.SetMinimum(0)
    max = h_data_func.GetMaximum()
    h_data_func.SetMaximum(max*1.3)
    h_data_func.Draw('E1')
    #legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E1 same')
    if side != 'none':
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
    for i,x in enumerate(overlay_signals_func):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func )
    #h_test = hs.GetHistogram().Clone()
    #h_div_func.Divide(h_test )
    h_div_func = make_stacked_histo_weight_pad2(h_div_func)
    h_div_func.SetXTitle(title)
    h_div_func.SetMarkerStyle(ROOT.kFullCircle);
    h_div_func.SetMarkerSize(0.9);
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



def make_stacked_histo_weight_pad2(h_div_func):
    h_div_func.SetMinimum(0.0)
    h_div_func.SetMaximum(2.0)
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    return h_div_func

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
    purity_func = overlay_sig*globale.scale[globale.overlay]*100/all_pass
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
    h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',1000,0,100)
    
    #calculating the mean weights
    globale.overlay_out.Draw(weight_name+'>>h_weight_func',cut+' && fidVol','0')
    overlay_pass_weight = h_weight_func.GetMean()  
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
    purity_func = overlay_sig*globale.scale[globale.overlay]*100/all_pass
    print '\nEfficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}\n'.format(purity_func*efficiency_func/100)
    
    #print amount of each sample
    print 'Amount of %7s'%'ext'+':\t{0:0.2f}%'.format( ext_pass*globale.scale_out[globale.ext_out]*100.0/all_pass )
    print 'Amount of %7s'%'dirt'+':\t{0:0.2f}%'.format( dirt_pass*globale.scale_out[globale.dirt_out]*100.0/all_pass )
    print 'Amount of %7s'%'overlay'+':\t{0:0.2f}%\n'.format( overlay_pass*globale.scale_out[globale.overlay_out]*100.0/all_pass )
    
    #print the amount per signal definition:
    h_weight_func2 = ROOT.TH1F("h_weight_func2",'h_weight_func2',1000,0,100)
    for x in globale.overlay_signals:
        globale.overlay_out.Draw(weight_name+'>>h_weight_func2',cut+' && '+x,'0')
        this_weight2 = h_weight_func2.GetMean()
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay_out.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0*this_weight2/all_pass)
    
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
    h_weight_func = ROOT.TH1F("h_weight_func",'h_weight_func',1000,0,100)
    #get weights of events with cut
    globale.overlay_out.Draw('EventWeight*TunedCentralValue_Genie'+'>>h_weight_func',cut+' && numu_signal','0')
    weight_cut = h_weight_func.GetMean()
    #get weights for all true events
    globale.overlay_out.Draw('EventWeight*TunedCentralValue_Genie'+'>>h_weight_func','numu_true','0')
    weight = h_weight_func.GetMean()
    #calculate efficiency
    efficiency = globale.overlay_out.GetEntries(cut+' && numu_signal')*100.0*weight_cut /(globale.overlay_out.GetEntries('numu_true')*weight+0.001)
    #purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum_out(cut)
    purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100*weight_cut /(getTotNum_out_w(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)

    xstart = 0
    xend = 2
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")

    xstart = -1
    xend = 1
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")

    xstart = -3.14159
    xend = 3.14159
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")

    xstart = -50
    xend = 300
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")

    xstart = -150
    xend = 150
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")

    xstart = -50
    xend = 1050
    xbin = 100
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
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")
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
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]
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
    
    h_data_func.Delete()
    h_ext_func.Delete()
    h_dirt_func.Delete()
    #h_overlay_func = {} # make an array of histograms for the different interactions
    for x in globale.overlay_signals:
        h_overlay_func[x].Delete()
    #sel.Delete()
    return normalization