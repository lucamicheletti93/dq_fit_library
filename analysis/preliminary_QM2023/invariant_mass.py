import yaml
import json
import sys
import argparse
from array import array
import os
import math 
from os import path
import numpy as np
import pandas as pd
import uncertainties
from uncertainties import ufloat, unumpy
import ROOT
from ROOT import TCanvas, TH1F, TH2F, TGraphErrors, TLegend
sys.path.append('../../utils')
from utils_library import LoadStyle, PropagateErrorsOnRatio, ToCArray
from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

ROOT.gROOT.ProcessLineSync(".x ../../fit_library/VWGPdf.cxx+")
ROOT.gROOT.ProcessLineSync(".x ../../fit_library/CB2Pdf.cxx+")

def main():
    LoadStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetHatchesSpacing(0.3)
    #gStyle.SetHatchesLineWidth(2)

    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.05)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)

    ptMin = "3"
    ptMax = "4"
    path = "/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/output/analysis/Run3_full_stat_apass4_data_run3_tails_matchedMchMid"
    histName = "hist_mass_pt_all_histo_PairsMuonSEPM_matchedMchMid_pt_" + ptMin + "_" + ptMax
    fInName = "CB2_CB2_VWG__2.1_4.9"
    
    fIn = ROOT.TFile("{}/{}__{}.root".format(path, histName, fInName), "READ")
    print("{}/{}__{}.root".format(path, histName, fInName))
    canvasIn = fIn.Get("fit_plot_{}".format(fInName))
    listOfPrimitives = canvasIn.GetListOfPrimitives()
    
    # Frame
    frame = listOfPrimitives.At(0)
    frame.GetXaxis().SetRangeUser(2.2, 4.5)
    frame.GetYaxis().SetRangeUser(1e2, 5e5)
    frame.SetTitle(" ")

    # Histograms
    histData = listOfPrimitives.At(1)

    # PDFs
    pdfSum = listOfPrimitives.At(4)
    pdfJpsi = listOfPrimitives.At(5)
    pdfPsi2s = listOfPrimitives.At(6)
    pdfBkg = listOfPrimitives.At(7)

    canvasOut = TCanvas("canvasOut", "canvasOut", 800, 600)
    ROOT.gPad.SetLogy(1)
    frame.Draw()
    histData.Draw("EP SAME")
    pdfSum.Draw("SAME")
    pdfBkg.Draw("SAME")
    pdfJpsi.Draw("SAME")
    pdfPsi2s.Draw("SAME")

    legend = TLegend(0.68, 0.50, 0.85, 0.85, " ", "brNDC")
    SetLegend(legend)
    legend.SetTextSize(0.045)
    legend.AddEntry(histData,"Data", "P")
    legend.AddEntry(pdfSum,"Fit", "L")
    legend.AddEntry(pdfJpsi,"J/#psi", "L")
    legend.AddEntry(pdfPsi2s,"#psi(2S)", "L")
    legend.AddEntry(pdfBkg,"Background", "L")
    legend.Draw()

    letexTitle.DrawLatex(0.20, 0.88, "ALICE Preliminary, #sqrt{#it{s}} = 13.6 TeV")
    letexTitle.DrawLatex(0.20, 0.81, "J/#psi, #psi(2S) #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4")
    letexTitle.DrawLatex(0.20, 0.74, ptMin + " < #it{p}_{T} < " + ptMax + " GeV/#it{c}")

    canvasOut.Update()

    input()
    canvasOut.SaveAs("invariantMass_pt_" + ptMin + "_" + ptMax + ".pdf")

if __name__ == '__main__':
    main()