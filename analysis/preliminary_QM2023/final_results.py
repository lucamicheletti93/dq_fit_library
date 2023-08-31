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



def main():
    LoadStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetHatchesSpacing(0.3)
    #gStyle.SetHatchesLineWidth(2)

    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.05)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)

    inputDir  = "systematics_full_stat_matchedMchMid"

    brJpsiToMuMu = 0.05961
    brPsi2sToMuMu = 0.008

    # Global systematics
    systBrJpsiToMuMu = 0.00033
    systBrPsi2sToMuMu = 0.0006
    systRelBrJpsiToMuMu = systBrJpsiToMuMu / brJpsiToMuMu
    systRelBrPsi2sToMuMu = systBrPsi2sToMuMu / brPsi2sToMuMu
    systRelBr = math.sqrt(systRelBrJpsiToMuMu * systRelBrJpsiToMuMu + systRelBrPsi2sToMuMu * systRelBrPsi2sToMuMu)
    systRelPsi2sWidth = 0.05
    relGlobal = math.sqrt(systRelBr * systRelBr + systRelPsi2sWidth * systRelPsi2sWidth)

    systRelAxeRealisticVsIdealY = [00.012, 0.017, 0.006, 0.000, 0.003, 0.015]
    systRelAxeRealisticVsIdealPt = [0.009, 0.014, 0.013, 0.002, 0.000, 0.000, 0.012, 0.023]
    systRelAxeRealisticVsIdealInt = [0.008]

    # pt-dependence
    dfYieldJpsiPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Jpsi_vs_pt.txt'.format(inputDir), sep=' ')
    ptMin = dfYieldJpsiPt["x_min"].to_numpy()
    ptMax = dfYieldJpsiPt["x_max"].to_numpy()
    ptCentr = (ptMin + ptMax) / 2.
    ptWidth = (ptMax - ptMin) / 2.
    ptArr = np.append(ptMin, ptMax[len(ptMin)-1],)
    yieldJpsiPt = dfYieldJpsiPt["val"].to_numpy()
    statYieldJpsiPt = dfYieldJpsiPt["stat"].to_numpy()
    systYieldJpsiPt = dfYieldJpsiPt["syst"].to_numpy()
    print("Sum J/psi vs pT: ", sum(yieldJpsiPt))

    dfYieldPsi2sPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Psi2s_vs_pt.txt'.format(inputDir), sep=' ')
    yieldPsi2sPt = dfYieldPsi2sPt["val"].to_numpy()
    statYieldPsi2sPt = dfYieldPsi2sPt["stat"].to_numpy()
    systYieldPsi2sPt = dfYieldPsi2sPt["syst"].to_numpy()
    print("Sum J/psi vs pT: ", sum(yieldPsi2sPt))

    dfRatioPsi2sOverJpsiPt = pd.read_csv('/Users/lucamicheletti/cernbox/run3_Psi2s_over_Jpsi/final_results_LHC22all_periods/Comparison_methods/RatioValuesCorentinPt.txt', sep=' ')
    ptMin = dfYieldJpsiPt["x_min"].to_numpy()
    ptMax = dfYieldJpsiPt["x_max"].to_numpy()
    ptCentr = (ptMin + ptMax) / 2.
    ptWidth = (ptMax - ptMin) / 2.
    ptArr = np.append(ptMin, ptMax[len(ptMin)-1],)
    ratioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["val"].to_numpy()
    statRatioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["stat"].to_numpy()
    systRatioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["syst"].to_numpy()

    dfAxeJpsiPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Jpsi_vs_pt.txt', sep=' ')
    axeJpsiPt = dfAxeJpsiPt["val"].to_numpy()
    statAxeJpsiPt = dfAxeJpsiPt["stat"].to_numpy()
    systAxeJpsiPt = dfAxeJpsiPt["syst"].to_numpy()

    dfAxePsi2sPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_vs_pt.txt', sep=' ')
    axePsi2sPt = dfAxePsi2sPt["val"].to_numpy()
    statAxePsi2sPt = dfAxePsi2sPt["stat"].to_numpy()
    systAxePsi2sPt = dfAxePsi2sPt["syst"].to_numpy()

    dfAxeRatioPsi2sOverJpsiPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_over_Jpsi_vs_pt.txt', sep=' ')
    axeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["val"].to_numpy()
    statAxeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["stat"].to_numpy()
    systAxeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["syst"].to_numpy()

    # CGC + NRQCD
    dfCsJpsiTheorCgcNrqcdPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_jpsi_cgc_nrqcd.txt', sep=' ')
    ptCentrTheorCgcNrqcd = dfCsJpsiTheorCgcNrqcdPt["x_centr"].to_numpy()
    minCsJpsiTheorCgcNrqcdPt = dfCsJpsiTheorCgcNrqcdPt["val_min"].to_numpy()
    maxCsJpsiTheorCgcNrqcdPt = dfCsJpsiTheorCgcNrqcdPt["val_max"].to_numpy()

    dfCsPsi2sTheorCgcNrqcdPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_psi2s_cgc_nrqcd.txt', sep=' ')
    minCsPsi2sTheorCgcNrqcdPt = dfCsPsi2sTheorCgcNrqcdPt["val_min"].to_numpy()
    maxCsPsi2sTheorCgcNrqcdPt = dfCsPsi2sTheorCgcNrqcdPt["val_max"].to_numpy()

    ptWidthTheorCgcNrqcd = np.full(len(ptCentrTheorCgcNrqcd), 0.05)
    minCsPsi2sOverJpsiTheorCgcNrqcdPt = minCsPsi2sTheorCgcNrqcdPt / minCsJpsiTheorCgcNrqcdPt
    maxCsPsi2sOverJpsiTheorCgcNrqcdPt = maxCsPsi2sTheorCgcNrqcdPt / maxCsJpsiTheorCgcNrqcdPt
    csPsi2sOverJpsiTheorCgcNrqcdPt = (maxCsPsi2sOverJpsiTheorCgcNrqcdPt + minCsPsi2sOverJpsiTheorCgcNrqcdPt) / 2.
    minCsPsi2sOverJpsiTheorCgcNrqcdPt = csPsi2sOverJpsiTheorCgcNrqcdPt - minCsPsi2sOverJpsiTheorCgcNrqcdPt
    maxCsPsi2sOverJpsiTheorCgcNrqcdPt = maxCsPsi2sOverJpsiTheorCgcNrqcdPt - csPsi2sOverJpsiTheorCgcNrqcdPt

    # NLO NRQCD
    dfCsJpsiTheorNloNrqcdPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_jpsi_nlo_nrqcd.txt', sep=' ')
    ptCentrTheorNloNrqcd = dfCsJpsiTheorNloNrqcdPt["x_centr"].to_numpy()
    minCsJpsiTheorNloNrqcdPt = dfCsJpsiTheorNloNrqcdPt["val_min"].to_numpy()
    maxCsJpsiTheorNloNrqcdPt = dfCsJpsiTheorNloNrqcdPt["val_max"].to_numpy()

    dfCsPsi2sTheorNloNrqcdPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_psi2s_nlo_nrqcd.txt', sep=' ')
    minCsPsi2sTheorNloNrqcdPt = dfCsPsi2sTheorNloNrqcdPt["val_min"].to_numpy()
    maxCsPsi2sTheorNloNrqcdPt = dfCsPsi2sTheorNloNrqcdPt["val_max"].to_numpy()

    ptWidthTheorNloNrqcd = np.full(len(ptCentrTheorNloNrqcd), 0.25)
    minCsPsi2sOverJpsiTheorNloNrqcdPt = minCsPsi2sTheorNloNrqcdPt / minCsJpsiTheorNloNrqcdPt
    maxCsPsi2sOverJpsiTheorNloNrqcdPt = maxCsPsi2sTheorNloNrqcdPt / maxCsJpsiTheorNloNrqcdPt
    csPsi2sOverJpsiTheorNloNrqcdPt = (maxCsPsi2sOverJpsiTheorNloNrqcdPt + minCsPsi2sOverJpsiTheorNloNrqcdPt) / 2.
    minCsPsi2sOverJpsiTheorNloNrqcdPt = csPsi2sOverJpsiTheorNloNrqcdPt - minCsPsi2sOverJpsiTheorNloNrqcdPt
    maxCsPsi2sOverJpsiTheorNloNrqcdPt = maxCsPsi2sOverJpsiTheorNloNrqcdPt - csPsi2sOverJpsiTheorNloNrqcdPt


    # CS NLO
    dfCsPsi2sOverJpsiTheorCsNloPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_nlo_pt.txt', sep=' ')
    ptMinTheorCsCoNlo = dfCsPsi2sOverJpsiTheorCsNloPt["x_min"].to_numpy()
    ptMaxTheorCsCoNlo = dfCsPsi2sOverJpsiTheorCsNloPt["x_max"].to_numpy()
    ptCentrTheorCsCoNlo = (ptMinTheorCsCoNlo + ptMaxTheorCsCoNlo) / 2.
    ptWidthTheorCsCoNlo = (ptMaxTheorCsCoNlo - ptMinTheorCsCoNlo) / 2.
    #csPsi2sOverJpsiTheorCsNloPt = dfCsPsi2sOverJpsiTheorCsNloPt["val"].to_numpy()
    #minCsPsi2sOverJpsiTheorCsNloPt = dfCsPsi2sOverJpsiTheorCsNloPt["val_min"].to_numpy()
    #maxCsPsi2sOverJpsiTheorCsNloPt = dfCsPsi2sOverJpsiTheorCsNloPt["val_max"].to_numpy()

    dfCsPsi2sOverJpsiTheorCsCoNloPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_co_nlo_pt.txt', sep=' ')
    csPsi2sOverJpsiTheorCsCoNloPt = dfCsPsi2sOverJpsiTheorCsCoNloPt["val"].to_numpy()
    minCsPsi2sOverJpsiTheorCsCoNloPt = dfCsPsi2sOverJpsiTheorCsCoNloPt["val_min"].to_numpy()
    maxCsPsi2sOverJpsiTheorCsCoNloPt = dfCsPsi2sOverJpsiTheorCsCoNloPt["val_max"].to_numpy()

    # ICEM + FONLL
    dfCsJpsiTheorIcemFonllPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/psi2s_over_jpsi_ICEM_FONLL_vs_pt.txt', sep=' ')
    ptMinTheorIcemFonll = dfYieldJpsiPt["x_min"].to_numpy()
    ptMaxTheorIcemFonll = dfYieldJpsiPt["x_max"].to_numpy()
    ptCentrTheorIcemFonll = (ptMinTheorIcemFonll + ptMaxTheorIcemFonll) / 2.
    ptWidthTheorIcemFonll = (ptMaxTheorIcemFonll - ptMinTheorIcemFonll) / 2.
    csPsi2sOverJpsiTheorIcemFonllPt = dfCsJpsiTheorIcemFonllPt["val"].to_numpy()
    minCsPsi2sOverJpsiTheorIcemFonllPt = dfCsJpsiTheorIcemFonllPt["err_min"].to_numpy()
    maxCsPsi2sOverJpsiTheorIcemFonllPt = dfCsJpsiTheorIcemFonllPt["err_max"].to_numpy()


    csRatioPsi2sOverJpsiPt = np.zeros((len(ptMin),), dtype=float)
    statCsRatioPsi2sOverJpsiPt = np.zeros((len(ptMin),), dtype=float)
    systCsRatioPsi2sOverJpsiPt = np.zeros((len(ptMin),), dtype=float)

    for iPt in range(0, len(ptMin)):
        relStatYield = statRatioPsi2sOverJpsiPt[iPt] / ratioPsi2sOverJpsiPt[iPt]
        relSystYield = systRatioPsi2sOverJpsiPt[iPt] / ratioPsi2sOverJpsiPt[iPt]
        relStatAxe = statAxeRatioPsi2sOverJpsiPt[iPt] / axeRatioPsi2sOverJpsiPt[iPt]
        relSystAxe = systAxeRatioPsi2sOverJpsiPt[iPt] / axeRatioPsi2sOverJpsiPt[iPt]
        value = ratioPsi2sOverJpsiPt[iPt] / axeRatioPsi2sOverJpsiPt[iPt]
        csRatioPsi2sOverJpsiPt[iPt] = (value * (brJpsiToMuMu / brPsi2sToMuMu))
        statCsRatioPsi2sOverJpsiPt[iPt] = ((value * math.sqrt(relStatYield * relStatYield + relStatAxe * relStatAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))
        systCsRatioPsi2sOverJpsiPt[iPt] = ((value * math.sqrt(relSystYield * relSystYield + relSystAxe * relSystAxe + systRelAxeRealisticVsIdealPt[iPt] * systRelAxeRealisticVsIdealPt[iPt])) * (brJpsiToMuMu / brPsi2sToMuMu))

    for i in range(len(ptMin)):
        print("%3.2f - %3.2f | %4.3f +/- %4.3f +/- %4.3f " % (ptMin[i], ptMax[i], csRatioPsi2sOverJpsiPt[i], statCsRatioPsi2sOverJpsiPt[i], systCsRatioPsi2sOverJpsiPt[i]))
    
    graStatCsRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, csRatioPsi2sOverJpsiPt, ptWidth, statCsRatioPsi2sOverJpsiPt)
    SetGraStat(graStatCsRatioPsi2sOverJpsiPt, 20, ROOT.kBlack)

    graSystCsRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, csRatioPsi2sOverJpsiPt, ptWidth, systCsRatioPsi2sOverJpsiPt)
    SetGraSyst(graSystCsRatioPsi2sOverJpsiPt, 20, ROOT.kBlack)

    graCsPsi2sOverJpsiTheorCgcNrqcdPt = ROOT.TGraphAsymmErrors(len(ptCentrTheorCgcNrqcd), ptCentrTheorCgcNrqcd, csPsi2sOverJpsiTheorCgcNrqcdPt, ptWidthTheorCgcNrqcd, ptWidthTheorCgcNrqcd, minCsPsi2sOverJpsiTheorCgcNrqcdPt, maxCsPsi2sOverJpsiTheorCgcNrqcdPt)
    graCsPsi2sOverJpsiTheorCgcNrqcdPt.SetFillStyle(3353)
    graCsPsi2sOverJpsiTheorCgcNrqcdPt.SetFillColorAlpha(ROOT.kOrange+7, 0.7)
    graCsPsi2sOverJpsiTheorCgcNrqcdPt.SetLineColor(ROOT.kOrange+7)

    graCsPsi2sOverJpsiTheorNloNrqcdPt = ROOT.TGraphAsymmErrors(len(ptCentrTheorNloNrqcd), ptCentrTheorNloNrqcd, csPsi2sOverJpsiTheorNloNrqcdPt, ptWidthTheorNloNrqcd, ptWidthTheorNloNrqcd, minCsPsi2sOverJpsiTheorNloNrqcdPt, maxCsPsi2sOverJpsiTheorNloNrqcdPt)
    graCsPsi2sOverJpsiTheorNloNrqcdPt.SetFillColorAlpha(ROOT.kGray+2, 0.5)
    graCsPsi2sOverJpsiTheorNloNrqcdPt.SetLineColor(ROOT.kGray+2)


    #graCsPsi2sOverJpsiTheorCsNloPt = ROOT.TGraphAsymmErrors(len(ptMinTheor), ptCentrTheor, csPsi2sOverJpsiTheorCsNloPt, ptWidthTheor, ptWidthTheor, minCsPsi2sOverJpsiTheorCsNloPt, maxCsPsi2sOverJpsiTheorCsNloPt)
    #graCsPsi2sOverJpsiTheorCsNloPt.SetFillStyle(3353)
    #graCsPsi2sOverJpsiTheorCsNloPt.SetFillColorAlpha(ROOT.kAzure+2, 0.5)
    #graCsPsi2sOverJpsiTheorCsNloPt.SetLineColor(ROOT.kAzure+2)

    graCsPsi2sOverJpsiTheorCsCoNloPt = ROOT.TGraphAsymmErrors(len(ptMinTheorCsCoNlo), ptCentrTheorCsCoNlo, csPsi2sOverJpsiTheorCsCoNloPt, ptWidthTheorCsCoNlo, ptWidthTheorCsCoNlo, minCsPsi2sOverJpsiTheorCsCoNloPt, maxCsPsi2sOverJpsiTheorCsCoNloPt)
    graCsPsi2sOverJpsiTheorCsCoNloPt.SetFillStyle(3353)
    graCsPsi2sOverJpsiTheorCsCoNloPt.SetFillColorAlpha(ROOT.kAzure+1, 0.5)
    graCsPsi2sOverJpsiTheorCsCoNloPt.SetLineColor(ROOT.kAzure+1)

    graCsPsi2sOverJpsiTheorIcemFonllPt = ROOT.TGraphAsymmErrors(len(ptCentrTheorIcemFonll), ptCentrTheorIcemFonll, csPsi2sOverJpsiTheorIcemFonllPt, ptWidthTheorIcemFonll, ptWidthTheorIcemFonll, minCsPsi2sOverJpsiTheorIcemFonllPt, maxCsPsi2sOverJpsiTheorIcemFonllPt)
    graCsPsi2sOverJpsiTheorIcemFonllPt.SetFillStyle(3335)
    graCsPsi2sOverJpsiTheorIcemFonllPt.SetFillColorAlpha(ROOT.kGreen+1, 0.5)
    graCsPsi2sOverJpsiTheorIcemFonllPt.SetLineColor(ROOT.kGreen+2)


    histStatYieldJpsiPt = TH1F("histStatYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldJpsiPt.SetBinContent(i+1, yieldJpsiPt[i]), histStatYieldJpsiPt.SetBinError(i+1, statYieldJpsiPt[i])
    SetGraStat(histStatYieldJpsiPt, 20, ROOT.kRed+1)
    histStatYieldJpsiPt.Scale(1, "WIDTH")

    histSystYieldJpsiPt = TH1F("histSystYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldJpsiPt.SetBinContent(i+1, yieldJpsiPt[i]), histSystYieldJpsiPt.SetBinError(i+1, systYieldJpsiPt[i])
    SetGraSyst(histSystYieldJpsiPt, 20, ROOT.kRed+1)
    histSystYieldJpsiPt.Scale(1, "WIDTH")

    histStatYieldPsi2sPt = TH1F("histStatYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldPsi2sPt.SetBinContent(i+1, yieldPsi2sPt[i]), histStatYieldPsi2sPt.SetBinError(i+1, statYieldPsi2sPt[i])
    SetGraStat(histStatYieldPsi2sPt, 20, ROOT.kAzure+4)
    histStatYieldPsi2sPt.Scale(1, "WIDTH")

    histSystYieldPsi2sPt = TH1F("histSystYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldPsi2sPt.SetBinContent(i+1, yieldPsi2sPt[i]), histSystYieldPsi2sPt.SetBinError(i+1, systYieldPsi2sPt[i])
    SetGraSyst(histSystYieldPsi2sPt, 20, ROOT.kAzure+4)
    histSystYieldPsi2sPt.Scale(1, "WIDTH")

    # Raw yield
    legendYieldJpsiVsPsi2sPt = TLegend(0.69, 0.59, 0.89, 0.79, " ", "brNDC")
    SetLegend(legendYieldJpsiVsPsi2sPt)
    legendYieldJpsiVsPsi2sPt.AddEntry(histSystYieldJpsiPt, "J/#psi", "FP")
    legendYieldJpsiVsPsi2sPt.AddEntry(histSystYieldPsi2sPt, "#psi(2S)", "FP")

    canvasYieldsPt = TCanvas("canvasYieldsPt", "canvasYieldsPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiPt  = TH2F("histGridYieldJpsiPt", "", 100, 0, 20, 100, 10, 1e7)
    histGridYieldJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (Gev/#it{c})")
    histGridYieldJpsiPt.GetYaxis().SetTitle("d^{2}#it{N} / d#it{p}_{T} d#it{y} (GeV/#it{c}^{-1})")
    histGridYieldJpsiPt.Draw()
    histSystYieldJpsiPt.Draw("E2 SAME")
    histStatYieldJpsiPt.Draw("EP SAME")
    histSystYieldPsi2sPt.Draw("E2 SAME")
    histStatYieldPsi2sPt.Draw("EP SAME")
    legendYieldJpsiVsPsi2sPt.Draw("SAME")
    letexTitle.DrawLatex(0.40, 0.86, "ALICE Preliminary, #sqrt{#it{s}} = 13.6 TeV")
    letexTitle.DrawLatex(0.40, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4")
    canvasYieldsPt.Update()

    # Corrected ratio
    legendCsPsi2sOverJpsiPt1 = TLegend(0.18, 0.5, 0.35, 0.82, " ", "brNDC")
    SetLegend(legendCsPsi2sOverJpsiPt1)
    #legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorCsNloPt,"CS, NLO","F")
    #legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorCsCoNloPt,"CS + CO, NLO","F")
    legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorCgcNrqcdPt,"CGC + NRQCD","F")
    legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorIcemFonllPt,"ICEM (V. Cheung #it{et al.}) + FONLL","F")
    legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorCsCoNloPt,"CS + CO, NLO (M. Butendsch#ddot{o}n #it{et al.})","F")
    legendCsPsi2sOverJpsiPt1.AddEntry(graCsPsi2sOverJpsiTheorNloNrqcdPt,"NLO NRQCD","F")
    legendCsPsi2sOverJpsiPt1.AddEntry(graSystCsRatioPsi2sOverJpsiPt,"Data","FP")

    #legendCsPsi2sOverJpsiPt2 = TLegend(0.43, 0.35, 0.60, 0.52, " ", "brNDC")
    #SetLegend(legendCsPsi2sOverJpsiPt2)
    #legendCsPsi2sOverJpsiPt2.AddEntry(graCsPsi2sOverJpsiTheorCsCoNloPt,"CS + CO, NLO","F")
    #legendCsPsi2sOverJpsiPt2.AddEntry(graCsPsi2sOverJpsiTheorIcemFonllPt,"ICEM + FONLL","F")

    canvasCsPsi2sOverJpsiPt = TCanvas("canvasCsPsi2sOverJpsiPt", "canvasCsPsi2sOverJpsiPt", 800, 600)
    #ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiPt  = TH2F("histGridCsPsi2sOverJpsiPt", "", 100, 0., 20., 100, 0.001, 0.9)
    histGridCsPsi2sOverJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    histGridCsPsi2sOverJpsiPt.GetYaxis().SetTitle("(d^{2}#sigma_{#psi(2S)} / d#it{p}_{T} d#it{y}) / (d^{2}#sigma_{J/#psi} / d#it{p}_{T} d#it{y})")
    histGridCsPsi2sOverJpsiPt.Draw()
    graCsPsi2sOverJpsiTheorIcemFonllPt.Draw("E2 SAME")
    #graCsPsi2sOverJpsiTheorCsNloPt.Draw("L SAME")
    graCsPsi2sOverJpsiTheorCsCoNloPt.Draw("E2 SAME")
    graCsPsi2sOverJpsiTheorCgcNrqcdPt.Draw("E2 SAME")
    graCsPsi2sOverJpsiTheorNloNrqcdPt.Draw("E2 SAME")
    graSystCsRatioPsi2sOverJpsiPt.Draw("E2 SAME")
    graStatCsRatioPsi2sOverJpsiPt.Draw("EP SAME")
    legendCsPsi2sOverJpsiPt1.Draw("EP SAME")
    #legendCsPsi2sOverJpsiPt2.Draw("EP SAME")
    letexTitle.DrawLatex(0.18, 0.86, "ALICE Preliminary, #sqrt{#it{s}} = 13.6 TeV")
    letexTitle.DrawLatex(0.18, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4")
    letexTitle.DrawLatex(0.63, 0.20, "Global unc. = %3.2f %%" % (relGlobal * 100))
    canvasCsPsi2sOverJpsiPt.Update()


    # y-dependence
    dfYieldJpsiY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Jpsi_vs_y.txt'.format(inputDir), sep=' ')
    yMin = dfYieldJpsiY["x_min"].to_numpy()
    yMax = dfYieldJpsiY["x_max"].to_numpy()
    yCentr = (yMin + yMax) / 2.
    yWidth = (yMax - yMin) / 2.
    yArr = np.append(yMin, yMax[len(yMin)-1],)
    yieldJpsiY = dfYieldJpsiY["val"].to_numpy()
    statYieldJpsiY = dfYieldJpsiY["stat"].to_numpy()
    systYieldJpsiY = dfYieldJpsiY["syst"].to_numpy()
    print("Sum J/psi vs pT: ", sum(yieldJpsiY))

    dfYieldPsi2sY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Psi2s_vs_y.txt'.format(inputDir), sep=' ')
    yieldPsi2sY = dfYieldPsi2sY["val"].to_numpy()
    statYieldPsi2sY = dfYieldPsi2sY["stat"].to_numpy()
    systYieldPsi2sY = dfYieldPsi2sY["syst"].to_numpy()
    print("Sum J/psi vs pT: ", sum(yieldPsi2sY))

    dfRatioPsi2sOverJpsiY = pd.read_csv('/Users/lucamicheletti/cernbox/run3_Psi2s_over_Jpsi/final_results_LHC22all_periods/Comparison_methods/RatioValuesCorentinY.txt', sep=' ')
    ratioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["val"].to_numpy()
    statRatioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["stat"].to_numpy()
    systRatioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["syst"].to_numpy()

    dfAxeJpsiY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Jpsi_vs_y.txt', sep=' ')
    axeJpsiY = dfAxeJpsiY["val"].to_numpy()
    statAxeJpsiY = dfAxeJpsiY["stat"].to_numpy()
    systAxeJpsiY = dfAxeJpsiY["syst"].to_numpy()

    dfAxePsi2sY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_vs_y.txt', sep=' ')
    axePsi2sY = dfAxePsi2sY["val"].to_numpy()
    statAxePsi2sY = dfAxePsi2sY["stat"].to_numpy()
    systAxePsi2sY = dfAxePsi2sY["syst"].to_numpy()

    dfAxeRatioPsi2sOverJpsiY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_over_Jpsi_vs_y.txt', sep=' ')
    axeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["val"].to_numpy()
    statAxeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["stat"].to_numpy()
    systAxeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["syst"].to_numpy()

    #dfCsPsi2sOverJpsiTheorCsNloY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_nlo_y.txt', sep=' ')
    #yMinTheor = dfCsPsi2sOverJpsiTheorCsNloY["x_min"].to_numpy()
    #yMaxTheor = dfCsPsi2sOverJpsiTheorCsNloY["x_max"].to_numpy()
    #yCentrTheor = (yMinTheor + yMaxTheor) / 2.
    #yWidthTheor = (yMaxTheor - yMinTheor) / 2.
    #csPsi2sOverJpsiTheorCsNloY = dfCsPsi2sOverJpsiTheorCsNloY["val"].to_numpy()
    #minCsPsi2sOverJpsiTheorCsNloY = dfCsPsi2sOverJpsiTheorCsNloY["val_min"].to_numpy()
    #maxCsPsi2sOverJpsiTheorCsNloY = dfCsPsi2sOverJpsiTheorCsNloY["val_max"].to_numpy()

    #dfCsPsi2sOverJpsiTheorCsCoNloY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_co_nlo_y.txt', sep=' ')
    #csPsi2sOverJpsiTheorCsCoNloY = dfCsPsi2sOverJpsiTheorCsCoNloY["val"].to_numpy()
    #minCsPsi2sOverJpsiTheorCsCoNloY = dfCsPsi2sOverJpsiTheorCsCoNloY["val_min"].to_numpy()
    #maxCsPsi2sOverJpsiTheorCsCoNloY = dfCsPsi2sOverJpsiTheorCsCoNloY["val_max"].to_numpy()

    # ICEM + FONLL
    dfCsJpsiTheorIcemFonllY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/psi2s_over_jpsi_ICEM_FONLL_vs_y.txt', sep=' ')
    yMinTheorIcemFonll = dfYieldJpsiY["x_min"].to_numpy()
    yMaxTheorIcemFonll = dfYieldJpsiY["x_max"].to_numpy()
    yCentrTheorIcemFonll = (yMinTheorIcemFonll + yMaxTheorIcemFonll) / 2.
    yWidthTheorIcemFonll = (yMaxTheorIcemFonll - yMinTheorIcemFonll) / 2.
    csPsi2sOverJpsiTheorIcemFonllY = dfCsJpsiTheorIcemFonllY["val"].to_numpy()
    minCsPsi2sOverJpsiTheorIcemFonllY = dfCsJpsiTheorIcemFonllY["err_min"].to_numpy()
    maxCsPsi2sOverJpsiTheorIcemFonllY = dfCsJpsiTheorIcemFonllY["err_max"].to_numpy()

    csRatioPsi2sOverJpsiY = np.zeros((len(yMin),), dtype=float)
    statCsRatioPsi2sOverJpsiY = np.zeros((len(yMin),), dtype=float)
    systCsRatioPsi2sOverJpsiY = np.zeros((len(yMin),), dtype=float)

    for iY in range(0, len(yMin)):
        relStatYield = statRatioPsi2sOverJpsiY[iY] / ratioPsi2sOverJpsiY[iY]
        relSystYield = systRatioPsi2sOverJpsiY[iY] / ratioPsi2sOverJpsiY[iY]
        relStatAxe = statAxeRatioPsi2sOverJpsiY[iY] / axeRatioPsi2sOverJpsiY[iY]
        relSystAxe = systAxeRatioPsi2sOverJpsiY[iY] / axeRatioPsi2sOverJpsiY[iY]
        value = ratioPsi2sOverJpsiY[iY] / axeRatioPsi2sOverJpsiY[iY]
        csRatioPsi2sOverJpsiY[iY] = (value * (brJpsiToMuMu / brPsi2sToMuMu))
        statCsRatioPsi2sOverJpsiY[iY] = ((value * math.sqrt(relStatYield * relStatYield + relStatAxe * relStatAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))
        systCsRatioPsi2sOverJpsiY[iY] = ((value * math.sqrt(relSystYield * relSystYield + relSystAxe * relSystAxe + systRelAxeRealisticVsIdealY[iY] * systRelAxeRealisticVsIdealY[iY])) * (brJpsiToMuMu / brPsi2sToMuMu))

    print("Mean of the ratio")
    for i in range(len(yMin)):
        print("%3.2f - %3.2f | %4.3f +/- %4.3f +/- %4.3f " % (yMin[i], yMax[i], csRatioPsi2sOverJpsiY[i], statCsRatioPsi2sOverJpsiY[i], systCsRatioPsi2sOverJpsiY[i]))

    graStatCsRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, csRatioPsi2sOverJpsiY, yWidth, statCsRatioPsi2sOverJpsiY)
    SetGraStat(graStatCsRatioPsi2sOverJpsiY, 20, ROOT.kBlack)

    graSystCsRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, csRatioPsi2sOverJpsiY, yWidth, systCsRatioPsi2sOverJpsiY)
    SetGraSyst(graSystCsRatioPsi2sOverJpsiY, 20, ROOT.kBlack)

    #graCsPsi2sOverJpsiTheorCsNloY = ROOT.TGraphAsymmErrors(len(yMinTheor), yCentrTheor, csPsi2sOverJpsiTheorCsNloY, yWidthTheor, yWidthTheor, minCsPsi2sOverJpsiTheorCsNloY, maxCsPsi2sOverJpsiTheorCsNloY)
    #graCsPsi2sOverJpsiTheorCsNloY.SetFillStyle(3353)
    #graCsPsi2sOverJpsiTheorCsNloY.SetFillColorAlpha(ROOT.kAzure+2, 0.5)
    #graCsPsi2sOverJpsiTheorCsNloY.SetLineColor(ROOT.kAzure+2)

    #graCsPsi2sOverJpsiTheorCsCoNloY = ROOT.TGraphAsymmErrors(len(yMinTheor), yCentrTheor, csPsi2sOverJpsiTheorCsCoNloY, yWidthTheor, yWidthTheor, minCsPsi2sOverJpsiTheorCsCoNloY, maxCsPsi2sOverJpsiTheorCsCoNloY)
    #graCsPsi2sOverJpsiTheorCsCoNloY.SetFillStyle(3353)
    #graCsPsi2sOverJpsiTheorCsCoNloY.SetFillColorAlpha(ROOT.kRed+1, 0.5)
    #graCsPsi2sOverJpsiTheorCsCoNloY.SetLineColor(ROOT.kRed+1)

    graCsPsi2sOverJpsiTheorIcemFonllY = ROOT.TGraphAsymmErrors(len(ptCentrTheorIcemFonll), ptCentrTheorIcemFonll, csPsi2sOverJpsiTheorIcemFonllY, ptWidthTheorIcemFonll, ptWidthTheorIcemFonll, minCsPsi2sOverJpsiTheorIcemFonllY, maxCsPsi2sOverJpsiTheorIcemFonllY)
    graCsPsi2sOverJpsiTheorIcemFonllY.SetFillColorAlpha(ROOT.kRed+1, 0.5)
    graCsPsi2sOverJpsiTheorIcemFonllY.SetLineColor(ROOT.kRed+1)


    histStatYieldJpsiY = TH1F("histStatYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldJpsiY.SetBinContent(i+1, yieldJpsiY[i]), histStatYieldJpsiY.SetBinError(i+1, statYieldJpsiY[i])
    SetGraStat(histStatYieldJpsiY, 20, ROOT.kRed+1)
    histStatYieldJpsiY.Scale(1, "WIDTH")

    histSystYieldJpsiY = TH1F("histSystYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldJpsiY.SetBinContent(i+1, yieldJpsiY[i]), histSystYieldJpsiY.SetBinError(i+1, systYieldJpsiY[i])
    SetGraSyst(histSystYieldJpsiY, 20, ROOT.kRed+1)
    histSystYieldJpsiY.Scale(1, "WIDTH")

    histStatYieldPsi2sY = TH1F("histStatYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldPsi2sY.SetBinContent(i+1, yieldPsi2sY[i]), histStatYieldPsi2sY.SetBinError(i+1, statYieldPsi2sY[i])
    SetGraStat(histStatYieldPsi2sY, 20, ROOT.kAzure+4)
    histStatYieldPsi2sY.Scale(1, "WIDTH")

    histSystYieldPsi2sY = TH1F("histSystYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldPsi2sY.SetBinContent(i+1, yieldPsi2sY[i]), histSystYieldPsi2sY.SetBinError(i+1, systYieldPsi2sY[i])
    SetGraSyst(histSystYieldPsi2sY, 20, ROOT.kAzure+4)
    histSystYieldPsi2sY.Scale(1, "WIDTH")


    # Raw yield
    legendYieldJpsiVsPsi2sY = TLegend(0.18, 0.2, 0.35, 0.42, " ", "brNDC")
    SetLegend(legendYieldJpsiVsPsi2sY)
    legendYieldJpsiVsPsi2sY.AddEntry(histSystYieldJpsiY, "J/#psi", "FP")
    legendYieldJpsiVsPsi2sY.AddEntry(histSystYieldPsi2sY, "#psi(2S)", "FP")

    canvasYieldsY = TCanvas("canvasYieldsY", "canvasYieldsY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiY  = TH2F("histGridYieldJpsiY", "", 100, 2.5, 4, 100, 10, 1e7)
    histGridYieldJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridYieldJpsiY.GetYaxis().SetTitle("d#it{N} / d#it{y}")
    histGridYieldJpsiY.Draw()
    histSystYieldJpsiY.Draw("E2 SAME")
    histStatYieldJpsiY.Draw("EP SAME")
    histSystYieldPsi2sY.Draw("E2 SAME")
    histStatYieldPsi2sY.Draw("EP SAME")
    legendYieldJpsiVsPsi2sY.Draw("SAME")
    letexTitle.DrawLatex(0.18, 0.42, "ALICE Preliminary, #sqrt{#it{s}} = 13.6 TeV")
    letexTitle.DrawLatex(0.18, 0.36, "J/#psi, #psi(2S) #rightarrow #mu^{+}#mu^{-}")
    canvasYieldsY.Update()

    # Corrected ratio
    legendCsPsi2sOverJpsiY = TLegend(0.18, 0.2, 0.35, 0.42, " ", "brNDC")
    SetLegend(legendCsPsi2sOverJpsiY)
    #legendCsPsi2sOverJpsiY.AddEntry(graCsPsi2sOverJpsiTheorCsNloY,"CS, NLO, 4 < #it{p}_{T} < 20 GeV/#it{c}", "F")
    #legendCsPsi2sOverJpsiY.AddEntry(graCsPsi2sOverJpsiTheorCsCoNloY,"CS + CO, NLO, 4 < #it{p}_{T} < 20 GeV/#it{c}", "F")
    legendCsPsi2sOverJpsiY.AddEntry(graSystCsRatioPsi2sOverJpsiY,"Data, #it{p}_{T} < 20 GeV/#it{c}","FP")

    canvasCsPsi2sOverJpsiY = TCanvas("canvasCsPsi2sOverJpsiY", "canvasCsPsi2sOverJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiY  = TH2F("histGridCsPsi2sOverJpsiY", "", 100, 2.5, 4, 100, 0.001, 1)
    histGridCsPsi2sOverJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridCsPsi2sOverJpsiY.GetYaxis().SetTitle("(d#sigma_{#psi(2S)} / d#it{y}) / (d#sigma_{J/#psi} / d#it{y})")
    histGridCsPsi2sOverJpsiY.Draw()
    #graCsPsi2sOverJpsiTheorCsNloY.Draw("L SAME")
    #graCsPsi2sOverJpsiTheorCsCoNloY.Draw("E2 SAME")
    graCsPsi2sOverJpsiTheorIcemFonllY.Draw("E2 SAME")
    graSystCsRatioPsi2sOverJpsiY.Draw("E2 SAME")
    graStatCsRatioPsi2sOverJpsiY.Draw("EP SAME")
    legendCsPsi2sOverJpsiY.Draw("EP SAME")
    letexTitle.DrawLatex(0.18, 0.46, "ALICE Preliminary, #sqrt{#it{s}} = 13.6 TeV")
    letexTitle.DrawLatex(0.18, 0.40, "J/#psi, #psi(2S) #rightarrow #mu^{+}#mu^{-}")
    letexTitle.DrawLatex(0.63, 0.20, "Global unc. = %3.2f %%" % (relGlobal * 100))
    canvasCsPsi2sOverJpsiY.Update()

    input()
    canvasYieldsPt.SaveAs("yieldsPt.pdf")
    canvasCsPsi2sOverJpsiPt.SaveAs("csPsi2sOverJpsiPt.pdf")
    canvasYieldsY.SaveAs("yieldsY.pdf")
    canvasCsPsi2sOverJpsiY.SaveAs("csPsi2sOverJpsiY.pdf")




if __name__ == '__main__':
    main()