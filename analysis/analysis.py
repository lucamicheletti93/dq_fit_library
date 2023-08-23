import yaml
import json
import sys
import argparse
import math
from array import array
import os
from os import path
import numpy as np
import pandas as pd
import uncertainties
from uncertainties import ufloat, unumpy
import ROOT
from ROOT import TCanvas, TH1F, TH2F, TGraphErrors, TLegend
sys.path.append('../utils')
from utils_library import LoadStyle, PropagateErrorsOnRatio, ToCArray
from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

inputDir  = "systematics_full_stat_matchedMchMid"
outputDir = "figures_full_stat_matchedMchMid"

brJpsiToMuMu = 0.05961
brPsi2sToMuMu = 0.008

################
################
def intResults():
    LoadStyle()
    ROOT.gStyle.SetOptStat(0)

    ###############
    # Load datasets
    ###############
    # integrated
    dfYieldJpsiInt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Jpsi_integrated.txt'.format(inputDir), sep=' ')
    intMin = np.array(0.00)
    intMax = np.array(20.00)
    intCentr = (intMin + intMax) / 2.
    intWidth = (intMax - intMin) / 2.
    yieldJpsiInt = dfYieldJpsiInt["val"].to_numpy()
    statYieldJpsiInt = dfYieldJpsiInt["stat"].to_numpy()
    systYieldJpsiInt = dfYieldJpsiInt["syst"].to_numpy()
    print("Integrated J/psi: ", sum(yieldJpsiInt))

    dfYieldPsi2sInt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Psi2s_integrated.txt'.format(inputDir), sep=' ')
    yieldPsi2sInt = dfYieldPsi2sInt["val"].to_numpy()
    statYieldPsi2sInt = dfYieldPsi2sInt["stat"].to_numpy()
    systYieldPsi2sInt = dfYieldPsi2sInt["syst"].to_numpy()
    print("Integrated Psi(2S): ", sum(yieldPsi2sInt))

    dfRatioPsi2sOverJpsiInt = pd.read_csv('/Users/lucamicheletti/cernbox/run3_Psi2s_over_Jpsi/final_results_LHC22all_periods/Comparison_methods/RatioValuesCorentinInt.txt', sep=' ')
    ratioPsi2sOverJpsiInt = dfRatioPsi2sOverJpsiInt["val"].to_numpy()
    statRatioPsi2sOverJpsiInt = dfRatioPsi2sOverJpsiInt["stat"].to_numpy()
    systRatioPsi2sOverJpsiInt = dfRatioPsi2sOverJpsiInt["syst"].to_numpy()

    dfAxeJpsiInt = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Jpsi_int.txt', sep=' ')
    axeJpsiInt = dfAxeJpsiInt["val"].to_numpy()
    statAxeJpsiInt = dfAxeJpsiInt["stat"].to_numpy()
    systAxeJpsiInt = dfAxeJpsiInt["syst"].to_numpy()

    dfAxePsi2sInt = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Psi2s_int.txt', sep=' ')
    axePsi2sInt = dfAxePsi2sInt["val"].to_numpy()
    statAxePsi2sInt = dfAxePsi2sInt["stat"].to_numpy()
    systAxePsi2sInt = dfAxePsi2sInt["syst"].to_numpy()

    dfAxeRatioPsi2sOverJpsiInt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_over_Jpsi_int.txt', sep=' ')
    axeRatioPsi2sOverJpsiInt = dfAxeRatioPsi2sOverJpsiInt["val"].to_numpy()
    statAxeRatioPsi2sOverJpsiInt = dfAxeRatioPsi2sOverJpsiInt["stat"].to_numpy()
    systAxeRatioPsi2sOverJpsiInt = dfAxeRatioPsi2sOverJpsiInt["syst"].to_numpy()

    ################################
    # Plot the Ratio Psi(2S) / J/psi
    ################################
    measStatRatioInt, statMeasStatRatioInt = PropagateErrorsOnRatio(yieldJpsiInt, statYieldJpsiInt, yieldPsi2sInt, statYieldPsi2sInt)
    measSystRatioInt, systMeasSystRatioInt = PropagateErrorsOnRatio(yieldJpsiInt, systYieldJpsiInt, yieldPsi2sInt, systYieldPsi2sInt)
    
    print(ToCArray(measStatRatioInt, ctype='double', name='ratioPsi2sOverJpsiInt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(statMeasStatRatioInt, ctype='double', name='statRatioPsi2sOverJpsiInt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(systMeasSystRatioInt, ctype='double', name='systRatioPsi2sOverJpsiInt', formatter=lambda x: '{:0.5f}'.format(x)))

    # Compute and plot the corrected yield
    corrYieldJpsiInt, statCorrYieldJpsiInt = PropagateErrorsOnRatio(axeJpsiInt, statAxeJpsiInt, yieldJpsiInt, statYieldJpsiInt)
    corrYieldJpsiInt, systCorrYieldJpsiInt = PropagateErrorsOnRatio(axeJpsiInt, systAxeJpsiInt, yieldJpsiInt, systYieldJpsiInt)
    corrYieldJpsiInt = corrYieldJpsiInt / brJpsiToMuMu
    statCorrYieldJpsiInt = statCorrYieldJpsiInt / brJpsiToMuMu
    systCorrYieldJpsiInt = systCorrYieldJpsiInt / brJpsiToMuMu

    corrYieldPsi2sInt, statCorrYieldPsi2sInt = PropagateErrorsOnRatio(axePsi2sInt, statAxePsi2sInt, yieldPsi2sInt, statYieldPsi2sInt)
    corrYieldPsi2sInt, systCorrYieldPsi2sInt = PropagateErrorsOnRatio(axePsi2sInt, systAxePsi2sInt, yieldPsi2sInt, systYieldPsi2sInt)
    corrYieldPsi2sInt = corrYieldPsi2sInt / brPsi2sToMuMu
    statCorrYieldPsi2sInt = statCorrYieldPsi2sInt / brPsi2sToMuMu
    systCorrYieldPsi2sInt = systCorrYieldPsi2sInt / brPsi2sToMuMu

    # Compute and plot the cross section ratio
    csRatioInt, statCsRatioInt = PropagateErrorsOnRatio(corrYieldJpsiInt, statCorrYieldJpsiInt, corrYieldPsi2sInt, statCorrYieldPsi2sInt)
    csRatioInt, systCsRatioInt = PropagateErrorsOnRatio(corrYieldJpsiInt, systCorrYieldJpsiInt, corrYieldPsi2sInt, systCorrYieldPsi2sInt)

    csRatioPsi2sOverJpsiInt = np.zeros((1,), dtype=float)
    statCsRatioPsi2sOverJpsiInt = np.zeros((1,), dtype=float)
    systCsRatioPsi2sOverJpsiInt = np.zeros((1,), dtype=float)

    for iInt in range(0, 1):
        relStatYield = statRatioPsi2sOverJpsiInt[iInt] / ratioPsi2sOverJpsiInt[iInt]
        relSystYield = systRatioPsi2sOverJpsiInt[iInt] / ratioPsi2sOverJpsiInt[iInt]
        relStatAxe = statAxeRatioPsi2sOverJpsiInt[iInt] / axeRatioPsi2sOverJpsiInt[iInt]
        relSystAxe = systAxeRatioPsi2sOverJpsiInt[iInt] / axeRatioPsi2sOverJpsiInt[iInt]
        value = ratioPsi2sOverJpsiInt[iInt] / axeRatioPsi2sOverJpsiInt[iInt]
        csRatioPsi2sOverJpsiInt[iInt] = (value * (brJpsiToMuMu / brPsi2sToMuMu))
        statCsRatioPsi2sOverJpsiInt[iInt] = ((value * math.sqrt(relStatYield * relStatYield + relStatAxe * relStatAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))
        systCsRatioPsi2sOverJpsiInt[iInt] = ((value * math.sqrt(relSystYield * relSystYield + relSystAxe * relSystAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))

    # Print all results for analysis    

    print("Ratio of the mean")
    for i in range(0, 1):
        relStatRatio = (statCsRatioInt[i] / csRatioInt[i]) * 100
        relSystRatio = (systCsRatioInt[i] / csRatioInt[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (intMin, intMax, csRatioInt[i], statCsRatioInt[i], relStatRatio, systCsRatioInt[i], relSystRatio))

    print("Mean of the ratio")
    for i in range(0, 1):
        relStatRatio = (statCsRatioPsi2sOverJpsiInt[i] / csRatioPsi2sOverJpsiInt[i]) * 100
        relSystRatio = (systCsRatioPsi2sOverJpsiInt[i] / csRatioPsi2sOverJpsiInt[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (intMin, intMax, csRatioPsi2sOverJpsiInt[i], statCsRatioPsi2sOverJpsiInt[i], relStatRatio, systCsRatioPsi2sOverJpsiInt[i], relSystRatio))

    intCentrSqrtsRun2 = np.array([5.0, 7.0, 8.0, 13.0])
    intWidthSqrtsRun2 = np.array([0.20, 0.20, 0.20, 0.20])
    csRatioIntSqrtsRun2 = np.array([0.128, 0.170, 0.138, 0.146])
    statCsRatioIntSqrtsRun2 = np.array([0.029, 0.011, 0.009, 0.0053])
    systCsRatioIntSqrtsRun2 = np.array([0.012, 0.013, 0.0153, 0.0087])

    graStatCsRatioIntSqrtsRun2 = TGraphErrors(4, intCentrSqrtsRun2, csRatioIntSqrtsRun2, 0, statCsRatioIntSqrtsRun2)
    SetGraStat(graStatCsRatioIntSqrtsRun2, 20, ROOT.kBlack)

    graSystCsRatioIntSqrtsRun2 = TGraphErrors(4, intCentrSqrtsRun2, csRatioIntSqrtsRun2, intWidthSqrtsRun2, systCsRatioIntSqrtsRun2)
    SetGraSyst(graSystCsRatioIntSqrtsRun2, 20, ROOT.kBlack)

    intCentrSqrts = np.array([13.6])
    intWidthSqrts = np.array([0.20])

    graStatCsRatioInt = TGraphErrors(1, intCentrSqrts, csRatioInt, 0, statCsRatioInt)
    SetGraStat(graStatCsRatioInt, 20, ROOT.kRed+1)

    graSystCsRatioInt = TGraphErrors(1, intCentrSqrts, csRatioInt, intWidthSqrts, systCsRatioInt)
    SetGraSyst(graSystCsRatioInt, 20, ROOT.kRed+1)

    graStatCsRatioPsi2sOverJpsiInt = TGraphErrors(1, intCentrSqrts, csRatioPsi2sOverJpsiInt, 0, statCsRatioPsi2sOverJpsiInt)
    SetGraStat(graStatCsRatioPsi2sOverJpsiInt, 20, ROOT.kAzure+2)

    graSystCsRatioPsi2sOverJpsiInt = TGraphErrors(1, intCentrSqrts, csRatioPsi2sOverJpsiInt, intWidthSqrts, systCsRatioPsi2sOverJpsiInt)
    SetGraSyst(graSystCsRatioPsi2sOverJpsiInt, 20, ROOT.kAzure+2)

    canvasCsRatioInt = TCanvas("canvasCsRatioInt", "canvasCsRatioInt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsRatioInt  = TH2F("histGridCsRatioInt", "", 100, 0, 20, 100, 0.01, 0.5)
    histGridCsRatioInt.GetXaxis().SetTitle("#sqrt{#it{s}} (TeV)")
    histGridCsRatioInt.GetYaxis().SetTitle("d#sigma_{#psi(2S)}/d#it{y} / d#sigma_{J/#psi}/d#it{y}")
    histGridCsRatioInt.Draw()
    graSystCsRatioIntSqrtsRun2.Draw("E2 SAME")
    graStatCsRatioIntSqrtsRun2.Draw("EP SAME")
    graSystCsRatioInt.Draw("E2 SAME")
    graStatCsRatioInt.Draw("EP SAME")
    graSystCsRatioPsi2sOverJpsiInt.Draw("E2 SAME")
    graStatCsRatioPsi2sOverJpsiInt.Draw("EP SAME")
    canvasCsRatioInt.SaveAs("{}/cross_section_psi2s_over_jpsi_int.pdf".format(outputDir))

################
################
def ptResults():
    LoadStyle()
    ROOT.gStyle.SetOptStat(0)

    ###############
    # Load datasets
    ###############
    dfYieldJpsiPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Jpsi_vs_pt.txt'.format(inputDir), sep=' ')
    ptMin = dfYieldJpsiPt["x_min"].to_numpy()
    ptMax = dfYieldJpsiPt["x_max"].to_numpy()
    ptArr = np.append(ptMin, ptMax[len(ptMin)-1],)
    ptCentr = (ptMin + ptMax) / 2.
    ptWidth = (ptMax - ptMin) / 2.
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
    ratioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["val"].to_numpy()
    statRatioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["stat"].to_numpy()
    systRatioPsi2sOverJpsiPt = dfRatioPsi2sOverJpsiPt["syst"].to_numpy()

    # Print all results for analysis
    for i in range(len(ptMin)):
        relStatJpsi = (statYieldJpsiPt[i] / yieldJpsiPt[i]) * 100
        relSystJpsi = (systYieldJpsiPt[i] / yieldJpsiPt[i]) * 100
        relStatPsi2s = (statYieldPsi2sPt[i] / yieldPsi2sPt[i]) * 100
        relSystPsi2s = (systYieldPsi2sPt[i] / yieldPsi2sPt[i]) * 100
        print("& %3.2f - %3.2f & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) ??" % (ptMin[i], ptMax[i], 
            yieldJpsiPt[i], statYieldJpsiPt[i], relStatJpsi, systYieldJpsiPt[i], relSystJpsi, 
            yieldPsi2sPt[i], statYieldPsi2sPt[i], relStatPsi2s, systYieldPsi2sPt[i], relSystPsi2s))

    # Run2 comparison
    dfYieldJpsiPt = pd.read_csv('run2_results/sig_Jpsi_vs_pt_run2.txt', sep=' ')
    ptMinRun2 = dfYieldJpsiPt["x_min"].to_numpy()
    ptMaxRun2 = dfYieldJpsiPt["x_max"].to_numpy()
    ptArrRun2 = np.append(ptMinRun2, ptMaxRun2[len(ptMinRun2)-1],)
    ptCentrRun2 = (ptMinRun2 + ptMaxRun2) / 2.
    ptWidthRun2 = (ptMaxRun2 - ptMinRun2) / 2.
    yieldJpsiPtRun2 = dfYieldJpsiPt["val"].to_numpy()
    statYieldJpsiPtRun2 = dfYieldJpsiPt["stat"].to_numpy()
    systYieldJpsiPtRun2 = dfYieldJpsiPt["syst"].to_numpy()

    dfYieldPsi2sPt = pd.read_csv('run2_results/sig_Psi2s_vs_pt_run2.txt', sep=' ')
    yieldPsi2sPtRun2 = dfYieldPsi2sPt["val"].to_numpy()
    statYieldPsi2sPtRun2 = dfYieldPsi2sPt["stat"].to_numpy()
    systYieldPsi2sPtRun2 = dfYieldPsi2sPt["syst"].to_numpy()

    dfCsRatioPt = pd.read_csv('run2_results/cs_ratio_Psi2s_Jpsi_vs_pt_run2.txt', sep=' ')
    csRatioPtRun2 = dfCsRatioPt["val"].to_numpy()
    statCsRatioPtRun2 = dfCsRatioPt["stat"].to_numpy()
    systCsRatioPtRun2 = dfCsRatioPt["syst"].to_numpy()

    # Axe
    dfAxeJpsiPt = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Jpsi_vs_pt.txt', sep=' ')
    axeJpsiPt = dfAxeJpsiPt["val"].to_numpy()
    statAxeJpsiPt = dfAxeJpsiPt["stat"].to_numpy()
    systAxeJpsiPt = dfAxeJpsiPt["syst"].to_numpy()

    dfAxePsi2sPt = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Psi2s_vs_pt.txt', sep=' ')
    axePsi2sPt = dfAxePsi2sPt["val"].to_numpy()
    statAxePsi2sPt = dfAxePsi2sPt["stat"].to_numpy()
    systAxePsi2sPt = dfAxePsi2sPt["syst"].to_numpy()

    dfAxeRatioPsi2sOverJpsiPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_over_Jpsi_vs_pt.txt', sep=' ')
    axeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["val"].to_numpy()
    statAxeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["stat"].to_numpy()
    systAxeRatioPsi2sOverJpsiPt = dfAxeRatioPsi2sOverJpsiPt["syst"].to_numpy()

    for i in range(len(ptMin)):
        relStatJpsi = (statAxeJpsiPt[i] / axeJpsiPt[i]) * 100
        relSystJpsi = (systAxeJpsiPt[i] / axeJpsiPt[i]) * 100
        relStatPsi2s = (statAxePsi2sPt[i] / axePsi2sPt[i]) * 100
        relSystPsi2s = (systAxePsi2sPt[i] / axePsi2sPt[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) & %4.3f $\pm$ %4.3f (%3.2f?%%) ??" % (ptMin[i], ptMax[i], 
            axeJpsiPt[i], statAxeJpsiPt[i], relStatJpsi, 
            axePsi2sPt[i], statAxePsi2sPt[i], relStatPsi2s))
        
    ############################
    # Create and fill histograms
    ############################
    histStatYieldJpsiPt = TH1F("histStatYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldJpsiPt.SetBinContent(i+1, yieldJpsiPt[i]), histStatYieldJpsiPt.SetBinError(i+1, statYieldJpsiPt[i])
    SetGraStat(histStatYieldJpsiPt, 20, ROOT.kRed+1)
    histStatYieldJpsiPt.Scale(1, "WIDTH")

    histSystYieldJpsiPt = TH1F("histSystYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldJpsiPt.SetBinContent(i+1, yieldJpsiPt[i]), histSystYieldJpsiPt.SetBinError(i+1, systYieldJpsiPt[i])
    SetGraSyst(histSystYieldJpsiPt, 20, ROOT.kRed+1)
    histSystYieldJpsiPt.Scale(1, "WIDTH")

    histStatYieldJpsiNormPt = TH1F("histStatYieldJpsiNormPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldJpsiNormPt.SetBinContent(i+1, yieldJpsiPt[i]), histStatYieldJpsiNormPt.SetBinError(i+1, statYieldJpsiPt[i])
    SetGraStat(histStatYieldJpsiNormPt, 20, ROOT.kRed+1)
    histStatYieldJpsiNormPt.Scale(1. / histStatYieldJpsiNormPt.Integral(), "WIDTH")

    histSystYieldJpsiNormPt = TH1F("histSystYieldJpsiNormPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldJpsiNormPt.SetBinContent(i+1, yieldJpsiPt[i]), histSystYieldJpsiNormPt.SetBinError(i+1, systYieldJpsiPt[i])
    SetGraSyst(histSystYieldJpsiNormPt, 20, ROOT.kRed+1)
    histSystYieldJpsiNormPt.Scale(1. / histSystYieldJpsiNormPt.Integral(), "WIDTH")

    histStatYieldPsi2sPt = TH1F("histStatYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldPsi2sPt.SetBinContent(i+1, yieldPsi2sPt[i]), histStatYieldPsi2sPt.SetBinError(i+1, statYieldPsi2sPt[i])
    SetGraStat(histStatYieldPsi2sPt, 20, ROOT.kAzure+4)
    histStatYieldPsi2sPt.Scale(1, "WIDTH")

    histSystYieldPsi2sPt = TH1F("histSystYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldPsi2sPt.SetBinContent(i+1, yieldPsi2sPt[i]), histSystYieldPsi2sPt.SetBinError(i+1, systYieldPsi2sPt[i])
    SetGraSyst(histSystYieldPsi2sPt, 20, ROOT.kAzure+4)
    histSystYieldPsi2sPt.Scale(1, "WIDTH")

    histStatYieldPsi2sNormPt = TH1F("histStatYieldPsi2sNormPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatYieldPsi2sNormPt.SetBinContent(i+1, yieldPsi2sPt[i]), histStatYieldPsi2sNormPt.SetBinError(i+1, statYieldPsi2sPt[i])
    SetGraStat(histStatYieldPsi2sNormPt, 20, ROOT.kAzure+4)
    histStatYieldPsi2sNormPt.Scale(1. / histStatYieldPsi2sNormPt.Integral(), "WIDTH")

    histSystYieldPsi2sNormPt = TH1F("histSystYieldPsi2sNormPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystYieldPsi2sNormPt.SetBinContent(i+1, yieldPsi2sPt[i]), histSystYieldPsi2sNormPt.SetBinError(i+1, systYieldPsi2sPt[i])
    SetGraSyst(histSystYieldPsi2sNormPt, 20, ROOT.kAzure+4)
    histSystYieldPsi2sNormPt.Scale(1. / histSystYieldPsi2sNormPt.Integral(), "WIDTH")

    # Run2
    histStatYieldJpsiNormPtRun2 = TH1F("histStatYieldJpsiNormPtRun2", "", len(ptArrRun2)-1, ptArrRun2)
    for i in range(0, len(ptArrRun2)-1) : histStatYieldJpsiNormPtRun2.SetBinContent(i+1, yieldJpsiPtRun2[i]), histStatYieldJpsiNormPtRun2.SetBinError(i+1, statYieldJpsiPtRun2[i])
    SetGraStat(histStatYieldJpsiNormPtRun2, 20, ROOT.kBlack)
    histStatYieldJpsiNormPtRun2.Scale(1. / histStatYieldJpsiNormPtRun2.Integral(), "WIDTH")

    histSystYieldJpsiNormPtRun2 = TH1F("histSystYieldJpsiNormPtRun2", "", len(ptArrRun2)-1, ptArrRun2)
    for i in range(0, len(ptArrRun2)-1) : histSystYieldJpsiNormPtRun2.SetBinContent(i+1, yieldJpsiPtRun2[i]), histSystYieldJpsiNormPtRun2.SetBinError(i+1, systYieldJpsiPtRun2[i])
    SetGraSyst(histSystYieldJpsiNormPtRun2, 20, ROOT.kBlack)
    histSystYieldJpsiNormPtRun2.Scale(1. / histSystYieldJpsiNormPtRun2.Integral(), "WIDTH")

    histStatYieldPsi2sNormPtRun2 = TH1F("histStatYieldPsi2sNormPtRun2", "", len(ptArrRun2)-1, ptArrRun2)
    for i in range(0, len(ptArrRun2)-1) : histStatYieldPsi2sNormPtRun2.SetBinContent(i+1, yieldPsi2sPtRun2[i]), histStatYieldPsi2sNormPtRun2.SetBinError(i+1, statYieldPsi2sPtRun2[i])
    SetGraStat(histStatYieldPsi2sNormPtRun2, 20, ROOT.kBlack)
    histStatYieldPsi2sNormPtRun2.Scale(1. / histStatYieldPsi2sNormPtRun2.Integral(), "WIDTH")

    histSystYieldPsi2sNormPtRun2 = TH1F("histSystYieldPsi2sNormPtRun2", "", len(ptArrRun2)-1, ptArrRun2)
    for i in range(0, len(ptArrRun2)-1) : histSystYieldPsi2sNormPtRun2.SetBinContent(i+1, yieldPsi2sPtRun2[i]), histSystYieldPsi2sNormPtRun2.SetBinError(i+1, systYieldPsi2sPtRun2[i])
    SetGraSyst(histSystYieldPsi2sNormPtRun2, 20, ROOT.kBlack)
    histSystYieldPsi2sNormPtRun2.Scale(1. / histSystYieldPsi2sNormPtRun2.Integral(), "WIDTH")
        
    # Acceptance-efficiency
    histAxeJpsiPt = TH1F("histAxeJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histAxeJpsiPt.SetBinContent(i+1, axeJpsiPt[i]), histAxeJpsiPt.SetBinError(i+1, statAxeJpsiPt[i])
    SetGraStat(histAxeJpsiPt, 20, ROOT.kRed+1)

    histAxePsi2sPt = TH1F("histAxePsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histAxePsi2sPt.SetBinContent(i+1, axePsi2sPt[i]), histAxePsi2sPt.SetBinError(i+1, statAxePsi2sPt[i])
    SetGraStat(histAxePsi2sPt, 20, ROOT.kAzure+4)

    #############
    # Plot Yields 
    #############
    legendYieldJpsiVsPsi2sPt = TLegend(0.69, 0.69, 0.89, 0.89, " ", "brNDC")
    SetLegend(legendYieldJpsiVsPsi2sPt)
    legendYieldJpsiVsPsi2sPt.AddEntry(histSystYieldJpsiPt, "J/#psi", "FP")
    legendYieldJpsiVsPsi2sPt.AddEntry(histSystYieldPsi2sPt, "#psi(2S)", "FP")

    canvasYieldJpsiVsPsi2sPt = TCanvas("canvasYieldJpsiVsPsi2sPt", "canvasYieldJpsiVsPsi2sPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiPt  = TH2F("histGridYieldJpsiPt", "", 100, 0, 20, 100, 10, 5e5)
    histGridYieldJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (Gev/#it{c})")
    histGridYieldJpsiPt.GetYaxis().SetTitle("Raw yield")
    histGridYieldJpsiPt.Draw()
    histSystYieldJpsiPt.Draw("E2 SAME")
    histStatYieldJpsiPt.Draw("EP SAME")
    histSystYieldPsi2sPt.Draw("E2 SAME")
    histStatYieldPsi2sPt.Draw("EP SAME")
    legendYieldJpsiVsPsi2sPt.Draw("SAME")
    canvasYieldJpsiVsPsi2sPt.SaveAs("{}/jpsi_vs_psi2s_yield_vs_pt.pdf".format(outputDir))

    ########################
    # Plot Yields normalized 
    ########################
    legendYieldJpsiPt = TLegend(0.69, 0.69, 0.89, 0.89, " ", "brNDC")
    SetLegend(legendYieldJpsiPt)
    legendYieldJpsiPt.AddEntry(histSystYieldJpsiNormPt, "Run3", "FP")
    legendYieldJpsiPt.AddEntry(histSystYieldJpsiNormPtRun2, "Run2", "FP")

    canvasYieldJpsiPt = TCanvas("canvasYieldJpsiPt", "canvasYieldJpsiPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiPt  = TH2F("histGridYieldJpsiPt", "", 100, 0, 20, 100, 0.001, 1)
    histGridYieldJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (Gev/#it{c})")
    histGridYieldJpsiPt.GetYaxis().SetTitle("Normalized yield")
    histGridYieldJpsiPt.Draw()
    histSystYieldJpsiNormPtRun2.Draw("E2 SAME")
    histStatYieldJpsiNormPtRun2.Draw("EP SAME")
    histSystYieldJpsiNormPt.Draw("E2 SAME")
    histStatYieldJpsiNormPt.Draw("EP SAME")
    legendYieldJpsiPt.Draw("SAME")
    canvasYieldJpsiPt.SaveAs("{}/jpsi_norm_yield_vs_pt.pdf".format(outputDir))

    legendYieldPsi2sPt = TLegend(0.69, 0.69, 0.89, 0.89, " ", "brNDC")
    SetLegend(legendYieldPsi2sPt)
    legendYieldPsi2sPt.AddEntry(histSystYieldPsi2sNormPt, "Run3", "FP")
    legendYieldPsi2sPt.AddEntry(histSystYieldPsi2sNormPtRun2, "Run2", "FP")

    canvasYieldPsi2sPt = TCanvas("canvasYieldPsi2sPt", "canvasYieldPsi2sPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldPsi2sPt  = TH2F("histGridYieldPsi2sPt", "", 100, 0, 20, 100, 0.001, 1)
    histGridYieldPsi2sPt.GetXaxis().SetTitle("#it{p}_{T} (Gev/#it{c})")
    histGridYieldPsi2sPt.GetYaxis().SetTitle("Normalized yield")
    histGridYieldPsi2sPt.Draw()
    histSystYieldPsi2sNormPtRun2.Draw("E2 SAME")
    histStatYieldPsi2sNormPtRun2.Draw("EP SAME")
    histSystYieldPsi2sNormPt.Draw("E2 SAME")
    histStatYieldPsi2sNormPt.Draw("EP SAME")
    legendYieldPsi2sPt.Draw("SAME")
    canvasYieldPsi2sPt.SaveAs("{}/psi2s_norm_yield_vs_pt.pdf".format(outputDir))

    ################################
    # Plot the Ratio Psi(2S) / J/psi
    ################################
    measStatRatioPt, statMeasStatRatioPt = PropagateErrorsOnRatio(yieldJpsiPt, statYieldJpsiPt, yieldPsi2sPt, statYieldPsi2sPt)
    measSystRatioPt, systMeasSystRatioPt = PropagateErrorsOnRatio(yieldJpsiPt, systYieldJpsiPt, yieldPsi2sPt, systYieldPsi2sPt)

    print(ToCArray(measStatRatioPt, ctype='double', name='ratioPsi2sOverJpsiPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(statMeasStatRatioPt, ctype='double', name='statRatioPsi2sOverJpsiPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(systMeasSystRatioPt, ctype='double', name='systRatioPsi2sOverJpsiPt', formatter=lambda x: '{:0.5f}'.format(x)))

    graStatRatioPt = TGraphErrors(len(ptMin), ptCentr, measStatRatioPt, ptWidth, statMeasStatRatioPt)
    SetGraStat(graStatRatioPt, 20, ROOT.kRed+1)

    graSystRatioPt = TGraphErrors(len(ptMin), ptCentr, measSystRatioPt, ptWidth, systMeasSystRatioPt)
    SetGraSyst(graSystRatioPt, 20, ROOT.kRed+1)

    measStatRatioPtRun2, statMeasStatRatioPtRun2 = PropagateErrorsOnRatio(yieldJpsiPtRun2, statYieldJpsiPtRun2, yieldPsi2sPtRun2, statYieldPsi2sPtRun2)
    measSystRatioPtRun2, systMeasSystRatioPtRun2 = PropagateErrorsOnRatio(yieldJpsiPtRun2, systYieldJpsiPtRun2, yieldPsi2sPtRun2, systYieldPsi2sPtRun2)

    graStatRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, ratioPsi2sOverJpsiPt, ptWidth, statRatioPsi2sOverJpsiPt)
    SetGraStat(graStatRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)

    graSystRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, ratioPsi2sOverJpsiPt, ptWidth, systRatioPsi2sOverJpsiPt)
    SetGraSyst(graSystRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)
    

    # Print all results for analysis
    for i in range(len(ptMin)):
        relStatJpsi = (statYieldJpsiPt[i] / yieldJpsiPt[i]) * 100
        relSystJpsi = (systYieldJpsiPt[i] / yieldJpsiPt[i]) * 100
        relStatPsi2s = (statYieldPsi2sPt[i] / yieldPsi2sPt[i]) * 100
        relSystPsi2s = (systYieldPsi2sPt[i] / yieldPsi2sPt[i]) * 100
        relStatRatio = (statMeasStatRatioPt[i] / measStatRatioPt[i]) * 100
        relSystRatio = (systMeasSystRatioPt[i] / measStatRatioPt[i]) * 100
        print("& %3.2f - %3.2f & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (ptMin[i], ptMax[i], 
            yieldJpsiPt[i], statYieldJpsiPt[i], relStatJpsi, systYieldJpsiPt[i], relSystJpsi, 
            yieldPsi2sPt[i], statYieldPsi2sPt[i], relStatPsi2s, systYieldPsi2sPt[i], relSystPsi2s,
            measStatRatioPt[i], statMeasStatRatioPt[i], relStatRatio, statMeasStatRatioPt[i], relSystRatio))

    graStatRatioPtRun2 = TGraphErrors(len(ptMinRun2), ptCentrRun2, measStatRatioPtRun2, ptWidthRun2, statMeasStatRatioPtRun2)
    SetGraStat(graStatRatioPtRun2, 20, ROOT.kBlack)

    graSystRatioPtRun2 = TGraphErrors(len(ptMinRun2), ptCentrRun2, measSystRatioPtRun2, ptWidthRun2, systMeasSystRatioPtRun2)
    SetGraSyst(graSystRatioPtRun2, 20, ROOT.kBlack)

    legendRatioPt = TLegend(0.69, 0.49, 0.89, 0.69, " ", "brNDC")
    SetLegend(legendRatioPt)
    legendRatioPt.AddEntry(graSystRatioPt, "Run3 (method 1)", "FP")
    legendRatioPt.AddEntry(graSystRatioPsi2sOverJpsiPt, "Run3 (method 2)", "FP")
    legendRatioPt.AddEntry(graSystRatioPtRun2, "Run2", "FP")

    canvasRatioPt = TCanvas("canvasRatioPt", "canvasRatioPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridRatioPt  = TH2F("histGridRatioPt", "", 100, 0, 20, 100, 0.001, 0.10)
    histGridRatioPt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c}")
    histGridRatioPt.GetYaxis().SetTitle("#psi(2S) / J/#psi")
    histGridRatioPt.Draw()
    graSystRatioPtRun2.Draw("E2 SAME")
    graStatRatioPtRun2.Draw("EP SAME")
    graSystRatioPt.Draw("E2 SAME")
    graStatRatioPt.Draw("EP SAME")
    graSystRatioPsi2sOverJpsiPt.Draw("E2 SAME")
    graStatRatioPsi2sOverJpsiPt.Draw("EP SAME")
    legendRatioPt.Draw("EP SAME")
    canvasRatioPt.SaveAs("{}/psi2s_over_jpsi_vs_pt.pdf".format(outputDir))

    ######################################
    # Plot corrected Ratio Psi(2S) / J/psi
    ######################################
    legendAxePt = TLegend(0.59, 0.29, 0.79, 0.59, " ", "brNDC")
    SetLegend(legendAxePt)
    legendAxePt.AddEntry(histAxeJpsiPt, "J/#psi", "L")
    legendAxePt.AddEntry(histAxePsi2sPt, "#psi(2S)", "L")

    canvasAxePt = TCanvas("canvasAxePt", "canvasAxePt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridAxePt  = TH2F("histGridAxePt", "", 100, 0, 20, 100, 0.1, 1)
    histGridAxePt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c}")
    histGridAxePt.GetYaxis().SetTitle("A#times#varepsilon")
    histGridAxePt.Draw()
    histAxeJpsiPt.Draw("HE SAME")
    histAxePsi2sPt.Draw("HE SAME")
    legendAxePt.Draw("EP SAME")
    canvasAxePt.SaveAs("{}/axe_jpsi_psi2s_vs_pt.pdf".format(outputDir))

    # Compute and plot the corrected yield
    corrYieldJpsiPt, statCorrYieldJpsiPt = PropagateErrorsOnRatio(axeJpsiPt, statAxeJpsiPt, yieldJpsiPt, statYieldJpsiPt)
    corrYieldJpsiPt, systCorrYieldJpsiPt = PropagateErrorsOnRatio(axeJpsiPt, systAxeJpsiPt, yieldJpsiPt, systYieldJpsiPt)
    corrYieldJpsiPt = corrYieldJpsiPt / brJpsiToMuMu
    statCorrYieldJpsiPt = statCorrYieldJpsiPt / brJpsiToMuMu
    systCorrYieldJpsiPt = systCorrYieldJpsiPt / brJpsiToMuMu

    corrYieldPsi2sPt, statCorrYieldPsi2sPt = PropagateErrorsOnRatio(axePsi2sPt, statAxePsi2sPt, yieldPsi2sPt, statYieldPsi2sPt)
    corrYieldPsi2sPt, systCorrYieldPsi2sPt = PropagateErrorsOnRatio(axePsi2sPt, systAxePsi2sPt, yieldPsi2sPt, systYieldPsi2sPt)
    corrYieldPsi2sPt = corrYieldPsi2sPt / brPsi2sToMuMu
    statCorrYieldPsi2sPt = statCorrYieldPsi2sPt / brPsi2sToMuMu
    systCorrYieldPsi2sPt = systCorrYieldPsi2sPt / brPsi2sToMuMu

    histStatCorrYieldJpsiPt = TH1F("histStatCorrYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatCorrYieldJpsiPt.SetBinContent(i+1, corrYieldJpsiPt[i]), histStatCorrYieldJpsiPt.SetBinError(i+1, statCorrYieldJpsiPt[i])
    SetGraStat(histStatCorrYieldJpsiPt, 20, ROOT.kRed+1)
    histStatCorrYieldJpsiPt.Scale(1., "WIDTH")

    histSystCorrYieldJpsiPt = TH1F("histSystCorrYieldJpsiPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystCorrYieldJpsiPt.SetBinContent(i+1, corrYieldJpsiPt[i]), histSystCorrYieldJpsiPt.SetBinError(i+1, systCorrYieldJpsiPt[i])
    SetGraSyst(histSystCorrYieldJpsiPt, 20, ROOT.kRed+1)
    histSystCorrYieldJpsiPt.Scale(1., "WIDTH")

    histStatCorrYieldPsi2sPt = TH1F("histStatCorrYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histStatCorrYieldPsi2sPt.SetBinContent(i+1, corrYieldPsi2sPt[i]), histStatCorrYieldPsi2sPt.SetBinError(i+1, statCorrYieldPsi2sPt[i])
    SetGraStat(histStatCorrYieldPsi2sPt, 20, ROOT.kAzure+4)
    histStatCorrYieldPsi2sPt.Scale(1., "WIDTH")

    histSystCorrYieldPsi2sPt = TH1F("histSystCorrYieldPsi2sPt", "", len(ptArr)-1, ptArr)
    for i in range(0, len(ptArr)-1) : histSystCorrYieldPsi2sPt.SetBinContent(i+1, corrYieldPsi2sPt[i]), histSystCorrYieldPsi2sPt.SetBinError(i+1, systCorrYieldPsi2sPt[i])
    SetGraSyst(histSystCorrYieldPsi2sPt, 20, ROOT.kAzure+4)
    histSystCorrYieldPsi2sPt.Scale(1., "WIDTH")

    canvasCorrYieldJpsiPt = TCanvas("canvasCorrYieldJpsiPt", "canvasCorrYieldJpsiPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCorrYieldJpsiPt  = TH2F("histGridCorrYieldJpsiPt", "", 100, 0, 20, 100, 1e3, 1e9)
    histGridCorrYieldJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (Gev/#it{c})")
    histGridCorrYieldJpsiPt.GetYaxis().SetTitle("d#it{N} / (d#it{p}_{T} * A#times#varepsilon * BR)")
    histGridCorrYieldJpsiPt.Draw()
    histSystCorrYieldJpsiPt.Draw("E2 SAME")
    histStatCorrYieldJpsiPt.Draw("EP SAME")
    histSystCorrYieldPsi2sPt.Draw("E2 SAME")
    histStatCorrYieldPsi2sPt.Draw("EP SAME")
    legendYieldJpsiVsPsi2sPt.Draw("SAME")
    canvasCorrYieldJpsiPt.SaveAs("{}/jpsi_vs_psi2s_corr_yield_vs_pt.pdf".format(outputDir))

    # Compute and plot the cross section ratio
    csRatioPt, statCsRatioPt = PropagateErrorsOnRatio(corrYieldJpsiPt, statCorrYieldJpsiPt, corrYieldPsi2sPt, statCorrYieldPsi2sPt)
    csRatioPt, systCsRatioPt = PropagateErrorsOnRatio(corrYieldJpsiPt, systCorrYieldJpsiPt, corrYieldPsi2sPt, systCorrYieldPsi2sPt)

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
        systCsRatioPsi2sOverJpsiPt[iPt] = ((value * math.sqrt(relSystYield * relSystYield + relSystAxe * relSystAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))


    # Print all results for analysis    
    print("Ratio of the mean")
    for i in range(len(ptMin)):
        relStatJpsi = (statCsRatioPt[i] / csRatioPt[i]) * 100
        relSystJpsi = (systCsRatioPt[i] / csRatioPt[i]) * 100
        relStatPsi2s = (statCsRatioPt[i] / csRatioPt[i]) * 100
        relSystPsi2s = (systCsRatioPt[i] / csRatioPt[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (ptMin[i], ptMax[i], csRatioPt[i], statCsRatioPt[i], relStatJpsi, systCsRatioPt[i], relSystJpsi))

    print("Mean of the ratio")
    for i in range(len(ptMin)):
        relStatJpsi = (statCsRatioPsi2sOverJpsiPt[i] / csRatioPsi2sOverJpsiPt[i]) * 100
        relSystJpsi = (systCsRatioPsi2sOverJpsiPt[i] / csRatioPsi2sOverJpsiPt[i]) * 100
        relStatPsi2s = (statCsRatioPsi2sOverJpsiPt[i] / csRatioPsi2sOverJpsiPt[i]) * 100
        relSystPsi2s = (systCsRatioPsi2sOverJpsiPt[i] / csRatioPsi2sOverJpsiPt[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (ptMin[i], ptMax[i], csRatioPsi2sOverJpsiPt[i], statCsRatioPsi2sOverJpsiPt[i], relStatJpsi, systCsRatioPsi2sOverJpsiPt[i], relSystJpsi))

    graStatCsRatioPt = TGraphErrors(len(ptMin), ptCentr, csRatioPt, ptWidth, statCsRatioPt)
    SetGraStat(graStatCsRatioPt, 20, ROOT.kRed+1)

    graSystCsRatioPt = TGraphErrors(len(ptMin), ptCentr, csRatioPt, ptWidth, systCsRatioPt)
    SetGraSyst(graSystCsRatioPt, 20, ROOT.kRed+1)

    graStatCsRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, csRatioPsi2sOverJpsiPt, ptWidth, statCsRatioPsi2sOverJpsiPt)
    SetGraStat(graStatCsRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)

    graSystCsRatioPsi2sOverJpsiPt = TGraphErrors(len(ptMin), ptCentr, csRatioPsi2sOverJpsiPt, ptWidth, systCsRatioPsi2sOverJpsiPt)
    SetGraSyst(graSystCsRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)

    graStatCsRatioPtRun2 = TGraphErrors(len(ptMinRun2), ptCentrRun2, csRatioPtRun2, ptWidthRun2, statCsRatioPtRun2)
    SetGraStat(graStatCsRatioPtRun2, 20, ROOT.kBlack)

    graSystCsRatioPtRun2 = TGraphErrors(len(ptMinRun2), ptCentrRun2, csRatioPtRun2, ptWidthRun2, systCsRatioPtRun2)
    SetGraSyst(graSystCsRatioPtRun2, 20, ROOT.kBlack)

    legendCsRatioPt = TLegend(0.69, 0.49, 0.89, 0.69, " ", "brNDC")
    SetLegend(legendCsRatioPt)
    legendCsRatioPt.AddEntry(graSystRatioPt, "Run3 (method 1)", "FP")
    legendCsRatioPt.AddEntry(graSystCsRatioPsi2sOverJpsiPt, "Run3 (method 2)", "FP")
    legendCsRatioPt.AddEntry(graSystRatioPtRun2, "Run2", "FP")

    canvasCsRatioPt = TCanvas("canvasCsRatioPt", "canvasCsRatioPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsRatioPt  = TH2F("histGridCsRatioPt", "", 100, 0, 20, 100, 0.001, 1)
    histGridCsRatioPt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c}")
    histGridCsRatioPt.GetYaxis().SetTitle("d^{2}#sigma_{#psi(2S)}/d#it{p}_{T}d#it{y} / d^{2}#sigma_{J/#psi}/d#it{p}_{T}d#it{y}")
    histGridCsRatioPt.Draw()
    graSystCsRatioPtRun2.Draw("E2 SAME")
    graStatCsRatioPtRun2.Draw("EP SAME")
    graSystCsRatioPt.Draw("E2 SAME")
    graStatCsRatioPt.Draw("EP SAME")
    graSystCsRatioPsi2sOverJpsiPt.Draw("E2 SAME")
    graStatCsRatioPsi2sOverJpsiPt.Draw("EP SAME")
    legendCsRatioPt.Draw("EP SAME")
    canvasCsRatioPt.SaveAs("{}/cross_section_psi2s_over_jpsi_vs_pt.pdf".format(outputDir))

################
################
def yResults():
    LoadStyle()
    ROOT.gStyle.SetOptStat(0)

    ###############
    # Load datasets
    ###############
    # y-dependence
    dfYieldJpsiY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Jpsi_vs_y.txt'.format(inputDir), sep=' ')
    yMin = dfYieldJpsiY["x_min"].to_numpy()
    yMax = dfYieldJpsiY["x_max"].to_numpy()
    yArr = np.append(yMin, yMax[len(yMin)-1],)
    yCentr = (yMin + yMax) / 2.
    yWidth = (yMax - yMin) / 2.
    yieldJpsiY = dfYieldJpsiY["val"].to_numpy()
    statYieldJpsiY = dfYieldJpsiY["stat"].to_numpy()
    systYieldJpsiY = dfYieldJpsiY["syst"].to_numpy()
    print("Sum J/psi vs y: ", sum(yieldJpsiY))

    dfYieldPsi2sY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/{}/sig_Psi2s_vs_y.txt'.format(inputDir), sep=' ')
    yieldPsi2sY = dfYieldPsi2sY["val"].to_numpy()
    statYieldPsi2sY = dfYieldPsi2sY["stat"].to_numpy()
    systYieldPsi2sY = dfYieldPsi2sY["syst"].to_numpy()
    print("Sum Psi(2S) vs y: ", sum(yieldPsi2sY))

    dfRatioPsi2sOverJpsiY = pd.read_csv('/Users/lucamicheletti/cernbox/run3_Psi2s_over_Jpsi/final_results_LHC22all_periods/Comparison_methods/RatioValuesCorentinY.txt', sep=' ')
    ratioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["val"].to_numpy()
    statRatioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["stat"].to_numpy()
    systRatioPsi2sOverJpsiY = dfRatioPsi2sOverJpsiY["syst"].to_numpy()

    # Print all results for analysis
    for i in range(len(yMin)):
        relStatJpsi = (statYieldJpsiY[i] / yieldJpsiY[i]) * 100
        relSystJpsi = (systYieldJpsiY[i] / yieldJpsiY[i]) * 100
        relStatPsi2s = (statYieldPsi2sY[i] / yieldPsi2sY[i]) * 100
        relSystPsi2s = (systYieldPsi2sY[i] / yieldPsi2sY[i]) * 100
        print("& %3.2f - %3.2f & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) & %1.0f $\pm$ %1.0f (%3.2f?%%) $\pm$ %1.0f (%3.2f?%%) ??" % (yMin[i], yMax[i], 
            yieldJpsiY[i], statYieldJpsiY[i], relStatJpsi, systYieldJpsiY[i], relSystJpsi, 
            yieldPsi2sY[i], statYieldPsi2sY[i], relStatPsi2s, systYieldPsi2sY[i], relSystPsi2s))

    dfYieldJpsiY = pd.read_csv('run2_results/sig_Jpsi_vs_y_run2.txt', sep=' ')
    yMinRun2 = dfYieldJpsiY["x_min"].to_numpy()
    yMaxRun2 = dfYieldJpsiY["x_max"].to_numpy()
    yArrRun2 = np.append(yMinRun2, yMaxRun2[len(yMinRun2)-1],)
    yCentrRun2 = (yMinRun2 + yMaxRun2) / 2.
    yWidthRun2 = (yMaxRun2 - yMinRun2) / 2.
    yieldJpsiYRun2 = dfYieldJpsiY["val"].to_numpy()
    statYieldJpsiYRun2 = dfYieldJpsiY["stat"].to_numpy()
    systYieldJpsiYRun2 = dfYieldJpsiY["syst"].to_numpy()

    dfYieldPsi2sY = pd.read_csv('run2_results/sig_Psi2s_vs_y_run2.txt', sep=' ')
    yieldPsi2sYRun2 = dfYieldPsi2sY["val"].to_numpy()
    statYieldPsi2sYRun2 = dfYieldPsi2sY["stat"].to_numpy()
    systYieldPsi2sYRun2 = dfYieldPsi2sY["syst"].to_numpy()

    dfCsRatioY = pd.read_csv('run2_results/cs_ratio_Psi2s_Jpsi_vs_y_run2.txt', sep=' ')
    csRatioYRun2 = dfCsRatioY["val"].to_numpy()
    statCsRatioYRun2 = dfCsRatioY["stat"].to_numpy()
    systCsRatioYRun2 = dfCsRatioY["syst"].to_numpy()

    dfAxeJpsiY = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Jpsi_vs_y.txt', sep=' ')
    axeJpsiY = dfAxeJpsiY["val"].to_numpy()
    statAxeJpsiY = dfAxeJpsiY["stat"].to_numpy()
    systAxeJpsiY = dfAxeJpsiY["syst"].to_numpy()

    dfAxePsi2sY = pd.read_csv('acceptance_efficiency/ideal_no_cut/axe_Psi2s_vs_y.txt', sep=' ')
    axePsi2sY = dfAxePsi2sY["val"].to_numpy()
    statAxePsi2sY = dfAxePsi2sY["stat"].to_numpy()
    systAxePsi2sY = dfAxePsi2sY["syst"].to_numpy()

    dfAxeRatioPsi2sOverJpsiY = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/acceptance_efficiency/ideal_no_cut/axe_Psi2s_over_Jpsi_vs_y.txt', sep=' ')
    axeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["val"].to_numpy()
    statAxeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["stat"].to_numpy()
    systAxeRatioPsi2sOverJpsiY = dfAxeRatioPsi2sOverJpsiY["syst"].to_numpy()

    # Print all results for analysis
    for i in range(len(yMin)):
        relStatJpsi = (statAxeJpsiY[i] / axeJpsiY[i]) * 100
        relSystJpsi = (systAxeJpsiY[i] / axeJpsiY[i]) * 100
        relStatPsi2s = (statAxePsi2sY[i] / axePsi2sY[i]) * 100
        relSystPsi2s = (systAxePsi2sY[i] / axePsi2sY[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) & %4.3f $\pm$ %4.3f (%3.2f?%%) ??" % (yMin[i], yMax[i], 
            axeJpsiY[i], statAxeJpsiY[i], relStatJpsi,
            axePsi2sY[i], statAxePsi2sY[i], relStatPsi2s))

    ############################
    # Create and fill histograms
    ############################
    histStatYieldJpsiY = TH1F("histStatYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldJpsiY.SetBinContent(i+1, yieldJpsiY[i]), histStatYieldJpsiY.SetBinError(i+1, statYieldJpsiY[i])
    SetGraStat(histStatYieldJpsiY, 20, ROOT.kRed+1)
    histStatYieldJpsiY.Scale(1, "WIDTH")

    histSystYieldJpsiY = TH1F("histSystYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldJpsiY.SetBinContent(i+1, yieldJpsiY[i]), histSystYieldJpsiY.SetBinError(i+1, systYieldJpsiY[i])
    SetGraSyst(histSystYieldJpsiY, 20, ROOT.kRed+1)
    histSystYieldJpsiY.Scale(1, "WIDTH")

    histStatYieldJpsiNormY = TH1F("histStatYieldJpsiNormY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldJpsiNormY.SetBinContent(i+1, yieldJpsiY[i]), histStatYieldJpsiNormY.SetBinError(i+1, statYieldJpsiY[i])
    SetGraStat(histStatYieldJpsiNormY, 20, ROOT.kRed+1)
    histStatYieldJpsiNormY.Scale(1. / histStatYieldJpsiNormY.Integral())

    histSystYieldJpsiNormY = TH1F("histSystYieldJpsiNormY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldJpsiNormY.SetBinContent(i+1, yieldJpsiY[i]), histSystYieldJpsiNormY.SetBinError(i+1, systYieldJpsiY[i])
    SetGraSyst(histSystYieldJpsiNormY, 20, ROOT.kRed+1)
    histSystYieldJpsiNormY.Scale(1. / histSystYieldJpsiNormY.Integral())

    histStatYieldPsi2sY = TH1F("histStatYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldPsi2sY.SetBinContent(i+1, yieldPsi2sY[i]), histStatYieldPsi2sY.SetBinError(i+1, statYieldPsi2sY[i])
    SetGraStat(histStatYieldPsi2sY, 20, ROOT.kAzure+4)
    histStatYieldPsi2sY.Scale(1, "WIDTH")

    histSystYieldPsi2sY = TH1F("histSystYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldPsi2sY.SetBinContent(i+1, yieldPsi2sY[i]), histSystYieldPsi2sY.SetBinError(i+1, systYieldPsi2sY[i])
    SetGraSyst(histSystYieldPsi2sY, 20, ROOT.kAzure+4)
    histSystYieldPsi2sY.Scale(1, "WIDTH")

    histStatYieldPsi2sNormY = TH1F("histStatYieldPsi2sNormY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatYieldPsi2sNormY.SetBinContent(i+1, yieldPsi2sY[i]), histStatYieldPsi2sNormY.SetBinError(i+1, statYieldPsi2sY[i])
    SetGraStat(histStatYieldPsi2sNormY, 20, ROOT.kAzure+4)
    histStatYieldPsi2sNormY.Scale(1. / histStatYieldPsi2sNormY.Integral())

    histSystYieldPsi2sNormY = TH1F("histSystYieldPsi2sNormY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystYieldPsi2sNormY.SetBinContent(i+1, yieldPsi2sY[i]), histSystYieldPsi2sNormY.SetBinError(i+1, systYieldPsi2sY[i])
    SetGraSyst(histSystYieldPsi2sNormY, 20, ROOT.kAzure+4)
    histSystYieldPsi2sNormY.Scale(1. / histSystYieldPsi2sNormY.Integral())

    # Run2
    histStatYieldJpsiNormYRun2 = TH1F("histStatYieldJpsiNormYRun2", "", len(yArrRun2)-1, yArrRun2)
    for i in range(0, len(yArrRun2)-1) : histStatYieldJpsiNormYRun2.SetBinContent(i+1, yieldJpsiYRun2[i]), histStatYieldJpsiNormYRun2.SetBinError(i+1, statYieldJpsiYRun2[i])
    SetGraStat(histStatYieldJpsiNormYRun2, 20, ROOT.kBlack)
    histStatYieldJpsiNormYRun2.Scale(1. / histStatYieldJpsiNormYRun2.Integral())

    histSystYieldJpsiNormYRun2 = TH1F("histSystYieldJpsiNormYRun2", "", len(yArrRun2)-1, yArrRun2)
    for i in range(0, len(yArrRun2)-1) : histSystYieldJpsiNormYRun2.SetBinContent(i+1, yieldJpsiYRun2[i]), histSystYieldJpsiNormYRun2.SetBinError(i+1, systYieldJpsiYRun2[i])
    SetGraSyst(histSystYieldJpsiNormYRun2, 20, ROOT.kBlack)
    histSystYieldJpsiNormYRun2.Scale(1. / histSystYieldJpsiNormYRun2.Integral())

    histStatYieldPsi2sNormYRun2 = TH1F("histStatYieldPsi2sNormYRun2", "", len(yArrRun2)-1, yArrRun2)
    for i in range(0, len(yArrRun2)-1) : histStatYieldPsi2sNormYRun2.SetBinContent(i+1, yieldPsi2sYRun2[i]), histStatYieldPsi2sNormYRun2.SetBinError(i+1, statYieldPsi2sYRun2[i])
    SetGraStat(histStatYieldPsi2sNormYRun2, 20, ROOT.kBlack)
    histStatYieldPsi2sNormYRun2.Scale(1. / histStatYieldPsi2sNormYRun2.Integral())

    histSystYieldPsi2sNormYRun2 = TH1F("histSystYieldPsi2sNormYRun2", "", len(yArrRun2)-1, yArrRun2)
    for i in range(0, len(yArrRun2)-1) : histSystYieldPsi2sNormYRun2.SetBinContent(i+1, yieldPsi2sYRun2[i]), histSystYieldPsi2sNormYRun2.SetBinError(i+1, systYieldPsi2sYRun2[i])
    SetGraSyst(histSystYieldPsi2sNormYRun2, 20, ROOT.kBlack)
    histSystYieldPsi2sNormYRun2.Scale(1. / histSystYieldPsi2sNormYRun2.Integral())

    # Acceptance-efficiency
    histAxeJpsiY = TH1F("histAxeJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histAxeJpsiY.SetBinContent(i+1, axeJpsiY[i]), histAxeJpsiY.SetBinError(i+1, statAxeJpsiY[i])
    SetGraStat(histAxeJpsiY, 20, ROOT.kRed+1)

    histAxePsi2sY = TH1F("histAxePsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histAxePsi2sY.SetBinContent(i+1, axePsi2sY[i]), histAxePsi2sY.SetBinError(i+1, statAxePsi2sY[i])
    SetGraStat(histAxePsi2sY, 20, ROOT.kAzure+4)

    #############
    # Plot Yields 
    #############
    legendYieldJpsiVsPsi2sY = TLegend(0.49, 0.29, 0.69, 0.49, " ", "brNDC")
    SetLegend(legendYieldJpsiVsPsi2sY)
    legendYieldJpsiVsPsi2sY.AddEntry(histSystYieldJpsiY, "J/#psi", "FP")
    legendYieldJpsiVsPsi2sY.AddEntry(histSystYieldPsi2sY, "#psi(2S)", "FP")

    canvasYieldJpsiY = TCanvas("canvasYieldJpsiY", "canvasYieldJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiY  = TH2F("histGridYieldJpsiY", "", 100, 2.5, 4, 100, 10, 5e6)
    histGridYieldJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridYieldJpsiY.GetYaxis().SetTitle("Raw yield")
    histGridYieldJpsiY.Draw()
    histSystYieldJpsiY.Draw("E2 SAME")
    histStatYieldJpsiY.Draw("EP SAME")
    histSystYieldPsi2sY.Draw("E2 SAME")
    histStatYieldPsi2sY.Draw("EP SAME")
    legendYieldJpsiVsPsi2sY.Draw("SAME")
    canvasYieldJpsiY.SaveAs("{}/jpsi_vs_psi2s_yield_vs_y.pdf".format(outputDir))

    ########################
    # Plot Yields normalized 
    ########################

    legendYieldJpsiY = TLegend(0.69, 0.69, 0.89, 0.89, " ", "brNDC")
    SetLegend(legendYieldJpsiY)
    legendYieldJpsiY.AddEntry(histSystYieldJpsiNormY, "Run3", "FP")
    legendYieldJpsiY.AddEntry(histSystYieldJpsiNormYRun2, "Run2", "FP")

    canvasYieldJpsiY = TCanvas("canvasYieldJpsiY", "canvasYieldJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldJpsiY  = TH2F("histGridYieldJpsiY", "", 100, 2.5, 4, 100, 0.03, 1)
    histGridYieldJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridYieldJpsiY.GetYaxis().SetTitle("Normalized yield")
    histGridYieldJpsiY.Draw()
    histSystYieldJpsiNormYRun2.Draw("E2 SAME")
    histStatYieldJpsiNormYRun2.Draw("EP SAME")
    histSystYieldJpsiNormY.Draw("E2 SAME")
    histStatYieldJpsiNormY.Draw("EP SAME")
    legendYieldJpsiY.Draw("SAME")
    canvasYieldJpsiY.SaveAs("{}/jpsi_norm_yield_vs_y.pdf".format(outputDir))

    legendYieldPsi2sY = TLegend(0.29, 0.29, 0.49, 0.49, " ", "brNDC")
    SetLegend(legendYieldPsi2sY)
    legendYieldPsi2sY.AddEntry(histSystYieldPsi2sNormY, "Run3", "FP")
    legendYieldPsi2sY.AddEntry(histSystYieldPsi2sNormYRun2, "Run2", "FP")

    canvasYieldPsi2sY = TCanvas("canvasYieldPsi2sY", "canvasYieldPsi2sY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridYieldPsi2sY  = TH2F("histGridYieldPsi2sY", "", 100, 2.5, 4, 100, 0.005, 1)
    histGridYieldPsi2sY.GetXaxis().SetTitle("#it{y}")
    histGridYieldPsi2sY.GetYaxis().SetTitle("Normalized yield")
    histGridYieldPsi2sY.Draw()
    histSystYieldPsi2sNormYRun2.Draw("E2 SAME")
    histStatYieldPsi2sNormYRun2.Draw("EP SAME")
    histSystYieldPsi2sNormY.Draw("E2 SAME")
    histStatYieldPsi2sNormY.Draw("EP SAME")
    legendYieldPsi2sY.Draw("SAME")
    canvasYieldPsi2sY.SaveAs("{}/psi2s_norm_yield_vs_y.pdf".format(outputDir))

    ################################
    # Plot the Ratio Psi(2S) / J/psi
    ################################
    measStatRatioY, statMeasStatRatioY = PropagateErrorsOnRatio(yieldJpsiY, statYieldJpsiY, yieldPsi2sY, statYieldPsi2sY)
    measSystRatioY, systMeasSystRatioY = PropagateErrorsOnRatio(yieldJpsiY, systYieldJpsiY, yieldPsi2sY, systYieldPsi2sY)

    print(ToCArray(measStatRatioY, ctype='double', name='ratioPsi2sOverJpsiY', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(statMeasStatRatioY, ctype='double', name='statRatioPsi2sOverJpsiY', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(systMeasSystRatioY, ctype='double', name='systRatioPsi2sOverJpsiY', formatter=lambda x: '{:0.5f}'.format(x)))

    graStatRatioY = TGraphErrors(len(yMin), yCentr, measStatRatioY, yWidth, statMeasStatRatioY)
    SetGraStat(graStatRatioY, 20, ROOT.kRed+1)

    graSystRatioY = TGraphErrors(len(yMin), yCentr, measSystRatioY, yWidth, systMeasSystRatioY)
    SetGraSyst(graSystRatioY, 20, ROOT.kRed+1)

    graStatRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, ratioPsi2sOverJpsiY, yWidth, statRatioPsi2sOverJpsiY)
    SetGraStat(graStatRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)

    graSystRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, ratioPsi2sOverJpsiY, yWidth, systRatioPsi2sOverJpsiY)
    SetGraSyst(graSystRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)

    measStatRatioYRun2, statMeasStatRatioYRun2 = PropagateErrorsOnRatio(yieldJpsiYRun2, statYieldJpsiYRun2, yieldPsi2sYRun2, statYieldPsi2sYRun2)
    measSystRatioYRun2, systMeasSystRatioYRun2 = PropagateErrorsOnRatio(yieldJpsiYRun2, systYieldJpsiYRun2, yieldPsi2sYRun2, systYieldPsi2sYRun2)

    # Print all results for analysis
    graStatRatioYRun2 = TGraphErrors(len(yMinRun2), yCentrRun2, measStatRatioYRun2, yWidthRun2, statMeasStatRatioYRun2)
    SetGraStat(graStatRatioYRun2, 20, ROOT.kBlack)

    graSystRatioYRun2 = TGraphErrors(len(yMinRun2), yCentrRun2, measSystRatioYRun2, yWidthRun2, systMeasSystRatioYRun2)
    SetGraSyst(graSystRatioYRun2, 20, ROOT.kBlack)

    legendRatioY = TLegend(0.29, 0.29, 0.49, 0.49, " ", "brNDC")
    SetLegend(legendRatioY)
    legendRatioY.AddEntry(graSystRatioY, "Run3 (method 1)", "FP")
    legendRatioY.AddEntry(graSystRatioPsi2sOverJpsiY, "Run3 (method 2)", "FP")
    legendRatioY.AddEntry(graSystRatioYRun2, "Run2", "FP")
        
    canvasRatioY = TCanvas("canvasRatioY", "canvasRatioY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridRatioY  = TH2F("histGridRatioY", "", 100, 2.5, 4, 100, 0, 0.10)
    histGridRatioY.GetXaxis().SetTitle("#it{y}")
    histGridRatioY.GetYaxis().SetTitle("#psi(2S) / J/#psi")
    histGridRatioY.Draw()
    graSystRatioYRun2.Draw("E2 SAME")
    graStatRatioYRun2.Draw("EP SAME")
    graSystRatioY.Draw("E2 SAME")
    graStatRatioY.Draw("EP SAME")
    graSystRatioPsi2sOverJpsiY.Draw("E2 SAME")
    graStatRatioPsi2sOverJpsiY.Draw("EP SAME")
    legendRatioY.Draw("SAME")
    canvasRatioY.SaveAs("{}/psi2s_over_jpsi_vs_y.pdf".format(outputDir))

    ######################################
    # Plot corrected Ratio Psi(2S) / J/psi
    ######################################
    legendAxeY = TLegend(0.29, 0.19, 0.49, 0.49, " ", "brNDC")
    SetLegend(legendAxeY)
    legendAxeY.AddEntry(histAxeJpsiY, "J/#psi (Run2 MC)", "L")
    legendAxeY.AddEntry(histAxePsi2sY, "#psi(2S) (Run2 MC)", "L")

    canvasAxeY = TCanvas("canvasAxeY", "canvasAxeY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridAxeY  = TH2F("histGridAxeY", "", 100, 2.5, 4, 100, 0.01, 1)
    histGridAxeY.GetXaxis().SetTitle("#it{y}")
    histGridAxeY.GetYaxis().SetTitle("A#times#varepsilon")
    histGridAxeY.Draw()
    histAxeJpsiY.Draw("HE SAME")
    histAxePsi2sY.Draw("HE SAME")
    legendAxeY.Draw("SAME")
    canvasAxeY.SaveAs("{}/axe_jpsi_psi2s_vs_y.pdf".format(outputDir))

    # Compute and plot the corrected yield
    corrYieldJpsiY, statCorrYieldJpsiY = PropagateErrorsOnRatio(axeJpsiY, statAxeJpsiY, yieldJpsiY, statYieldJpsiY)
    corrYieldJpsiY, systCorrYieldJpsiY = PropagateErrorsOnRatio(axeJpsiY, systAxeJpsiY, yieldJpsiY, systYieldJpsiY)
    corrYieldJpsiY = corrYieldJpsiY / brJpsiToMuMu
    statCorrYieldJpsiY = statCorrYieldJpsiY / brJpsiToMuMu
    systCorrYieldJpsiY = systCorrYieldJpsiY / brJpsiToMuMu

    corrYieldPsi2sY, statCorrYieldPsi2sY = PropagateErrorsOnRatio(axePsi2sY, statAxePsi2sY, yieldPsi2sY, statYieldPsi2sY)
    corrYieldPsi2sY, systCorrYieldPsi2sY = PropagateErrorsOnRatio(axePsi2sY, systAxePsi2sY, yieldPsi2sY, systYieldPsi2sY)
    corrYieldPsi2sY = corrYieldPsi2sY / brPsi2sToMuMu
    statCorrYieldPsi2sY = statCorrYieldPsi2sY / brPsi2sToMuMu
    systCorrYieldPsi2sY = systCorrYieldPsi2sY / brPsi2sToMuMu

    histStatCorrYieldJpsiY = TH1F("histStatCorrYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatCorrYieldJpsiY.SetBinContent(i+1, corrYieldJpsiY[i]), histStatCorrYieldJpsiY.SetBinError(i+1, statCorrYieldJpsiY[i])
    SetGraStat(histStatCorrYieldJpsiY, 20, ROOT.kRed+1)
    histStatCorrYieldJpsiY.Scale(1., "WIDTH")

    histSystCorrYieldJpsiY = TH1F("histSystCorrYieldJpsiY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystCorrYieldJpsiY.SetBinContent(i+1, corrYieldJpsiY[i]), histSystCorrYieldJpsiY.SetBinError(i+1, systCorrYieldJpsiY[i])
    SetGraSyst(histSystCorrYieldJpsiY, 20, ROOT.kRed+1)
    histSystCorrYieldJpsiY.Scale(1., "WIDTH")

    histStatCorrYieldPsi2sY = TH1F("histStatCorrYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histStatCorrYieldPsi2sY.SetBinContent(i+1, corrYieldPsi2sY[i]), histStatCorrYieldPsi2sY.SetBinError(i+1, statCorrYieldPsi2sY[i])
    SetGraStat(histStatCorrYieldPsi2sY, 20, ROOT.kAzure+4)
    histStatCorrYieldPsi2sY.Scale(1., "WIDTH")

    histSystCorrYieldPsi2sY = TH1F("histSystCorrYieldPsi2sY", "", len(yArr)-1, yArr)
    for i in range(0, len(yArr)-1) : histSystCorrYieldPsi2sY.SetBinContent(i+1, corrYieldPsi2sY[i]), histSystCorrYieldPsi2sY.SetBinError(i+1, systCorrYieldPsi2sY[i])
    SetGraSyst(histSystCorrYieldPsi2sY, 20, ROOT.kAzure+4)
    histSystCorrYieldPsi2sY.Scale(1., "WIDTH")

    canvasCorrYieldJpsiY = TCanvas("canvasCorrYieldJpsiY", "canvasCorrYieldJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCorrYieldJpsiY  = TH2F("histGridCorrYieldJpsiY", "", 100, 2.5, 4, 100, 1e3, 1e9)
    histGridCorrYieldJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridCorrYieldJpsiY.GetYaxis().SetTitle("d#it{N} / (d#it{p}_{T} * A#times#varepsilon * BR)")
    histGridCorrYieldJpsiY.Draw()
    histSystCorrYieldJpsiY.Draw("E2 SAME")
    histStatCorrYieldJpsiY.Draw("EP SAME")
    histSystCorrYieldPsi2sY.Draw("E2 SAME")
    histStatCorrYieldPsi2sY.Draw("EP SAME")
    legendYieldJpsiVsPsi2sY.Draw("SAME")
    canvasCorrYieldJpsiY.SaveAs("{}/jpsi_vs_psi2s_corr_yield_vs_y.pdf".format(outputDir))

    # Compute and plot the cross section ratio
    csRatioY, statCsRatioY = PropagateErrorsOnRatio(corrYieldJpsiY, statCorrYieldJpsiY, corrYieldPsi2sY, statCorrYieldPsi2sY)
    csRatioY, systCsRatioY = PropagateErrorsOnRatio(corrYieldJpsiY, systCorrYieldJpsiY, corrYieldPsi2sY, systCorrYieldPsi2sY)

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
        systCsRatioPsi2sOverJpsiY[iY] = ((value * math.sqrt(relSystYield * relSystYield + relSystAxe * relSystAxe)) * (brJpsiToMuMu / brPsi2sToMuMu))

    # Print all results for analysis    
    print("Ratio of the mean")
    for i in range(len(yMin)):
        relStatJpsi = (statCsRatioY[i] / csRatioY[i]) * 100
        relSystJpsi = (systCsRatioY[i] / csRatioY[i]) * 100
        relStatPsi2s = (statCsRatioY[i] / csRatioY[i]) * 100
        relSystPsi2s = (systCsRatioY[i] / csRatioY[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (yMin[i], yMax[i], csRatioY[i], statCsRatioY[i], relStatJpsi, systCsRatioY[i], relSystJpsi))

    print("Mean of the ratio")
    for i in range(len(yMin)):
        relStatJpsi = (statCsRatioPsi2sOverJpsiY[i] / csRatioPsi2sOverJpsiY[i]) * 100
        relSystJpsi = (systCsRatioPsi2sOverJpsiY[i] / csRatioPsi2sOverJpsiY[i]) * 100
        relStatPsi2s = (statCsRatioPsi2sOverJpsiY[i] / csRatioPsi2sOverJpsiY[i]) * 100
        relSystPsi2s = (systCsRatioPsi2sOverJpsiY[i] / csRatioPsi2sOverJpsiY[i]) * 100
        print("& %3.2f - %3.2f & %4.3f $\pm$ %4.3f (%3.2f?%%) $\pm$ %4.3f (%3.2f?%%) ??" % (yMin[i], yMax[i], csRatioPsi2sOverJpsiY[i], statCsRatioPsi2sOverJpsiY[i], relStatJpsi, systCsRatioPsi2sOverJpsiY[i], relSystJpsi))

    graStatCsRatioY = TGraphErrors(len(yMin), yCentr, csRatioY, yWidth, statCsRatioY)
    SetGraStat(graStatCsRatioY, 20, ROOT.kRed+1)

    graSystCsRatioY = TGraphErrors(len(yMin), yCentr, csRatioY, yWidth, systCsRatioY)
    SetGraSyst(graSystCsRatioY, 20, ROOT.kRed+1)

    graStatCsRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, csRatioPsi2sOverJpsiY, yWidth, statCsRatioPsi2sOverJpsiY)
    SetGraStat(graStatCsRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)

    graSystCsRatioPsi2sOverJpsiY = TGraphErrors(len(yMin), yCentr, csRatioPsi2sOverJpsiY, yWidth, systCsRatioPsi2sOverJpsiY)
    SetGraSyst(graSystCsRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)

    graStatCsRatioYRun2 = TGraphErrors(len(yMinRun2), yCentrRun2, csRatioYRun2, yWidthRun2, statCsRatioYRun2)
    SetGraStat(graStatCsRatioYRun2, 20, ROOT.kBlack)

    graSystCsRatioYRun2 = TGraphErrors(len(yMinRun2), yCentrRun2, csRatioYRun2, yWidthRun2, systCsRatioYRun2)
    SetGraSyst(graSystCsRatioYRun2, 20, ROOT.kBlack)

    legendCsRatioY = TLegend(0.69, 0.39, 0.89, 0.59, " ", "brNDC")
    SetLegend(legendCsRatioY)
    legendCsRatioY.AddEntry(graSystRatioY, "Run3 (method 1)", "FP")
    legendCsRatioY.AddEntry(graSystCsRatioPsi2sOverJpsiY, "Run3 (method 2)", "FP")
    legendCsRatioY.AddEntry(graSystRatioYRun2, "Run2", "FP")

    canvasCsRatioY = TCanvas("canvasCsRatioY", "canvasCsRatioY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsRatioY  = TH2F("histGridCsRatioY", "", 100, 2.5, 4, 100, 0, 0.35)
    histGridCsRatioY.GetXaxis().SetTitle("#it{y}")
    histGridCsRatioY.GetYaxis().SetTitle("d#sigma_{#psi(2S)}/d#it{y} / d#sigma_{J/#psi}/d#it{y}")
    histGridCsRatioY.Draw()
    graSystCsRatioYRun2.Draw("E2 SAME")
    graStatCsRatioYRun2.Draw("EP SAME")
    graSystCsRatioY.Draw("E2 SAME")
    graStatCsRatioY.Draw("EP SAME")
    graSystCsRatioPsi2sOverJpsiY.Draw("E2 SAME")
    graStatCsRatioPsi2sOverJpsiY.Draw("EP SAME")
    legendCsRatioY.Draw("EP SAME")
    canvasCsRatioY.SaveAs("{}/cross_section_psi2s_over_jpsi_vs_y.pdf".format(outputDir))


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--int_res", help="plot results integrated", action="store_true")
    parser.add_argument("--pt_res", help="plot results vs pT", action="store_true")
    parser.add_argument("--y_res", help="plot results vs y", action="store_true")
    args = parser.parse_args()
    print(args)

    if args.int_res:
        intResults()

    if args.pt_res:
        ptResults()

    if args.y_res:
        yResults()


if __name__ == '__main__':
    main()