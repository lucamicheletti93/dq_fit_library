import yaml
import json
import sys
import csv
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

    dirInNames = [
        "forward-01",
        "forward-02",
        "forward-03",
        "forward-04",
        "forward-05",
        "forward-06"
    ]

    fInNames = [
        "central",
        "double_renorm",
        "half_renorm",
        "high_mass",
        "low_mass",
    ]


    dfCsJpsiIcemPt = pd.DataFrame()
    dfCsPsi2sIcemPt = pd.DataFrame()
    dfCsJpsiIcemY = pd.DataFrame()
    dfCsPsi2sIcemY = pd.DataFrame()

    for iDir, dirInName in enumerate(dirInNames, start=1):
        ptCentr = []
        valCentr = []

        for iFile, fInName in enumerate(fInNames, start=1):

            with open("ICEM/1-prompt Jpsi pt-dist/13.6 TeV/{}/{}.txt".format(dirInName, fInName), 'r') as fIn:
                for line in fIn:
                    values = line.strip().split()
                    if len(values) >= 2:
                        ptCentr.append(float(values[0]))
                        valCentr.append(float(values[1]))

            if dirInName == "forward-01" and fInName == "central":
                dfCsJpsiIcemPt['ptCentr'] = ptCentr
                dfCsJpsiIcemPt["{}{}".format(fInName, iDir)] = valCentr
                dfCsJpsiIcemY['ptCentr'] = ptCentr
                dfCsJpsiIcemY["{}{}".format(fInName, iDir)] = valCentr
            else:
                dfCsJpsiIcemPt["{}{}".format(fInName, iDir)] = valCentr
                dfCsJpsiIcemY["{}{}".format(fInName, iDir)] = valCentr
            
            ptCentr = []
            valCentr = []
            with open("ICEM/2-psi2S pt-dist/13.6 TeV/{}/{}.txt".format(dirInName, fInName), 'r') as fIn:
                for line in fIn:
                    values = line.strip().split()
                    if len(values) >= 2:
                        ptCentr.append(float(values[0]))
                        valCentr.append(float(values[1]))

            if dirInName == "forward-01" and fInName == "central":
                dfCsPsi2sIcemPt['ptCentr'] = ptCentr
                dfCsPsi2sIcemPt["{}{}".format(fInName, iDir)] = valCentr
                dfCsPsi2sIcemY['ptCentr'] = ptCentr
                dfCsPsi2sIcemY["{}{}".format(fInName, iDir)] = valCentr
            else:
                dfCsPsi2sIcemPt["{}{}".format(fInName, iDir)] = valCentr
                dfCsPsi2sIcemY["{}{}".format(fInName, iDir)] = valCentr
            
            ptCentr = []
            valCentr = []
            
    centralCols = ['central1', 'central2', 'central3', 'central4', 'central5', 'central6']
    doubleRenormCols = ['double_renorm1', 'double_renorm2', 'double_renorm3', 'double_renorm4', 'double_renorm5', 'double_renorm6']
    halfRenormCols = ['half_renorm1', 'half_renorm2', 'half_renorm3', 'half_renorm4', 'half_renorm5', 'half_renorm6']
    highMassCols = ['high_mass1', 'high_mass2', 'high_mass3', 'high_mass4', 'high_mass5', 'high_mass6']
    lowMassCols = ['low_mass1', 'low_mass2', 'low_mass3', 'low_mass4', 'low_mass5', 'low_mass6']

    ptMins = [0, 1, 2, 3, 4, 5, 7, 10]
    ptMaxs = [1, 2, 3, 4, 5, 7, 10, 20]
    yMins = [2.50, 2.75, 3.00, 3.25, 3.50, 3.75]
    yMaxs = [2.75, 3.00, 3.25, 3.50, 3.75, 4.00]

    yMinArr = np.asarray(yMins)
    yMaxArr = np.asarray(yMaxs)
    yCentrArr = (yMinArr + yMaxArr) / 2.
    yWidthArr = (yMaxArr - yMinArr) / 2.

    ptMinArr = np.asarray(ptMins)
    ptMaxArr = np.asarray(ptMaxs)
    ptCentrArr = (ptMinArr + ptMaxArr) / 2.
    ptWidthArr = (ptMaxArr - ptMinArr) / 2.

    # creating a new column to store the pT binning used in data
    for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        if iPt == 0:
            dfCsJpsiIcemPt.loc[(dfCsJpsiIcemPt.ptCentr >= ptMin) & (dfCsJpsiIcemPt.ptCentr <= ptMax), 'ptGroup'] = iPt
            dfCsPsi2sIcemPt.loc[(dfCsPsi2sIcemPt.ptCentr >= ptMin) & (dfCsPsi2sIcemPt.ptCentr <= ptMax), 'ptGroup'] = iPt
        else:
            dfCsJpsiIcemPt.loc[(dfCsJpsiIcemPt.ptCentr > ptMin) & (dfCsJpsiIcemPt.ptCentr <= ptMax), 'ptGroup'] = iPt
            dfCsPsi2sIcemPt.loc[(dfCsPsi2sIcemPt.ptCentr > ptMin) & (dfCsPsi2sIcemPt.ptCentr <= ptMax), 'ptGroup'] = iPt

    #######################################
    # Pt dependence
    #######################################
    dfCsJpsiIcemPt['ptWeights'] = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.25]
    dfCsPsi2sIcemPt['ptWeights'] = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.25]
    dfCsJpsiIcemPt['centralSum'] = dfCsJpsiIcemPt[centralCols].sum(axis=1)
    dfCsPsi2sIcemPt['centralSum'] = dfCsPsi2sIcemPt[centralCols].sum(axis=1)

    #dfMergedCsPsi2s = dfCsPsi2sIcemPt.groupby('ptGroup').agg({'centralSumPtWeight': 'sum'}).reset_index()
    #csCentralJpsiIcemPt = (dfMergedCsJpsi['centralSumPtWeight'].to_numpy()) * 1e-3
    #csCentralPsi2sIcemPt = (dfMergedCsPsi2s['centralSumPtWeight'].to_numpy()) * 1e-3
    #csCentralPsi2sOverJpsiIcemPt = csCentralPsi2sIcemPt / csCentralJpsiIcemPt
    #print(ToCArray(csCentralJpsiIcemPt, ctype='double', name='Ramona_136TeV_Central_jpsi', formatter=lambda x: '{:0.5f}'.format(x)))
    #print(ToCArray(csCentralPsi2sIcemPt, ctype='double', name='Ramona_136TeV_Central_psi2S', formatter=lambda x: '{:0.5f}'.format(x)))
    
    # double renorm
    dfCsJpsiIcemPt['doubleRenormSum'] = dfCsJpsiIcemPt[doubleRenormCols].sum(axis=1)
    dfCsPsi2sIcemPt['doubleRenormSum'] = dfCsPsi2sIcemPt[doubleRenormCols].sum(axis=1)
    #dfMergedCsJpsi = dfCsJpsiIcemPt.groupby('ptGroup').agg({'doubleRenormSum': 'sum'}).reset_index()
    #dfMergedCsPsi2s = dfCsPsi2sIcemPt.groupby('ptGroup').agg({'doubleRenormSum': 'sum'}).reset_index()
    #csDoubleRenormJpsiIcemPt = (dfMergedCsJpsi['doubleRenormSum'].to_numpy()) * 1e-3
    #csDoubleRenormPsi2sIcemPt = (dfMergedCsPsi2s['doubleRenormSum'].to_numpy()) * 1e-3
    #csDoubleRenormPsi2sOverJpsiIcemPt = csDoubleRenormPsi2sIcemPt / csDoubleRenormJpsiIcemPt
    #print(ToCArray(csDoubleRenormJpsiIcemPt, ctype='double', name='Ramona_136TeV_Scale_High_jpsi', formatter=lambda x: '{:0.5f}'.format(x)))
    #print(ToCArray(csDoubleRenormPsi2sIcemPt, ctype='double', name='Ramona_136TeV_Scale_High_psi2S', formatter=lambda x: '{:0.5f}'.format(x)))

    #deltaCsDoubleRenormJpsiIcemPt = np.abs(csDoubleRenormJpsiIcemPt - csCentralJpsiIcemPt)
    #deltaCsDoubleRenormPsi2sIcemPt = np.abs(csDoubleRenormPsi2sIcemPt - csCentralPsi2sIcemPt)
    #relErrCsDoubleRenormPsi2sOverJpsiIcemPt = (deltaCsDoubleRenormJpsiIcemPt / csCentralJpsiIcemPt) - (deltaCsDoubleRenormPsi2sIcemPt / csCentralPsi2sIcemPt)
    
    # half renorm
    dfCsJpsiIcemPt['halfRenormSum'] = dfCsJpsiIcemPt[halfRenormCols].sum(axis=1)
    dfCsPsi2sIcemPt['halfRenormSum'] = dfCsPsi2sIcemPt[halfRenormCols].sum(axis=1)
    #dfMergedCsJpsi = dfCsJpsiIcemPt.groupby('ptGroup').agg({'halfRenormSum': 'sum'}).reset_index()
    #dfMergedCsPsi2s = dfCsPsi2sIcemPt.groupby('ptGroup').agg({'halfRenormSum': 'sum'}).reset_index()
    #csHalfRenormJpsiIcemPt = (dfMergedCsJpsi['halfRenormSum'].to_numpy()) * 1e-3
    #csHalfRenormPsi2sIcemPt = (dfMergedCsPsi2s['halfRenormSum'].to_numpy()) * 1e-3
    #csHalfRenormPsi2sOverJpsiIcemPt = csHalfRenormPsi2sIcemPt / csHalfRenormJpsiIcemPt
    #print(ToCArray(csHalfRenormJpsiIcemPt, ctype='double', name='Ramona_136TeV_Scale_Low_jpsi', formatter=lambda x: '{:0.5f}'.format(x)))
    #print(ToCArray(csHalfRenormPsi2sIcemPt, ctype='double', name='Ramona_136TeV_Scale_Low_psi2S', formatter=lambda x: '{:0.5f}'.format(x)))

    #deltaCsHalfRenormJpsiIcemPt = np.abs(csHalfRenormJpsiIcemPt - csCentralJpsiIcemPt)
    #deltaCsHalfRenormPsi2sIcemPt = np.abs(csHalfRenormPsi2sIcemPt - csCentralPsi2sIcemPt)
    #relErrCsHalfRenormPsi2sOverJpsiIcemPt = (deltaCsHalfRenormJpsiIcemPt / csCentralJpsiIcemPt) - (deltaCsHalfRenormPsi2sIcemPt / csCentralPsi2sIcemPt)

    # high mass
    dfCsJpsiIcemPt['highMassSum'] = dfCsJpsiIcemPt[highMassCols].sum(axis=1)
    dfCsPsi2sIcemPt['highMassSum'] = dfCsPsi2sIcemPt[highMassCols].sum(axis=1)
    #dfMergedCsJpsi = dfCsJpsiIcemPt.groupby('ptGroup').agg({'highMassSum': 'sum'}).reset_index()
    #dfMergedCsPsi2s = dfCsPsi2sIcemPt.groupby('ptGroup').agg({'highMassSum': 'sum'}).reset_index()
    #csHighMassJpsiIcemPt = (dfMergedCsJpsi['highMassSum'].to_numpy()) * 1e-3
    #csHighMassPsi2sIcemPt = (dfMergedCsPsi2s['highMassSum'].to_numpy()) * 1e-3
    #csHighMassPsi2sOverJpsiIcemPt = csHighMassPsi2sIcemPt / csHighMassJpsiIcemPt
    #print(ToCArray(csHighMassJpsiIcemPt, ctype='double', name='Ramona_136TeV_Mass_High_jpsi', formatter=lambda x: '{:0.5f}'.format(x)))
    #print(ToCArray(csHighMassPsi2sIcemPt, ctype='double', name='Ramona_136TeV_Mass_High_psi2S', formatter=lambda x: '{:0.5f}'.format(x)))

    #deltaCsHighMassJpsiIcemPt = np.abs(csHighMassJpsiIcemPt - csCentralJpsiIcemPt)
    #deltaCsHighMassPsi2sIcemPt = np.abs(csHighMassPsi2sIcemPt - csCentralPsi2sIcemPt)
    #relErrCsHighMassPsi2sOverJpsiIcemPt = (deltaCsHighMassJpsiIcemPt / csCentralJpsiIcemPt) - (deltaCsHighMassPsi2sIcemPt / csCentralPsi2sIcemPt)

    # low mass
    dfCsJpsiIcemPt['lowMassSum'] = dfCsJpsiIcemPt[lowMassCols].sum(axis=1)
    dfCsPsi2sIcemPt['lowMassSum'] = dfCsPsi2sIcemPt[lowMassCols].sum(axis=1)
    #dfMergedCsJpsi = dfCsJpsiIcemPt.groupby('ptGroup').agg({'lowMassSum': 'sum'}).reset_index()
    #dfMergedCsPsi2s = dfCsPsi2sIcemPt.groupby('ptGroup').agg({'lowMassSum': 'sum'}).reset_index()
    #csLowMassJpsiIcemPt = (dfMergedCsJpsi['lowMassSum'].to_numpy()) * 1e-3
    #csLowMassPsi2sIcemPt = (dfMergedCsPsi2s['lowMassSum'].to_numpy()) * 1e-3
    #csLowMassPsi2sOverJpsiIcemPt = csLowMassPsi2sIcemPt / csLowMassJpsiIcemPt
    #print(ToCArray(csLowMassJpsiIcemPt, ctype='double', name='Ramona_136TeV_Mass_Low_jpsi', formatter=lambda x: '{:0.5f}'.format(x)))
    #print(ToCArray(csLowMassPsi2sIcemPt, ctype='double', name='Ramona_136TeV_Mass_Low_psi2S', formatter=lambda x: '{:0.5f}'.format(x)))

    #deltaCsLowMassJpsiIcemPt = np.abs(csLowMassJpsiIcemPt - csCentralJpsiIcemPt)
    #deltaCsLowMassPsi2sIcemPt = np.abs(csLowMassPsi2sIcemPt - csCentralPsi2sIcemPt)
    #relErrCsLowMassPsi2sOverJpsiIcemPt = (deltaCsLowMassJpsiIcemPt / csCentralJpsiIcemPt) - (deltaCsLowMassPsi2sIcemPt / csCentralPsi2sIcemPt)

    dfCsJpsiIcemPt.to_csv('jpsi_icem_aggregated.csv', index=False)
    dfCsPsi2sIcemPt.to_csv('psi2s_icem_aggregated.csv', index=False)

    exit()


    # J/psi
    dfCsJpsiFonll = pd.read_csv('FONLL/cs_jpsi_fonll_pt.txt', sep=' ')
    csCentralJpsiFonllPt = (dfCsJpsiFonll["central"].to_numpy()) * 2e-6
    csMinJpsiFonllPt = (dfCsJpsiFonll["min"].to_numpy()) * 2e-6
    csMaxJpsiFonllPt = (dfCsJpsiFonll["max"].to_numpy()) * 2e-6
    csMinScJpsiFonllPt = (dfCsJpsiFonll["min_sc"].to_numpy()) * 2e-6
    csMaxScJpsiFonllPt = (dfCsJpsiFonll["max_sc"].to_numpy()) * 2e-6
    csMinMassJpsiFonllPt = (dfCsJpsiFonll["min_mass"].to_numpy()) * 2e-6
    csMaxMassJpsiFonllPt = (dfCsJpsiFonll["max_mass"].to_numpy()) * 2e-6
    csMinPdfJpsiFonllPt = (dfCsJpsiFonll["min_pdf"].to_numpy()) * 2e-6
    csMaxPdfJpsiFonllPt = (dfCsJpsiFonll["max_pdf"].to_numpy()) * 2e-6

    # Psi(2S)
    dfCsPsi2sFonll = pd.read_csv('FONLL/cs_psi2s_fonll_pt.txt', sep=' ')
    csCentralPsi2sFonllPt = (dfCsPsi2sFonll["central"].to_numpy()) * 2e-6
    csMinPsi2sFonllPt = (dfCsPsi2sFonll["min"].to_numpy()) * 2e-6
    csMaxPsi2sFonllPt = (dfCsPsi2sFonll["max"].to_numpy()) * 2e-6
    csMinScPsi2sFonllPt = (dfCsPsi2sFonll["min_sc"].to_numpy()) * 2e-6
    csMaxScPsi2sFonllPt = (dfCsPsi2sFonll["max_sc"].to_numpy()) * 2e-6
    csMinMassPsi2sFonllPt = (dfCsPsi2sFonll["min_mass"].to_numpy()) * 2e-6
    csMaxMassPsi2sFonllPt = (dfCsPsi2sFonll["max_mass"].to_numpy()) * 2e-6
    csMinPdfPsi2sFonllPt = (dfCsPsi2sFonll["min_pdf"].to_numpy()) * 2e-6
    csMaxPdfPsi2sFonllPt = (dfCsPsi2sFonll["max_pdf"].to_numpy()) * 2e-6

    deltaCsMinScJpsiFonllPt = np.abs(csMinScJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMinScPsi2sFonllPt = np.abs(csMinScPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMinScPsi2sOverJpsiFonllPt = (deltaCsMinScJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMinScPsi2sFonllPt / csCentralPsi2sFonllPt)

    deltaCsMaxScJpsiFonllPt = np.abs(csMaxScJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMaxScPsi2sFonllPt = np.abs(csMaxScPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMaxScPsi2sOverJpsiFonllPt = (deltaCsMaxScJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMaxScPsi2sFonllPt / csCentralPsi2sFonllPt)

    deltaCsMinMassJpsiFonllPt = np.abs(csMinMassJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMinMassPsi2sFonllPt = np.abs(csMinMassPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMinMassPsi2sOverJpsiFonllPt = (deltaCsMinMassJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMinMassPsi2sFonllPt / csCentralPsi2sFonllPt)

    deltaCsMaxMassJpsiFonllPt = np.abs(csMaxMassJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMaxMassPsi2sFonllPt = np.abs(csMaxMassPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMaxMassPsi2sOverJpsiFonllPt = (deltaCsMaxMassJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMaxMassPsi2sFonllPt / csCentralPsi2sFonllPt)

    deltaCsMinPdfJpsiFonllPt = np.abs(csMinPdfJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMinPdfPsi2sFonllPt = np.abs(csMinPdfPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMinPdfPsi2sOverJpsiFonllPt = (deltaCsMinPdfJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMinPdfPsi2sFonllPt / csCentralPsi2sFonllPt)

    deltaCsMaxPdfJpsiFonllPt = np.abs(csMaxPdfJpsiFonllPt - csCentralJpsiFonllPt)
    deltaCsMaxPdfPsi2sFonllPt = np.abs(csMaxPdfPsi2sFonllPt - csCentralPsi2sFonllPt)
    relErrCsMaxPdfPsi2sOverJpsiFonllPt = (deltaCsMaxPdfJpsiFonllPt / csCentralJpsiFonllPt) - (deltaCsMaxPdfPsi2sFonllPt / csCentralPsi2sFonllPt)


    # Ratio ICEM + FONLL
    # Central value
    csCentralJpsiIcemFonllPt = csCentralJpsiIcemPt + csCentralJpsiFonllPt
    csCentralPsi2sIcemFonllPt = csCentralPsi2sIcemPt + csCentralPsi2sFonllPt
    csCentralPsi2sOverJpsiIcemFonllPt = csCentralPsi2sIcemFonllPt / csCentralJpsiIcemFonllPt

    # Full Error
    #errFullMinPsi2sOverJpsiPt = csCentralPsi2sOverJpsiIcemFonllPt * np.sqrt(
        #relErrCsHalfRenormPsi2sOverJpsiIcemPt * relErrCsHalfRenormPsi2sOverJpsiIcemPt +
        #relErrCsLowMassPsi2sOverJpsiIcemPt * relErrCsLowMassPsi2sOverJpsiIcemPt +
        #relErrCsMinScPsi2sOverJpsiFonllPt * relErrCsMinScPsi2sOverJpsiFonllPt +
        #relErrCsMinMassPsi2sOverJpsiFonllPt * relErrCsMinMassPsi2sOverJpsiFonllPt +
        #relErrCsMinPdfPsi2sOverJpsiFonllPt * relErrCsMinPdfPsi2sOverJpsiFonllPt
    #)

    #errFullMaxPsi2sOverJpsiPt = csCentralPsi2sOverJpsiIcemFonllPt * np.sqrt(
        #relErrCsDoubleRenormPsi2sOverJpsiIcemPt * relErrCsDoubleRenormPsi2sOverJpsiIcemPt +
        #relErrCsHighMassPsi2sOverJpsiIcemPt * relErrCsHighMassPsi2sOverJpsiIcemPt +
        #relErrCsMaxScPsi2sOverJpsiFonllPt * relErrCsMaxScPsi2sOverJpsiFonllPt +
        #relErrCsMaxMassPsi2sOverJpsiFonllPt * relErrCsMaxMassPsi2sOverJpsiFonllPt +
        #relErrCsMaxPdfPsi2sOverJpsiFonllPt * relErrCsMaxPdfPsi2sOverJpsiFonllPt
    #)

    errFullMinPsi2sOverJpsiIcemFonllPt = csCentralPsi2sOverJpsiIcemFonllPt * np.sqrt(
        relErrCsDoubleRenormPsi2sOverJpsiIcemPt * relErrCsDoubleRenormPsi2sOverJpsiIcemPt +
        relErrCsHighMassPsi2sOverJpsiIcemPt * relErrCsHighMassPsi2sOverJpsiIcemPt +
        relErrCsMinScPsi2sOverJpsiFonllPt * relErrCsMinScPsi2sOverJpsiFonllPt +
        relErrCsMinMassPsi2sOverJpsiFonllPt * relErrCsMinMassPsi2sOverJpsiFonllPt +
        relErrCsMinPdfPsi2sOverJpsiFonllPt * relErrCsMinPdfPsi2sOverJpsiFonllPt
    )

    errFullMaxPsi2sOverJpsiIcemFonllPt = csCentralPsi2sOverJpsiIcemFonllPt * np.sqrt(
        relErrCsHalfRenormPsi2sOverJpsiIcemPt * relErrCsHalfRenormPsi2sOverJpsiIcemPt +
        relErrCsLowMassPsi2sOverJpsiIcemPt * relErrCsLowMassPsi2sOverJpsiIcemPt +
        relErrCsMaxScPsi2sOverJpsiFonllPt * relErrCsMaxScPsi2sOverJpsiFonllPt +
        relErrCsMaxMassPsi2sOverJpsiFonllPt * relErrCsMaxMassPsi2sOverJpsiFonllPt +
        relErrCsMaxPdfPsi2sOverJpsiFonllPt * relErrCsMaxPdfPsi2sOverJpsiFonllPt
    )


    # CS + CO NLO
    # Come propagare al ratio FONLL ?
    #dfCsPsi2sOverJpsiCsCoNloPt = pd.read_csv('/Users/lucamicheletti/GITHUB/dq_fit_library/analysis/theory_predictions/cs_co_nlo_pt.txt', sep=' ')
    #ptMinCsCoNlo = dfCsPsi2sOverJpsiCsCoNloPt["x_min"].to_numpy()
    #ptMaxCsCoNlo = dfCsPsi2sOverJpsiCsCoNloPt["x_max"].to_numpy()
    #csPsi2sOverJpsiCsCoNloPt = (dfCsPsi2sOverJpsiCsCoNloPt["val"].to_numpy())
    #csMinPsi2sOverJpsiCsCoNloPt = (dfCsPsi2sOverJpsiCsCoNloPt["val_min"].to_numpy())
    #csMaxPsi2sOverJpsiCsCoNloPt = (dfCsPsi2sOverJpsiCsCoNloPt["val_max"].to_numpy())

    #deltaCsMinPsi2sOverJpsiCsCoNloPt = np.abs(csMinPsi2sOverJpsiCsCoNloPt - csPsi2sOverJpsiCsCoNloPt)
    #relErrCsMinPsi2sOverJpsiCsCoNloPt = (deltaCsMinPsi2sOverJpsiCsCoNloPt / csPsi2sOverJpsiCsCoNloPt)

    #deltaCsMaxPsi2sOverJpsiCsCoNloPt = np.abs(csMaxPsi2sOverJpsiCsCoNloPt - csPsi2sOverJpsiCsCoNloPt)
    #relErrCsMaxPsi2sOverJpsiCsCoNloPt = (deltaCsMaxPsi2sOverJpsiCsCoNloPt / csPsi2sOverJpsiCsCoNloPt)

    #errFullMinPsi2sOverJpsiCsCoNloFonllPt = csCentralPsi2sOverJpsiCsCoNloFonllPt * np.sqrt(
        #relErrCsDoubleRenormPsi2sOverJpsiIcemPt * relErrCsDoubleRenormPsi2sOverJpsiIcemPt +
        #relErrCsHighMassPsi2sOverJpsiIcemPt * relErrCsHighMassPsi2sOverJpsiIcemPt +
        #relErrCsMinScPsi2sOverJpsiFonllPt * relErrCsMinScPsi2sOverJpsiFonllPt +
        #relErrCsMinMassPsi2sOverJpsiFonllPt * relErrCsMinMassPsi2sOverJpsiFonllPt +
        #relErrCsMinPdfPsi2sOverJpsiFonllPt * relErrCsMinPdfPsi2sOverJpsiFonllPt
    #)

    #errFullMaxPsi2sOverJpsiCsCoNloFonllPt = csCentralPsi2sOverJpsiCsCoNloFonllPt * np.sqrt(
        #relErrCsHalfRenormPsi2sOverJpsiIcemPt * relErrCsHalfRenormPsi2sOverJpsiIcemPt +
        #relErrCsLowMassPsi2sOverJpsiIcemPt * relErrCsLowMassPsi2sOverJpsiIcemPt +
        #relErrCsMaxScPsi2sOverJpsiFonllPt * relErrCsMaxScPsi2sOverJpsiFonllPt +
        #relErrCsMaxMassPsi2sOverJpsiFonllPt * relErrCsMaxMassPsi2sOverJpsiFonllPt +
        #relErrCsMaxPdfPsi2sOverJpsiFonllPt * relErrCsMaxPdfPsi2sOverJpsiFonllPt
    #)

    #print(deltaCsMinPdfJpsiFonllPt)
    #print(deltaCsMinPdfJpsiFonllPt[3:])


    graErrCsRatioPsi2sOverJpsiPt = ROOT.TGraphAsymmErrors(len(ptCentrArr), ptCentrArr, csCentralPsi2sOverJpsiIcemFonllPt, ptWidthArr, ptWidthArr, errFullMinPsi2sOverJpsiIcemFonllPt, errFullMaxPsi2sOverJpsiIcemFonllPt)
    graErrCsRatioPsi2sOverJpsiPt.SetFillColorAlpha(ROOT.kOrange+7, 0.8)
    graErrCsRatioPsi2sOverJpsiPt.SetLineColor(ROOT.kOrange+7)

    graCsRatioPsi2sOverJpsiPt = ROOT.TGraphErrors(len(ptCentrArr), ptCentrArr, csCentralPsi2sOverJpsiIcemFonllPt, ptWidthArr, 0)
    SetGraStat(graCsRatioPsi2sOverJpsiPt, 20, ROOT.kBlack)

    canvasCsPsi2sOverJpsiPt = TCanvas("canvasCsPsi2sOverJpsiPt", "canvasCsPsi2sOverJpsiPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiPt  = TH2F("histGridCsPsi2sOverJpsiPt", "", 100, 0., 20., 100, 0.1, 0.7)
    histGridCsPsi2sOverJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    histGridCsPsi2sOverJpsiPt.GetYaxis().SetTitle("(d^{2}#sigma_{#psi(2S)} / d#it{p}_{T} d#it{y}) / (d^{2}#sigma_{J/#psi} / d#it{p}_{T} d#it{y})")
    histGridCsPsi2sOverJpsiPt.Draw()
    graErrCsRatioPsi2sOverJpsiPt.Draw("E2 SAME")
    graCsRatioPsi2sOverJpsiPt.Draw("EP SAME")
    canvasCsPsi2sOverJpsiPt.Update()
    canvasCsPsi2sOverJpsiPt.SaveAs("psi2s_over_jpsi_ICEM_FONLL_vs_pt.pdf")

    with open("psi2s_over_jpsi_ICEM_FONLL_vs_pt.txt", 'w') as fOut:
        fOut.write("x_min x_max val err_min err_max \n")
        for iPt, (ptMin, ptMax) in enumerate(zip(ptMinArr, ptMaxArr)):
            fOut.write(f"{ptMin} {ptMax} {csCentralPsi2sOverJpsiIcemFonllPt[iPt]} {errFullMinPsi2sOverJpsiIcemFonllPt[iPt]} {errFullMaxPsi2sOverJpsiIcemFonllPt[iPt]} \n")

    #######################################
    # Rapidity dependence
    #######################################
    csCentralJpsiIcemY = np.zeros(len(yMinArr))
    csCentralPsi2sIcemY = np.zeros(len(yMinArr))
    csDoubleRenormJpsiIcemY = np.zeros(len(yMinArr))
    csDoubleRenormPsi2sIcemY = np.zeros(len(yMinArr))
    csHalfRenormJpsiIcemY = np.zeros(len(yMinArr))
    csHalfRenormPsi2sIcemY = np.zeros(len(yMinArr))
    csHighMassJpsiIcemY = np.zeros(len(yMinArr))
    csHighMassPsi2sIcemY = np.zeros(len(yMinArr))
    csLowMassJpsiIcemY = np.zeros(len(yMinArr))
    csLowMassPsi2sIcemY = np.zeros(len(yMinArr))

    for iY, (centralCol, doubleRenormCols, halfRenormCols, highMassCols, lowMassCols) in enumerate(zip(centralCols, doubleRenormCols, halfRenormCols, highMassCols, lowMassCols)):
        csCentralJpsiIcemY[iY] = (dfCsJpsiIcemY[centralCol].sum()) * 1e-3
        csCentralPsi2sIcemY[iY] = (dfCsPsi2sIcemY[centralCol].sum()) * 1e-3
        csDoubleRenormJpsiIcemY[iY] = (dfCsJpsiIcemY[doubleRenormCols].sum()) * 1e-3
        csDoubleRenormPsi2sIcemY[iY] = (dfCsPsi2sIcemY[doubleRenormCols].sum()) * 1e-3
        csHalfRenormJpsiIcemY[iY] = (dfCsJpsiIcemY[halfRenormCols].sum()) * 1e-3
        csHalfRenormPsi2sIcemY[iY] = (dfCsPsi2sIcemY[halfRenormCols].sum()) * 1e-3
        csHighMassJpsiIcemY[iY] = (dfCsJpsiIcemY[highMassCols].sum()) * 1e-3
        csHighMassPsi2sIcemY[iY] = (dfCsPsi2sIcemY[highMassCols].sum()) * 1e-3
        csLowMassJpsiIcemY[iY] = (dfCsJpsiIcemY[lowMassCols].sum()) * 1e-3
        csLowMassPsi2sIcemY[iY] = (dfCsPsi2sIcemY[lowMassCols].sum()) * 1e-3


    deltaCsDoubleRenormJpsiIcemY = np.abs(csDoubleRenormJpsiIcemY - csCentralJpsiIcemY)
    deltaCsDoubleRenormPsi2sIcemY = np.abs(csDoubleRenormPsi2sIcemY - csCentralPsi2sIcemY)
    relErrCsDoubleRenormPsi2sOverJpsiIcemY = (deltaCsDoubleRenormJpsiIcemY / csCentralJpsiIcemY) - (deltaCsDoubleRenormPsi2sIcemY / csCentralPsi2sIcemY)

    deltaCsHalfRenormJpsiIcemY = np.abs(csHalfRenormJpsiIcemY - csCentralJpsiIcemY)
    deltaCsHalfRenormPsi2sIcemY = np.abs(csHalfRenormPsi2sIcemY - csCentralPsi2sIcemY)
    relErrCsHalfRenormPsi2sOverJpsiIcemY = (deltaCsHalfRenormJpsiIcemY / csCentralJpsiIcemY) - (deltaCsHalfRenormPsi2sIcemY / csCentralPsi2sIcemY)

    deltaCsHighMassJpsiIcemY = np.abs(csHighMassJpsiIcemY - csCentralJpsiIcemY)
    deltaCsHighMassPsi2sIcemY = np.abs(csHighMassPsi2sIcemY - csCentralPsi2sIcemY)
    relErrCsHighMassPsi2sOverJpsiIcemY = (deltaCsHighMassJpsiIcemY / csCentralJpsiIcemY) - (deltaCsHighMassPsi2sIcemY / csCentralPsi2sIcemY)

    deltaCsLowMassJpsiIcemY = np.abs(csLowMassJpsiIcemY - csCentralJpsiIcemY)
    deltaCsLowMassPsi2sIcemY = np.abs(csLowMassPsi2sIcemY - csCentralPsi2sIcemY)
    relErrCsLowMassPsi2sOverJpsiIcemY = (deltaCsLowMassJpsiIcemY / csCentralJpsiIcemY) - (deltaCsLowMassPsi2sIcemY / csCentralPsi2sIcemY)


    # J/psi
    dfCsJpsiFonll = pd.read_csv('FONLL/cs_jpsi_fonll_y.txt', sep=' ')
    csCentralJpsiFonllY = (dfCsJpsiFonll["central"].to_numpy()) * 2e-6
    csMinScJpsiFonllY = (dfCsJpsiFonll["min_sc"].to_numpy()) * 2e-6
    csMaxScJpsiFonllY = (dfCsJpsiFonll["max_sc"].to_numpy()) * 2e-6
    csMinMassJpsiFonllY = (dfCsJpsiFonll["min_mass"].to_numpy()) * 2e-6
    csMaxMassJpsiFonllY = (dfCsJpsiFonll["max_mass"].to_numpy()) * 2e-6
    csMinPdfJpsiFonllY = (dfCsJpsiFonll["min_pdf"].to_numpy()) * 2e-6
    csMaxPdfJpsiFonllY = (dfCsJpsiFonll["max_pdf"].to_numpy()) * 2e-6

    # Psi(2S)
    dfCsPsi2sFonll = pd.read_csv('FONLL/cs_psi2s_fonll_y.txt', sep=' ')
    csCentralPsi2sFonllY = (dfCsPsi2sFonll["central"].to_numpy()) * 2e-6
    csMinScPsi2sFonllY = (dfCsPsi2sFonll["min_sc"].to_numpy()) * 2e-6
    csMaxScPsi2sFonllY = (dfCsPsi2sFonll["max_sc"].to_numpy()) * 2e-6
    csMinMassPsi2sFonllY = (dfCsPsi2sFonll["min_mass"].to_numpy()) * 2e-6
    csMaxMassPsi2sFonllY = (dfCsPsi2sFonll["max_mass"].to_numpy()) * 2e-6
    csMinPdfPsi2sFonllY = (dfCsPsi2sFonll["min_pdf"].to_numpy()) * 2e-6
    csMaxPdfPsi2sFonllY = (dfCsPsi2sFonll["max_pdf"].to_numpy()) * 2e-6

    deltaCsMinScJpsiFonllY = np.abs(csMinScJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMinScPsi2sFonllY = np.abs(csMinScPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMinScPsi2sOverJpsiFonllY = (deltaCsMinScJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMinScPsi2sFonllY / csCentralPsi2sFonllY)

    deltaCsMaxScJpsiFonllY = np.abs(csMaxScJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMaxScPsi2sFonllY = np.abs(csMaxScPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMaxScPsi2sOverJpsiFonllY = (deltaCsMaxScJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMaxScPsi2sFonllY / csCentralPsi2sFonllY)


    deltaCsMinMassJpsiFonllY = np.abs(csMinMassJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMinMassPsi2sFonllY = np.abs(csMinMassPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMinMassPsi2sOverJpsiFonllY = (deltaCsMinMassJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMinMassPsi2sFonllY / csCentralPsi2sFonllY)

    deltaCsMaxMassJpsiFonllY = np.abs(csMaxMassJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMaxMassPsi2sFonllY = np.abs(csMaxMassPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMaxMassPsi2sOverJpsiFonllY = (deltaCsMaxMassJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMaxMassPsi2sFonllY / csCentralPsi2sFonllY)

    deltaCsMinPdfJpsiFonllY = np.abs(csMinPdfJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMinPdfPsi2sFonllY = np.abs(csMinPdfPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMinPdfPsi2sOverJpsiFonllY = (deltaCsMinPdfJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMinPdfPsi2sFonllY / csCentralPsi2sFonllY)

    deltaCsMaxPdfJpsiFonllY = np.abs(csMaxPdfJpsiFonllY - csCentralJpsiFonllY)
    deltaCsMaxPdfPsi2sFonllY = np.abs(csMaxPdfPsi2sFonllY - csCentralPsi2sFonllY)
    relErrCsMaxPdfPsi2sOverJpsiFonllY = (deltaCsMaxPdfJpsiFonllY / csCentralJpsiFonllY) - (deltaCsMaxPdfPsi2sFonllY / csCentralPsi2sFonllY)


    # Ratio ICEM + FONLL
    # Central value
    csCentralJpsiIcemFonllY = csCentralJpsiIcemY + csCentralJpsiFonllY
    csCentralPsi2sIcemFonllY = csCentralPsi2sIcemY + csCentralPsi2sFonllY
    csCentralPsi2sOverJpsiIcemFonllY = csCentralPsi2sIcemFonllY / csCentralJpsiIcemFonllY


    # Full Error
    #errFullMinPsi2sOverJpsiIcemFonllY = csCentralPsi2sOverJpsiIcemFonllY * np.sqrt(
        #relErrCsHalfRenormPsi2sOverJpsiIcemY * relErrCsHalfRenormPsi2sOverJpsiIcemY +
        #relErrCsLowMassPsi2sOverJpsiIcemY * relErrCsLowMassPsi2sOverJpsiIcemY +
        #relErrCsMinScPsi2sOverJpsiFonllY * relErrCsMinScPsi2sOverJpsiFonllY +
        #relErrCsMinMassPsi2sOverJpsiFonllY * relErrCsMinMassPsi2sOverJpsiFonllY +
        #relErrCsMinPdfPsi2sOverJpsiFonllY * relErrCsMinPdfPsi2sOverJpsiFonllY
    #)

    #errFullMaxPsi2sOverJpsiIcemFonllY = csCentralPsi2sOverJpsiIcemFonllY * np.sqrt(
        #relErrCsDoubleRenormPsi2sOverJpsiIcemY * relErrCsDoubleRenormPsi2sOverJpsiIcemY +
        #relErrCsHighMassPsi2sOverJpsiIcemY * relErrCsHighMassPsi2sOverJpsiIcemY +
        #relErrCsMaxScPsi2sOverJpsiFonllY * relErrCsMaxScPsi2sOverJpsiFonllY +
        #relErrCsMaxMassPsi2sOverJpsiFonllY * relErrCsMaxMassPsi2sOverJpsiFonllY +
        #relErrCsMaxPdfPsi2sOverJpsiFonllY * relErrCsMaxPdfPsi2sOverJpsiFonllY
    #)

    errFullMinPsi2sOverJpsiIcemFonllY = csCentralPsi2sOverJpsiIcemFonllY * np.sqrt(
        relErrCsDoubleRenormPsi2sOverJpsiIcemY * relErrCsDoubleRenormPsi2sOverJpsiIcemY +
        relErrCsHighMassPsi2sOverJpsiIcemY * relErrCsHighMassPsi2sOverJpsiIcemY +
        relErrCsMinScPsi2sOverJpsiFonllY * relErrCsMinScPsi2sOverJpsiFonllY +
        relErrCsMinMassPsi2sOverJpsiFonllY * relErrCsMinMassPsi2sOverJpsiFonllY +
        relErrCsMinPdfPsi2sOverJpsiFonllY * relErrCsMinPdfPsi2sOverJpsiFonllY
    )

    errFullMaxPsi2sOverJpsiIcemFonllY = csCentralPsi2sOverJpsiIcemFonllY * np.sqrt(
        relErrCsHalfRenormPsi2sOverJpsiIcemY * relErrCsHalfRenormPsi2sOverJpsiIcemY +
        relErrCsLowMassPsi2sOverJpsiIcemY * relErrCsLowMassPsi2sOverJpsiIcemY +
        relErrCsMaxScPsi2sOverJpsiFonllY * relErrCsMaxScPsi2sOverJpsiFonllY +
        relErrCsMaxMassPsi2sOverJpsiFonllY * relErrCsMaxMassPsi2sOverJpsiFonllY +
        relErrCsMaxPdfPsi2sOverJpsiFonllY * relErrCsMaxPdfPsi2sOverJpsiFonllY
    )

    graErrCsRatioPsi2sOverJpsiY = ROOT.TGraphAsymmErrors(len(yCentrArr), yCentrArr, csCentralPsi2sOverJpsiIcemFonllY, yWidthArr, yWidthArr, errFullMinPsi2sOverJpsiIcemFonllY, errFullMaxPsi2sOverJpsiIcemFonllY)
    graErrCsRatioPsi2sOverJpsiY.SetFillColorAlpha(ROOT.kOrange+7, 0.8)
    graErrCsRatioPsi2sOverJpsiY.SetLineColor(ROOT.kOrange+7)

    graCsRatioPsi2sOverJpsiY = ROOT.TGraphErrors(len(yCentrArr), yCentrArr, csCentralPsi2sOverJpsiIcemFonllY, yWidthArr, 0)
    SetGraStat(graCsRatioPsi2sOverJpsiY, 20, ROOT.kBlack)

    canvasCsPsi2sOverJpsiY = TCanvas("canvasCsPsi2sOverJpsiY", "canvasCsPsi2sOverJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiY  = TH2F("histGridCsPsi2sOverJpsiY", "", 100, 2.5, 4., 100, 0.1, 0.7)
    histGridCsPsi2sOverJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridCsPsi2sOverJpsiY.GetYaxis().SetTitle("(d#sigma_{#psi(2S)} / d#it{y}) / (d#sigma_{J/#psi} / d#it{y})")
    histGridCsPsi2sOverJpsiY.Draw()
    graErrCsRatioPsi2sOverJpsiY.Draw("E2 SAME")
    graCsRatioPsi2sOverJpsiY.Draw("EP SAME")
    canvasCsPsi2sOverJpsiY.Update()
    canvasCsPsi2sOverJpsiY.SaveAs("psi2s_over_jpsi_ICEM_FONLL_vs_y.pdf")


    with open("psi2s_over_jpsi_ICEM_FONLL_vs_y.txt", 'w') as fOut:
        fOut.write("x_min x_max val err_min err_max \n")
        for iY, (yMin, yMax) in enumerate(zip(yMinArr, yMaxArr)):
            fOut.write(f"{yMin} {yMax} {csCentralPsi2sOverJpsiIcemFonllY[iY]} {errFullMinPsi2sOverJpsiIcemFonllY[iY]} {errFullMaxPsi2sOverJpsiIcemFonllY[iY]} \n")




    input()

    exit()
    graStatCsCentralRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csCentralPsi2sOverJpsiIcemPt, ptWidthArr, 0)
    SetGraStat(graStatCsCentralRatioPsi2sOverJpsiPt, 20, ROOT.kBlack)

    
    graStatCsDoubleRenormRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csDoubleRenormPsi2sOverJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsDoubleRenormRatioPsi2sOverJpsiPt, 20, ROOT.kRed+1)

    
    graStatCsHalfRenormRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHalfRenormPsi2sOverJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsHalfRenormRatioPsi2sOverJpsiPt, 20, ROOT.kRed+1)
    graStatCsHalfRenormRatioPsi2sOverJpsiPt.SetLineStyle(ROOT.kDashed)

    
    graStatCsHighMassRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHighMassPsi2sOverJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsHighMassRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)

    graStatCsLowMassRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csLowMassPsi2sOverJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsLowMassRatioPsi2sOverJpsiPt, 20, ROOT.kAzure+2)
    graStatCsLowMassRatioPsi2sOverJpsiPt.SetLineStyle(ROOT.kDashed)

    
    graSystCsCentralRatioPsi2sOverJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csCentralPsi2sOverJpsiIcemPt, ptWidthArr, relErrFullPsi2sOverJpsiPt)
    graSystCsCentralRatioPsi2sOverJpsiPt.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
    graSystCsCentralRatioPsi2sOverJpsiPt.SetLineColor(ROOT.kOrange+7)

    canvasCsPsi2sOverJpsiPt = TCanvas("canvasCsPsi2sOverJpsiPt", "canvasCsPsi2sOverJpsiPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiPt  = TH2F("histGridCsPsi2sOverJpsiPt", "", 100, 0., 20., 100, 0.1, 0.7)
    histGridCsPsi2sOverJpsiPt.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    histGridCsPsi2sOverJpsiPt.GetYaxis().SetTitle("(d^{2}#sigma_{#psi(2S)} / d#it{p}_{T} d#it{y}) / (d^{2}#sigma_{J/#psi} / d#it{p}_{T} d#it{y})")
    histGridCsPsi2sOverJpsiPt.Draw()
    graSystCsCentralRatioPsi2sOverJpsiPt.Draw("E2 SAME")
    graStatCsCentralRatioPsi2sOverJpsiPt.Draw("P SAME")
    graStatCsDoubleRenormRatioPsi2sOverJpsiPt.Draw("SAME")
    graStatCsHalfRenormRatioPsi2sOverJpsiPt.Draw("SAME")
    graStatCsHighMassRatioPsi2sOverJpsiPt.Draw("SAME")
    graStatCsLowMassRatioPsi2sOverJpsiPt.Draw("SAME")
    canvasCsPsi2sOverJpsiPt.Update()

    # Psi(2S) vs Pt
    graStatCsCentralPsi2sPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csCentralPsi2sIcemPt, ptWidthArr, 0)
    SetGraStat(graStatCsCentralPsi2sPt, 20, ROOT.kBlack)

    graStatCsDoubleRenormPsi2sPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csDoubleRenormPsi2sIcemPt, 0, 0)
    SetGraStat(graStatCsDoubleRenormPsi2sPt, 20, ROOT.kRed+1)

    graStatCsHalfRenormPsi2sPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHalfRenormPsi2sIcemPt, 0, 0)
    SetGraStat(graStatCsHalfRenormPsi2sPt, 20, ROOT.kRed+1)

    graStatCsHighMassPsi2sPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHighMassPsi2sIcemPt, 0, 0)
    SetGraStat(graStatCsHighMassPsi2sPt, 20, ROOT.kAzure+2)

    graStatCsLowMassPsi2sPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csLowMassPsi2sIcemPt, 0, 0)
    SetGraStat(graStatCsLowMassPsi2sPt, 20, ROOT.kAzure+2)

    # J/psi vs Pt
    graStatCsCentralJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csCentralJpsiIcemPt, ptWidthArr, 0)
    SetGraStat(graStatCsCentralJpsiPt, 20, ROOT.kBlack)
    graStatCsCentralJpsiPt.SetLineStyle(ROOT.kDashed)

    graStatCsDoubleRenormJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csDoubleRenormJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsDoubleRenormJpsiPt, 20, ROOT.kRed+1)
    graStatCsDoubleRenormJpsiPt.SetLineStyle(ROOT.kDashed)

    graStatCsHalfRenormJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHalfRenormJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsHalfRenormJpsiPt, 20, ROOT.kRed+1)
    graStatCsHalfRenormJpsiPt.SetLineStyle(ROOT.kDashed)

    graStatCsHighMassJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csHighMassJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsHighMassJpsiPt, 20, ROOT.kAzure+2)
    graStatCsHighMassJpsiPt.SetLineStyle(ROOT.kDashed)

    graStatCsLowMassJpsiPt = TGraphErrors(len(ptCentrArr), ptCentrArr, csLowMassJpsiIcemPt, 0, 0)
    SetGraStat(graStatCsLowMassJpsiPt, 20, ROOT.kAzure+2)
    graStatCsLowMassJpsiPt.SetLineStyle(ROOT.kDashed)

    canvasCsPt = TCanvas("canvasCsPt", "canvasCsPt", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPt  = TH2F("histGridCsPt", "", 100, 0., 20., 100, 1, 1e4)
    histGridCsPt.GetXaxis().SetTitle("#it{y}")
    histGridCsPt.GetYaxis().SetTitle("d#sigma / d#it{p}_{T}")
    histGridCsPt.Draw()
    graStatCsCentralPsi2sPt.Draw("SAME")
    graStatCsDoubleRenormPsi2sPt.Draw("SAME")
    graStatCsHalfRenormPsi2sPt.Draw("SAME")
    graStatCsHighMassPsi2sPt.Draw("SAME")
    graStatCsLowMassPsi2sPt.Draw("SAME")
    graStatCsCentralJpsiPt.Draw("SAME")
    graStatCsDoubleRenormJpsiPt.Draw("SAME")
    graStatCsHalfRenormJpsiPt.Draw("SAME")
    graStatCsHighMassJpsiPt.Draw("SAME")
    graStatCsLowMassJpsiPt.Draw("SAME")
    canvasCsPt.Update()



    #######################################
    # Rapidity dependence
    #######################################
    csCentralJpsiY = np.zeros(len(yMinArr))
    csCentralPsi2sY = np.zeros(len(yMinArr))
    csDoubleRenormJpsiY = np.zeros(len(yMinArr))
    csDoubleRenormPsi2sY = np.zeros(len(yMinArr))
    csHalfRenormJpsiY = np.zeros(len(yMinArr))
    csHalfRenormPsi2sY = np.zeros(len(yMinArr))
    csHighMassJpsiY = np.zeros(len(yMinArr))
    csHighMassPsi2sY = np.zeros(len(yMinArr))
    csLowMassJpsiY = np.zeros(len(yMinArr))
    csLowMassPsi2sY = np.zeros(len(yMinArr))

    for iY, (centralCol, doubleRenormCols, halfRenormCols, highMassCols, lowMassCols) in enumerate(zip(centralCols, doubleRenormCols, halfRenormCols, highMassCols, lowMassCols)):
        csCentralJpsiY[iY]= dfCsJpsiIcemY[centralCol].sum()
        csCentralPsi2sY[iY]= dfCsPsi2sIcemY[centralCol].sum()
        csDoubleRenormJpsiY[iY]= dfCsJpsiIcemY[doubleRenormCols].sum()
        csDoubleRenormPsi2sY[iY]= dfCsPsi2sIcemY[doubleRenormCols].sum()
        csHalfRenormJpsiY[iY]= dfCsJpsiIcemY[halfRenormCols].sum()
        csHalfRenormPsi2sY[iY]= dfCsPsi2sIcemY[halfRenormCols].sum()
        csHighMassJpsiY[iY]= dfCsJpsiIcemY[highMassCols].sum()
        csHighMassPsi2sY[iY]= dfCsPsi2sIcemY[highMassCols].sum()
        csLowMassJpsiY[iY]= dfCsJpsiIcemY[lowMassCols].sum()
        csLowMassPsi2sY[iY]= dfCsPsi2sIcemY[lowMassCols].sum()

    csCentralPsi2sOverJpsiY = csCentralPsi2sY / csCentralJpsiY
    csDoubleRenormPsi2sOverJpsiY = csDoubleRenormPsi2sY / csDoubleRenormJpsiY
    csHalfRenormPsi2sOverJpsiY = csHalfRenormPsi2sY / csHalfRenormJpsiY
    csHighMassPsi2sOverJpsiY = csHighMassPsi2sY / csHighMassJpsiY
    csLowMassPsi2sOverJpsiY = csLowMassPsi2sY / csLowMassJpsiY


    deltaCsDoubleRenormJpsiY = np.abs(csDoubleRenormJpsiY - csCentralJpsiY)
    deltaCsDoubleRenormPsi2sY = np.abs(csDoubleRenormPsi2sY - csCentralPsi2sY)
    relErrCsDoubleRenormPsi2sOverJpsiY = (deltaCsDoubleRenormJpsiY / csCentralJpsiY) - (deltaCsDoubleRenormPsi2sY / csCentralPsi2sY)

    deltaCsHalfRenormJpsiY = np.abs(csHalfRenormJpsiY - csCentralJpsiY)
    deltaCsHalfRenormPsi2sY = np.abs(csHalfRenormPsi2sY - csCentralPsi2sY)
    relErrCsHalfRenormPsi2sOverJpsiY = (deltaCsHalfRenormJpsiY / csCentralJpsiY) - (deltaCsHalfRenormPsi2sY / csCentralPsi2sY)

    deltaCsHighMassJpsiY = np.abs(csHighMassJpsiY - csCentralJpsiY)
    deltaCsHighMassPsi2sY = np.abs(csHighMassPsi2sY - csCentralPsi2sY)
    relErrCsHighMassPsi2sOverJpsiY = (deltaCsHighMassJpsiY / csCentralJpsiY) - (deltaCsHighMassPsi2sY / csCentralPsi2sY)

    deltaCsLowMassJpsiY = np.abs(csLowMassJpsiY - csCentralJpsiY)
    deltaCsLowMassPsi2sY = np.abs(csLowMassPsi2sY - csCentralPsi2sY)
    relErrCsLowMassPsi2sOverJpsiY = (deltaCsLowMassJpsiY / csCentralJpsiY) - (deltaCsLowMassPsi2sY / csCentralPsi2sY)

    # Full Error
    relErrFullPsi2sOverJpsiY = csCentralPsi2sOverJpsiY * np.sqrt(
        relErrCsDoubleRenormPsi2sOverJpsiY * relErrCsDoubleRenormPsi2sOverJpsiY +
        relErrCsHalfRenormPsi2sOverJpsiY * relErrCsHalfRenormPsi2sOverJpsiY +
        relErrCsHighMassPsi2sOverJpsiY * relErrCsHighMassPsi2sOverJpsiY +
        relErrCsLowMassPsi2sOverJpsiY * relErrCsLowMassPsi2sOverJpsiY
    )

    graStatCsCentralRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csCentralPsi2sOverJpsiY, yWidthArr, 0)
    SetGraStat(graStatCsCentralRatioPsi2sOverJpsiY, 20, ROOT.kBlack)

    graStatCsDoubleRenormRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csDoubleRenormPsi2sOverJpsiY, 0, 0)
    SetGraStat(graStatCsDoubleRenormRatioPsi2sOverJpsiY, 20, ROOT.kRed+1)

    graStatCsHalfRenormRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csHalfRenormPsi2sOverJpsiY, 0, 0)
    SetGraStat(graStatCsHalfRenormRatioPsi2sOverJpsiY, 20, ROOT.kRed+1)
    graStatCsHalfRenormRatioPsi2sOverJpsiY.SetLineStyle(ROOT.kDashed)

    graStatCsHighMassRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csHighMassPsi2sOverJpsiY, 0, 0)
    SetGraStat(graStatCsHighMassRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)

    graStatCsLowMassRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csLowMassPsi2sOverJpsiY, 0, 0)
    SetGraStat(graStatCsLowMassRatioPsi2sOverJpsiY, 20, ROOT.kAzure+2)
    graStatCsLowMassRatioPsi2sOverJpsiY.SetLineStyle(ROOT.kDashed)

    graSystCsCentralRatioPsi2sOverJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csCentralPsi2sOverJpsiY, yWidthArr, relErrFullPsi2sOverJpsiY)
    graSystCsCentralRatioPsi2sOverJpsiY.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
    graSystCsCentralRatioPsi2sOverJpsiY.SetLineColor(ROOT.kOrange+7)

    

    canvasCsPsi2sOverJpsiY = TCanvas("canvasCsPsi2sOverJpsiY", "canvasCsPsi2sOverJpsiY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsPsi2sOverJpsiY  = TH2F("histGridCsPsi2sOverJpsiY", "", 100, 2.5, 4.0, 100, 0.01, 1.)
    histGridCsPsi2sOverJpsiY.GetXaxis().SetTitle("#it{y}")
    histGridCsPsi2sOverJpsiY.GetYaxis().SetTitle("(d#sigma_{#psi(2S)} / d#it{y}) / (d#sigma_{J/#psi} / d#it{y})")
    histGridCsPsi2sOverJpsiY.Draw()
    graSystCsCentralRatioPsi2sOverJpsiY.Draw("E2 SAME")
    graStatCsCentralRatioPsi2sOverJpsiY.Draw("P SAME")
    graStatCsDoubleRenormRatioPsi2sOverJpsiY.Draw("SAME")
    graStatCsHalfRenormRatioPsi2sOverJpsiY.Draw("SAME")
    graStatCsHighMassRatioPsi2sOverJpsiY.Draw("SAME")
    graStatCsLowMassRatioPsi2sOverJpsiY.Draw("SAME")
    canvasCsPsi2sOverJpsiY.Update()

    # Psi(2S) vs Y
    graStatCsCentralPsi2sY = TGraphErrors(len(yCentrArr), yCentrArr, csCentralPsi2sY, yWidthArr, 0)
    SetGraStat(graStatCsCentralPsi2sY, 20, ROOT.kBlack)

    graStatCsDoubleRenormPsi2sY = TGraphErrors(len(yCentrArr), yCentrArr, csDoubleRenormPsi2sY, 0, 0)
    SetGraStat(graStatCsDoubleRenormPsi2sY, 20, ROOT.kRed+1)

    graStatCsHalfRenormPsi2sY = TGraphErrors(len(yCentrArr), yCentrArr, csHalfRenormPsi2sY, 0, 0)
    SetGraStat(graStatCsHalfRenormPsi2sY, 20, ROOT.kRed+1)

    graStatCsHighMassPsi2sY = TGraphErrors(len(yCentrArr), yCentrArr, csHighMassPsi2sY, 0, 0)
    SetGraStat(graStatCsHighMassPsi2sY, 20, ROOT.kAzure+2)

    graStatCsLowMassPsi2sY = TGraphErrors(len(yCentrArr), yCentrArr, csLowMassPsi2sY, 0, 0)
    SetGraStat(graStatCsLowMassPsi2sY, 20, ROOT.kAzure+2)

    # J/psi vs Y
    graStatCsCentralJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csCentralJpsiY, yWidthArr, 0)
    SetGraStat(graStatCsCentralJpsiY, 20, ROOT.kBlack)
    graStatCsCentralJpsiY.SetLineStyle(ROOT.kDashed)

    graStatCsDoubleRenormJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csDoubleRenormJpsiY, 0, 0)
    SetGraStat(graStatCsDoubleRenormJpsiY, 20, ROOT.kRed+1)
    graStatCsDoubleRenormJpsiY.SetLineStyle(ROOT.kDashed)

    graStatCsHalfRenormJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csHalfRenormJpsiY, 0, 0)
    SetGraStat(graStatCsHalfRenormJpsiY, 20, ROOT.kRed+1)
    graStatCsHalfRenormJpsiY.SetLineStyle(ROOT.kDashed)

    graStatCsHighMassJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csHighMassJpsiY, 0, 0)
    SetGraStat(graStatCsHighMassJpsiY, 20, ROOT.kAzure+2)
    graStatCsHighMassJpsiY.SetLineStyle(ROOT.kDashed)

    graStatCsLowMassJpsiY = TGraphErrors(len(yCentrArr), yCentrArr, csLowMassJpsiY, 0, 0)
    SetGraStat(graStatCsLowMassJpsiY, 20, ROOT.kAzure+2)
    graStatCsLowMassJpsiY.SetLineStyle(ROOT.kDashed)

    canvasCsY = TCanvas("canvasCsY", "canvasCsY", 800, 600)
    ROOT.gPad.SetLogy(1)
    histGridCsY  = TH2F("histGridCsY", "", 100, 2.5, 4.0, 100, 1e2, 1e4)
    histGridCsY.GetXaxis().SetTitle("#it{y}")
    histGridCsY.GetYaxis().SetTitle("d#sigma / d#it{y}")
    histGridCsY.Draw()
    graStatCsCentralPsi2sY.Draw("SAME")
    graStatCsDoubleRenormPsi2sY.Draw("SAME")
    graStatCsHalfRenormPsi2sY.Draw("SAME")
    graStatCsHighMassPsi2sY.Draw("SAME")
    graStatCsLowMassPsi2sY.Draw("SAME")
    graStatCsCentralJpsiY.Draw("SAME")
    graStatCsDoubleRenormJpsiY.Draw("SAME")
    graStatCsHalfRenormJpsiY.Draw("SAME")
    graStatCsHighMassJpsiY.Draw("SAME")
    graStatCsLowMassJpsiY.Draw("SAME")
    canvasCsY.Update()





    dfCsJpsiFonll = pd.read_csv('cs_jpsi_fonll.txt', sep=' ')
    ptMinFonll = dfCsJpsiFonll["ptmin"].to_numpy()
    ptMaxFonll = dfCsJpsiFonll["ptmax"].to_numpy()
    ptCentrFonll = (ptMinFonll + ptMaxFonll) / 2.
    ptWidthFonll = (ptMaxFonll - ptMinFonll) / 2.
    csCentralJpsiFonllPt = dfCsJpsiFonll["central"].to_numpy()
    csMinJpsiFonllPt = dfCsJpsiFonll["min"].to_numpy()
    csMaxJpsiFonllPt = dfCsJpsiFonll["max"].to_numpy()
    csMinScJpsiFonllPt = dfCsJpsiFonll["min_sc"].to_numpy()
    csMaxScJpsiFonllPt = dfCsJpsiFonll["max_sc"].to_numpy()
    csMinMassJpsiFonllPt = dfCsJpsiFonll["min_mass"].to_numpy()
    csMaxMassJpsiFonllPt = dfCsJpsiFonll["max_mass"].to_numpy()
    csMinPdfJpsiFonllPt = dfCsJpsiFonll["min_pdf"].to_numpy()
    csMaxPdfJpsiFonllPt = dfCsJpsiFonll["max_pdf"].to_numpy()

    print(ToCArray(csCentralJpsiFonllPt, ctype='double', name='csCentralJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinJpsiFonllPt, ctype='double', name='csMinJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxJpsiFonllPt, ctype='double', name='csMaxJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinScJpsiFonllPt, ctype='double', name='csMinScJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxScJpsiFonllPt, ctype='double', name='csMaxScJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinMassJpsiFonllPt, ctype='double', name='csMinMassJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxMassJpsiFonllPt, ctype='double', name='csMaxMassJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinPdfJpsiFonllPt, ctype='double', name='csMinPdfJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxPdfJpsiFonllPt, ctype='double', name='csMaxPdfJpsiFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print("------------------------------------------")

    dfCsPsi2sFonll = pd.read_csv('cs_psi2s_fonll.txt', sep=' ')
    ptMinFonll = dfCsPsi2sFonll["ptmin"].to_numpy()
    ptMaxFonll = dfCsPsi2sFonll["ptmax"].to_numpy()
    csCentralPsi2sFonllPt = dfCsPsi2sFonll["central"].to_numpy()
    csMinPsi2sFonllPt = dfCsPsi2sFonll["min"].to_numpy()
    csMaxPsi2sFonllPt = dfCsPsi2sFonll["max"].to_numpy()
    csMinScPsi2sFonllPt = dfCsPsi2sFonll["min_sc"].to_numpy()
    csMaxScPsi2sFonllPt = dfCsPsi2sFonll["max_sc"].to_numpy()
    csMinMassPsi2sFonllPt = dfCsPsi2sFonll["min_mass"].to_numpy()
    csMaxMassPsi2sFonllPt = dfCsPsi2sFonll["max_mass"].to_numpy()
    csMinPdfPsi2sFonllPt = dfCsPsi2sFonll["min_pdf"].to_numpy()
    csMaxPdfPsi2sFonllPt = dfCsPsi2sFonll["max_pdf"].to_numpy()

    print(ToCArray(csCentralPsi2sFonllPt, ctype='double', name='csCentralPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinPsi2sFonllPt, ctype='double', name='csMinPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxPsi2sFonllPt, ctype='double', name='csMaxPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinScPsi2sFonllPt, ctype='double', name='csMinScPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxScPsi2sFonllPt, ctype='double', name='csMaxScPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinMassPsi2sFonllPt, ctype='double', name='csMinMassPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxMassPsi2sFonllPt, ctype='double', name='csMaxMassPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMinPdfPsi2sFonllPt, ctype='double', name='csMinPdfPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))
    print(ToCArray(csMaxPdfPsi2sFonllPt, ctype='double', name='csMaxPdfPsi2sFonllPt', formatter=lambda x: '{:0.5f}'.format(x)))



    input()
    exit()












    # Stampa il DataFrame con la nuova colonna
    #print(dfCsJpsiIcemPt)



    #sum_columns = dfCsJpsiIcemPt.iloc[:, 1:7].sum()

    #print("Somma delle colonne dalla 1 alla 6:")
    #print(sum_columns)





        



if __name__ == '__main__':
    main()