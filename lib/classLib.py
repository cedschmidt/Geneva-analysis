# coding: utf-8

from __future__ import absolute_import, division, print_function

import os
import glob
import numpy as np
from scipy.misc import imread
from scipy.optimize import curve_fit


#from progressbar import ProgressBar, Percentage, Bar
from lib.funLib import *


class scanData(object):
    def __init__(self, path = os.getcwd()):
        self.foldername = path
        self.noImage = int
        self.noPeaks = int
        
        self.imSizeInit = np.zeros(2, int)
        self.imSizeCurrent = np.zeros(2, int)
        
        self.currentImNum = 0
        self.rect = np.zeros(4, int)
        self.rectSpec = np.zeros(4, int)
        
        self.wavelRef = np.zeros(3, float)
        self.specRef = np.array([0, 20, 50], int)
        
        self.t0Ref = 35
        self.deltat = -20

        
        self.openInit()        
         
    def openInit(self):
        self.nameList = np.array(sorted(glob.glob(self.foldername + "/" + "*.png")))
        self.temp = imread(self.nameList[0])
        self.imSizeInit = np.shape(self.temp)
        self.imSizeCurrent = np.shape(self.temp)
             
        self.noImage = self.nameList.size
        self.refState = 0
        
        self.createStackInit()
        self.createCalibration()
        
        self.rect = [0, self.imSizeInit[0]-1, 0, self.imSizeInit[1]-1]
        self.rectPlotSpec = [0, self.imSizeInit[0]-1, 0, self.imSizeInit[1]-1]
        self.currentIm = self.imStackscanTot[self.rect[0]:self.rect[1], \
        self.rect[2]:self.rect[3], 0]
        self.Im1dSpec =self.imStackscanTot[:, :, 0]        
        
    def changeCurrentIm(self, value):
        self.currentIm = self.imStackscanTot[self.rect[0]:self.rect[1], \
        self.rect[2]:self.rect[3], value]
        self.Im1dSpec = self.imStackscanTot[self.rectSpec[0]:self.rectSpec[1]\
        , self.rectSpec[2]:self.rectSpec[3], value] 
                
    def cropIm(self):
        self.currentIm = self.imStackscanTot[self.rect[0]:self.rect[1], \
        self.rect[2]:self.rect[3], self.currentImNum]
        
        self.imStackscanTotCrop = np.zeros((self.rect[1]-self.rect[0], \
        self.rect[3]-self.rect[2], self.noImage), int)
        
        imStackscanCrop = np.zeros((self.rect[1]-self.rect[0], \
        self.rect[3]-self.rect[2], int(self.noImage/2)), int)
        
        imStackrefCrop = np.zeros((self.rect[1]-self.rect[0], \
        self.rect[3]-self.rect[2], int(self.noImage/2)), int)
        
        self.xaxisInit = np.arange(self.rect[2],self.rect[3], 1)
        self.xaxisCalib = self.xaxisInit
        
        for idx in range(self.noImage):
            self.imStackscanTotCrop[:, :, idx] = self.imStackscanTot[\
            self.rect[0]:self.rect[1], self.rect[2]:self.rect[3], idx]
            if idx%2 == 0:
                imStackscanCrop[:, :, int(idx/2)] = self.imStackscanTotCrop[:, :, idx]
            else:
                imStackrefCrop[:, :, int((idx-1)/2)] = self.imStackscanTotCrop[:, :, idx]
        
        
        self.spectrumMean = np.mean(self.imStackscanTotCrop, axis = (0,2))
        self.imScanCrop = np.mean(imStackscanCrop, axis = 0)
        
        self.imRefCrop = np.mean(imStackrefCrop, axis = 0)
        self.imNoRefCrop = 0.5*(self.imScanCrop + self.imRefCrop)
        self.imDiffCrop = 0.5*(self.imScanCrop - self.imRefCrop)/self.imRefCrop
        
        if self.refState == 0:
            self.imScanCurrent = self.imNoRefCrop
        elif self.refState == 2:
            self.imScanCurrent = self.imScanCrop
        
        
        if (self.specRef[0]<self.rect[2]) | (self.specRef[0]>self.rect[3]):
            self.specRef[0] = self.rect[2] + 10
            
        if (self.specRef[2]<self.rect[2]) | (self.specRef[2]>self.rect[3]):
            self.specRef[2] = self.rect[3] - 10
            
        if (self.specRef[1]<self.rect[2]) | (self.specRef[1]>self.rect[3]):
            self.specRef[1] = int(0.5*(self.specRef[0] + self.specRef[2]))
            
        
    def cropIm1dSpec(self):
        self.Im1dSpec = self.imStackscanTot[self.rectSpec[0]:self.rectSpec[1],\
        self.rectSpec[2]:self.rectSpec[3], self.currentImNum]
        
    def createStackInit(self):
        self.imStackscanTot = np.zeros((self.imSizeInit[0], self.imSizeInit[1], self.noImage), int)
        self.imStackscanTotCrop = np.zeros((self.imSizeInit[0], self.imSizeInit[1], self.noImage), int)
        imStackscanInit = np.zeros((self.imSizeInit[0], self.imSizeInit[1], int(self.noImage/2)), int)
        imStackrefInit = np.zeros((self.imSizeInit[0], self.imSizeInit[1], int(self.noImage/2)), int)
        
        #pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=self.noImage).start()
        print('Loading images in stack')
                
        for idx in range(self.noImage):
            #pbar.update(idx+1)
            self.imStackscanTot[:, :, idx] = imread(self.nameList[idx])
            if idx%2 == 0:
                imStackscanInit[:, :, int(idx/2)] = self.imStackscanTot[:, :, idx]
            else:
                imStackrefInit[:, :, int((idx-1)/2)] = self.imStackscanTot[:, :, idx]
                
        print("\n")
        
        self.imStackscanTotCrop = self.imStackscanTot
        self.imScanInit = np.mean(imStackscanInit, axis = 0)
        self.imRefInit = np.mean(imStackrefInit, axis = 0)
        self.imScanCrop = self.imScanInit
        self.imrefCrop = self.imRefInit
        self.imNoRefInit = 0.5*(self.imScanInit + self.imRefInit)
        self.imDiffInit = 0.5*(self.imScanInit - self.imRefInit)/self.imRefInit
        self.imNoRefCrop = self.imNoRefInit
        self.imDiffCrop = self.imDiffInit
        
        if self.refState == 0:
            self.imScanCurrent = self.imNoRefCrop
        elif self.refState == 2:
            self.imScanCurrent = self.imScanCrop
        
    def createCalibration(self):
        self.spectrumMean = np.mean(self.imStackscanTotCrop, axis = (0,2))
        self.spectrumScanPump = np.mean(self.imStackscanTotCrop, axis = (0,2))
        self.xaxisInit = np.arange(self.imSizeInit[1])
        self.xaxisCalib = self.xaxisInit
        self.taxis = np.arange(int(self.noImage/2))
        self.specRef[0] = np.amin(self.xaxisCalib + 10)
        self.specRef[1] = np.amax(self.xaxisCalib - 10)
        self.specRef[2] = int(0.5*(self.specRef[0] + self.specRef[2]))
        
    def specCalibrate(self):
        F = lambda x, D, alpha0: D/(x + alpha0)
        poptinit, pcov = curve_fit(F, self.specRef, self.wavelRef)

        F2 = lambda x, A, D, alpha0: D/(x/(np.sqrt(1 + A*x**2.)) + alpha0)
        popt, pcov = curve_fit(F2, self.specRef, self.wavelRef, p0 = [0., poptinit[0], poptinit[1]])
        
        print("Calibration results:\n")
        print("Spectrum calibrated with the following function : A/(x/(np.sqrt(1 + K*x**2.)) + alpha0)\n")
        print("A = ", popt[1])
        print("K = ", popt[0])
        print("alpha0 = ", popt[2])
        
        tempaxis = self.xaxisCalib
        self.xaxisCalib = F2(tempaxis, popt[0], popt[1], popt[2])

            
        
    
class analysisDatas(object):
    def __init__(self, *args):
        self.taxis = args[0]
        self.energy = args[1]
        
        try:
            if len(args) == 3:
                self.scanInit = args[2]
                self.scan = args[2]
                self.refInit = 0.
                self.ref = 0.
                self.normInit = 1.
                self.norm = 1.
                self.background = args[2]
    
            elif len(args) == 4:
                self.scanInit = args[2]
                self.scan = args[2]
                self.refInit = args[3]
                self.ref = args[3]
                self.normInit = args[2]
                self.norm = args[2]
                self.background = args[2]
        except:
            print("Analysis parameters are wrong")
            
        self.itemDict = {"Raw Scan":self.scanInit, "Raw ref":self.refInit, \
        "Raw norm":self.normInit, "1":1., "0":0.}
        
        self.updateDisplay()
        
        self.noProcessedScan = 0
        self.noProcessedRef = 0
        self.noProcessedNorm = 0  
        self.noBackground = 0
        
        self.lowerBoundNorm = np.amin(self.energy) + 10
        self.upperBoundNorm = np.amax(self.energy) - 10
        
        self.lowerBoundRef = np.amin(self.taxis) + 30
        self.upperBoundRef = np.amax(self.taxis) - 30
        
        self.axisSmoothList = {"Time":0, "Energy":1}
        self.axisSmoothCurrent = 0
        self.windowLen = 3
        self.windowType = 'hamming'
        
    def updateDisplay(self):
        self.display = (self.scan - self.ref)/self.norm
        self.displaytemp = (self.scan - self.ref)/self.norm
        
    def movingAverage(self):
        axis = self.axisSmoothCurrent
        windowLen = self.windowLen
        windowType = self.windowType
        
        if axis == 0:
            self.displaytemp = movingSmooth(self.display, axis, windowLen, windowType)
        elif axis == 1:
            self.background = movingSmooth(self.display, axis, windowLen, windowType)
            self.displaytemp = -(self.display - self.background)
        
    def normalization(self):
        lmin = np.argmin(np.abs(self.energy - self.lowerBoundNorm))
        lmax = np.argmin(np.abs(self.energy - self.upperBoundNorm))
        sumRef = np.sum(self.display[:, np.amin([lmin, lmax]):np.amax([lmin, lmax])], axis = 1)
        
        for idx in range(len(sumRef)):
            self.displaytemp[idx, :] = self.display[idx, :]/sumRef[idx]
            
    def createRef(self):
        lmin = np.argmin(np.abs(self.taxis - self.lowerBoundRef))
        lmax = np.argmin(np.abs(self.taxis - self.upperBoundRef))
        refmean = np.mean(self.display[np.amin([lmin, lmax]):np.amax([lmin, lmax]), :], axis = 0)
        
        self.ref = np.zeros_like(self.display)
        
        for i in range(len(self.taxis)):
            self.ref[i, :] = refmean
        
        self.displaytemp = self.display - self.ref  
        
    def getPeaks(self, lam, p, niter=10):
        self.background = peakExtract(-self.display, lam, p, niter)
        self.displaytemp = -(self.display - self.background)
        
    def validateChange(self):
        self.display = self.displaytemp
        
    def storeProcessedToScan(self):
        scan = self.display
        self.itemDict.update({"Processed Scan " + str(self.noProcessedScan):scan})
        self.noProcessedScan += 1
        
    def storeProcessedToRef(self):
        self.itemDict.update({"Processed ref " + str(self.noProcessedRef):self.ref})
        self.noProcessedRef += 1
        
    def storeProcessedToNorm(self):
        norm = self.display
        self.itemDict.update({"Processed norm " + str(self.noProcessedNorm):norm})
        self.noProcessedNorm += 1
        
    def storeToBack(self):
        back = self.background
        self.itemDict.update({"Background " + str(self.noBackground):back})
        self.noBackground += 1

        
    def saveAnalysisData(self, scanDate, scanName):
        
        path = os.getcwd() + "/Datas/"             
        pathname2 = "current/"
        pathname3 = "other"
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        if not os.path.exists(path + scanDate + "/" + scanName + "/" + pathname2):
            os.makedirs(path + scanDate + "/" + scanName + "/" + pathname2)
            
        if not os.path.exists(path + scanDate + "/" + scanName + "/" + pathname3):
            os.makedirs(path + scanDate + "/" + scanName + "/" + pathname3)
            
        savingPath1 = path + scanDate + "/" + scanName + "/" + pathname2 + "/"
        savingPath2 = path + scanDate + "/" + scanName + "/" + pathname3 + "/"
        
        np.savetxt(savingPath1 + "energy.txt", self.energy)
        np.savetxt(savingPath1 + "time.txt", self.taxis)
        np.savetxt(savingPath2 + "energy.txt", self.energy)
        np.savetxt(savingPath2 + "time.txt", self.taxis)
        
        np.savetxt(savingPath1 + "currentScan.txt", self.display)
        
        for key, value in self.itemDict.iteritems():
            if np.size(value)>1:
                np.savetxt(savingPath2 + str(key) + ".txt", value)
            
        
        
        

        
        

            
        
            
        
        



