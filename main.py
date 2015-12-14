# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:37:48 2015

@author: pertoty
"""


#  Import libraries and classes
from __future__ import absolute_import, division, print_function
import sys, gc
    
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.uic import loadUiType

import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from lib.classLib import scanData, analysisDatas
from lib.funLib import  DraggableVLine, DraggableHLine, RectangleSelector


Ui_MainWindow, QMainWindow = loadUiType('gui_layout/lh_gui.ui')


# Creation of classes for analysis


class Main(QMainWindow, Ui_MainWindow):

###########################    Initlization part    ###########################
    def __init__(self, foldername = None):
        '''
        Initialization of main window
        '''
        if foldername == None:
            #Initiate main window
            super(Main, self).__init__()
            self.setupUi(self)
            self.isDirectlyClose = False
            
            #Menu creation
            self.create_menu()
            
            #Creation of figures
                #Images for area selection and image scrolling
            self.fig2d = Figure(facecolor='1')
            self.canvas2d = FigureCanvas(self.fig2d)
                
                # Spectrum for wavelength calibration
            self.fig1dSpecCalib = Figure(facecolor='1')
            self.canvas1dSpecCalib = FigureCanvas(self.fig1dSpecCalib)
                #Image for 2d scan (spectrum vs time)
            self.fig2dScan = Figure(facecolor='1')
            self.canvas2dScan = FigureCanvas(self.fig2dScan)
                #Image for 2d analysis
            self.figAnalysis = Figure(facecolor='1')
            self.canvasAnalysis = FigureCanvas(self.figAnalysis)
                #Figure for 1d background baseline
            self.figBackBaseline = Figure(facecolor='1')
            self.canvasBackBaseline = FigureCanvas(self.figBackBaseline)
        else:
            #Initiate main window
            self.setupUi(self)
            
            gc.collect()
            self.isDirectlyClose = False
            
            #Creation of figures
                #Images for area selection and image scrolling
            self.fig2d = Figure(facecolor='1')
            self.canvas2d = FigureCanvas(self.fig2d)
                
                # Spectrum for wavelength calibration
            self.fig1dSpecCalib = Figure(facecolor='1')
            self.canvas1dSpecCalib = FigureCanvas(self.fig1dSpecCalib)
                #Image for 2d scan (spectrum vs time)
            self.fig2dScan = Figure(facecolor='1')
            self.canvas2dScan = FigureCanvas(self.fig2dScan)
                #Image for 2d analysis
            self.figAnalysis = Figure(facecolor='1')
            self.canvasAnalysis = FigureCanvas(self.figAnalysis)
                #Figure for 1d background baseline
            self.figBackBaseline = Figure(facecolor='1')
            self.canvasBackBaseline = FigureCanvas(self.figBackBaseline)
            
            # Creation of scan object and range for image display (specific to each scan)                  
            self.initScanObject(foldername)
            
            # Initilization of widgets related to the Calibration scan object 
            self.initCalibrationTab()
            self.tabWidget.setCurrentIndex(0)
            
        
    def close (self):
        for childQWidget in self.findChildren(QWidget):
            childQWidget.close()
        self.isDirectlyClose = True
        return QMainWindow.close(self)

    def closeEvent (self, eventQCloseEvent):
        if self.isDirectlyClose:
            eventQCloseEvent.accept()
        else:
            answer = QMessageBox.question (
                self,
                'Quit',
                'Are you sure you want to quit ?',
                QMessageBox.Yes,
                QMessageBox.No)
            if (answer == QMessageBox.Yes) or (self.isDirectlyClose == True):
                eventQCloseEvent.accept()
            else:
                eventQCloseEvent.ignore()

 
    def create_menu(self):
        '''
        This function will create all the menu entrance and connect them to the
        corresponding function
        '''
        exitAction = QAction(QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qApp.quit)

        openAction = QAction('Select folder', self)
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Select a folder')
        openAction.triggered.connect(self.openFile)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openAction)
        fileMenu.addAction(exitAction)

    def openFile(self):
        '''
        This function is connected to the 'Select folder' menu entrance. It will
        create a scan object which will be used for all the calibration.
        '''
        foldername = str(QFileDialog.getExistingDirectory(
                              self, 'Select folder'))
                              
        self.__init__(foldername)

#        # Creation of scan object and range for image display (specific to each scan)                  
#        self.initScanObject()
#        
#        # Initilization of widgets related to the Calibration scan object 
#        self.initCalibrationTab()
    
        
################################   Calibration tab ############################
        
    def initCalibrationTab(self):
        print("test2")
        self.initSliderH()
        self.initSliderV()
        self.plot2d()
        self.plotCalibSpec()
        self.plot2dScan()
        self.initButton()
        self.initRefText()
        self.initCheckRef()
        self.initCheckSub()
        self.initCheckScanUpdate()
        self.firstCheck = 1
        self.initdeltatText()
        self.initSaveScan()
        
    def initScanObject(self, folderName):
        self.lhObj = scanData(folderName)
        self.vmin = np.amin(self.lhObj.imStackscanTot)
        self.vmax = np.amax(self.lhObj.imStackscanTot)
        self.imNumber.setText(str(self.lhObj.currentImNum)+'/'+ str(self.lhObj.noImage))


######################################################## Overview and selection

    def plot2d(self):
        '''
        Once the folder has been selected and scan object created, this function 
        initiate the display of the 2d images containing all the images of the
        scan
        '''
        self.fig2d.clear()
        self.axes2d = self.fig2d.add_subplot(111)
        
        self.im2d = self.axes2d.imshow(self.lhObj.currentIm, aspect='auto',
                                       interpolation='nearest',
                                       vmin=self.vmin,
                                       vmax=self.vmax)

        self.rect = RectangleSelector(self.axes2d, self._on_mouse_release, useblit = True)
                                      

        self.axes2d.set_xlabel('Pixel number')
        self.axes2d.set_ylabel('Pixel number')
        self.canvas2d.draw()
        self.mpl2dIm.addWidget(self.canvas2d)
        
    def changeImage(self, value):
        self.lhObj.currentImNum = value
        self.lhObj.changeCurrentIm(value)
        self.updateImage2d()
        self.imNumber.setText(str(self.lhObj.currentImNum)+'/'+ str(self.lhObj.noImage))

    def updateImage2d(self):
        self.axes2d.cla()
        self.im2d = self.axes2d.imshow(self.lhObj.currentIm,  aspect='auto',
                                       interpolation='nearest',
                                       vmin=self.vmin,
                                       vmax=self.vmax)
        self.canvas2d.draw()
        self.rect.update()
        
    def changeColor(self, value):
        self.vmax = value
        self.im2d.set_clim(vmax=value)
        self.canvas2d.draw()
        self.rect.update()
        
    def initButton(self):
        self.imCrop.clicked.connect(self._onpressButtonCrop)
        self.imInit.clicked.connect(self._onpressButtonInit)
        self.push_calib.clicked.connect(self._onpressCalib)

    def initSliderH(self):
        self.imSlider.setMaximum(self.lhObj.noImage-1)
        self.imSlider.valueChanged[int].connect(self.changeImage)

    def initSliderV(self):
        self.colorSlider.setMinimum(np.amin(self.lhObj.currentIm))
        self.colorSlider.setMaximum(np.amax(self.lhObj.currentIm))
        self.colorSlider.valueChanged[int].connect(self.changeColor)  
        
    def _on_mouse_release(self, epress=None, erelease=None):
        x0, y0, width, height = self.rect._rect_bbox
        self.lhObj.rectSpec = [y0 + self.lhObj.rect[0], 
                                   self.lhObj.rect[0] + y0 + height,
                                   self.lhObj.rect[2] + x0,
                                   self.lhObj.rect[2] + x0 + width]
        
    def _onpressButtonCrop(self):
        x0, y0, width, height = np.int64(self.rect._rect_bbox)
        self.lhObj.rect = [y0 + self.lhObj.rect[0], self.lhObj.rect[0] + y0 +
                           height, self.lhObj.rect[2] + x0,
                           self.lhObj.rect[2] + x0 + width]

        self.lhObj.cropIm()
        self.updateImage2d() 
        self.updatePlotCalibSpec()
        self.updatePlot2dScan()

    def _onpressButtonInit(self):
        self.lhObj.rect = [0, self.lhObj.imSizeInit[0]-1, 0, self.lhObj.imSizeInit[1]-1]
        self.lhObj.rectSpec = [0, self.lhObj.imSizeInit[0]-1, 0, self.lhObj.imSizeInit[1]-1]
        
        self.lhObj.changeCurrentIm(self.lhObj.currentImNum)
        self.updateImage2d()
        
        
    
######################################################## Wavelength calib box
        
        
    def plotCalibSpec(self):
        '''
        Initialization of 1d plot of calibration spectrum with cropped images
        '''
        self.fig1dSpecCalib.clear()
        self.axes1dSpecCalib = self.fig1dSpecCalib.add_subplot(111)

        self.im1dSpecCalib = self.axes1dSpecCalib.plot(self.lhObj.xaxisCalib, \
        self.lhObj.spectrumMean)
        self.line1 = DraggableVLine(self.moveRef1, self.axes1dSpecCalib, self.canvas1dSpecCalib, self.lhObj.specRef[0], colorV = 'red')
        self.line2 = DraggableVLine(self.moveRef2,self.axes1dSpecCalib, self.canvas1dSpecCalib, self.lhObj.specRef[1], colorV = 'blue')
        self.line3 = DraggableVLine(self.moveRef3,self.axes1dSpecCalib, self.canvas1dSpecCalib, self.lhObj.specRef[2], colorV = 'green')
        
        self.canvas1dSpecCalib.draw()
        self.mplSpecCalib.addWidget(self.canvas1dSpecCalib)
        
    ########## 1D spectrum for calibration
        
    def updatePlotCalibSpec(self):
        
        self.axes1dSpecCalib.cla()

        self.im1dSpecCalib = self.axes1dSpecCalib.plot(self.lhObj.xaxisCalib, \
        self.lhObj.spectrumMean)
        self.line1.update(self.lhObj.specRef[0])
        self.line2.update(self.lhObj.specRef[1])
        self.line3.update(self.lhObj.specRef[2])
        
        self.axes2d.set_xlabel('Photon energy in eV')
        self.axes2d.set_ylabel('Intensity')
        
        self.canvas1dSpecCalib.draw()
        
    def moveRef1(self, value):
        self.lhObj.specRef[0] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
        
    def moveRef2(self, value):
        self.lhObj.specRef[1] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
        
    def moveRef3(self, value):
        self.lhObj.specRef[2] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
           
    def initRefText(self):
        self.ref1.setText("196.2")
        self.ref1.returnPressed.connect(self.readRef1)
        
        self.ref2.setText("184")
        self.ref2.returnPressed.connect(self.readRef2)
        
        self.ref3.setText("173")
        self.ref3.returnPressed.connect(self.readRef3)
        
    def readRef1(self):
        self.lhObj.wavelRef[0] = float(self.ref1.text())
        
    def readRef2(self):
        self.lhObj.wavelRef[1] = float(self.ref2.text())
        
    def readRef3(self):
        self.lhObj.wavelRef[2] = float(self.ref3.text())
        
    def _onpressCalib(self):
        self.readRef1()
        self.readRef2()
        self.readRef3()
        self.lhObj.specCalibrate()
        self.lhObj.specRef[0] = float(self.ref1.text())
        self.lhObj.specRef[1] = float(self.ref2.text())
        self.lhObj.specRef[2] = float(self.ref3.text())
        
        self.updatePlotCalibSpec()
        self.updatePlot2dScan()
        
        
        
        
######################################################## time calibration

        
    def plot2dScan(self):
        '''
        Once the folder has been selected and scan object created, this function 
        initiate the display of the 2d images containing all the images of the
        scan
        '''
        self.fig2dScan.clear()
        self.axes2dScan = self.fig2dScan.add_subplot(111)
        
        self.im2dScan = self.axes2dScan.pcolormesh(self.lhObj.xaxisCalib,
                                               self.lhObj.taxis,                                                                              
                                               np.transpose(self.lhObj.imNoRefCrop))
                                               
        self.linet0 = DraggableHLine(self.movet0Ref, self.axes2dScan, self.canvas2dScan, self.lhObj.t0Ref, colorH = 'white')
                                              
        self.axes2dScan.set_xlabel('Photon Energy')
        self.axes2dScan.set_ylabel('Time')
        self.canvas2dScan.draw()
        self.scanRaw.addWidget(self.canvas2dScan)


    def initdeltatText(self):
        self.deltat.setText("-20")
        self.deltat.returnPressed.connect(self.readDeltat)
        
    def initCheckRef(self):
        self.ref_onoff.stateChanged.connect(self.updateRef)
        
    def initCheckSub(self):
        self.ref_substract.stateChanged.connect(self.substractRef)
        
    def initCheckScanUpdate(self):
        self.updateGraph.clicked.connect(self._onpressScanUpdate)
        
    def initSaveScan(self):
        self.saveScan.clicked.connect(self._onpressSaveScan)
        
    def _onpressSaveScan(self):
        self.initAnalysis()      
        
    ######### 2D image for Scan
    def updatePlot2dScan(self):
        '''
        Once the folder has been selected and scan object created, this function 
        initiate the display of the 2d images containing all the images of the
        scan
        '''
        self.axes2dScan.cla()
        
        self.im2dScan = self.axes2dScan.pcolormesh(self.lhObj.xaxisCalib,
                                               self.lhObj.taxis,                                                                              
                                               np.transpose(self.lhObj.imScanCurrent))
                                               
        self.linet0.update(self.lhObj.t0Ref)
                                               
        self.axes2dScan.set_xlabel('Photon Energy')
        self.axes2dScan.set_ylabel('Time')
        self.canvas2dScan.draw()
        
    
    def updateRef(self, state): 
        self.lhObj.refState = state
        
        if state == 2:
            self.lhObj.imScanCurrent = self.lhObj.imScanCrop
        elif state == 0:
            self.lhObj.imScanCurrent = self.lhObj.imNoRefCrop
            
        self.updatePlot2dScan()
        
    def substractRef(self, state):
        if self.lhObj.refState == 2:
            if state == 2:
                self.lhObj.imScanCurrent = self.lhObj.imDiffCrop
            elif state == 0:
                    self.lhObj.imScanCurrent = self.lhObj.imScanCrop
                    
        self.updatePlot2dScan()
        
    def movet0Ref(self, value):
        self.lhObj.t0Ref = value
        
    def readDeltat(self):
        self.lhObj.deltat = float(self.deltat.text())
        
    def _onpressScanUpdate(self):
        self.readDeltat()
        
        if self.firstCheck == 1:
            self.lhObj.taxis = (self.lhObj.taxis - self.lhObj.t0Ref)*self.lhObj.deltat
        elif self.firstCheck == 0:
            self.lhObj.taxis = (self.lhObj.taxis - self.lhObj.t0Ref)
            
        self.firstCheck = 0
        self.lhObj.t0Ref = 0
        self.updatePlot2dScan()
      
        
        
        


        
     

###########################      Analysis tab      ############################
    #### Initialization of analysis tab
    def initAnalysis(self):
        
        if self.lhObj.refState == 0:
            try:
                self.AObj = analysisDatas(self.lhObj.taxis, self.lhObj.xaxisCalib,\
                (self.lhObj.imNoRefCrop).T)
            except:
                print("Cannot initiate analysis class")
        elif self.lhObj.refState == 2:
            print("Begin initiating", self.lhObj.refState)
            try:
                self.AObj = analysisDatas(self.lhObj.taxis, self.lhObj.xaxisCalib,\
                (self.lhObj.imScanCrop).T, (self.lhObj.imRefCrop).T)
            except:
                print("Cannot initiate analysis class") 
        else:
            print("refState is wrong")
            
        try:    
            self.initComboBoxes()
        except:
            print("Cannot initiate combo boxes")
          
        try:
            self.initStoreButton()
        except:
            print("Cannot initiate Store button")
            
        try:   
            self.initAnalysisGraph()
        except:
            print("Cannot initiate Analysis Graph")            
            
        try:
            self.initBackBaseline()
        except:
            print("Cannot initiate Baseline graph")
            
        try:
            self.initNormPreviewButton()
            self.initNormCancelButton()
            self.initNormApplyButton()
        except:
            print("Cannot initiate Normalization buttons")
        
        try:
            self.initRefPreviewButton()
            self.initRefAppButton()
            self.initRefCancelButton()
        except:
            print("Cannot initiate Reference button")
        
        try:
            self.initAxisBox()
            self.initWindowBox()
            self.initSmoothWindowLen()
            self.initSmoothPreview()
            self.initMovingCancelButton()
            self.initMovingApplyButton()
        except:
            print("Cannot initiate smooting box")
         
        try:
            self.initBaselineDetect()
        except:
            print("Baseline detection not initiated")
        
        try:
            self.initSavingEdit()
            self.initSaveButton()
        except:
            print("Cannot initiate saving box")
            
            
    ##### Graph initialization
            
    def initAnalysisGraph(self):

        self.figAnalysis.clear()
        self.axesAnalysis = self.figAnalysis.add_subplot(111)
        
        self.axesAnalysis.pcolormesh(self.AObj.energy, 
                                     self.AObj.taxis,                                                                           
                                     self.AObj.displaytemp)
                                     
        self.lineVlow = DraggableVLine(self.onmoveVlow, self.axesAnalysis, self.canvasAnalysis, self.AObj.lowerBoundNorm, colorV = 'white')
        self.lineVup = DraggableVLine(self.onmoveVup, self.axesAnalysis, self.canvasAnalysis, self.AObj.upperBoundNorm, colorV = 'white')
        
        self.lineHlow = DraggableHLine(self.onmoveHlow, self.axesAnalysis, self.canvasAnalysis, self.AObj.lowerBoundRef, colorH = 'white')
        self.lineHup = DraggableHLine(self.onmoveHup, self.axesAnalysis, self.canvasAnalysis, self.AObj.upperBoundRef, colorH = 'white')
                
        self.axesAnalysis.set_xlabel('Photon energy in eV')
        self.axesAnalysis.set_ylabel('Time delay (fs)')
        self.canvasAnalysis.draw()
        self.analysisGraph.addWidget(self.canvasAnalysis)
        
        
    def updateAnalysisGraph(self):
        self.axesAnalysis.cla()
        self.axesAnalysis.pcolormesh(self.AObj.energy, 
                                     self.AObj.taxis,                                                                           
                                     self.AObj.displaytemp)
                                     
        self.lineVlow.update(self.AObj.lowerBoundNorm)
        self.lineVup.update(self.AObj.upperBoundNorm)
        self.lineHlow.update(self.AObj.lowerBoundRef)
        self.lineHup.update(self.AObj.upperBoundRef)
        

        self.axesAnalysis.set_xlabel('Photon energy in eV')
        self.axesAnalysis.set_ylabel('Time delay (fs)')
        self.canvasAnalysis.draw()
        
        
    def initBackBaseline(self):

        self.figBackBaseline.clear()
        self.axesBackBaseline = self.figBackBaseline.add_subplot(111)
        
        self.axesBackBaseline.plot(self.AObj.energy, np.mean(self.AObj.display, axis = 0))
        self.axesBackBaseline.plot(self.AObj.energy, np.mean(self.AObj.background, axis = 0))

        self.axesBackBaseline.set_xlabel('Photon energy in eV')
        self.axesBackBaseline.set_ylabel('Time delay (fs)')
        self.canvasBackBaseline.draw()
        self.backBaseline.addWidget(self.canvasBackBaseline)
        
    def updateBackBaseline(self):

        self.axesBackBaseline.cla()
        
        self.axesBackBaseline.plot(self.AObj.energy, np.mean(self.AObj.display, axis = 0))
        self.axesBackBaseline.plot(self.AObj.energy, np.mean(self.AObj.background, axis = 0))


        self.axesBackBaseline.set_xlabel('Photon energy in eV')
        self.axesBackBaseline.set_ylabel('Time delay (fs)')
        self.canvasBackBaseline.draw()
        
    #####Graph display box   
        
    def initComboBoxes(self):
        self.scanbox.clear()
        self.refbox.clear()
        self.normbox.clear()
        
        if self.lhObj.refState == 0:
            self.scanbox.addItem("Raw Scan")            
            self.refbox.addItem("0")
            self.normbox.addItem("1")
            
        elif self.lhObj.refState == 2:
            self.scanbox.addItem("Raw Scan")
            self.refbox.addItem("Raw ref")
            self.refbox.addItem("0")
            self.normbox.addItem("Raw Scan")
            self.normbox.addItem("1")
            
        self.scanbox.activated[str].connect(self.onScanBoxActivated)
        self.refbox.activated[str].connect(self.onRefBoxActivated)
        self.normbox.activated[str].connect(self.onNormBoxActivated)
        
    def initStoreButton(self):
        self.scanStore.clicked.connect(self._onpressscanStore)
        self.refStore.clicked.connect(self._onpressrefStore)
        self.normStore.clicked.connect(self._onpressnormStore)
        
    def onScanBoxActivated(self, text):
        self.AObj.scan = self.AObj.itemDict[str(text)]
        self.AObj.updateDisplay()
        self.updateAnalysisGraph()
        
    def onRefBoxActivated(self, text):
        self.AObj.ref = self.AObj.itemDict[str(text)]
        self.AObj.updateDisplay()
        self.updateAnalysisGraph()
        
    def onNormBoxActivated(self, text):
        self.AObj.norm = self.AObj.itemDict[str(text)]
        self.AObj.updateDisplay()
        self.updateAnalysisGraph()
        
    def _onpressscanStore(self):
        self.scanbox.addItem("Processed Scan " + str(self.AObj.noProcessedScan))
        self.AObj.storeProcessedToScan()
        
    def _onpressrefStore(self):
        self.refbox.addItem("Processed ref " + str(self.AObj.noProcessedRef))
        self.AObj.storeProcessedToRef()
        
    def _onpressnormStore(self):
        self.normbox.addItem("Processed norm " + str(self.AObj.noProcessedNorm))
        self.AObj.storeProcessedToNorm()
        
        
    ##### Normalization box
    
        
    def initNormPreviewButton(self):
        self.normalizationPre.clicked.connect(self.updateNormPreview)
        
    def initNormCancelButton(self):
        self.normalizationCan.clicked.connect(self.updateNormCan)
        
    def initNormApplyButton(self):
        self.normalizationApp.clicked.connect(self.updateNormApp)
        
    def updateNormPreview(self):
        self.AObj.normalization()
        self.updateAnalysisGraph() 
        
    def updateNormCan(self):
        self.AObj.displaytemp = self.AObj.display
        self.updateAnalysisGraph()
        
    def updateNormApp(self):
        self.scanbox.addItem("Processed Scan " + str(self.AObj.noProcessedScan))
        index = self.scanbox.findText("Processed Scan " + str(self.AObj.noProcessedScan))
        self.AObj.validateChange()
        self.AObj.storeProcessedToScan()        
        self.scanbox.setCurrentIndex(index)    
        
    def onmoveVlow(self, value):
        self.AObj.lowerBoundNorm = value
     
    def onmoveVup(self, value):
        self.AObj.upperBoundNorm = value
        
        
        
    ##### Reference box    
      
        
    def initRefPreviewButton(self):
        self.referencePre.clicked.connect(self.updateRefPreview)
        
    def initRefCancelButton(self):
        self.referenceCan.clicked.connect(self.updateRefCan)
        
    def initRefAppButton(self):
        self.referenceApp.clicked.connect(self.updateRefApp)  
        
    def updateRefCan(self):
        self.AObj.displaytemp = self.AObj.display
        self.updateAnalysisGraph()
        
    def updateRefPreview(self):
        self.AObj.createRef()
        self.updateAnalysisGraph() 
        
    def updateRefApp(self):
        self.AObj.display = self.AObj.displaytemp
        self.refbox.addItem("Processed ref " + str(self.AObj.noProcessedRef))
        index = self.refbox.findText("Processed ref " + str(self.AObj.noProcessedRef))
        self.AObj.storeProcessedToRef()
        self.refbox.setCurrentIndex(index)
        
    def onmoveHlow(self, value):
        self.AObj.lowerBoundRef = value
        
    def onmoveHup(self, value):
        self.AObj.upperBoundRef = value
        
        
        
    ##### Moving average box
        
    def initAxisBox(self):
        self.axisBox.clear()
        self.axisBox.addItem("Time")
        self.axisBox.addItem("Energy")        
        self.axisBox.activated[str].connect(self.setSmoothAxis)
        
    def initWindowBox(self):
        self.windowBox.clear()
        self.windowBox.addItems(['flat', 'hanning', 'hamming', 'bartlett', 'blackman'])       
        self.windowBox.activated[str].connect(self.setSmoothWindow)
        
    def initSmoothWindowLen(self):
        self.windowLenEdit.setValue(int(self.AObj.windowLen))
        self.windowLenEdit.valueChanged.connect(self.setSmoothWindowLen)
        
    def initSmoothPreview(self):
        self.movingPre.clicked.connect(self.updateSmoothGraph)   
        
    def initMovingCancelButton(self):
        self.movingCan.clicked.connect(self.updateNormCan)
        
    def initMovingApplyButton(self):
        if self.AObj.axisSmoothCurrent == 1:
            self.movingApp.clicked.connect(self.updateNormApp)
        elif self.AObj.axisSmoothCurrent == 0:
            self.movingApp.clicked.connect(self.updateBackApp)
          
    def setSmoothAxis(self, text):
        self.AObj.axisSmoothCurrent = self.AObj.axisSmoothList[str(text)]
        
    def setSmoothWindow(self, text):
        self.AObj.windowType = str(text)
        
    def setSmoothWindowLen(self):
        self.AObj.windowLen = int(self.windowLenEdit.value())
        
    def updateSmoothGraph(self):
        self.AObj.movingAverage()
        if self.AObj.axisSmoothCurrent == 0:
            self.updateAnalysisGraph() 
        elif self.AObj.axisSmoothCurrent == 1:
            self.updateAnalysisGraph() 
            self.updateBackBaseline()
            
    def updateBackApp(self):
        self.scanbox.addItem("Processed Scan " + str(self.AObj.noProcessedScan))
        index1 = self.scanbox.findText("Processed Scan " + str(self.AObj.noProcessedScan))
        
        self.refbox.addItem("Background " + str(self.AObj.noProcessedScan))
        index2 = self.refbox.findText("0")
        
        self.AObj.validateChange()
        self.AObj.storeProcessedToScan()
        self.AObj.storeToBack()
        self.scanbox.setCurrentIndex(index1) 
        self.refbox.setCurrentIndex(index2) 
            
            
    ##### baseline detection
    def initBaselineDetect(self):
        self.initsmoothness()
        self.initAsymmetry()
        self.initNiter()
        self.initbasePre()
        self.initbaseApp()
        self.initBaseCan()
            
    def initsmoothness(self):
        self.smoothness.setText("10")
        
    def initAsymmetry(self):
        self.asymmetry.setText("0.001")
        
    def initNiter(self):
        self.Niter.setValue(10)
        
    def initbasePre(self):
        self.basePre.clicked.connect(self.updateBasePre)
        
    def initbaseApp(self):
        self.baseApp.clicked.connect(self.updateBaseApp)
        
    def initBaseCan(self):
        self.baseCan.clicked.connect(self.updateNormCan)
              
    def updateBasePre(self):
        lam = float(self.smoothness.text())
        p = float(self.asymmetry.text())
        niter = self.Niter.value()
        self.AObj.getPeaks(lam, p, niter)
        self.updateAnalysisGraph()
        self.updateBackBaseline()
        
    def updateBaseApp(self):
        self.scanbox.addItem("Processed Scan " + str(self.AObj.noProcessedScan))
        index1 = self.scanbox.findText("Processed Scan " + str(self.AObj.noProcessedScan))
        
        self.refbox.addItem("Background " + str(self.AObj.noBackground))
        index2 = self.refbox.findText("0")
        
        self.AObj.validateChange()
        self.AObj.storeProcessedToScan()
        self.AObj.storeToBack()
        self.scanbox.setCurrentIndex(index1) 
        self.refbox.setCurrentIndex(index2) 

        
        
        
    ##### Saving box
        
    def initSavingEdit(self):
        self.dateEdit.setText("Enter date")
        self.nameEdit.setText("Enter scan name")
        
    def initSaveButton(self):   
        self.saveButton.clicked.connect(self.onpressSaveButton) 
        
    def onpressSaveButton(self):
        scanDate = str(self.dateEdit.text())
        scanName = str(self.nameEdit.text())
        self.AObj.saveAnalysisData(scanDate, scanName)
        
        

            
        




if __name__ == '__main__':

    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())
