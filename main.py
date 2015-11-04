# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:37:48 2015

@author: pertoty
"""


#  Import libraries and classes

from __future__ import division, print_function
import numpy as np

import sys
from PyQt4.uic import loadUiType
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    

from lib.classLib import scanData, analysisDatas
from lib.funLib import RectangleSelector

Ui_MainWindow, QMainWindow = loadUiType('gui_layout/lh_gui.ui')


# Creation of classes for analysis


class Main(QMainWindow, Ui_MainWindow):

###########################    Initlization part    ###########################
    def __init__(self):
        '''
        Initialization of main window
        '''
        #Initiate main window
        super(Main, self).__init__()
        self.setupUi(self)
        
        #Menu creation
        self.create_menu()
        
        #Creation of figures
            #Images for area selection and image scrolling
        self.fig2d = Figure(facecolor='1')
        self.canvas2d = FigureCanvas(self.fig2d)
            #Spectra over selection
        self.fig1dSpec = Figure(facecolor='1')
        self.canvas1dSpec = FigureCanvas(self.fig1dSpec)
            # Spectrum for wavelength calibration
        self.fig1dSpecCalib = Figure(facecolor='1')
        self.canvas1dSpecCalib = FigureCanvas(self.fig1dSpecCalib)
            #Image for 2d scan (spectrum vs time)
        self.fig2dScan = Figure(facecolor='1')
        self.canvas2dScan = FigureCanvas(self.fig2dScan)
            #Image for 2d analysis
        self.figAnalysis = Figure(facecolor='1')
        self.canvasAnalysis = FigureCanvas(self.figAnalysis)


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
        self.folderName = str(QFileDialog.getExistingDirectory(
                              self, 'Select folder'))

        # Creation of scan object and range for image display (specific to each scan)                  
        self.initScanObject()
        
        # Initilization of widgets related to the Calibration scan object 
        self.initCalibrationTab()
    
        
################################   Calibration tab ############################
        
    def initCalibrationTab(self):
        self.initSliderH()
        self.initSliderV()
        self.initSliderCalib()
        self.initt0Slider()
        self.plot2d()
        self.plot1dSpec()
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
        
    def initScanObject(self):
        self.lhObj = scanData(self.folderName)
        self.vmin = np.amin(self.lhObj.imStackscanTot)
        self.vmax = np.amax(self.lhObj.imStackscanTot)
        self.imNumber.setText(str(self.lhObj.currentImNum)+'/'+ str(self.lhObj.noImage))

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

        self.rect = RectangleSelector(self.axes2d, self._on_mouse_release,
                                      [0, 100, 0, 100], useblit=True)
        self.axes2d.set_xlabel('Pixel number')
        self.axes2d.set_ylabel('Pixel number')
        self.canvas2d.draw()
        self.mpl2dIm.addWidget(self.canvas2d)

    def plot1dSpec(self):
        '''
        Initialization of 1d plot of instantaneous spectrum according rectangle
        selection on image graph.
        '''
        self.fig1dSpec.clear()
        self.axes1dSpec = self.fig1dSpec.add_subplot(111)

        self.im1dSpec = self.axes1dSpec.plot(np.mean(self.lhObj.Im1dSpec,
                                                     axis=0))
        self.canvas1dSpec.draw()
        self.mpl1dSpec.addWidget(self.canvas1dSpec)
        
    def plotCalibSpec(self):
        '''
        Initialization of 1d plot of calibration spectrum with cropped images
        '''
        self.fig1dSpecCalib.clear()
        self.axes1dSpecCalib = self.fig1dSpecCalib.add_subplot(111)

        self.im1dSpecCalib = self.axes1dSpecCalib.plot(self.lhObj.xaxisCalib, \
        self.lhObj.spectrumMean)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[0], color = 'red', linewidth = 2)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[1], color = 'blue', linewidth = 2)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[2], color = 'green', linewidth = 2)
        
        self.canvas1dSpecCalib.draw()
        self.mplSpecCalib.addWidget(self.canvas1dSpecCalib)
        
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
                                               
        self.axes2dScan.axhline(self.lhObj.t0Ref, color = 'white', linewidth = 2)

        self.axes2dScan.set_xlabel('Photon Energy')
        self.axes2dScan.set_ylabel('Time')
        self.canvas2dScan.draw()
        self.scanRaw.addWidget(self.canvas2dScan)

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
        
    def initSliderCalib(self):
        self.setSliderLim()
        self.ref1_slider.valueChanged[int].connect(self.moveRef1)        
        self.ref2_slider.valueChanged[int].connect(self.moveRef2)
        self.ref3_slider.valueChanged[int].connect(self.moveRef3)
        
    def initt0Slider(self):
        self.sett0SliderLim()
        self.t0_slider.valueChanged[int].connect(self.movet0Ref)
        
    def initRefText(self):
        self.ref1.setText("196.2")
        self.ref1.returnPressed.connect(self.readRef1)
        
        self.ref2.setText("184")
        self.ref2.returnPressed.connect(self.readRef2)
        
        self.ref3.setText("173")
        self.ref3.returnPressed.connect(self.readRef3)
        
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

    ################     Updating graph   

        ########### 2d image
    def changeImage(self, value):
        self.lhObj.currentImNum = value
        self.lhObj.changeCurrentIm(value)
        self.updateImage2d()
        self.updateImage1dSpec()
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
        
        ######### 1D spectrum for scrolling

    def updateImage1dSpec(self):
        self.axes1dSpec.cla()
        self.im1dSpec = self.axes1dSpec.plot(np.mean(self.lhObj.Im1dSpec,
                                                     axis=0))
        self.canvas1dSpec.draw()
        
        ########## 1D spectrum for calibration
        
    def updatePlotCalibSpec(self):
        
        self.axes1dSpecCalib.cla()

        self.im1dSpecCalib = self.axes1dSpecCalib.plot(self.lhObj.xaxisCalib, \
        self.lhObj.spectrumMean)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[0], color = 'red', linewidth = 2)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[1], color = 'blue', linewidth = 2)
        self.axes1dSpecCalib.axvline(self.lhObj.specRef[2], color = 'green', linewidth = 2)
        
        self.axes2d.set_xlabel('Photon energy in eV')
        self.axes2d.set_ylabel('Intensity')
        
        self.canvas1dSpecCalib.draw()
        
    def moveRef1(self, value):
        self.lhObj.specRef[0] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
        self.updatePlotCalibSpec()
        
    def moveRef2(self, value):
        self.lhObj.specRef[1] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
        self.updatePlotCalibSpec()
        
    def moveRef3(self, value):
        self.lhObj.specRef[2] = value*(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0])
        self.updatePlotCalibSpec()
        
    def setSliderLim(self):
        xmin = int(np.amin(self.lhObj.xaxisCalib)/(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0]))
        xmax = int(np.amax(self.lhObj.xaxisCalib)/(self.lhObj.xaxisCalib[1] - self.lhObj.xaxisCalib[0]))
        self.ref1_slider.setMinimum(xmin)
        self.ref1_slider.setMaximum(xmax)
        
        self.ref2_slider.setMinimum(xmin)
        self.ref2_slider.setMaximum(xmax)
        
        self.ref3_slider.setMinimum(xmin)
        self.ref3_slider.setMaximum(xmax)
        
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
                                               
        self.axes2dScan.axhline(self.lhObj.t0Ref, color = 'white', linewidth = 2)

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
        
    def sett0SliderLim(self):
        tmin = int(np.amin(self.lhObj.taxis))
        tmax = int(np.amax(self.lhObj.taxis))

        
        self.t0_slider.setMinimum(tmin)
        self.t0_slider.setMaximum(tmax)
        
    def movet0Ref(self, value):
        self.lhObj.t0Ref = value
        self.updatePlot2dScan()

    ####### Widget interaction

    def _on_mouse_release(self, epress=None, erelease=None):
        x0, y0, width, height = self.rect._rect_bbox
        self.lhObj.rectSpec = [y0 + self.lhObj.rect[0],
                                   self.lhObj.rect[0] + y0 + height,
                                   self.lhObj.rect[2] + x0,
                                   self.lhObj.rect[2] + x0 + width]
        self.lhObj.cropIm1dSpec()
        self.updateImage1dSpec()

    def _onpressButtonCrop(self):
        x0, y0, width, height = np.int64(self.rect._rect_bbox)
        self.lhObj.rect = [y0 + self.lhObj.rect[0], self.lhObj.rect[0] + y0 +
                           height, self.lhObj.rect[2] + x0,
                           self.lhObj.rect[2] + x0 + width]

        self.lhObj.cropIm()
        self.updateImage2d() 
        self.setSliderLim()       
        self.updatePlotCalibSpec()
        self.updatePlot2dScan()

    def _onpressButtonInit(self):
        self.lhObj.rect = [0, self.lhObj.imSizeInit[0]-1, 0, self.lhObj.imSizeInit[1]-1]
        self.lhObj.rectSpec = [0, self.lhObj.imSizeInit[0]-1, 0, self.lhObj.imSizeInit[1]-1]
        
        self.lhObj.changeCurrentIm(self.lhObj.currentImNum)
        self.updateImage2d()
        self.updateImage1dSpec()
        
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
        
    def _onpressScanUpdate(self):
        self.readDeltat()
        
        if self.firstCheck == 1:
            self.lhObj.taxis = (self.lhObj.taxis - self.lhObj.t0Ref)*self.lhObj.deltat
        elif self.firstCheck == 0:
            self.lhObj.taxis = (self.lhObj.taxis - self.lhObj.t0Ref)
            
        self.firstCheck = 0
        self.lhObj.t0Ref = 0
        self.sett0SliderLim()
        self.updatePlot2dScan()
        
    def _onpressSaveScan(self):
        self.initAnalysis()
        
 
              
    def readRef1(self):
        self.lhObj.wavelRef[0] = float(self.ref1.text())
        
    def readRef2(self):
        self.lhObj.wavelRef[1] = float(self.ref2.text())
        
    def readRef3(self):
        self.lhObj.wavelRef[2] = float(self.ref3.text())
        
    def readDeltat(self):
        self.lhObj.deltat = float(self.deltat.text())
        

###########################      Analysis tab      ############################
    #### Initialization of analysis tab
    def initAnalysis(self):
        
        if self.lhObj.refState == 0:
            self.AObj = analysisDatas(self.lhObj.taxis, self.lhObj.xaxisCalib,\
            (self.lhObj.imNoRefCrop).T)
        elif self.lhObj.refState == 2:
            self.AObj = analysisDatas(self.lhObj.taxis, self.lhObj.xaxisCalib,\
            (self.lhObj.imScanCrop).T, (self.lhObj.imRefCrop).T)
            
        self.initComboBoxes()
        self.initStoreButton()
        self.initAnalysisGraph()
        self.initSliderNormalization()
        self.initSliderReference()
        self.initNormPreviewButton()
        self.initNormCancelButton()
        self.initNormApplyButton()
        self.initRefPreviewButton()
        self.initRefAppButton()
        self.initAxisBox()
        self.initWindowBox()
        self.initSmoothWindowLen()
        self.initSmoothPreview()
        self.initMovingCancelButton()
        self.initMovingApplyButton()
        self.initSavingEdit()
        self.initSaveButton()
            
    ##### Graph initialization   
    def initAnalysisGraph(self):

        self.figAnalysis.clear()
        self.axesAnalysis = self.figAnalysis.add_subplot(111)
        
        self.axesAnalysis.pcolormesh(self.AObj.energy, 
                                     self.AObj.taxis,                                                                           
                                     self.AObj.displaytemp)
                                     
        self.axesAnalysis.axvline(self.AObj.energy[self.AObj.lowerBoundNorm], color = 'white', linewidth = 2)
        self.axesAnalysis.axvline(self.AObj.energy[self.AObj.upperBoundNorm], color = 'white', linewidth = 2)

        self.axesAnalysis.axhline(self.AObj.taxis[self.AObj.lowerBoundRef], color = 'white', linewidth = 2)
        self.axesAnalysis.axhline(self.AObj.taxis[self.AObj.upperBoundRef], color = 'white', linewidth = 2)
        
        self.mpl_toolbar = NavigationToolbar(self.canvasAnalysis, self)

        self.axesAnalysis.set_xlabel('Photon energy in eV')
        self.axesAnalysis.set_ylabel('Time delay (fs)')
        self.canvasAnalysis.draw()
        self.analysisGraph.addWidget(self.canvasAnalysis)
        self.analysisGraph.addWidget(self.mpl_toolbar)
        
        
    def updateAnalysisGraph(self):
        self.axesAnalysis.cla()
        self.axesAnalysis.pcolormesh(self.AObj.energy, 
                                     self.AObj.taxis,                                                                           
                                     self.AObj.displaytemp)
                                     
        self.axesAnalysis.axvline(self.AObj.energy[-self.AObj.lowerBoundNorm], color = 'white', linewidth = 2)
        self.axesAnalysis.axvline(self.AObj.energy[-self.AObj.upperBoundNorm], color = 'white', linewidth = 2)

        self.axesAnalysis.axhline(self.AObj.taxis[self.AObj.lowerBoundRef], color = 'white', linewidth = 2)
        self.axesAnalysis.axhline(self.AObj.taxis[self.AObj.upperBoundRef], color = 'white', linewidth = 2)

        self.axesAnalysis.set_xlabel('Photon energy in eV')
        self.axesAnalysis.set_ylabel('Time delay (fs)')
        self.canvasAnalysis.draw()
        
    ##### Widgets initialization   
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
        
    def initSliderNormalization(self):
        self.setSliderAnalysisLim()
        self.lowBoundNorm.valueChanged[int].connect(self.setNormLowBound)
        self.upBoundNorm.valueChanged[int].connect(self.setNormUpBound)

    def initSliderReference(self):
        self.setSliderAnalysisLim()
        self.lowBoundRef.valueChanged[int].connect(self.setRefLowBound)
        self.upBoundRef.valueChanged[int].connect(self.setRefUpBound)
        
    def setSliderAnalysisLim(self):
        self.lowBoundNorm.setMinimum(0)
        self.lowBoundNorm.setMaximum(len(self.AObj.energy)-1)
        
        self.upBoundNorm.setMinimum(0)
        self.upBoundNorm.setMaximum(len(self.AObj.energy)-1)
        
        self.lowBoundRef.setMinimum(0)
        self.lowBoundRef.setMaximum(len(self.AObj.taxis)-1)
        
        self.upBoundRef.setMinimum(0)
        self.upBoundRef.setMaximum(len(self.AObj.taxis)-1)
        
    def initNormPreviewButton(self):
        self.normalizationPre.clicked.connect(self.updateNormPreview)
        
    def initNormCancelButton(self):
        self.normalizationCan.clicked.connect(self.updateNormCan)
        
    def initNormApplyButton(self):
        self.normalizationApp.clicked.connect(self.updateNormApp)
        
    def initRefPreviewButton(self):
        self.referencePre.clicked.connect(self.updateRefPreview)
        
    def initRefCancelButton(self):
        self.referenceCan.clicked.connect(self.updateRefCan)
        
    def initRefAppButton(self):
        self.referenceApp.clicked.connect(self.updateRefApp)   
        
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
        self.windowLenEdit.setText(str(self.AObj.windowLen))
        self.windowLenEdit.returnPressed.connect(self.setSmoothWindowLen)
        
    def initSmoothPreview(self):
        self.movingPre.clicked.connect(self.updateSmoothGraph)   
        
    def initMovingCancelButton(self):
        self.movingCan.clicked.connect(self.updateNormCan)
        
    def initMovingApplyButton(self):
        self.movingApp.clicked.connect(self.updateNormApp)   
        
    def initSavingEdit(self):
        self.dateEdit.setText("Enter date")
        self.nameEdit.setText("Enter scan name")
        
    def initSaveButton(self):   
        self.saveButton.clicked.connect(self.onpressSaveButton) 
        
        
        
    
    ##### callback functions for widgets    
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
        self.refbox.addItem("Processed ref " + str(self.AObj.noProcessedScan))
        self.AObj.storeProcessedToRef()
        
    def _onpressnormStore(self):
        self.normbox.addItem("Processed norm " + str(self.AObj.noProcessedScan))
        self.AObj.storeProcessedToNorm()
        
    def updateNormPreview(self):
        self.AObj.normalization()
        self.updateAnalysisGraph()
        
    def updateRefPreview(self):
        self.AObj.createRef()
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
        
    def updateRefCan(self):
        self.AObj.displaytemp = self.AObj.display
        self.updateAnalysisGraph()
        
    def updateRefApp(self):
        self.AObj.display = self.AObj.displaytemp
        self.refbox.addItem("Processed ref " + str(self.AObj.noProcessedRef))
        index = self.refbox.findText("Processed ref " + str(self.AObj.noProcessedRef))
        self.AObj.storeProcessedToRef()
        self.refbox.setCurrentIndex(index)
        
    def setNormLowBound(self, value):
        self.AObj.lowerBoundNorm = int(value)
        self.updateAnalysisGraph()
        
    def setNormUpBound(self, value):
        self.AObj.upperBoundNorm = int(value)
        self.updateAnalysisGraph()
        
    def setRefLowBound(self, value):
        self.AObj.lowerBoundRef = int(value)
        self.updateAnalysisGraph()
        
    def setRefUpBound(self, value):
        self.AObj.upperBoundRef = int(value)
        self.updateAnalysisGraph()
        
    def setSmoothAxis(self, text):
        self.AObj.axisSmoothCurrent = self.AObj.axisSmoothList[str(text)]
        
    def setSmoothWindow(self, text):
        self.AObj.windowType = str(text)
        
    def setSmoothWindowLen(self):
        self.AObj.windowLen = float(self.windowLenEdit.text())
        
    def updateSmoothGraph(self):
        self.AObj.windowLen = float(self.windowLenEdit.text())
        
        self.AObj.movingAverage()
        self.updateAnalysisGraph()
        
    def onpressSaveButton(self):
        scanDate = str(self.dateEdit.text())
        scanName = str(self.nameEdit.text())
        self.AObj.saveAnalysisData(scanDate, scanName)
        
        
        
        
        

            
        




if __name__ == '__main__':

    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())
