import sys,os
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication,QWidget,QPlainTextEdit, QMessageBox,QAction,QMenu, QActionGroup,QToolButton,QToolBar
from PyQt5.QtWidgets import QScrollArea,QVBoxLayout, QAbstractItemView, QListWidget,QListWidgetItem, QPushButton, QFileDialog,QDialog
from PyQt5.QtCore import pyqtSlot
from PyQt5 import uic,QtGui
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import shapely
import time
from copy import copy


from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import pickle #import/export text files
import skimage as ski

import json



#set proper current directory
current_dir=os.path.abspath(os.path.dirname(sys.argv[0]))
os.chdir(current_dir)

uifile_1 = "plot_builder.ui" # Main window
uifile_2 = "Block_props.ui" # Block properties
uifile_3 = "physics.ui" # Enter file here.

form_1, base_1 = uic.loadUiType(uifile_1)
form_2, base_2 = uic.loadUiType(uifile_2)
form_3, base_3 = uic.loadUiType(uifile_3)

#/Users/geolog/Documents/Geology/myCDF/intrusion_emplacement_gerya/
# py/25_12_2022_Parallel/test_joblib_speedup/КРАУНЦ23/py/github_for_jit/cython/
# C_numpy_wrapper/model_data_process_numpy_c/new2/Intrusion_emplacement2.pyx 
#

class Example(base_1, form_1):
    def __init__(self):


        super().__init__()
        self.setupUi(self)
        self.layout = QVBoxLayout(self)
                
        #block parameters window
        self.block_params=[]
        
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # set the layout
        #layout = QVBoxLayout()
        #layout.addWidget(self.toolbar)
        #layout.addWidget(self.canvas)
        self.scrollArea.setWidget(self.canvas)
        #self.setLayout(layout)
         
        #block list management
        self.list1.setSelectionMode(QAbstractItemView.SingleSelection)
        self.list1.itemSelectionChanged.connect(self.on_SelectFromList)
        
        ###
        self.modelWidth=int(self.mWtxt.text());
        self.modelHeight=int(self.mHtxt.text());

        #physics defaults
        self.modelPhysics_ini = {
            "heightpx":50,
            "wigthpx":150,
            "heightkm":50,
            "widthkm":150,
            "gx":0,
            "gy":9.81,
            "ttop":273, #K
            "tbottom":1600, #K
            "tintrus":1700, #Ks
            "melting":True,
            "tmode":'stationary', # stationary, extension collision subduction
            "name":'',
            'horcolvel':2*3.168808781402895e-10,
            'horextvel':2*3.168808781402895e-10,
            'horextstrainpow':-14
            }

        self.settings_fname = 'settings.p'
        if not os.path.isfile(self.settings_fname):
            self.modelPhysics = copy(self.modelPhysics_ini)
            self.saveSettings(self.modelPhysics,self.settings_fname)
        else:
            self.modelPhysics = self.loadSettings(self.settings_fname)
            #TODO implement setting of the model name when model is loaded
            if 'name' in self.modelPhysics: self.leModelName.setText(self.modelPhysics['name'])

        print(self.modelPhysics)

        #pixel model and block array
        self.model=np.zeros([self.modelHeight,self.modelWidth])
        self.bm=[] #model dictionary
        #self.timeout = 0.05 #timeout when drawing
        
        ###show initial picture
        self.figure.clear()
        
        #self.figure.tight_layout();
        
        self.ax = self.figure.add_subplot(111)
        self.ax.cla() #or picture will not be updated
        
        #self.figure.tight_layout() #must be called twice
        
        #set label size to avoid setting it twice
        self.ax.tick_params(axis='both', which='major', labelsize=4)
        self.ax.imshow(self.model, interpolation='nearest');
        #self.drawCanvas()
        self.canvas.draw() 
       
        #click on Screen
        self.event_listener1=self.canvas.mpl_connect('button_press_event', self.onClick)
        self.event_listener2=self.canvas.mpl_connect('button_release_event', self.onRelease)
        self.event_listener3=[]
        self.event_listener4=self.canvas.mpl_connect('motion_notify_event', self.onHover)
        #TODO make fileOpen and saveModelBut buttons functionality
        #self.fileOpen.clicked.connect(self.file_open_dialogue)
        self.removeSelBtn.clicked.connect(self.removeSel)
        self.moveUpBtn.clicked.connect(self.moveUp)
        self.moveDnBtn.clicked.connect(self.moveDn)
        self.clearModelBut.clicked.connect(self.clearModel)
        self.applyBtn.clicked.connect(self.applyNewDimensions)
        self.saveModelBut.clicked.connect(self.exportModel)
        self.addFullSizeBtn.clicked.connect(self.addFullsizeBlock)
        #self.buildBut.clicked.connect(self.build_plot)
        self.df=pd.DataFrame() #global dataframe
        
        
        # getting slider moved signal
        self.sldCurSize.sliderMoved.connect(lambda: do_action()) #process slider moving
        self.sldCurSize.setRange(1, 30) #set qslider minmax
        self.sldCurFreq.sliderMoved.connect(lambda: do_action2()) #process slider moving
        self.sldCurFreq.setRange(5, 500) #set qslider minmax
 
        # method called when signal is emitted
        def do_action(): #for qslider
             # setting text to the label
            self.leCurSize.setText(str(self.sldCurSize.value()))
        def do_action2(): #for qslider
             # setting text to the label
            self.leCurFreq.setText(str(self.sldCurFreq.value()))
        
        self.leCurSize.setText(str(self.sldCurSize.value()))
        self.leCurFreq.setText(str(self.sldCurFreq.value()))
        
        #update slider value when line edit text is changed
        self.leCurSize.textChanged.connect(lambda:onChange())
        self.leCurFreq.textChanged.connect(lambda:onChange2())
        self.physicsBtn.clicked.connect(lambda:setPhysics())

        def onChange():
            try:
                self.sldCurSize.setValue(int(self.leCurSize.text()))
            except ValueError:
                self.sldCurSize.setValue(1)
        def onChange2():
            try:
                self.sldCurFreq.setValue(int(self.leCurFreq.text()))
            except ValueError:
                self.sldCurFreq.setValue(1)

        def setPhysics():
            #TODO show up physics window
            self.physicsWindow = EnterPhysics()
            self.physicsWindow.show()
            #TODO set widget values from 
            #self.modelPhysics
            '''
            "heightpx"
            "wigthpx"
            "heightkm":50,
            "widthkm":150,
            "gx":0,
            "gy":9.81,
            "ttop":273, #K
            "tbottom":1600, #K
            "tintrus":1700, #Ks
            "melting":True,
            "tmode":'stationary' 
            '''
            mf = self.modelPhysics
            if 'gx' in mf: self.physicsWindow.sb_gx.setValue(mf['gx'])
            if 'gy' in mf: self.physicsWindow.sb_gy.setValue(mf['gy'])
            if 'heightkm' in mf: self.physicsWindow.sb_heightkm.setValue(mf['heightkm'])
            if 'widthkm' in mf: self.physicsWindow.sb_widthkm.setValue(mf['widthkm'])
            if 'seadepth' in mf: self.physicsWindow.sb_seadepth.setValue(mf['seadepth'])
            if 'ttop' in mf: self.physicsWindow.sb_Ttop.setValue(mf['ttop'])
            if 'tbottom' in mf: self.physicsWindow.sb_Tbottom.setValue(mf['tbottom'])
            if 'tintrus' in mf: self.physicsWindow.sb_TIntrusion.setValue(mf['tintrus'])
            if 'melting' in mf: self.physicsWindow.chbMelt.setChecked(bool(mf['melting']))

            if 'tmoho' in mf: self.physicsWindow.sb_TMoho.setValue(mf['tmoho'])
            if 'ymohokm' in mf: self.physicsWindow.sb_YMoho.setValue(mf['ymohokm'])
            if 'lithbottom' in mf: self.physicsWindow.sb_YLithBottom.setValue(mf['lithbottom'])

            if 'horcolvel' in mf: self.physicsWindow.sb_horColVel.setValue(mf['horcolvel']/3.168808781402895e-10)
            if 'horextvel' in mf: self.physicsWindow.sb_horExtVel.setValue(mf['horextvel']/3.168808781402895e-10)
            if 'lithbottom' in mf: self.physicsWindow.sb_YLithBottom.setValue(mf['lithbottom'])

            if 'tmode' in mf:
                if mf['tmode']=='stationary':
                    self.physicsWindow.rbStable.setChecked(True)
                elif mf['tmode']=='extension':
                    self.physicsWindow.rbExtension.setChecked(True)
                elif mf['tmode']=='collision':
                    self.physicsWindow.rbShort.setChecked(True)  
                elif mf['tmode']=='subduction':
                    self.physicsWindow.rbSubduction.setChecked(True)        
            del mf #remove pointer
        
        #TODO add menu bar
        self._createActions()
        self._createMenuBar()
        self._createToolBars()




    def restorePhysics(self):
            
            mf = self.modelPhysics_ini
            if 'gx' in mf: self.physicsWindow.sb_gx.setValue(mf['gx'])
            if 'gy' in mf: self.physicsWindow.sb_gy.setValue(mf['gy'])
            if 'heightkm' in mf: self.physicsWindow.sb_heightkm.setValue(mf['heightkm'])
            if 'widthkm' in mf: self.physicsWindow.sb_widthkm.setValue(mf['widthkm'])
            if 'seadepth' in mf: self.physicsWindow.sb_seadepth.setValue(mf['seadepth'])
            if 'ttop' in mf: self.physicsWindow.sb_Ttop.setValue(mf['ttop'])
            if 'tbottom' in mf: self.physicsWindow.sb_Tbottom.setValue(mf['tbottom'])
            if 'tintrus' in mf: self.physicsWindow.sb_TIntrusion.setValue(mf['tintrus'])
            if 'melting' in mf: self.physicsWindow.chbMelt.setChecked(bool(mf['melting']))

            if 'tmoho' in mf: self.physicsWindow.sb_TMoho.setValue(mf['tmoho'])
            if 'ymohokm' in mf: self.physicsWindow.sb_YMoho.setValue(mf['ymohokm'])
            if 'lithbottom' in mf: self.physicsWindow.sb_YLithBottom.setValue(mf['lithbottom'])

            if 'tmode' in mf:
                if mf['tmode']=='stationary':
                    self.physicsWindow.rbStable.setChecked(True)
                elif mf['tmode']=='extension':
                    self.physicsWindow.rbExtension.setChecked(True)
                elif mf['tmode']=='collision':
                    self.physicsWindow.rbShort.setChecked(True)        
            del mf #remove pointer

        #self.poly_list = [] #list of polygons for graphical input
    #some flags
    is_drawing = False #user is not drawing currently  
    
    #TODO saveSettings
    def saveSettings(self,mf,fname='defaults.json'):
        info = {}
        with open(fname,mode='w',encoding='utf-8') as f:
            f.write(json.dumps(mf,ensure_ascii=False,indent=4))

    #TODO loadSettings
    def loadSettings(self,fname='defaults.json'):
        with open(fname,encoding='utf-8') as f:
                mf = json.loads(f.read())
        return mf
    

    def exportModel(self):
        fileName = QFileDialog.getSaveFileName(self,("Save model into file..."),'Model',("(*.cfg)"))
        #update some values from text fields
        self.modelPhysics.update({"name":self.leModelName.text()})        
        myvars=(self.bm,self.model)
        #try:
            #with open(fileName[0],"wb") as f:
            #    pickle.dump( myvars, f)
        self.exportModelJson(self.bm,fileName[0])
        #except:
        #    print('unknown error')

    def exportModel_old(self):
        fileName = QFileDialog.getSaveFileName(self,("Save model into file..."),'Model',("(*.cfg)"))
                
        myvars=(self.bm,self.model)
        try:
            with open(fileName[0],"wb") as f:
                pickle.dump(myvars, f)
        except:
            print('no file name provided')

    #TODO exportModelJson 
    def exportModelJson(self,bm,filename):
        #TODO 1) iterate over blocks in the model 2) create list of block dictionaries 3) arrange model_out of blocks and physics 
        block_dict_list = []
        for block in bm:
            block_dict = dict()
            block_dict.update({'blockname':block[0]})
            block_dict.update({'rectangle_ullr':block[1]})
            block_dict.update({'points':block[2]})
            block_dict_list.append(block_dict)


        #TODO JSON export that is the model file for model export into json
        model_out = {
            "name":self.modelPhysics['name'],
            "description":"model situation",
            "sizepx":[self.modelHeight,self.modelWidth],
            "sizekm":[self.modelPhysics['heightkm'],self.modelPhysics['widthkm']],
            "blockunits":"percents",
            "blocks": block_dict_list,
            "gravity":[self.modelPhysics['gx'],self.modelPhysics['gy']],
            "ttop":self.modelPhysics['ttop'],
            "tbottom":self.modelPhysics['tbottom'],
            'tintrus': self.modelPhysics['tintrus'],
            'melting': self.modelPhysics['melting'],
            'tmode':self.modelPhysics['tmode'],
            'horcolvel':self.modelPhysics['horcolvel'],
            'horextvel':self.modelPhysics['horextvel'],
            'horextstrainpow':self.modelPhysics['horextstrainpow'],
            'seadepth':self.modelPhysics['seadepth'],

            'tmoho':self.modelPhysics['tmoho'],
            'ymohokm':self.modelPhysics['ymohokm'],
            'lithbottom':self.modelPhysics['lithbottom']   #km
        }
        #add all remaining keys 
        for key in self.modelPhysics:
            if key not in model_out:
                model_out.update({key:self.modelPhysics[key]})

        #запись в json
        try:
            with open(filename,mode='w',encoding='utf-8') as f:
                f.write(json.dumps(model_out,ensure_ascii=False,indent=4))
        except FileNotFoundError:
            print('file not found')




    def importModelJson(self,filename):
        #TODO importModelJsons
        #читаем из json
        info_read = []
        with open(filename,encoding='utf-8') as f:
            info_read = json.loads(f.read())

        self.modelPhysics.update({'name':info_read['name']})
        self.modelPhysics.update({'description':info_read['description']})
        self.modelPhysics.update({'sizepx':[info_read['sizepx'][0],info_read['sizepx'][1]]})
        self.modelPhysics.update({'heightkm':info_read['sizekm'][0]})
        self.modelPhysics.update({'widthkm':info_read['sizekm'][1]})
        self.modelPhysics.update({'blockunits':info_read['blockunits']})
        self.modelPhysics.update({'gx':info_read['gravity'][0]})
        self.modelPhysics.update({'gy':info_read['gravity'][1]})
        self.modelPhysics.update({'ttop':info_read['ttop']})
        self.modelPhysics.update({'tbottom':info_read['tbottom']})
        self.modelPhysics.update({'tmode':info_read['tmode']})

        self.bm = []
        for block in info_read['blocks']:
            self.bm.append([block['blockname'],block['rectangle_ullr'],block['points']])

        QMessageBox.information(self,'Model is loaded','The model file has been loaded successfully!',QMessageBox.Ok)        

        self.mHtxt.setText(str(self.modelPhysics['sizepx'][0]))
        self.mWtxt.setText(str(self.modelPhysics['sizepx'][1]))
        self.leModelName.setText(self.modelPhysics['name']) 
        self.modelResize()
        #self.drawCanvas()
        #self.drawModel()


    def modelResize(self):
        self.modelWidth=int(self.mWtxt.text());
        self.modelHeight=int(self.mHtxt.text());
        #self.model=np.zeros([self.modelHeight,self.modelWidth])
        #self.bm=[] #set initial values
        #show initial picture
        self.drawModel()
        self.drawCanvas()

    def applyNewDimensions(self):
        msgBox=QMessageBox.question(self, 'Info','Model should be REDRAWN to \
                                    apply new size. Are you agree?',
                                    QMessageBox.Yes|QMessageBox.No, QMessageBox.No)
        if msgBox==QMessageBox.Yes:
            self.modelResize()

            
    def clearModel(self):
        msgBox=QMessageBox.question(self, 'Info','Model should be RESETTED \
                                    Are you agree?',
                                    QMessageBox.Yes|QMessageBox.No, QMessageBox.No)
        if msgBox==QMessageBox.Yes:
            self.modelWidth=int(self.mWtxt.text());
            self.modelHeight=int(self.mHtxt.text());
            self.model=np.zeros([self.modelHeight,self.modelWidth])
            self.bm=[] #set initial values
            #show initial picture
            self.drawCanvas()
            #clear list
            self.list1.clear()
            self.leModelName.setText('')
    
    def exitApp(self):
        msgBox=QMessageBox.question(self, 'Info','Are you sure \
        want program to terminate?',
                                    QMessageBox.Yes|QMessageBox.No, QMessageBox.No)
        if msgBox==QMessageBox.Yes:
            app.quit()
    
    #returns patch object for overlay
    def get_cursor(self,ixr:int,iyr:int,size:int, shape:str='undef')->patches:
        if shape=='brush':
            rect = patches.Circle((ixr,iyr),size,linewidth=0.25,edgecolor='r',facecolor='none')
        elif shape=='pencil':
            rect = patches.Rectangle((ixr-size,iyr-size),size*2,size*2,linewidth=0.25,edgecolor='r',facecolor='none')
        else:
            rect = patches.Circle((ixr,iyr),size,linewidth=0.25,edgecolor='r',facecolor='r')
        
        return rect

    def removeSel(self):
        idx=self.list1.currentRow()
        del self.bm[idx]
        
        self.drawCanvas()
        self.drawModel()
        

    def moveUp(self):
        idx=self.list1.currentRow()
        self.changeOrder(self.bm,idx,idx-1)

    def moveDn(self):
        idx=self.list1.currentRow()
        self.changeOrder(self.bm,idx,idx+1)

    def changeOrder(self,somelist,oldidx,newidx):
        #TODO swap objects with indeces
        item2swap = somelist[oldidx]
        somelist.pop(oldidx)
        somelist.insert(newidx,item2swap)
        self.drawModel()
        self.drawCanvas()
    
    def addFullsizeBlock(self):
        global ixp, iyp, ixr, iyr
        self.poly_list = [] #create list of cursor polygons 
        #self.modelWidth*0.005 and np.abs(iyr-iyp)>self.modelHeight
        ixp,iyp = 0,0
        ixr, iyr = self.modelWidth,self.modelHeight
        self.block_params=EnterSource()
        self.block_params.show();
    
    def drawCanvas(self):
        self.figure.clear()
        
        #self.figure.tight_layout();

        #show initial picture
        self.ax = self.figure.add_subplot(111)
        self.ax.cla() #or picture will not be updated
        self.figure.tight_layout(); #must be called twice
        self.ax.tick_params(axis='both', which='major', labelsize=4)
        self.ax.imshow(self.model, interpolation='nearest')
        self.canvas.draw()  
    
    def onClick(self,event): #click on app window
        global ixp, iyp
        self.poly_list = [] #create list of cursor polygons 
        self.event_listener3=self.canvas.mpl_connect('motion_notify_event', self.onMove)
        ixp, iyp = int(event.xdata), int(event.ydata)
        if self.rbRegion.isChecked():  #если включен режим добавления блоков
            #glue start edges to model boundaries
            if ixp<self.modelWidth*0.1:
                ixp=0
            if iyp<self.modelHeight*0.1:
                iyp=0
            print ('x_p ='+ str(ixp)+'y_p='+str(iyp))   
        else:
            self.is_drawing = True    
    
    def getCursorPoly(self,size,ixr, iyr):
        
        if self.rbBrush.isChecked():
            poly = shapely.Point(ixr,iyr).buffer(size)
        #elif self.rbPencil.isChecked():
        else: #self.rbPencil.isChecked()
            x = [ixr-size, ixr-size, ixr+size, ixr+size, ixr-size]
            y = [iyr-size, iyr+size, iyr+size, iyr-size, iyr-size]
            poly = shapely.Polygon(zip(x,y))
        # else:
        #     x = [ixr-size[0], ixr-size[0], ixr+size[0], ixr+size[0], ixr-size[0]]
        #     y = [iyr-size[1], iyr+size[1], iyr+size[1], iyr-size[1], iyr-size[1]]
        #     poly = shapely.Polygon(zip(x,y))

        return poly

    #return a patch to show on screen
    def getCursorPatch(self,size)->patches:
        # Create a Cursor patch
        if self.rbBrush.isChecked():
            rect = self.get_cursor(ixr=ixr,iyr=iyr,size=size, shape='brush')
        elif self.rbPencil.isChecked():
            rect = self.get_cursor(ixr=ixr,iyr=iyr,size=size, shape='pencil')
        else:
            rect = self.get_cursor(ixr=ixr,iyr=iyr,size=size, shape='cross')
        return rect


    def graphicInput(self,event):
        global ixr, iyr
        
        time.sleep(1/int(self.leCurFreq.text())) #sleep 0.1 s to avoid hanging

        if event.xdata and event.ydata:
            ixr, iyr = int(event.xdata), int(event.ydata)
            #glue edges to model boundaries
            if np.abs(ixr-self.modelWidth)<self.modelWidth*0.05:
                ixr=self.modelWidth-1
            if np.abs(iyr-self.modelHeight)<self.modelHeight*0.05:
                iyr=self.modelHeight-1
            print ('x_r ='+ str(ixr)+'y_r='+str(iyr))
        
            self.drawCanvas()
            
            if self.rbRegion.isChecked():  #если включен режим добавления блоков
                # Create a Rectangle patch
                rect = patches.Rectangle((ixp,iyp),ixr-ixp,iyr-iyp,linewidth=1,edgecolor='r',facecolor='none')
                # Add the patch to the Axes
                self.ax.add_patch(rect)
                self.canvas.draw()

            else:
                # Create a Cursor patch
                try:
                    size = int(self.leCurSize.text())/2
                except:
                    self.leCurSize.setText(str(1))
                rect = self.getCursorPatch(size)
                
                # Add the patch to the Axes (when MB is pressed)
                self.ax.add_patch(rect)
                self.canvas.draw()    

                #return poly
                poly =  self.getCursorPoly(size,ixr,iyr)
                self.poly_list.append(poly)
                print('new len of poly list is',len(self.poly_list))
                print('poly',poly)

            if self.is_drawing == True:
                print('ixr, iyr=',ixr, iyr) #points log

        
    def onMove(self,event): #move after click on app window with LMB pressed
        self.graphicInput(event)



    def onHover(self,event): #click on app window
        global ixr, iyr
        
        if event.xdata and event.ydata: #inside figure
            ixr, iyr = int(event.xdata), int(event.ydata)
            self.drawCanvas()
            #TODO create a canvas patch for pan and pencil
            if not self.rbRegion.isChecked():
                try:
                    size = int(self.leCurSize.text())/2
                except:
                    self.leCurSize.setText(str(1))

                # Create a Cursor patch
                rect = self.getCursorPatch(size)
                
                # Add the patch to the Axes (when MB is released)
                self.ax.add_patch(rect)
                self.canvas.draw()  
        
        else:   #outside figure (to remove patch)
            self.drawCanvas()
  
            
    def onRelease(self,event): #click on app window
        global ixp, iyp, ixr, iyr
        
        #call graphic input at least once
        self.graphicInput(event)

        ixr, iyr = int(event.xdata), int(event.ydata)
        self.canvas.mpl_disconnect(self.event_listener3)

        print ('x_p ='+ str(ixp)+'y_p='+str(iyp))
        print ('x_r ='+ str(ixr)+'y_r='+str(iyr))
        
        if self.rbRegion.isChecked(): #если включен режим добавления блоков
            #если размер блока больше 1.5% по одному из направлений, выводим окно добавления блока
            if np.abs(ixr-ixp)>self.modelWidth*0.015 and np.abs(iyr-iyp)>self.modelHeight*0.015: 
                
                #glue
                if np.abs(ixr-self.modelWidth)<self.modelWidth*0.05:
                    ixr=self.modelWidth-1
                if np.abs(iyr-self.modelHeight)<self.modelHeight*0.05:
                    iyr=self.modelHeight-1

                self.block_params=EnterSource()
                self.block_params.show();
        #check if pencil or brush drawing is ended
        if self.rbPencil.isChecked() or self.rbBrush.isChecked():
            print('End of brush or pencil output is detected')
            self.block_params=EnterSource()
            self.block_params.show();
            

        print(self.bm)
        #turn off drawing
        self.is_drawing = False

        #TODO add adding sum of polygons to bm element to EnterSource()
        #self.poly_list.clear() #clear poly_list 
        
    def file_open_dialogue(self):
        fileName = QFileDialog.getOpenFileName(self,("Open File"),'',("(*.cfg)"))
        
        if fileName[0] != '':
            
            # #info_read = []
            # #with open(fileName[0],encoding='utf-8') as f:
            #     info_read = json.loads(f.read())


            # #(self.bm,self.model)=pickle.load(open( fileName[0], "rb" ))
            # #print(self.bm,np.shape(self.model));
            # self.mWtxt.setText(str(np.shape(self.model)[1]));
            # self.mHtxt.setText(str(np.shape(self.model)[0]));
            
            # self.modelWidth=np.shape(self.model)[1];
            # self.modelHeight=np.shape(self.model)[0];
            #self.filePathText.setText(fileName[0]);
            self.setWindowTitle(fileName[0])
            self.importModelJson(fileName[0])
            self.drawModel()
            self.drawCanvas()    
        
            
    def on_SelectFromList(self):
        #print([item.text() for item in self.list1.selectedItems()])
        print(self.bm)
        print(self.list1.currentRow()) #index of the selected row
        self.drawCanvas()
        
        # Create a Rectangle patch      
        #yStartPer,xStartPer,yEndPer,xEndPer
        idx=self.list1.currentRow()
        try:
            xc=self.bm[idx][1][1]*(self.modelWidth*0.01)
            yc=self.bm[idx][1][0]*(self.modelHeight*0.01)
            width=(self.bm[idx][1][3]-self.bm[idx][1][1])*(self.modelWidth*0.01);
            height=(self.bm[idx][1][2]-self.bm[idx][1][0])*(self.modelHeight*0.01)
            #block model structure 
            #[[['transition zone'], [36, 51, 72, 77]], [['upper mantle'], [10, 5, 64, 42]]]

            rect = patches.Rectangle((xc,yc),width,height,linewidth=0.5,
                                    edgecolor='r',facecolor='none',
                                    hatch='/',linestyle='--')
            # Add the patch to the Axes
            self.ax.add_patch(rect)
            self.canvas.draw()   
        except:
            print('List items obj does not exist. Already deleted?')

    def getIndexReversed(self,idx:int,mlist:list)->int:
        idx_rev = len(mlist)-1-idx
        return idx_rev    
        
    def drawModel(self):
        print("start drawing a model")
        #reinitialize model
        self.model=np.zeros([self.modelHeight,self.modelWidth])
        xproc=self.modelWidth/100; yproc=self.modelHeight/100;  #percent values

        labels=['Air/water',
                'Sediments',
                'Sediments(m)',
                'Sediments(rm)',
                'Sediments(f)',
                'Basalts',
                'Gabbro',
                'Lithospheric mantle',
                'Asthenospheric mantle',
                'Hydrated mantle',
                'Upper continental crust (granodiorite)',
                'Lower continental crust (diorite)',
                'Basalt(melt)',
                'Diorite(melt)',
                'Granodiorite(melt)']

        for block in self.bm[::-1]:               #self.bm blocks of model REVERSE ORDER
            print('looking for:')
            print(block[0][0])
            try:
                value=labels.index(block[0][0])+1
            except ValueError:
                value='unknown label'
             #idx=[block[1][0]*yproc,block[1][1]*xproc,block[1][2]*yproc,block[1][3]*xproc]
             #idx=np.uint16(np.round(idx))
            for obj in block[2]:
                yy = [y*yproc for y in obj[1]] #convertto pixels
                xx = [x*xproc for x in obj[0]] 
                #add polygon to model 
                rr, cc = ski.draw.polygon(yy, xx)
                try:
                    self.model[rr,cc]=value
                except:
                    print('can not draw, coordinate network error')

             #paint model pixels into colors of block
             #TODO draw polygons instead (by points) like patches
             #self.model[idx[0]:idx[2]+1,idx[1]:idx[3]+1]=value
        
        print(self.model)
        
        self.drawCanvas()
        # #fill list on Window
        
        self.list1.clear() #clear items before 
        
        for block in self.bm:
            item1 = QListWidgetItem() #need to copy theese items twice
            item1.setText(str(block[0][0]))
            self.list1.addItem(item1)

    def drawModel_old(self):
        print("start drawing a model")
        #reinitialize model
        self.model=np.zeros([self.modelHeight,self.modelWidth])
        xproc=self.modelWidth/100; yproc=self.modelHeight/100;  #percent values

        #TODO need to [['transition zone'], [36, 51, 72, 77]]
        #[['transition zone'], [[36, 72], [51, 77]]]  #block name pairs of points

        for block in self.bm:               #self.bm blocks of model
             if block[0][0]=='air':
                 value=0 
             elif block[0][0]=='water':
                 value=1
             elif block[0][0]=='oceanic litho':
                 value=2
             elif block[0][0]=='continental litho':
                 value=3
             elif block[0][0]=='upper mantle':
                 value=4
             elif block[0][0]=='transition zone':
                 value=5
             elif block[0][0]=='lower mantle':
                 value=6
             elif block[0][0]=='fault':
                 value=7
             idx=[block[1][0]*yproc,block[1][1]*xproc,block[1][2]*yproc,block[1][3]*xproc]
             idx=np.uint16(np.round(idx))

             #paint model pixels into colors of block
             #TODO draw polygons instead (by points) like patches
             self.model[idx[0]:idx[2]+1,idx[1]:idx[3]+1]=value
        
        print(self.model)
        
        self.drawCanvas()
        # #fill list on Window
        
        self.list1.clear() #clear items before 
        
        for block in self.bm:
            item1 = QListWidgetItem() #need to copy theese items twice
            item1.setText(str(block[0][0]))
            self.list1.addItem(item1)    
    
    def _createActions(self):
        self.newAction = QAction(self)
        self.newAction.setText('New model')
        self.openAction = QAction('Open model', self)
        #self.saveAction = QAction('Save model', self)
        self.saveAsAction = QAction('Save model as', self)
        self.exitAction = QAction('Exit', self)
        self.viewPreferences = QAction('Settings', self)
        self.aboutAction = QAction('About', self)

        #TODO connect actions and processing functions
        self.newAction.triggered.connect(self.clearModel)
        self.openAction.triggered.connect(self.file_open_dialogue)
        #self.saveAction.triggered.connect(self.exportModel)
        self.saveAsAction.triggered.connect(self.exportModel)
        # self.viewProjectFiles.triggered.connect(self.showLayersPane)
        # self.viewPreferences.triggered.connect(self.showPreferences)
        self.aboutAction.triggered.connect(self.about_dialogue)
        self.exitAction.triggered.connect(self.exitApp)
        

    def _createMenuBar(self):
        # menu bar on top of the window
        self.menuBar = self.menuBar()
        self.fileMenu = QMenu('&File', self)
        self.menuBar.addMenu(self.fileMenu)
        self.fileMenu.addAction(self.newAction)
        self.fileMenu.addAction(self.openAction)
        #self.fileMenu.addAction(self.saveAction)
        self.fileMenu.addAction(self.saveAsAction)
        self.fileMenu.addAction(self.exitAction)
        # self.viewMenu = self.menuBar.addMenu("&View")
        # self.viewMenu.addAction(self.viewProjectFiles)
        # self.viewMenu.addAction(self.viewPreferences)
        # self.menuLanguage = self.menuBar.addMenu('&Language')
        self.helpMenu = self.menuBar.addMenu('&Help')
        self.helpMenu.addAction(self.aboutAction)
        self.menuBar.setNativeMenuBar(False)
    
    def _createToolBars(self):
        # toolbars
        fileToolBar = self.addToolBar("File")
        editToolBar = QToolBar("Edit", self)

        # buttons
        new_project_button = QToolButton(self)
        new_project_button_icon = QtGui.QIcon()
        new_project_button_icon_path = os.path.join('resources', 'b_new.png')
        new_project_button_icon.addPixmap(QtGui.QPixmap(new_project_button_icon_path), QtGui.QIcon.Normal,
                                          QtGui.QIcon.Off)
        new_project_button_tooltip = 'Create new project'

        save_project_button = QToolButton(self)
        save_project_button_icon = QtGui.QIcon()
        save_project_button_icon_path = os.path.join('resources', 'b_save.png')
        save_project_button_icon.addPixmap(QtGui.QPixmap(save_project_button_icon_path), QtGui.QIcon.Normal,
                                           QtGui.QIcon.Off)
        save_project_button_tooltip = 'Save project'

        open_project_button = QToolButton(self)
        open_project_button_icon = QtGui.QIcon()
        open_project_button_icon_path = os.path.join('resources', 'b_open.png')
        open_project_button_icon.addPixmap(QtGui.QPixmap(open_project_button_icon_path), QtGui.QIcon.Normal,
                                           QtGui.QIcon.Off)
        open_project_button_tooltip = 'Open existing project...'


        # adding icon to the toolbuttonss
        new_project_button.setIcon(new_project_button_icon)
        new_project_button.clicked.connect(self.clearModel)
        new_project_button.setToolTip(new_project_button_tooltip)

        save_project_button.setIcon(save_project_button_icon)
        save_project_button.clicked.connect(self.exportModel)
        save_project_button.setToolTip(save_project_button_tooltip)

        open_project_button.setIcon(open_project_button_icon)
        open_project_button.clicked.connect(self.file_open_dialogue)
        open_project_button.setToolTip(open_project_button_tooltip)
        

        # add buttons to toolbars
        fileToolBar.addWidget(new_project_button)
        fileToolBar.addWidget(save_project_button)
        fileToolBar.addWidget(open_project_button)
   
        self.addToolBar(editToolBar)

    def about_dialogue(self):
        txt_label = """
        GMD (geodynamic model designer) by Shevyrev Sergei
        Software was developed in the FEGI FEB FAS
        It is used for development of the geodynamic models
        for further processing with GML (geodynamic model loader) 

        GML impements computation techniques of MIC method,
        which geodynamic implementation was proposed by Gerya T.V.(2019)
        """
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setTextFormat(Qt.RichText)
        msgBox.setText(
            txt_label + '\n <br><a href=\'http://fegi.ru\'>fegi.ru</a> <br> <a href=\'http://lefa.geologov.net\'>lefa.geologov.net</a>')
        msgBox.exec_()

        
class EnterSource(base_2, form_2):
    def __init__(self):
        super(base_2,self).__init__()
        self.setupUi(self)
        self.layout = QVBoxLayout(self)
        
        self.okBtn.clicked.connect(self.okFcn)
        self.closeBtn.clicked.connect(self.closeFcn)
        
        labels=['Air/water',
                'Sediments',
                'Sediments(m)',
                'Sediments(rm)',
                'Sediments(f)',
                'Basalts',
                'Gabbro',
                'Lithospheric mantle',
                'Asthenospheric mantle',
                'Hydrated mantle',
                'Upper continental crust (granodiorite)',
                'Lower continental crust (diorite)',
                'Basalt(melt)',
                'Diorite(melt)',
                'Granodiorite(melt)']
        for i in labels:
                item1 = QListWidgetItem() #need to copy theese items twice
                item1.setText(str(i))
                self.selObjectList.addItem(item1)
        self.selObjectList.setSelectionMode(QAbstractItemView.SingleSelection)
        self.selObjectList.itemSelectionChanged.connect(self.on_change1)

    def getXPercValues(self,pix:int)->int:
        return int((pix/ex.modelWidth)*100)

    def getYPercValues(self,pix:int)->int:
        return int((pix/ex.modelHeight)*100)
    
    def okFcn(self):
        global ixp, iyp, ixr, iyr, pnts_list
        #if ex.sourcetype=='area':
        #    if(self.threshMax)
        #msgBox=QMessageBox.question(self, 'Warning','Data may not be entered properly. Are you sure?',
        #                              QMessageBox.Yes|QMessageBox.No, QMessageBox.No)
        #if msgBox==QMessageBox.Yes:
        if ixp<ixr:
            xStartPer,xEndPer = self.getXPercValues(ixp),self.getXPercValues(ixr)
            #xStartPer=int((ixp/ex.modelWidth)*100); xEndPer=int((ixr/ex.modelWidth)*100); 
        else:
            #xStartPer=int((ixr/ex.modelWidth)*100); xEndPer=int((ixp/ex.modelWidth)*100);
            xStartPer,xEndPer = self.getXPercValues(ixr),self.getXPercValues(ixp) 
        if iyp<iyr:    
            #yStartPer=int((iyp/ex.modelHeight)*100); yEndPer=int((iyr/ex.modelHeight)*100)
            yStartPer,yEndPer = self.getYPercValues(iyp),self.getYPercValues(iyr)
        else:
            #yStartPer=int((iyr/ex.modelHeight)*100); yEndPer=int((iyp/ex.modelHeight)*100)
            yStartPer,yEndPer = self.getYPercValues(iyr),self.getYPercValues(iyp)

        #join polygons by cicle
        print('join polygons by cicle')
        print('len(ex.poly_list)=',len(ex.poly_list))
        if len(ex.poly_list)>0:
            print('len(ex.poly_list)>0')
            #convert xx and yy to percent xxp yyp
            xyp = [] #pair x and y lists 
            #TODO merge all polygons in poly_list
            #if polygons are separated MultiPolygon is formed => no exterior points
            #ex.poly_list = [shapely.union_all(ex.poly_list)]

            for p in ex.poly_list:
                #xyp = [] #pair x and y lists 
                xx, yy = p.exterior.coords.xy
                print('xx=',xx)
                print('yy=',yy)
                xyp.append([[self.getXPercValues(x) for x in xx],[self.getYPercValues(y) for y in yy]])


        else:
            print('len(ex.poly_list)==0')
            #TODO взять точки из границ, не забыть перобразовать в проценты
            #xxp = [xStartPer,xStartPer,xEndPer,xStartPer,xStartPer]
            #yyp = [yStartPer,yEndPer,yEndPer,yEndPer,yStartPer]
            # dx = ixp - ixr
            # dy = iyp - iyr
            #xyp = [[[ixp,ixp+dx,ixr,ixr+dx,ixp],[iyp,iyp,iyr+dy,iyr,iyp]]]
            xyp = [[[xStartPer,xStartPer,xEndPer,xEndPer,xStartPer],[yStartPer,yEndPer,yEndPer,yStartPer,yStartPer]]]

        xmin=min([min([x for x in xy[0]]) for xy in xyp]);ymin=min([min([y for y in xy[1]]) for xy in xyp])
        xmax=max([max([x for x in xy[0]]) for xy in xyp]);ymax=max([max([y for y in xy[1]]) for xy in xyp])

        #ex.bm.append([[ self.selObjectList.selectedItems()[0].text()],[yStartPer,xStartPer, yEndPer,xEndPer],[xxp,yyp]])
        
        #add to bottom
        #ex.bm.append([[ self.selObjectList.selectedItems()[0].text()],[ymin,xmin, ymax,xmax],xyp])
        #add to top
        ex.bm.insert(0,[[ self.selObjectList.selectedItems()[0].text()],[ymin,xmin, ymax,xmax],xyp])

        print([[ self.selObjectList.selectedItems()[0].text()],[yStartPer,xStartPer, yEndPer,xEndPer],xyp])
        ex.drawModel(); #recalculate and output the model
        ex.block_params.close();
        #clear list of polygons
        ex.poly_list.clear()
        
    
    def closeFcn(self):
        #if ex.sourcetype=='area':
        #    if(self.threshMax)
        #msgBox=QMessageBox.question(self, 'Warning','Data may not be entered properly. Are you sure?',
        #                              QMessageBox.Yes|QMessageBox.No, QMessageBox.No)
        #if msgBox==QMessageBox.Yes:
        ex.block_params.close();
    
    def on_change1(self):
        if self.selObjectList.selectedItems():
            self.okBtn.setEnabled(True)
               
class EnterPhysics(base_3, form_3):
    def __init__(self):
        super(base_3,self).__init__()
        self.setupUi(self)
        self.layout = QVBoxLayout(self)
        
        self.btnApply.clicked.connect(self.applyFcn)
        self.btnReset.clicked.connect(self.resetFcn)
    
    def applyFcn(self):
        #set from widgets to settings dictionary

        msgBox=QMessageBox.question(self, 'Info','Do you wish to apply these \
                                    values? Window will be closed',
                                    QMessageBox.Yes|QMessageBox.No, QMessageBox.No)

        if msgBox == QMessageBox.Yes:
            mf = ex.modelPhysics

            mf.update({'gx':self.sb_gx.value()})
            mf.update({'gy':self.sb_gy.value()})
            mf.update({'heightkm':self.sb_heightkm.value()})
            mf.update({'widthkm':self.sb_widthkm.value()}) 
            mf.update({'seadepth':self.sb_seadepth.value()}) 
            mf.update({'ttop':self.sb_Ttop.value()})
            mf.update({'tbottom':self.sb_Tbottom.value()})  
            mf.update({'tintrus':self.sb_TIntrusion.value()}) 
            mf.update({'tmoho':self.sb_TMoho.value()})     #temperature at Moho boundary (crust bottom)
            mf.update({'ymohokm':self.sb_YMoho.value()})   #depth at Moho boundary (km)
            mf.update({'lithbottom':self.sb_YLithBottom.value()})   #depth at Moho boundary (km)
            mf.update({'melting':self.chbMelt.isChecked()})
            if self.rbStable.isChecked(): mf.update({'tmode':'stationary'}) 
            if self.rbExtension.isChecked(): mf.update({'tmode':'extension'}) 
            if self.rbShort.isChecked(): mf.update({'tmode':'collision'})  
            if self.rbSubduction.isChecked(): mf.update({'tmode':'subduction'}) 
            mf.update({'horcolvel':self.sb_horColVel.value()*3.168808781402895e-10}) #cm/s #TODO convers to m/s in a model
            mf.update({'horextvel':self.sb_horExtVel.value()*3.168808781402895e-10}) #cm/s #TODO convers to m/s in a model
            mf.update({'horextstrainpow':self.sb_horExtStr.value()}) #cm/s #TODO convers to m/s in a model

            del mf #remove link 
            #write to file settings.p
            ex.saveSettings(ex.modelPhysics,ex.settings_fname)
            self.close()
    
    
    def resetFcn(self):
        #restore from ini
        msgBox=QMessageBox.question(self, 'Info','Do you really want to \
                                    restore default values?',
                                    QMessageBox.Yes|QMessageBox.No, QMessageBox.No)

        if msgBox == QMessageBox.Yes:
            ex.restorePhysics()
        

            
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    ex.show()
    sys.exit(app.exec_())