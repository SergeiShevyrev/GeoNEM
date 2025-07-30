#cython: language_level=2
"""
Implements intrusion_emplacement.m model in Python


Created on Thu Jun 16 16:21:38 2022

@model and Matlab (tm) software author: T.Gerya
@cython/python porter and GUI developer: S.Shevyrev
"""


import numpy as np
cimport numpy as cnp
cnp.import_array()


from cython.parallel import prange

#get core renredering into the order
from cython import boundscheck, wraparound


from numpy.random import rand
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix 
from scipy.sparse.linalg import spsolve
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from skimage.transform import resize #to rescale demo images
import pickle
import time
from shapely.geometry import Point, Polygon,multipolygon,LineString,box
from shapely import intersection
import random
import os,sys

import json
import skimage as ski


from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double,c_int,c_void_p,Structure,byref,cast,c_longdouble,c_long
#set cython datatype
#adopt datetime for cython
from cpython.datetime cimport datetime 

import uilib as ui


# set proper current directory
current_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
os.chdir(current_dir)

#start_params = {} 
cdef dict start_params={}; #rewrite or update

def setup_defaults():
    global start_params
    start_params.update({'modelfile':'test.cfg',
                        'outdirname':'out',
                        'csv':0,'is_load':'n',
                        'dislocation_type':0,
                        'melt_type':0,
                        'log_file':'log.out',
                        'stepmax':10000,
                        'savestep':100,
                        'savepicstep':10
                        })

if os.path.isfile('startup.p'):
    start_params = ui.parse_startup('startup.p')


    if 'modelfile' not in start_params:
        start_params.update({'modelfile':'test.cfg'})
    if 'outdirname' not in start_params:
        start_params['outdirname'] = 'out'
    if 'csv' not in start_params:
        start_params['csv'] = 0
    if 'is_load' not in start_params:
        start_params['is_load'] = 'n'
    if 'dislocation_type' not in start_params:
        start_params['dislocation_type'] = 0
    if 'melt_type' not in start_params:
        start_params['melt_type'] = 0
    if 'log_file' not in start_params:
        start_params['log_file'] = 'log.out'
    if 'stepmax' not in start_params:
        start_params['stepmax'] = 10000
    if 'savestep' not in start_params:
        start_params['savestep'] = 100
    if 'savepicstep' not in start_params:
        start_params['savepicstep'] = 10
    if 'savepicstep' not in start_params:
        start_params['savepicstep'] = 10
else:
    setup_defaults()

###############functions block

#set cython datatype
DTYPE2 = np.int64
DTYPE = np.double #on gcc float64 is double  float128 is double
ctypedef cnp.int64_t DTYPE_t #define c datatype


#returns datetime for logging
def get_now():
    cdef datetime ts=datetime.now();
    cdef str out='{}-{}-{} {}.{}.{}: '.format(ts.year,ts.month,ts.day,ts.hour,ts.minute,ts.second);
    return out;

#TODO returns numpy array of given shape for C pointer
def pntd2d_toarray(pntr,rows:int,cols:int,DTYPE):
    out = []
    for i in range(rows):
        row = [] 
        for j in range(cols):
            row.append(pntr[i][j])
        out.append(row)
    return np.array(out,dtype=DTYPE)

# Function Melt_fraction()
# This function compute melt fraction (xmelt) and latent heat (hlat)
# at given pressure (ppa), temperature (mtk)and rock type (rock)
# Function returns solution for new temperature (vx,vy,pr)
# and distribution of residuals (resx,resy,resc)

cdef Melt_fraction(DTYPE_t ppa, double mtk, int rock):
    
    cdef double tl,ts,P 
    cdef int HL
    # Calculate melt fraction using marker type
    P=ppa*1e-6; # Pressure, MPa
    tl=0; # Liquidus temperature
    
    #switch rock
    
    # 1 = Sticky air/water, no melting
    if rock==1:
        tl=0;
    
        # 2 = Sediments: latent heat 300 kJ/kg
    elif rock==2:
        # Solidus Temperature
        if (P<1200): 
            ts=889+17900/(P+54)+20200/(P+54)**2; 
        else:
            ts=831+0.06*P;
        
        # Liquidus temperature
        tl=1262+0.09*P;
        # Latent heat
        HL=300000;
    
        # 3 = Basalt: latent heat 380 kJ/kg
    elif rock==3:
        # Solidus Temperature
        if (P<1600): 
            ts=973-70400/(P+354)+77800000/(P+354)**2; 
        else:
            ts=935+0.0035*P+0.0000062*P**2;
        
        # Liquidus temperature
        tl=1423+0.105*P;
        # Latent heat
        HL=380000;

    # 4 = Gabbro: latent heat 380 kJ/kg
    elif rock==4:
        # Solidus Temperature
        if (P<1600): 
            ts=973-70400/(P+354)+77800000/(P+354)**2; 
        else:
            ts=935+0.0035*P+0.0000062*P**2;
        
        # Liquidus temperature
        tl=1423+0.105*P;
        # Latent heat
        HL=380000;
    
    # 5 = Lithospheric mantle (dry): latent heat 400 kJ/kg
    elif rock==5:
        # Solidus Temperature
        if (P<10000): 
            ts=1394+0.132899*P-0.000005104*P**2; 
        else:
            ts=2212+0.030819*(P-10000);
        
        # Liquidus temperature
        tl=2073+0.114*P;
        # Latent heat
        HL=400000;

    # 6 = Asthenospheric mantle (dry): latent heat 400 kJ/kg
    elif rock==6:
        # Solidus Temperature
        if (P<10000): 
            ts=1394+0.132899*P-0.000005104*P**2; 
        else:
            ts=2212+0.030819*(P-10000);
    
        # Liquidus temperature
        tl=2073+0.114*P;
        # Latent heat
        HL=400000;
    
    # 7 = Hydrated mantle (wet): latent heat 400 kJ/kg
    elif rock==7:
        # Solidus Temperature
        if (P<2400): 
            ts=1240+49800/(P+323); 
        else:
            ts=1266-0.0118*P+0.0000035*P**2;
        # Liquidus temperature
        tl=2073+0.114*P;
        # Latent heat
        HL=400000;

    # 8 = Upper continental crust: latent heat 300 kJ/kg
    elif rock==8:
    # Solidus Temperature
        if (P<1200): 
            ts=889+17900/(P+54)+20200/(P+54)**2; 
        else:
            ts=831+0.06*P;
        
        # Liquidus temperature
        tl=1262+0.09*P;
        # Latent heat
        HL=300000;

    # 9 = Lower continental crust: latent heat 380 kJ/kg
    elif rock==9:
    # Solidus Temperature
        if (P<1600): 
            ts=973-70400/(P+354)+77800000/(P+354)**2; 
        else:
            ts=935+0.0035*P+0.0000062*P**2;
        # Liquidus temperature
        tl=1423+0.105*P;
        # Latent heat
        HL=380000;
    else:
        print('error - unknown rock type {}'.format(rock));
    # Melt fraction and latent heat calc, check
    xmelt=0;
    hlat=0;
    if (tl>0):
    	# Solidus and liquidus must not entersect
        # in the extrapolation region
        if (ts>tl-100): 
            ts=tl-100;
        # Melt fraction
        xmelt=(mtk-ts)/(tl-ts);
        if (xmelt<0): 
            xmelt=0;
        if (xmelt>1):
            xmelt=1;
        	# Latent heat calc 
        hlat=HL*xmelt;
 	
    return xmelt,hlat


#some global variable
labels  =      ['Air/water',
                'Sediments',
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

def importModelJson(filename):
    #TODO importModelJsons
    #читаем из json
    info_read = []
    with open(filename,encoding='utf-8') as f:
        info_read = json.loads(f.read())

    modelPhysics = {}
    modelPhysics.update({'modelname':info_read['name']})
    modelPhysics.update({'description':info_read['description']})
    modelPhysics.update({'sizepx':[info_read['sizepx'][0],info_read['sizepx'][1]]})
    modelPhysics.update({'heightkm':info_read['sizekm'][0]})
    modelPhysics.update({'widthkm':info_read['sizekm'][1]})
    modelPhysics.update({'seadepth':info_read['seadepth']})
    modelPhysics.update({'blockunits':info_read['blockunits']})
    modelPhysics.update({'gx':info_read['gravity'][0]})
    modelPhysics.update({'gy':info_read['gravity'][1]})
    modelPhysics.update({'ttop':info_read['ttop']})
    modelPhysics.update({'tbottom':info_read['tbottom']})
    modelPhysics.update({'tintrus':info_read['tintrus']})
    modelPhysics.update({'tmode':info_read['tmode']})

    modelPhysics.update({'tmoho':info_read['tmoho']})
    modelPhysics.update({'ymohokm':info_read['ymohokm']})

    #if key is unlisted above
    for key in info_read:
        if key not in [*modelPhysics]:
            modelPhysics.update({key:info_read[key]})

    bm = []
    for block in info_read['blocks']:
        bm.append([block['blockname'],block['rectangle_ullr'],block['points']])

    print('Model is loaded','The model file has been loaded successfully!')
    return modelPhysics, bm    

def percent2pixels(percval:float,pxdimtot:int)->int:
    """
    converts percent value into pixel value, takes percent value and dimension 
    """
    return int(pxdimtot*percval/100)

def percent2m(percval:float,kmdimtot:int)->float:
    """
    converts percent value into pixel value, takes percent value and dimension 
    """
    return int(kmdimtot*percval/100) #m

#draws raster model for blocks loaded from JSON
def drawModel3(xnum,ynum,bm):   #this function is really applied
    """
    draw model into matrix by polygons with the artifacts
    """
    print("start drawing a model")
    #sedimentary layer thickness
    ltpx = 3 #pixels
    #reinitialize model
    model=np.ones([ynum,xnum])*6 #use MI 6 for background pixels
    xproc=xnum/100; yproc=ynum/100;  #percent conversion indeces

    for block in bm[::-1]:               #self.bm blocks of model REVERSE ORDER
        print('looking for:')
        print(block[0][0])
        try:
            value=labels.index(block[0][0])+1
        except ValueError:
            value=-1
            #idx=[block[1][0]*yproc,block[1][1]*xproc,block[1][2]*yproc,block[1][3]*xproc]
            #idx=np.uint16(np.round(idx))
        for obj in block[2]:
            yy = [y*yproc for y in obj[1]] #convertto pixels
            xx = [x*xproc for x in obj[0]] 
            #add polygon to model 
            rr, cc = ski.draw.polygon(yy, xx)
            #if these are sediments
            try:
                if block[0][0] == 'Sediments': #sediments MI 2 are interlaced with MI 8 - upper continental crust 
                    #get min and max row indeces
                    rmin = min(rr); rmax = max(rr)
                    for n,l in enumerate(range(rmin,rmax+1,ltpx)):
                        for r,c in zip(rr,cc):
                            if l<=r<=(l+ltpx):
                                if n%2 == 0:
                                    model[r,c] = 2
                                else:
                                    model[r,c] = 8
                else:
                    model[rr,cc]=value
            except:
                print('can not draw, coordinate network error')
    return model

def generate_sedimentary_dislocaitons(minxpx,minypx,maxxpx,maxypx,ltpx=3,dislocation_type=1)->list:
    """
    generates list of polygons of shapely to apply Marker-in-Cell method
    """
    if dislocation_type==1 or dislocation_type==2: #monocline
        if dislocation_type==1:
            dip_angle = 35          #straight
        elif dislocation_type==2:
            dip_angle = 125         #reversed
        h=ltpx;      #monocline layer thickness in px
        hv=h/np.sin(np.deg2rad(dip_angle));
        #dx=(32e3-7e3)/np.tan(np.deg2rad(dip_angle));
        dx=(maxypx-minypx)/np.tan(np.deg2rad(dip_angle));
    
        pnts_poly=[]
        #for n in range(0,int(maxxpx/h*2)):
        for n in range(int(minxpx-(maxxpx-minxpx)/4),int(2*maxxpx)):
            if not n%2==0:
                hv_r=int(np.random.randint(5,10)*0.1*hv) #randomized hv thickness
                poly_margins=box(minxpx,minypx,maxxpx,maxypx)
                try:
                    poly_layer=Polygon([(-8*h+n*h,minypx),(-8*h+hv_r+n*h,minypx),(-8*h+dx+hv_r+n*h,maxypx),(-8*h+dx+n*h,maxypx)])
                    poly_inter = intersection(poly_margins,poly_layer)
                    
                    yy, xx = poly_inter.exterior.coords.xy
                    #yy, xx = poly_layer.exterior.coords.xy
                    
                    #xx=list(xx)
                    #yy=list(yy)
                    if len(xx)>0 and len(yy)>0:
                        pnts_poly.append([xx,yy])
                except:
                    print('can not get intersect polygon')

    elif dislocation_type==3: #folds
        #dislocation description
       
        folds = 6 # how many folds cycles
        resolution = 100 # how many datapoints to generate
        start_depth=minypx-ltpx*10  #interval start
        end_depth= maxypx #interval end
        thickness=ltpx #layer thickness

        height= thickness * 2 #vertical amplitude of folds
        num_layers=((end_depth-start_depth)//thickness)*3
        
        length = np.pi * 2 * folds
        sine=np.sin(np.arange(0, length, length / resolution))*height  #wave
        
        #shapely polygon for layer
        poly_margins=box(minxpx,minypx,maxxpx,maxypx)

        layers=[] #two lines, lower and upper margin
        
        for n in range(num_layers):
            if n%2 !=0:
                my_wave = start_depth+n*thickness + sine
                my_wave2=my_wave+thickness*np.random.randint(6,10)*0.1
                x_points=np.linspace(0, xnum,len(my_wave))
                layers.append((np.concatenate((x_points,np.flip(x_points))),\
                                np.concatenate((my_wave,np.flip(my_wave2)))))

        #create polygons for layers polylines
        pnts_poly=[]
        for i in range(len(layers)):
            poly=[]
            for x,y in zip(layers[i][0],layers[i][1]):
                poly.append((x,y))

            #try:
            poly_layer=Polygon(poly)
            poly_inter = intersection(poly_margins,poly_layer)
            
            if poly_inter.geom_type == 'MultiPolygon':
                for polygon in list(poly_inter.geoms):
                    yy, xx = polygon.exterior.coords.xy
                    xx=list(xx)
                    yy=list(yy)
                    if len(xx)>0 and len(yy)>0:
                        pnts_poly.append([xx,yy])
            elif poly_inter.geom_type == 'Polygon':
                yy, xx = poly_inter.exterior.coords.xy
                xx=list(xx)
                yy=list(yy)
                if len(xx)>0 and len(yy)>0:
                    pnts_poly.append([xx,yy])
    else:
        print('Dislocation type is unknown')
        #shapely polygon for layer
        pnts_poly=[]
        poly_margins=box(minxpx,minypx,maxxpx,maxypx)
        yy, xx = poly_margins.exterior.coords.xy
        pnts_poly.append([xx,yy])
    
    return pnts_poly


def drawModel4(xnum,ynum,bm):
        """
        draw model into matrix by polygons with the artifacts
        """
        print("start drawing a model")
        #sedimentary layer thickness
        ltpx = 3 #pixels
        #reinitialize model
        model=np.ones([ynum,xnum])*6
        xproc=xnum/100; yproc=ynum/100;  #percent conversion indeces

        for block in bm[::-1]:               #self.bm blocks of model REVERSE ORDER
            print('looking for:')
            print(block[0][0])
            try:
                value=labels.index(block[0][0])+1
            except ValueError:
                value=-1
             #idx=[block[1][0]*yproc,block[1][1]*xproc,block[1][2]*yproc,block[1][3]*xproc]
             #idx=np.uint16(np.round(idx))
            for obj in block[2]:
                yy = [y*yproc for y in obj[1]] #convertto pixels
                xx = [x*xproc for x in obj[0]] 
                #add polygon to model 
                rr, cc = ski.draw.polygon(yy, xx)
                rr = [r for r in rr if r<=ynum]
                cc = [c for c in cc if c<=xnum]
                # print('rr=',rr)
                # print('cc=',cc)
                #if these are sediments
                #try:
                if block[0][0] == 'Sediments': #sediments MI 2 are interlaced with MI 8 - upper continental crust 
                    #get min and max row indeces
                    rmin = min(rr); rmax = max(rr)
                    for n,l in enumerate(range(rmin,rmax+1,ltpx)):
                        for r,c in zip(rr,cc):
                            if l<=r<=(l+ltpx):
                                if n%2 == 0:
                                    model[r,c] = 2
                                else:
                                    model[r,c] = 8
                elif block[0][0] == 'Sediments(m)' or block[0][0] == 'Sediments(rm)': #sediments MONOCLINE MI 2 are interlaced with MI 8 - upper continental crust 
                    dislocation_type = 1 if block[0][0] == 'Sediments(m)' else 2
                    #get min and max row indeces
                    rmin = min(rr); rmax = max(rr)
                    cmin = min(cc); cmax = max(cc)
                    pnts_poly = generate_sedimentary_dislocaitons(cmin,rmin,cmax,rmax,ltpx=3,dislocation_type=dislocation_type)

                    model[rr,cc] = 8
                    for pnts in pnts_poly:
                        xpoly,ypoly = pnts[0],pnts[1] 
                        rows, cols = ski.draw.polygon(xpoly, ypoly)
                        model[rows,cols] = 2

                elif block[0][0] == 'Sediments(f)': #sediments FOLDS MI 2 are interlaced with MI 8 - upper continental crust 
                    print('folds detected')
                    time.sleep(1)
                    #get min and max row indeces
                    rmin = min(rr); rmax = max(rr)
                    cmin = min(cc); cmax = max(cc)
                    pnts_poly = generate_sedimentary_dislocaitons(cmin,rmin,cmax,rmax,ltpx=3,dislocation_type=3)

                    model[rr,cc] = 8
                    for pnts in pnts_poly:
                        xpoly,ypoly = pnts[0],pnts[1] 
                        rows, cols = ski.draw.polygon(xpoly, ypoly)
                        model[rows,cols] = 2

                else:
                    try:
                        model[rr,cc]=value
                    except:
                        print('can not draw, coordinate network error')
        return model

#def Stokes_Continuity_solver_sandbox(cnp.ndarray L, cnp.ndarray R,cnp.ndarray prfirst,
#                                     cnp.ndarray etas,cnp.ndarray etan,int xnum,int ynum,
#                                     cnp.ndarray gridx, cnp.ndarray gridy,cnp.ndarray RX,cnp.ndarray RY,
#                                     cnp.ndarray RC,cnp.ndarray bleft,cnp.ndarray bright,
#                                     cnp.ndarray btop,cnp.ndarray bbottom,cnp.ndarray bintern):
#@boundscheck(False)
#@wraparound(False)
#cdef Stokes_Continuity_solver_sandbox(double [:,:] L, double [:] R, double [:] prfirst,
#                                     double [:,:] etas, double [:,:] etan, int xnum, int ynum,
#                                    double [:] gridx, double [:] gridy, double [:,:] RX, double [:,:] RY,
#                                     double [:,:] RC, double [:,:] bleft, double [:,:] bright,
#                                     double [:,:] btop, double [:,:] bbottom, double [:,:] bintern):

    
#
#@boundscheck(False)
#@wraparound(False)
#cdef Stokes_Continuity_solver_sandbox(L, double[:] R, prfirst, etas, etan, Py_ssize_t xnum,Py_ssize_t ynum,
#                                     double[:] gridx, double[:] gridy, double[:,:]RX, double[:,:]RY, double[:,:]RC,long[:,:]bleft, long[:,:]bright,
#                                     long[:,:]btop, long[:,:]bbottom, long[:]bintern):    

#paper how to speed up cython loops
#https://nealhughes.net/cython1/
    

@boundscheck(False)
@wraparound(False)
cdef Stokes_Continuity_solver_sandbox(L, double[:] R, prfirst, etas, etan, Py_ssize_t xnum,Py_ssize_t ynum,
                                     double[:] gridx, double[:] gridy, double[:,:]RX, double[:,:]RY,  double[:,:]RC,   double[:,:]bleft,   double[:,:]bright,
                                       double[:,:]btop,   double[:,:]bbottom, bintern):     
    
    # 
    # Staggered Grid for Multigrid
    # 
    #     vx       vx       vx    
    #long
    # vy  +---vy---+---vy---+   vy
    #     |        |        |
    #     vx   P   vx   P   vx    
    #     |        |        |
    # vy  +---vy---+---vy---+   vy
    #     |        |        |
    #     vx   P   vx   P   vx    
    #     |        |        |
    # vy  +---vy---+---vy---+   vy
    #
    #     vx       vx       vx    
    # 
    # Lines show basic grid
    # Basic (density) nodes are shown with +
    # Ghost nodes shown outside the basic grid
    # are used for boundary conditions
    cdef int i;
    cdef int j;
    # Boundary conditions
    # Pressure boundary condition
    # Pressure in first cell
    cdef int bpres=0;
    cdef int prnorm=prfirst[1];
    # Channel flow top->bottom
    if (prfirst[0]==1):
        bpres=1;
        prnorm=prfirst[1];
    
    
    # Velocity boundary conditions
    cdef double[:,:] btopx=np.zeros([btop.shape[0],2],dtype=DTYPE); 
    cdef double[:,:] btopy=np.zeros([btop.shape[0],2],dtype=DTYPE);
    cdef double[:,:] bbottomx=np.zeros([bbottom.shape[0],2],dtype=DTYPE);
    cdef double[:,:] bbottomy=np.zeros([bbottom.shape[0],2],dtype=DTYPE);
    cdef double[:,:] bleftx=np.zeros([bleft.shape[0],2],dtype=DTYPE);
    cdef double[:,:] blefty=np.zeros([bleft.shape[0],2],dtype=DTYPE);
    cdef double[:,:] brightx=np.zeros([bright.shape[0],2],dtype=DTYPE);
    cdef double[:,:] brighty=np.zeros([bright.shape[0],2],dtype=DTYPE);
    
    #weights 
    cdef double [:,:] dsxyn=np.zeros([ynum,xnum],dtype=DTYPE);
    cdef double [:,:] dsxxn=np.zeros([ynum-1,xnum-1],dtype=DTYPE);
    cdef double [:,:] wtetas=np.zeros([ynum,xnum],dtype=DTYPE);
    cdef double [:,:] wtetan=np.zeros([ynum,xnum],dtype=DTYPE);
        
    #cdef int N = btop.shape[0]  #переполнение буфера
    # for i in prange(N, nogil=True):
    #     btopx[i,0]=btop[i,0];
    #     btopx[i,1]=btop[i,1];
    #     btopy[i,0]=btop[i,2];
    #     btopy[i,1]=btop[i,3];
    #     bbottomx[i,0]=bbottom[i,0];
    #     bbottomx[i,1]=bbottom[i,1];
    #     bbottomy[i,0]=bbottom[i,2];
    #     bbottomy[i,1]=bbottom[i,3];
    btopx[:,0]=btop[:,0];            #TypeError: only size-1 arrays can be converted to Python scalars
    btopx[:,1]=btop[:,1];
    btopy[:,0]=btop[:,2];
    btopy[:,1]=btop[:,3];
    bbottomx[:,0]=bbottom[:,0];
    bbottomx[:,1]=bbottom[:,1];
    bbottomy[:,0]=bbottom[:,2];
    bbottomy[:,1]=bbottom[:,3];
    btopx[:,0]=btop[:,0];
    btopx[:,1]=btop[:,1];
    btopy[:,0]=btop[:,2];
    btopy[:,1]=btop[:,3];
    bbottomx[:,0]=bbottom[:,0];
    bbottomx[:,1]=bbottom[:,1];
    bbottomy[:,0]=bbottom[:,2];
    bbottomy[:,1]=bbottom[:,3];
    # cdef int M = btop.shape[0]  #переполнение буфера
    # for i in prange(M, nogil=True):
    #     bleftx[i,0]=bleft[i,0];
    #     bleftx[i,1]=bleft[i,1];
    #     blefty[i,0]=bleft[i,2];
    #     blefty[i,1]=bleft[i,3];
    #     brightx[i,0]=bright[i,0];
    #     brightx[i,1]=bright[i,1];
    #     brighty[i,0]=bright[i,2];
    #     brighty[i,1]=bright[i,3];
    
    # Prescribed internal velocity condition
    
    
    # Computing grid steps for basic nodes
    #cdef double[:] Y = np.zeros(N)
    #cdef cnp.ndarray xstp=np.zeros([xnum-1,1]); #datatype for numpy arrays
    #cdef cnp.ndarray ystp=np.zeros([ynum-1,1]);
    cdef  double[:] xstp=np.zeros([xnum-1],dtype=DTYPE); #datatype for numpy arrays
    cdef  double[:] ystp=np.zeros([ynum-1],dtype=DTYPE);
    
    
    
    
    for i in prange(0,xnum-1,nogil=True):
        xstp[i]=gridx[i+1]-gridx[i];
    
    for i in prange(0,ynum-1,nogil=True):
        ystp[i]=gridy[i+1]-gridy[i];
    #TODO - errors when compose equation
    #xstp=gridx[1:]-gridx[0:-1];
    #ystp=gridy[1:]-gridy[0:-1];
    
    # Computing grid steps for vx and vy nodes
    cdef  double[:] xstpc=np.zeros([xnum],dtype=DTYPE);
    cdef  double[:] ystpc=np.zeros([ynum],dtype=DTYPE);
    
    # First and last steps (for ghost nodes)
    xstpc[0]=xstp[0];
    xstpc[len(xstpc)-1]=xstp[len(xstp)-1];
    ystpc[0]=ystp[0];
    ystpc[len(ystpc)-1]=ystp[len(ystp)-1];
    
    for i in range(1,xnum-1):
        xstpc[i]=(gridx[i+1]-gridx[i-1])/2;
    
    for i in range(1,ynum-1):
        ystpc[i]=(gridy[i+1]-gridy[i-1])/2;
    #xstpc=(gridx[2:]-gridx[0:-2])/2;    
    #ystpc=(gridy[2:]-gridy[0:-2])/2;
    
    # Average x and y steps
    xstpavr=(gridx[len(gridx)-1]-gridx[0])/len(gridx);
    ystpavr=(gridy[len(gridy)-1]-gridy[0])/len(gridy);
    
    # Coefficient for scaling pressure
    pscale=2*etan[0,0]/(xstpavr+ystpavr); #select top left node as for Gerya etan[0,0]
    
    # Horizontal shift index
    ynum3=(ynum-1)*3;
    
    # Solving of Stokes and continuity equations on nodes
    
    cdef int ivx,ivy,ipr;
    
    for i in range(0,ynum-1):
        for j in range(0,xnum-1):
            # Indexes for P,vx,vy
            #ivx=((j-1)*(ynum-1)+(i-1))*3+1; #тут считает неверно!!!!!!!!!!!!! 
            #print('i={}'.format(i));
            #print('j={}'.format(j));
            
            #ind=np.ravel_multi_index([[i],[j]], (ynum-1,xnum-1), order='F')
            #ivx=ind*3-2; #
            #ivy=ind*3-1;
            #ipr=ind*3;
            
            
            #УБЕДИТЬСЯ ЧТО ЛЕВАЯ И ПРАВАЯ ЧАСТЬ L и R совпадают. Возможно, они потерялись еще раньше.
            #проверить как они собираются ниже. В солверах, рассмотренных ранее R работает с одним индексом
            
            #linear and their unknowns are described in Gerya, 2nd ed. 2019 page 99-100 (108-109)
            
            ivx=(((j)*(ynum-1)+(i))*3);
            ivy=ivx+1;
            ipr=ivx+2;
            
            #ivx=((j-1)*(ynum-1)+(i-1))*3;
            #ivy=ivx+1;
            #ipr=ivx+2;
            
            #ivx=ind*3;
            #ivy=ivx+1;
            #ipr=ivx+2;
            
            
            # x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
            if (j<xnum-2 and (j!=bintern[0] or i<bintern[1] or i>bintern[2])): #ind changed
                # x-Stokes equation stensil
                #     +----------------------+----------------------+   
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx(i-1,j)                 ------------    
                #     |                      |                      |
                #     |                      |                      |
                #     +-----vy(i-1,j)---etas(i,j+1)---vy(i-1,j+1)----    ystpc(i)--------
                #     |                      |                      |
                #     |                      |                      |
                # vx(i,j-1)  pr(i,j)      vx(i,j)     P(i,j+1)   vx(i,j+1)-----   ystp(i)
                #     |     etan(i,j)        |       etan(i,j+1)    |
                #     |                      |                      |
                #     +------vy(i,j)---etas(i+1,j+1)---vy(i,j+1)-----    ystpc(i+1)------
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx(i+1,j)                 -----------    
                #     |                      |                      |
                #     |                      |                      |
                #     +----------------------+----------------------+   
                #     |         xstp(j)      |      xstp(j+1)       |   
                #                 |      xstpc(j+1)     |   
                # Right Part
                R[ivx]=RX[i+1,j+1];
                # Computing Current x-Stokes coefficients
                # Central Vx node
                #print('ivx={}'.format(ivx))
                L[ivx,ivx]=-2*(etan[i,j+1]/xstp[j+1]+etan[i,j]/xstp[j])/xstpc[j+1]-(etas[i+1,j+1]/ystpc[i+1]+etas[i,j+1]/ystpc[i])/ystp[i];
                #L[ivx,ivx]=-2*(etan[i,j+1]/xstp[j+1]+etan[i,j]/xstp[j])/xstpc[j+1]-(etas[i+1,j+1]/ystpc[i+1]+etas[i,j+1]/ystpc[i])/ystp[i];
                # Left Vx node
                if (j>0):
                    ivxleft=ivx-ynum3;
                    #print('ivxleft=');
                    #print(ivxleft);
                    #time.sleep(5);
                    L[ivx,ivxleft]=2*etan[i,j]/xstp[j]/xstpc[j+1];
                else:
                    L[ivx,ivx]=L[ivx,ivx]+bleftx[i+1,1]*2*etan[i,j]/xstp[j]/xstpc[j+1];
                    R[ivx]=R[ivx]-bleftx[i+1,0]*2*etan[i,j]/xstp[j]/xstpc[j+1];
                
                # Right Vx node
                if (j<xnum-2): #ind changed
                    ivxright=ivx+ynum3;
                    L[ivx,ivxright]=2*etan[i,j+1]/xstp[j+1]/xstpc[j+1];
                else:
                    L[ivx,ivx]=L[ivx,ivx]+brightx[i+1,1]*2*etan[i,j+1]/xstp[j+1]/xstpc[j+1];
                    R[ivx]=R[ivx]-brightx[i+1,0]*2*etan[i,j+1]/xstp[j+1]/xstpc[j+1];
                
                # Top Vx node
                if (i>0):
                    ivxtop=ivx-3;
                    L[ivx,ivxtop]=etas[i,j+1]/ystpc[i]/ystp[i];
                else:
                    L[ivx,ivx]=L[ivx,ivx]+btopx[j+1,1]*etas[i,j+1]/ystpc[i]/ystp[i];
                    R[ivx]=R[ivx]-btopx[j+1,0]*etas[i,j+1]/ystpc[i]/ystp[i];
                
                # Bottom Vx node
                if (i<ynum-2):
                    ivxbottom=ivx+3;
                    L[ivx,ivxbottom]=etas[i+1,j+1]/ystpc[i+1]/ystp[i];
                else:
                    L[ivx,ivx]=L[ivx,ivx]+bbottomx[j+1,1]*etas[i+1,j+1]/ystpc[i+1]/ystp[i];
                    R[ivx]=R[ivx]-bbottomx[j+1,0]*etas[i+1,j+1]/ystpc[i+1]/ystp[i];
                
                # Top Left Vy node
                if (i>0):
                    ivytopleft=ivx-3+1;
                    L[ivx,ivytopleft]=etas[i,j+1]/xstpc[j+1]/ystp[i];
                else:
                    ivybottomleft=ivx+1;
                    L[ivx,ivybottomleft]=btopy[j+1,1]*etas[i,j+1]/xstpc[j+1]/ystp[i];
                    R[ivx]=R[ivx]-btopy[j+1,0]*etas[i,j+1]/xstpc[j+1]/ystp[i];
                # Top Right Vy node
                if (i>0):
                    ivytopright=ivx-3+1+ynum3;
                    L[ivx,ivytopright]=-etas[i,j+1]/xstpc[j+1]/ystp[i];
                else:
                    ivybottomright=ivx+1+ynum3;
                    L[ivx,ivybottomright]=-btopy[j+2,1]*etas[i,j+1]/xstpc[j+1]/ystp[i];
                    R[ivx]=R[ivx]+btopy[j+2,1]*etas[i,j+1]/xstpc[j+1]/ystp[i];
                
                # Bottom Left Vy node
                if (i<ynum-2):
                    ivybottomleft=ivx+1;
                    if (i>0):
                        L[ivx,ivybottomleft]=-etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                    else:
                        L[ivx,ivybottomleft]=L[ivx,ivybottomleft]-etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                else:
                    ivytopleft=ivx-3+1;
                    L[ivx,ivytopleft]=L[ivx,ivytopleft]-bbottomy[j+1,1]*etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                    R[ivx]=R[ivx]+bbottomy[j+1,0]*etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                # Bottom Right Vy node
                if (i<ynum-2):
                    ivybottomright=ivx+1+ynum3;
                    if (i>0):
                        L[ivx,ivybottomright]=etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                    else:
                        L[ivx,ivybottomright]=L[ivx,ivybottomright]+etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                    
                else:
                    ivytopright=ivx-3+1+ynum3;
                    L[ivx,ivytopright]=L[ivx,ivytopright]+bbottomy[j+2,1]*etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                    R[ivx]=R[ivx]-bbottomy[j+2,0]*etas[i+1,j+1]/xstpc[j+1]/ystp[i];
                
                # Left P node
                iprleft=ivx+2;
                #print('iprleft={}'.format(iprleft));
                #print('pscale={}'.format(pscale));
                #print('xstpc[j+1]={}'.format(xstpc[j+1]))
                
                L[ivx,iprleft]=pscale/xstpc[j+1];
                # Right P node
                iprright=ivx+2+ynum3;
                L[ivx,iprright]=-pscale/xstpc[j+1];
                
            # Ghost Vx_parameter=0 used for numbering, internal prescribed horizontal velocity
            else:
                L[ivx,ivx]=2*pscale/(xstpavr+ystpavr);
                if (j!=bintern[0] or i<bintern[1] or i>bintern[2]):
                    R[ivx]=0;
                else:
                    # Internal prescribed horizontal velocity
                    R[ivx]=2*pscale/(xstpavr+ystpavr)*bintern[3];
                
                
            # y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
            if (i<ynum-2 and (j!=bintern[4] or i<bintern[5] or i>bintern[6])): #ind changed
                # y-Stokes equation stensil
                #     +-------------------- -+-------vy(i-1,j)------+----------------------+-----    
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx(i,j-1)     P(i,j)    vx(i,j)                   |ystp(i)-------    
                #     |                      |        etan(i,j)     |                      |
                #     |                      |                      |                      |
                #     +-----vy(i,j-1)---etas(i+1,j)---vy(i,j)--etas(i+1,j+1)---vy(i,j+1)---+----- ystpc(i+1)
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx(i+1,j-1)  P(i+1,j)   vx(i+1,j)                 |ystp(i+1)------    
                #     |                      |      etan(i+1,j)     |                      |
                #     |                      |                      |                      |
                #     +----------------------+-------vy(i+1,j)------+----------------------+-----
                #               |          xstpc(j)      |         xstpc(j+1)      |   
                #                            |          xstp(j)     |
                #
                # Right Part
                R[ivy]=RY[i+1,j+1];
                # Computing Current y-Stokes coefficients
                # Central Vy node
                L[ivy,ivy]=-2*(etan[i+1,j]/ystp[i+1]+etan[i,j]/ystp[i])/ystpc[i+1]-(etas[i+1,j+1]/xstpc[j+1]+etas[i+1,j]/xstpc[j])/xstp[j];
                # Top Vy node
                if(i>0):
                    ivytop=ivy-3;
                    L[ivy,ivytop]=2*etan[i,j]/ystp[i]/ystpc[i+1];
                else:
                    # Add boundary condition for the top Vy node
                    L[ivy,ivy]=L[ivy,ivy]+btopy[j+1,1]*2*etan[i,j]/ystp[i]/ystpc[i+1];
                    R[ivy]=R[ivy]-btopy[j+1,0]*2*etan[i,j]/ystp[i]/ystpc[i+1];
                
                # Bottom Vy node
                if(i<ynum-3): #ind changed
                    ivybottom=ivy+3;
                    L[ivy,ivybottom]=2*etan[i+1,j]/ystp[i+1]/ystpc[i+1];
                else:
                    # Add boundary condition for the bottom Vy node
                    L[ivy,ivy]=L[ivy,ivy]+bbottomy[j+1,1]*2*etan[i+1,j]/ystp[i+1]/ystpc[i+1];
                    R[ivy]=R[ivy]-bbottomy[j+1,0]*2*etan[i+1,j]/ystp[i+1]/ystpc[i+1];
                
                # Left Vy node
                if(j>0):
                    ivyleft=ivy-ynum3;
                    L[ivy,ivyleft]=etas[i+1,j]/xstpc[j]/xstp[j];
                else:
                    # Add boundary condition for the left Vy node
                    L[ivy,ivy]=L[ivy,ivy]+blefty[i+1,1]*etas[i+1,j]/xstpc[j]/xstp[j];
                    R[ivy]=R[ivy]-blefty[i+1,0]*etas[i+1,j]/xstpc[j]/xstp[j];
                
                # Right Vy node
                if(j<xnum-2):
                    ivyright=ivy+ynum3;
                    L[ivy,ivyright]=etas[i+1,j+1]/xstpc[j+1]/xstp[j];
                else:
                    # Add boundary condition for the right Vy node
                    L[ivy,ivy]=L[ivy,ivy]+brighty[i+1,1]*etas[i+1,j+1]/xstpc[j+1]/xstp[j];
                    R[ivy]=R[ivy]-brighty[i+1,0]*etas[i+1,j+1]/xstpc[j+1]/xstp[j];
                
                # Top left Vx node
                if (j>0):
                    ivxtopleft=ivy-1-ynum3;
                    L[ivy,ivxtopleft]=etas[i+1,j]/ystpc[i+1]/xstp[j];
                else:
                    ivxtopright=ivy-1;
                    L[ivy,ivxtopright]=bleftx[i+1,1]*etas[i+1,j]/ystpc[i+1]/xstp[j];
                    R[ivy]=R[ivy]-bleftx[i+1,0]*etas[i+1,j]/ystpc[i+1]/xstp[j];
                
                # Bottom left Vx node
                if (j>0):
                    ivxbottomleft=ivy-1+3-ynum3;
                    L[ivy,ivxbottomleft]=-etas[i+1,j]/ystpc[i+1]/xstp[j];
                else:
                    ivxbottomright=ivy-1+3;
                    L[ivy,ivxbottomright]=-bleftx[i+2,1]*etas[i+1,j]/ystpc[i+1]/xstp[j];
                    R[ivy]=R[ivy]+bleftx[i+2,0]*etas[i+1,j]/ystpc[i+1]/xstp[j];

                # Top right Vx node
                if (j<xnum-2): #ind changed
                    ivxtopright=ivy-1;
                    if(j>0):
                        L[ivy,ivxtopright]=-etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    else:
                        L[ivy,ivxtopright]=L[ivy,ivxtopright]-etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    
                else:
                    ivxtopleft=ivy-1-ynum3;
                    L[ivy,ivxtopleft]=L[ivy,ivxtopleft]-brightx[i+1,1]*etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    R[ivy]=R[ivy]+brightx[i+1,0]*etas[i+1,j+1]/ystpc[i+1]/xstp[j];
               
                # Bottom right Vx node
                if (j<xnum-2): #ind changed
                    ivxbottomright=ivy-1+3;
                    if(j>0):
                        L[ivy,ivxbottomright]=etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    else:
                        L[ivy,ivxbottomright]=L[ivy,ivxbottomright]+etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    
                else:
                    ivxbottomleft=ivy-1+3-ynum3;
                    L[ivy,ivxbottomleft]=L[ivy,ivxbottomleft]+brightx[i+2,1]*etas[i+1,j+1]/ystpc[i+1]/xstp[j];
                    R[ivy]=R[ivy]-brightx[i+2,0]*etas[i+1,j+1]/ystpc[i+1]/xstp[j];
               
                # Top P node
                iprtop=ivy+1;
                L[ivy,iprtop]=pscale/ystpc[i+1];
                # Bottom P node
                iprbottom=ivy+1+3;
                L[ivy,iprbottom]=-pscale/ystpc[i+1];
                
            # Ghost Vy_parameter=0 used for numbering
            else:
                L[ivy,ivy]=2*pscale/(xstpavr+ystpavr);
                if (j!=bintern[4] or i<bintern[5] or i>bintern[6]):
                    R[ivy]=0;
                else:
                    # Internal prescribed horizontal velocity
                    R[ivy]=2*pscale/(xstpavr+ystpavr)*bintern[7];
            
            # Continuity equation dvx/dx+dvy/dy=RC
            if ( ((j>0 or i>0) and bpres==0) or (i>0 and i<ynum-2 and bpres==1) or (j>0 and j<xnum-2 and bpres==1) ): 
                # Continuity equation stensil
                #     +-----vy(i-1,j)--------+--------
                #     |                      |
                #     |                      |
                # vx(i,j-1)  pr(i,j)      vx(i,j) ystp(i)
                #     |                      |
                #     |                      |
                #     +------vy(i,j)---------+--------
                #     |        xstp(j)       |
                #
                # Right Part
                R[ipr]=RC[i,j];
                # Computing Current Continuity coefficients
                # Left Vx node
                if (j>0):
                    ivxleft=ipr-2-ynum3;
                    L[ipr,ivxleft]=-pscale/xstp[j];
                    # Add boundary condition for the right Vx node
                    if (j==xnum-1):
                        L[ipr,ivxleft]=L[ipr,ivxleft]+brightx[i+1,1]*pscale/xstp[j];
                        R[ipr]=R[ipr]-brightx[i+1,0]*pscale/xstp[j];

                # Right Vx node
                if (j<xnum-2): #ind changed
                    ivxright=ipr-2;
                    L[ipr,ivxright]=pscale/xstp[j];
                    # Add boundary condition for the left Vx node
                    if (j==0): #ind changed
                        L[ipr,ivxright]=L[ipr,ivxright]-bleftx[i+1,1]*pscale/xstp[j];
                        R[ipr]=R[ipr]+bleftx[i+1,0]*pscale/xstp[j];
                    
                # Top Vy node
                if (i>0):
                    ivytop=ipr-1-3;
                    L[ipr,ivytop]=-pscale/ystp[i];
                    # Add boundary condition for the bottom Vy node
                    if (i==ynum-2): #ind changed
                        L[ipr,ivytop]=L[ipr,ivytop]+bbottomy[j+1,1]*pscale/ystp[i];
                        R[ipr]=R[ipr]-bbottomy[j+1,0]*pscale/ystp[i];
                    
                # Bottom Vy node
                if (i<ynum-2): #ind changed
                    ivybottom=ipr-1;
                    L[ipr,ivybottom]=pscale/ystp[i];
                    # Add boundary condition for the top Vy node
                    if (i==0):
                        L[ipr,ivybottom]=L[ipr,ivybottom]-btopy[j+1,1]*pscale/ystp[i];
                        R[ipr]=R[ipr]+btopy[j+1,0]*pscale/ystp[i];
                    
                
            # Pressure definition for the boundary condition regions
            else:
                # Pressure definition in one cell
                if (bpres==0):
                    L[ipr,ipr]=2*pscale/(xstpavr+ystpavr);
                    R[ipr]=2*prnorm/(xstpavr+ystpavr);
                
                # Pressure definition at the top and bottom
                if (bpres==1):
                    L[ipr,ipr]=2*pscale/(xstpavr+ystpavr);
                    if (i==0):
                        R[ipr]=2*prnorm/(xstpavr+ystpavr);
                    else:
                        R[ipr]=0;

                # Pressure definition at the left and right
                if (bpres==2):
                    L[ipr,ipr]=2*pscale/(xstpavr+ystpavr);
                    if (j==0):
                        R[ipr]=2*prnorm/(xstpavr+ystpavr);
                    else:
                        R[ipr]=0;
                        
    # Solve matrix
    S=spsolve(L,R);
    #S=np.linalg.solve(L,R);
    #print('S=')
    #print(S)
    
    # Reload solution
    vx=np.zeros([ynum+1,xnum]);
    vy=np.zeros([ynum,xnum+1]);
    pr=np.zeros([ynum-1,xnum-1]);
    for i in range(0,ynum-1):
        for j in range(0,xnum-1):
            # Indexes for P,vx,vy
            #ivx=((j-1)*(ynum-1)+(i-1))*3+1;
            #ivy=ivx+1;
            #ipr=ivx+2;
            ivx=(((j)*(ynum-1)+(i))*3);
            ivy=ivx+1;
            ipr=ivx+2;
            
            # Reload Vx
            if (j<xnum-2):
                vx[i+1,j+1]=S[ivx];
            
            # Reload Vy
            if (i<ynum-2):
                vy[i+1,j+1]=S[ivy];
            
            # Reload P
            pr[i,j]=S[ipr]*pscale;
    
    # Apply vx boundary conditions
    # Left,Right Boundary
    for i in range(0,ynum):
        vx[i,0]=bleftx[i,0]+bleftx[i,1]*vx[i,1];
        vx[i,xnum-1]=brightx[i,0]+brightx[i,1]*vx[i,xnum-2];
        
    # Top, Bottom Boundary
    for j in range(0,xnum):
        vx[0,j]=btopx[j,0]+btopx[j,1]*vx[1,j];
        vx[ynum,j]=bbottomx[j,0]+bbottomx[j,1]*vx[ynum-1,j];
    
    # Apply vy boundary conditions
    # Left,Right Boundary
    for i in range(0,ynum):
        vy[i,0]=blefty[i,0]+blefty[i,1]*vy[i,1];
        vy[i,xnum]=brighty[i,0]+brighty[i,1]*vy[i,xnum];
        
    # Top, Bottom Boundary
    for j in range(0,xnum+1):
        vy[0,j]=btopy[j,0]+btopy[j,1]*vy[1,j];
        vy[ynum-1,j]=bbottomy[j,0]+bbottomy[j,1]*vy[ynum-2,j];
        
    # Initialize residual arrays
    resx=np.zeros([ynum+1,xnum]);
    resy=np.zeros([ynum,xnum+1]);
    resc=np.zeros([ynum-1,xnum-1]);
    # Computing residuals
    for i in range(0,ynum+1):
        for j in range(0,xnum+1):
    
            # x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
            if (j<xnum and (j!=bintern[0] or i<bintern[1] or i>bintern[2])):
                # vx-Boundrary conditions 
                if (i==0 or i==ynum or j==0 or j==xnum-1):
                    resx[i,j]=0;
                # x-Stokes equation
                # x-Stokes equation stensil
                #     +----------------------+----------------------+   
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx(i-1,j)                 ------------    
                #     |                      |                      |
                #     |                      |                      |
                #     +-----vy(i-1,j)---etas(i-1,j)---vy(i-1,j+1)----    ystpc(i-1)----
                #     |                      |                      |
                #     |                      |                      |
                # vx(i,j-1)  pr(i-1,j-1)  vx(i,j)     P(i-1,j)   vx(i,j+1)-----   ystp(i-1)
                #     |    etan(i-1,j-1)     |       etan(i-1,j)    |
                #     |                      |                      |
                #     +-------vy(i,j)----etas(i,j)-----vy(i,j+1)-----    ystpc(i)------
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx(i+1,j)                 -----------    
                #     |                      |                      |
                #     |                      |                      |
                #     +----------------------+----------------------+   
                #     |         xstp(j-1)    |      xstp(j)         |   
                #                 |      xstpc(j)        |   
                else:
                    # Computing Current x-Stokes residual
                    # dSIGMAxx/dx-dP/dx
                    resx[i,j]=RX[i,j]-(2*(etan[i-1,j]*(vx[i,j+1]-vx[i,j])/xstp[j]-etan[i-1,j-1]*(vx[i,j]-vx[i,j-1])/xstp[j-1])-(pr[i-1,j]-pr[i-1,j-1]))/xstpc[j];
                    # dSIGMAxy/dy
                    resx[i,j]=resx[i,j]-(etas[i,j]*((vx[i+1,j]-vx[i,j])/ystpc[i]+(vy[i,j+1]-vy[i,j])/xstpc[j])-etas[i-1,j]*((vx[i,j]-vx[i-1,j])/ystpc[i-1]+(vy[i-1,j+1]-vy[i-1,j])/xstpc[j]))/ystp[i-1];

                
            # y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
            if (i<ynum and (j!=bintern[4] or i<bintern[5] or i>bintern[6])):
                # vy-Boundrary conditions 
                if (i==0 or i==ynum-1 or j==0 or j==xnum):
                    resy[i,j]=0;
                #y-Stokes equation
                # y-Stokes equation stensil
                #     +-------------------- -+-------vy(i-1,j)------+----------------------+-----    
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx(i,j-1)  P(i-1,j-1)   vx(i,j)                   |ystp(i-1)-------    
                #     |                      |      etan(i-1,j-1)   |                      |
                #     |                      |                      |                      |
                #     +-----vy(i,j-1)---etas(i,j-1)----vy(i,j)----etas(i,j)---vy(i,j+1)----+----- ystpc(i)
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx(i+1,j-1)  P(i,j-1)   vx(i+1,j)                 |ystp(i)------    
                #     |                      |      etan(i,j-1)     |                      |
                #     |                      |                      |                      |
                #     +----------------------+-------vy(i+1,j)------+----------------------+-----
                #               |        xstpc(j-1)      |        xstpc(j)      |   
                #                            |       xstp(j-1)      |
                #
                else:
                    # Computing current residual
                    # dSIGMAyy/dy-dP/dy
                    resy[i,j]=RY[i,j]-(2*(etan[i,j-1]*(vy[i+1,j]-vy[i,j])/ystp[i]-etan[i-1,j-1]*(vy[i,j]-vy[i-1,j])/ystp[i-1])-(pr[i,j-1]-pr[i-1,j-1]))/ystpc[i];
                    # dSIGMAxy/dx
                    resy[i,j]=resy[i,j]-(etas[i,j]*((vy[i,j+1]-vy[i,j])/xstpc[j]+(vx[i+1,j]-vx[i,j])/ystpc[i])-etas[i,j-1]*((vy[i,j]-vy[i,j-1])/xstpc[j-1]+(vx[i+1,j-1]-vx[i,j-1])/ystpc[i]))/xstp[j-1];

                
            # Continuity equation dvx/dx+dvy/dy=RC
            if (i<ynum-1 and j<xnum-1):
                # Continuity equation stensil
                #     +------vy(i,j+1)-------+--------
                #     |                      |
                #     |                      |
                # vx(i+1,j)   pr(i,j)   vx(i+1,j+1) ystp(i)
                #     |                      |
                #     |                      |
                #     +-----vy(i+1,j+1)------+--------
                #     |        xstp(j)       |
                #
                # Computing current residual
                resc[i,j]=RC[i,j]-((vx[i+1,j+1]-vx[i+1,j])/xstp[j]+(vy[i+1,j+1]-vy[i,j+1])/ystp[i]);

    return L,R,vx,resx,vy,resy,pr,resc


# Function Temperature_solver_grid()
# This function formulates and solves  
# Heat conservation equation defined on 2D irregularly spaced grid
# with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy)
# given distribution of right parts for all equations (RT) on the grid 
# and given RHO*CP and k values
#
# Thermal Boundary condition specified by arrays bleft(),bright(),btop(),bbot() 
#
# Function returns solution for new temperature (vx,vy,pr)
# and distribution of residuals (resx,resy,resc)

@boundscheck(False)
@wraparound(False)
cdef Temperature_solver_grid(L,R,timestep,xnum,ynum,gridx,gridy,kt,rhocp,tk,RT,bleft,bright,btop,bbottom):
    #print('Thermal Boundary condition specified by arrays bleft(),bright(),btop(),bbot() ');
    #print('Function returns solution for new temperature (vx,vy,pr)');
    # 
    # Staggered Grid for Multigrid
    #
    #     T--------T--------T 
    #     |        |        |
    #     |        |        |
    #     |        |        |
    #     T--------T--------T 
    #     |        |        |
    #     |        |        |
    #     |        |        |
    #     T--------T--------T 
    # 
    # Lines show basic grid
    
    
    # Computing grid steps for basic nodes
    xstp=np.zeros([xnum-1]);
    ystp=np.zeros([ynum-1]);
    for i in range(0,xnum-1):
        xstp[i]=gridx[i+1]-gridx[i];
        
    for i in range(0,ynum-1):
        ystp[i]=gridy[i+1]-gridy[i];

    # Computing grid steps for qx and qy nodes
    xstpc=np.zeros([xnum-2]);
    ystpc=np.zeros([ynum-2]);
    
    for i in range(0,xnum-2):
        xstpc[i]=(gridx[i+2]-gridx[i])/2;
    
    for i in range(0,ynum-2):
        ystpc[i]=(gridy[i+2]-gridy[i])/2;
      
    
    # Solving of Temperature equation 
    for i in range(0,ynum):
        for j in range(0,xnum):
            # Index for T
            itk=np.ravel_multi_index([[i],[j]], (ynum,xnum), order='F');
            #itk=(j-1)*ynum+i; # когда было (j)*ynum+i - выдавало ошибку о сингулярности матрицы
                              # когда было (j-1)*ynum+i - выдавало ошибку о сингулярности матрицы
            # Boundary conditions
            if (j==0 or j==xnum-1 or i==0 or i==ynum-1):
                # Upper boundary: tk(i,j)=btop(1)+btop(2)*tk(i+1,j)
                if (i==0 and j>0 and j<xnum-1):
                    # Right part
                    R[itk,0]=btop[j,0];
                    # Left part: 1*tk(i,j)-btop(2)*tk(i+1,j)
                    L[itk,itk]=1;
                    L[itk,itk+1]=-btop[j,1];
                
                # Lower boundary: tk(i,j)=bbottom(1)+bbottom(2)*tk(i-1,j)
                if (i==ynum-1 and j>0 and j<xnum-1):
                    # Right part
                    R[itk,0]=bbottom[j,0];
                    # Left part: 1*tk(i,j)-bbottom(2)*tk(i-1,j)
                    L[itk,itk]=1;
                    L[itk,itk-1]=-bbottom[j,1];
                
                # Left boundary: tk(i,j)=bleft(1)+bleft(2)*tk(i,j+1)
                if (j==0):
                    # Right part
                    R[itk,0]=bleft[i,0];
                    # Left part: 1*tk(i,j)-bleft(2)*tk(i,j+1)
                    L[itk,itk]=1;
                    L[itk,itk+ynum]=-bleft[i,1];
                
                # Right boundary: tk(i,j)=bright(1)+bright(2)*tk(i,j-1)
                if (j==xnum-1):
                    # Right part
                    R[itk,0]=bright[i,0];
                    # Left part: 1*tk(i,j)-bright(2)*tk(i,j-1)
                    L[itk,itk]=1;
                    L[itk,itk-ynum]=-bright[i,1];
            
            # Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
            else:
                # Temperature equation stensil
                #
                #           | xstpc(j-1) |  
                #     | xstp(j-1) |  xstp(j)  |
                #     +-------tk(i-1,j)-------+--------- 
                #     |       kt(i-1,j)       |
                #     |           |           |
                #     |           |           |    ystp(i-1)-----
                #     |           |           |
                #     |           |           |
                # tk(i,j-1)----tk(i,j)-----tk(i,j+1)------  ystpc(i-1) 
                # kt(i,j-1)    kt(i,j)     kt(i,j+1) 
                #     |       rhocp(i,j)      |
                #     |           |           |    ystp(i)-------
                #     |           |           |
                #     |           |           |
                #     +-------tk(i+1,j)-------+---------
                #             kt(i+1,j)
                #
                # Right Part
                R[itk,0]=RT[i,j]+tk[i,j]*rhocp[i,j]/timestep;
                # Computing coefficients for the left part
                # Central T node
                #print('Central T node=')
                #print(rhocp[i,j]/timestep+((kt[i,j-1]+kt[i,j])/xstp[j-1]+(kt[i,j]+kt[i,j+1])/xstp[j])/2/xstpc[j-1]+((kt[i-1,j]+kt[i,j])/ystp[i-1]+(kt[i,j]+kt[i+1,j])/ystp[i])/2/ystpc[i-1])
                L[itk,itk]=rhocp[i,j]/timestep+((kt[i,j-1]+kt[i,j])/xstp[j-1]+(kt[i,j]+kt[i,j+1])/xstp[j])/2/xstpc[j-1]+((kt[i-1,j]+kt[i,j])/ystp[i-1]+(kt[i,j]+kt[i+1,j])/ystp[i])/2/ystpc[i-1];
                # Left T node
                L[itk,itk-ynum]=-(kt[i,j-1]+kt[i,j])/2/xstp[j-1]/xstpc[j-1];
                # Right T node
                L[itk,itk+ynum]=-(kt[i,j]+kt[i,j+1])/2/xstp[j]/xstpc[j-1];
                # Upper T node
                L[itk,itk-1]=-(kt[i-1,j]+kt[i,j])/2/ystp[i-1]/ystpc[i-1];
                # Lower T node
                L[itk,itk+1]=-(kt[i,j]+kt[i+1,j])/2/ystp[i]/ystpc[i-1];
            
    # Solve matrix
    S=spsolve(L,R);
    tknew=np.zeros([ynum,xnum]);
    
    
    # Reload solution
    for i in range(0,ynum):
        for j in range(0,xnum):
            # Index for T
            #itk=j*ynum+i;
            itk=np.ravel_multi_index([[i],[j]], (ynum,xnum), order='F');
            # Reload T
            tknew[i,j]=S[itk];
    
    
    rest=np.zeros([ynum,xnum]);
    # Computing residuals
    for i in range(0,ynum):
        for j in range(0,xnum):
    
            # Boundary conditions
            if (j==0 or j==xnum-1 or i==0 or i==ynum-1):
                rest[i,j]=0;
            # Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
            else:
                # Computing Current Temperature equation residual
                # Ht-DT/dt
                rest[i,j]=RT[i,j]-rhocp[i,j]*(tknew[i,j]-tk[i,j])/timestep;
                # -dqx/dx
                rest[i,j]=rest[i,j]+((kt[i,j]+kt[i,j+1])*(tknew[i,j+1]-tknew[i,j])/xstp[j]-(kt[i,j-1]+kt[i,j])*(tknew[i,j]-tknew[i,j-1])/xstp[j-1])/2/xstpc[j-1];
                # -dqy/dy
                rest[i,j]=rest[i,j]+((kt[i,j]+kt[i+1,j])*(tknew[i+1,j]-tknew[i,j])/ystp[i]-(kt[i-1,j]+kt[i,j])*(tknew[i,j]-tknew[i-1,j])/ystp[i-1])/2/ystpc[i-1];

    
    return L,R,tknew,rest 

#cnp.int64_t instead of long
#https://github.com/autonomousvision/differentiable_volumetric_rendering/pull/22

#use library function Interpolate_markers_nodes


# load the library (however it is not used yet)
libpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'mylib.so')
mylib = CDLL(libpath)   #already loaded


# C-type corresponding to numpy array 
ND_FLOAT_PNTR1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")

ND_INT_PNTR1 = np.ctypeslib.ndpointer(dtype=np.int64, 
                                      ndim=1,
                                      flags="C")

ND_FLOAT_PNTR2 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=2,
                                      flags="C")

ND_INT_PNTR2 = np.ctypeslib.ndpointer(dtype=np.int64, 
                                      ndim=2,
                                      flags="C")

#markup output
class ReVal(Structure):
    _fields_ = [("rho1", POINTER(POINTER(c_double))),
                ("tk1",POINTER(POINTER(c_double))),
                ("kt1",POINTER(POINTER(c_double))),
                ("rhocp1",POINTER(POINTER(c_double))),
                ("hr1",POINTER(POINTER(c_double))),
                ("ha1",POINTER(POINTER(c_double))),
                ("wtnodes",POINTER(POINTER(c_double))),
                ("etas1",POINTER(POINTER(c_double))),
                ("mus1",POINTER(POINTER(c_double))),
                ("sxy1",POINTER(POINTER(c_double))),
                ("wtetas",POINTER(POINTER(c_double))),
                ("etan1",POINTER(POINTER(c_double))),
                ("mun1",POINTER(POINTER(c_double))),
                ("sxx1",POINTER(POINTER(c_double))),
                ("wtetan",POINTER(POINTER(c_double))),
                ("timesum", c_double)]

#rho1,tk1,kt1,rhocp1,hr1,ha1,wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum




# define prototypes
mylib.interpolate_markers_nodes.argtypes = [
                                            c_int,          #const int xnum, //same as xnum NEW input parameter
                                            c_int,          #const int ynum, //same as xnum NEW input parameter
                                            #c_int,          #const int mark_num, //same as mark_num NEW input parameter
                                            c_int,          #const int layer_num, //same as layer_num NEW input parameter
                                            c_int,          #const int marknum
                                            ND_FLOAT_PNTR1, #gridx[xnum],
                                            ND_FLOAT_PNTR1,  #gridx[ynum],
                                            
                                            #c_size_t,       #size_t xnum
                                            #c_size_t,       #size_t ynum
                                            
                                            ND_FLOAT_PNTR2,  #MX[mark_num][1],
                                            ND_FLOAT_PNTR2,  #MY[mark_num][1],
                                            
                                            ND_INT_PNTR2,    #int MI[mark_num][1]
                                            ND_FLOAT_PNTR2,  #double MRAT[mark_num][1], 
                                            ND_FLOAT_PNTR2,  #double MGII[mark_num][1],
                                            ND_FLOAT_PNTR2,   #double MXN[mark_num][1],
                                            ND_FLOAT_PNTR2,  #double MYN[mark_num][1],
                                            ND_FLOAT_PNTR2,  #double MFLOW[][6],
                                            ND_FLOAT_PNTR2,    #double MPL[][7],
                                            ND_FLOAT_PNTR2,   #double MXM[mark_num][1]
                                            c_double,        #double tstp,
                                            c_int,          #long int tnum,
                                            ND_FLOAT_PNTR2,   #double gridt[][tnum], #какой тип для матриц? 
                                            c_int,           #int waterlev,
                                            c_double,        #double stressmin,
                                            c_int,           #int ttop,
                                            c_double,       #double etamelt,
                                            ND_FLOAT_PNTR2,  #double rho1[][xnum],
                                            ND_FLOAT_PNTR2,  #double tk1[][xnum],
                                            ND_FLOAT_PNTR2,   #double kt1[][xnum],
                                            ND_FLOAT_PNTR2,   #double rhocp1[][xnum],
                                            ND_FLOAT_PNTR2,   #double hr1[][xnum],
                                            ND_FLOAT_PNTR2,    #double ha1[][xnum],
                                            ND_FLOAT_PNTR2,    #wtnodes[][xnum],
                                            
                                            ND_FLOAT_PNTR2,    #double etas1[][xnum],
                                            ND_FLOAT_PNTR2,    #double mus1[][xnum],
                                            ND_FLOAT_PNTR2,   #double sxy1[][xnum],
                                            ND_FLOAT_PNTR2,    #double wtetas[][xnum],
                                            ND_FLOAT_PNTR2,    #double etan1[][xnum-1],
                                            ND_FLOAT_PNTR2,    #double mun1[][xnum-1],
                                            ND_FLOAT_PNTR2,    #double sxx1[][xnum-1],
                                            ND_FLOAT_PNTR2,    #double wtetan[][xnum-1],
                                            c_double,           #double timesum,
                                            ND_FLOAT_PNTR2,     #double xstp1[xnum-1][1],
                                            ND_FLOAT_PNTR2,    #double ystp1[xnum-1][1],
                                            
                                            ND_FLOAT_PNTR2,  #double MRHO[][5],
                                            ND_FLOAT_PNTR2,  #double MKT[][3],
                                            ND_FLOAT_PNTR2,  #double MTK[mark_num][1],
                                            ND_FLOAT_PNTR2,  #double MPR[mark_num][1],
                                            ND_FLOAT_PNTR1,  #double MCP[layer_num],
                                            ND_FLOAT_PNTR1,  #double MHR[layer_num],
                                            
                                            ND_FLOAT_PNTR2,  #double MSXX[mark_num][1],
                                            ND_FLOAT_PNTR2,   #double MSXY[mark_num][1],
                                            ND_FLOAT_PNTR2,   #double MEXX[mark_num][1],
                                            ND_FLOAT_PNTR2,   #double MEXY[mark_num][1],
                                            ND_FLOAT_PNTR2,   #double META[mark_num][1],
                                            ND_FLOAT_PNTR1,   #double MMU[layer_num],
                                            c_double,        #double timestep,
                                            c_double,        #double etamin,
                                            c_double,        #double etamax,
                                            c_int,           #int ntimestep,
                                            c_double,           #double etawt,
                                            c_int]               #int plastyn
                                            
mylib.interpolate_markers_nodes.restype = ReVal



#rho1,tk1,kt1,rhocp1,hr1,ha1,wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum

#convertion 2d tuple returned into tuple
def ptr2d_to_mat(ptr, rows, cols):
    return tuple(tuple(ptr[i][j] for j in range(cols)) for i in range(rows))





@boundscheck(False)
@wraparound(False)
cdef Interpolate_markers_nodes(int marknum,double[:] gridx,double[:] gridy,Py_ssize_t xnum,
                               Py_ssize_t ynum, double [:,:] MX, double [:,:] MY,
                               DTYPE_t[:,:] MI, double[:,:] MRAT, 
                               double[:,:] MGII,double[:,:] MXN, double[:,:] MYN,
                               double[:,:]MFLOW, double[:,:] MPL,double[:,:]MXM,double tstp,
                               long int tnum,double [:,:] gridt,int waterlev,double stressmin,int ttop,
                               double etamelt,double[:,:] rho1,double[:,:] tk1,double[:,:] kt1,
                               double[:,:] rhocp1,double[:,:]hr1, double[:,:]ha1,
                               double[:,:]wtnodes,double[:,:] etas1,double[:,:] mus1,
                               double[:,:] sxy1,double[:,:] wtetas,double[:,:] etan1,
                               double[:,:] mun1,double[:,:] sxx1, double[:,:]wtetan,long long timesum):
# Interpolating parameters from markers to nodes
    cdef long int mm1;
    cdef long int xn,yn,dp;
    cdef double dx, dy,mwt;
    cdef long int xnmin,xnmax;
    cdef long plawiter,grid_ind;
    cdef double MRHOCUR,MRHOCPCUR,MKTCUR,METACUR,sii0,plawexp,\
            xelvis, sxxnew,sxynew,siinewe,sii1,eii1,siicur,eiicur,MCOHES,MFRICT,\
            xmelt, hlat,hlat1,hlat0;
 
    #time.sleep(100);
    for mm1 in range(1,marknum):
        #custom_print('mm1={}'.format(mm1));
        # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and \
            MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 
            # Erosion-sedimentation
            # Find topography node for the marker
            xn=int((MX[mm1,0]-gridt[0,0])/tstp-0.5)+1;
            #custom_print('computed xn');
            #custom_print(xn);
            
            if (xn<0):
                xn=0;
            
            if (xn>tnum-2):
                xn=tnum-2;
            
            # Compute relative distance to topography node
            dx=(MX[mm1,0]-gridt[0,xn])/tstp;
            # Compute topograhy elevation above the marker
            dy=gridt[1,xn]*(1-dx)+gridt[1,xn+1]*dx;
            
            
            # water/air to sediments transformation
            if (MI[mm1,0]==1 and MY[mm1,0]>dy):
                MI[mm1,0]=2; # Change marker type
                MRAT[mm1,0]=1; # Reset strain rate ratio                  
                MGII[mm1,0]=0; # Reset strain
            
            # Rocks to water/air transformation
            if (MI[mm1,0]>1 and MY[mm1,0]<dy):
                MI[mm1,0]=1; # Change marker type
                MRAT[mm1,0]=1; # Reset strain rate ratio
                MGII[mm1,0]=0; # Reset strain
           
            #  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o Mrho(xm,ym)       ?
            #           ?                              ?
            #           ?                              ?
            #  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            #
            # Define indexes for upper left node in the cell where the marker is
            # using bisection
            # Find horizontal index
            #custom_print('Find horizontal index')
            xnmin=0;
            xnmax=xnum;
            while ((xnmax-xnmin)>1):
                # !!! SUBTRACT 0.5 since int16(0.5)=1
                xn=int((xnmax+xnmin)/2);
                #custom_print(xnmax-xnmin)
                #time.sleep(3)
                grid_ind = int(xn);
                if(gridx[grid_ind]>MX[mm1,0]):
                    xnmax=xn;
                else:
                    xnmin=xn;
            
            xn=xnmin;
            # Check horizontal index
            if (xn<0):
                xn=0;
            
            if (xn>xnum-2):
                xn=xnum-2;
            
            # Save horizontal index
            MXN[mm1,0]=xn;
            
            #custom_print('Find vertical index')
            # Find vertical index
            ynmin=0;
            ynmax=ynum;
            while ((ynmax-ynmin)>1):
                # !!! SUBTRACT 0.5 since int16(0.5)=1
                yn=int((ynmax+ynmin)/2);
                if(gridy[yn]>MY[mm1,0]):
                    ynmax=yn;
                else:
                    ynmin=yn;

            #custom_print('vertical yn index was found');  
            #custom_print(yn)
            yn=ynmin;
            # Check vertical index
            if (yn<0):
                yn=0;
            if (yn>ynum-2):
                yn=ynum-2;
            # Save Vertical index
            MYN[mm1,0]=yn;

            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
            dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

            # Compute marker weight koefficient from cell dimensions
            # Number of markers in a cell is in invert proportion to the cell volume
            mwt=1;#/xstp1(xn)/ystp1(yn);                                                        #TODO stopped types explication HERE

            # Compute density from marker temperature
            MRHOCUR=MRHO[MI[mm1,0],1]*(1-MRHO[MI[mm1,0],2]*(MTK[mm1,0]-273))*(1+MRHO[MI[mm1,0],3]*(MPR[mm1,0]-1e+5));

            # Compute rho*Cp for marker 
            MRHOCPCUR=MRHOCUR*MCP[MI[mm1,0]];
            
            #custom_print('Change density for \"air\"')
            # Change density for "air"
            if (MI[mm1,0]==1 and MY[mm1,0]<waterlev):
                MRHOCUR=1;
            
            # Compute thermal conductivity from marker temperature
            # Rock thermal conductivity (W/m/K): k=k0+a/(T+77)
            MKTCUR=MKT[MI[mm1,0],1]+MKT[MI[mm1,0],2]/(MTK[mm1,0]+77);
            
            # Compute adiabatic heating term (alp*T*DP/Dt)
            MHACUR=MRHO[MI[mm1,0],2]*MTK[mm1,0];
            
            #custom_print('Computing Marker Viscosity')
            # Computing Marker Viscosity
            if(MFLOW[MI[mm1,0],1]==0):
                # Constant viscosity
                METACUR=MFLOW[MI[mm1,0],2];
            else:
                # Power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
                # Iterate for viscosity
                # First viscosity value
                # Compute and check old marker stress invariant in Pa
                sii0=(MSXX[mm1,0]**2+MSXY[mm1,0]**2)**0.5;
                # Check old marker stress invariant (should be allways positive to be used in power law)
                if(sii0<stressmin):
                    sii0=stressmin;
                
                # Check marker temperature
                plawexp=MTK[mm1,0];
                if(plawexp<ttop):
                    plawexp=ttop;
                
                # Compute exponential term: 
                # Ea is in J/mol(=1000*kJ/mol)
                # Va is in J/Pa (=1e-6*cm^3) 
                # Cut if too big (at cold temperature);
                plawexp=(MFLOW[MI[mm1,0],4]*1000+MFLOW[MI[mm1,0],5]*1e-6*MPR[mm1,0])/RGAS/plawexp;
                if(plawexp>150):
                    plawexp=150;
                
                # Compute AD*exp[-Ea/RT)
                plawexp=MFLOW[MI[mm1,0],2]*np.exp(-plawexp);
                # Compute strain rate invariant from power law
                eii0=plawexp*(1e-6*sii0)**MFLOW[MI[mm1,0],3];
                # Compute effective viscosity
                eta0=sii0/2/eii0;
                # Forcasting second invariant of future marker stress for given viscoelastic timestep
                xelvis=eta0/(MMU[MI[mm1,0]]*timestep+eta0);
                sxxnew=MSXX[mm1,0]*xelvis+2*eta0*MEXX[mm1,0]*MRAT[mm1,0]*(1-xelvis);
                sxynew=MSXY[mm1,0]*xelvis+2*eta0*MEXY[mm1,0]*MRAT[mm1,0]*(1-xelvis);
                sii1=(sxxnew**2+sxynew**2)**0.5;
                # Check new marker stress invariant (should be allways positive to be used in power law)
                if(sii1<stressmin):
                    sii1=stressmin;
                
                # Compute strain rate invariant from power law
                eii1=plawexp*(1e-6*sii1)**MFLOW[MI[mm1,0],3];
                # Compute effective viscosity
                METACUR=sii1/2/eii1;
                # Iterate for viscosity which corresponds to future stress invariant using bisection
                # Iteration counter
                plawiter=0;
                while(plawiter<20 and abs(sii1-sii0)>1):
                    # Add iteration counter
                    plawiter=plawiter+1;
                    # Compute middle stress
                    siicur=(sii0+sii1)/2;
                    # Compute strain rate invariant from power law
                    eiicur=plawexp*(1e-6*siicur)**MFLOW[MI[mm1,0],3];
                    # Compute effective viscosity
                    METACUR=siicur/2/eiicur;
                    # Forcasting second invariant of future marker stress for given viscoelastic timestep
                    xelvis=METACUR/(MMU[MI[mm1,0]]*timestep+METACUR);
                    sxxnew=MSXX[mm1,0]*xelvis+2*METACUR*MEXX[mm1,0]*MRAT[mm1,0]*(1-xelvis);
                    sxynew=MSXY[mm1,0]*xelvis+2*METACUR*MEXY[mm1,0]*MRAT[mm1,0]*(1-xelvis);
                    siinew=(sxxnew**2+sxynew**2)**0.5;
                    # Changing bisection limits
                    if((sii0<sii1 and siicur<siinew) or (sii0>sii1 and siicur>siinew)):
                        sii0=siicur;
                    else:
                        sii1=siicur;
                   
                # Limiting viscosity for the power law
                if (METACUR<etamin): 
                    METACUR=etamin;
                
                if (METACUR>etamax): 
                    METACUR=etamax;
            
            #custom_print(get_now()+'Check if any plastic yeiding condition is present.');     
            # Check if any plastic yeiding condition is present 
            if (ntimestep>1 and (MPL[MI[mm1,0],1]>0 or MPL[MI[mm1,0],3]>0)):
                # Checking for plastic yeilding
                # Compute second invariant for a purely elastic stress build-up
                sxxnewe=MSXX[mm1,0]+2*MMU[MI[mm1,0]]*timestep*MEXX[mm1,0]*MRAT[mm1,0];
                sxynewe=MSXY[mm1,0]+2*MMU[MI[mm1,0]]*timestep*MEXY[mm1,0]*MRAT[mm1,0];
                siinewe=(sxxnewe**2+sxynewe**2)**0.5;
                # Checking yeilding criterion for strain weakening/hardening
                # C=C0, FI=FI0 for strain<=GAM0
                # C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
                # C=C1, FI=FI1 for strain>=GAM0
                MCOHES=MPL[MI[mm1,0],1];
                MFRICT=MPL[MI[mm1,0],3];
                if (MGII[mm1,0]>=MPL[MI[mm1,0],6]):
                    MCOHES=MPL[MI[mm1,0],2];
                    MFRICT=MPL[MI[mm1,0],4];
                
                if (MGII[mm1,0]>MPL[MI[mm1,0],5] and MGII[mm1,0]<MPL[MI[mm1,0],6]):
                    MCOHES=MPL[MI[mm1,0],1]+(MPL[MI[mm1,0],2]-MPL[MI[mm1,0],1])/(MPL[MI[mm1,0],6]-MPL[MI[mm1,0],5])*(MGII[mm1,0]-MPL[MI[mm1,0],5]);
                    MFRICT=MPL[MI[mm1,0],3]+(MPL[MI[mm1,0],4]-MPL[MI[mm1,0],3])/(MPL[MI[mm1,0],6]-MPL[MI[mm1,0],5])*(MGII[mm1,0]-MPL[MI[mm1,0],5]);
                
                # Computing yelding stress for the marker
                siiyeld=MCOHES+MFRICT*MPR[mm1,0];
                if (siiyeld<0): 
                    siiyeld=0;
                #custom_print(get_now()+'Correcting rock properties for yeilding.');  
                # Correcting rock properties for yeilding 
                if (siiyeld<siinewe):
                    # Bringing marker viscosity to yeilding stress
                    METACURP=MMU[MI[mm1,0]]*timestep*siiyeld/(siinewe-siiyeld);
                    METACURP=METACURP**(1-etawt)*META[mm1,0]**etawt;
                    if(METACURP<METACUR):
                        METACUR=METACURP;
                        # Limiting viscosity for the yeilding
                        if (METACUR<etamin): 
                            METACUR=etamin;
                        
                        if (METACUR>etamax): 
                            METACUR=etamax;
                        
                        # Mark that plastic yeildind occur
                        plastyn=1;
                        # Mark that plastic strain needs to be accumulated
                        MGII[mm1,0]=abs(MGII[mm1,0]);
                    else:
                        # Reset plastic strain if no yelding
                        MGII[mm1,0]=-1e-20;

            # Compute 1/MU values (MU is shear modulus) 
            MMUCUR=1/MMU[MI[mm1,0]];
            
            #custom_print(get_now()+'Molten rocks.'); 
            # Molten rocks
            xmelt, hlat=Melt_fraction(MPR[mm1,0],MTK[mm1,0],MI[mm1,0]);
            # Save marker melting
            #custom_print('timesum={}'.format(timesum))
            if(timesum>0):
                #custom_print('Save marker melt ratio');
                MXM[mm1,0]=xmelt;
            #custom_print('xmelt={}'.format(xmelt));
            if(xmelt>0 and timesum>0):
                # Reset creep parameters for molten rocks                   
                MRAT[mm1,0]=1; # Reset strain rate ratio
                MGII[mm1,0]=0; # Reset strain
                # Viscosity of partially molten rocks
                if(xmelt>0.1):
                    METACUR=etamelt;
                
                # Density
                MRHOCUR=MRHOCUR*((1-xmelt)+MRHO[MI[mm1,0],4]/MRHO[MI[mm1,0],1]*xmelt);
                # RHO*CP
                MRHOCPCUR=MRHOCPCUR*((1-xmelt)+MRHO[MI[mm1,0],4]/MRHO[MI[mm1,0],1]*xmelt);
                # Compute adiabatic heating term (alp*T*DP/Dt)
                MHACUR=MHACUR*((1-xmelt)+MRHO[MI[mm1,0],4]/MRHO[MI[mm1,0],1]*xmelt);
                # Latent heating: effective adiabatic term, RHOCP
                if(xmelt<1):
                    # Melting adiabatic term: alpham=-rho*(dHlat/dP)/T
                    # Numerical differentiation
                    dp=1000; # Pressure increment, Pa
                    xmelt, hlat0=Melt_fraction(MPR[mm1,0]-dp,MTK[mm1,0],MI[mm1,0]);
                    xmelt, hlat1=Melt_fraction(MPR[mm1,0]+dp,MTK[mm1,0],MI[mm1,0]);
                    MHACUR=MHACUR-(hlat1-hlat0)/(2.0*dp);
                    # Melting heat capacity term: cpm=dHlat/dT 
                    # Numerical differentiation 
                    dt=1.0; # Temperature increment, Pa
                    xmelt, hlat0=Melt_fraction(MPR[mm1,0],MTK[mm1,0]-dt,MI[mm1,0]);
                    xmelt, hlat1=Melt_fraction(MPR[mm1,0],MTK[mm1,0]+dt,MI[mm1,0]);
                    MRHOCPCUR=MRHOCPCUR+MRHOCUR*(hlat1-hlat0)/(2.0*dt);
            
            # Save marker viscosity
            META[mm1,0]=copy.copy(METACUR);

            # Add properties to 4 surrounding nodes
            rho1[yn,xn]=rho1[yn,xn]+(1.0-dx)*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn,xn]=tk1[yn,xn]+(1.0-dx)*(1.0-dy)*MTK[mm1,0]*mwt;
            kt1[yn,xn]=kt1[yn,xn]+(1.0-dx)*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn,xn]=rhocp1[yn,xn]+(1.0-dx)*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn,xn]=hr1[yn,xn]+(1.0-dx)*(1.0-dy)*MHR[MI[mm1,0]]*mwt;
            ha1[yn,xn]=ha1[yn,xn]+(1.0-dx)*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn,xn]=wtnodes[yn,xn]+(1.0-dx)*(1.0-dy)*mwt;

            rho1[yn+1,xn]=rho1[yn+1,xn]+(1.0-dx)*dy*MRHOCUR*mwt;
            tk1[yn+1,xn]=tk1[yn+1,xn]+(1.0-dx)*dy*MTK[mm1,0]*mwt;
            kt1[yn+1,xn]=kt1[yn+1,xn]+(1.0-dx)*dy*MKTCUR*mwt;
            rhocp1[yn+1,xn]=rhocp1[yn+1,xn]+(1.0-dx)*dy*MRHOCPCUR*mwt;
            hr1[yn+1,xn]=hr1[yn+1,xn]+(1.0-dx)*dy*MHR[MI[mm1,0]]*mwt;
            ha1[yn+1,xn]=ha1[yn+1,xn]+(1.0-dx)*dy*MHACUR*mwt;
            wtnodes[yn+1,xn]=wtnodes[yn+1,xn]+(1.0-dx)*dy*mwt;

            rho1[yn,xn+1]=rho1[yn,xn+1]+dx*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn,xn+1]=tk1[yn,xn+1]+dx*(1.0-dy)*MTK[mm1,0]*mwt;
            kt1[yn,xn+1]=kt1[yn,xn+1]+dx*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn,xn+1]=rhocp1[yn,xn+1]+dx*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn,xn+1]=hr1[yn,xn+1]+dx*(1.0-dy)*MHR[MI[mm1,0]]*mwt;
            ha1[yn,xn+1]=ha1[yn,xn+1]+dx*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn,xn+1]=wtnodes[yn,xn+1]+dx*(1.0-dy)*mwt;

            rho1[yn+1,xn+1]=rho1[yn+1,xn+1]+dx*dy*MRHOCUR*mwt;
            tk1[yn+1,xn+1]=tk1[yn+1,xn+1]+dx*dy*MTK[mm1,0]*mwt;
            kt1[yn+1,xn+1]=kt1[yn+1,xn+1]+dx*dy*MKTCUR*mwt;
            rhocp1[yn+1,xn+1]=rhocp1[yn+1,xn+1]+dx*dy*MRHOCPCUR*mwt;
            hr1[yn+1,xn+1]=hr1[yn+1,xn+1]+dx*dy*MHR[MI[mm1,0]]*mwt;
            ha1[yn+1,xn+1]=ha1[yn+1,xn+1]+dx*dy*MHACUR*mwt;
            wtnodes[yn+1,xn+1]=wtnodes[yn+1,xn+1]+dx*dy*mwt;

            # Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
            # only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 and dy<=0.5):
                etas1[yn,xn]=etas1[yn,xn]+(1.0-dx)*(1.0-dy)*METACUR*mwt;
                mus1[yn,xn]=mus1[yn,xn]+(1.0-dx)*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn,xn]=sxy1[yn,xn]+(1.0-dx)*(1.0-dy)*MSXY[mm1,0]*mwt;
                wtetas[yn,xn]=wtetas[yn,xn]+(1.0-dx)*(1.0-dy)*mwt;
            
            if(dx<=0.5 and dy>=0.5):
                etas1[yn+1,xn]=etas1[yn+1,xn]+(1.0-dx)*dy*METACUR*mwt;
                mus1[yn+1,xn]=mus1[yn+1,xn]+(1.0-dx)*dy*MMUCUR*mwt;
                sxy1[yn+1,xn]=sxy1[yn+1,xn]+(1.0-dx)*dy*MSXY[mm1,0]*mwt;
                wtetas[yn+1,xn]=wtetas[yn+1,xn]+(1.0-dx)*dy*mwt;
            
            if(dx>=0.5 and dy<=0.5):
                etas1[yn,xn+1]=etas1[yn,xn+1]+dx*(1.0-dy)*METACUR*mwt;
                mus1[yn,xn+1]=mus1[yn,xn+1]+dx*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn,xn+1]=sxy1[yn,xn+1]+dx*(1.0-dy)*MSXY[mm1,0]*mwt;
                wtetas[yn,xn+1]=wtetas[yn,xn+1]+dx*(1.0-dy)*mwt;
            
            if(dx>=0.5 and dy>=0.5):
                etas1[yn+1,xn+1]=etas1[yn+1,xn+1]+dx*dy*METACUR*mwt;
                mus1[yn+1,xn+1]=mus1[yn+1,xn+1]+dx*dy*MMUCUR*mwt;
                sxy1[yn+1,xn+1]=sxy1[yn+1,xn+1]+dx*dy*MSXY[mm1,0]*mwt;
                wtetas[yn+1,xn+1]=wtetas[yn+1,xn+1]+dx*dy*mwt;
            

            # Add viscosity etan(), normal stress sxx() and shear modulus mun() to the center of current cell
            etan1[yn,xn]=etan1[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy))*METACUR*mwt;
            mun1[yn,xn]=mun1[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy))*MMUCUR*mwt;
            sxx1[yn,xn]=sxx1[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy))*MSXX[mm1,0]*mwt;
            wtetan[yn,xn]=wtetan[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy))*mwt;
    
    #end of interpolating parameters from markers to nodes
    return rho1,tk1,kt1,rhocp1,hr1, ha1, wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum


@boundscheck(False)
@wraparound(False)
cdef Compute_subgrid_diffusion(int marknum,double[:] gridx,double[:] gridy,Py_ssize_t xnum,
                               Py_ssize_t ynum, double [:,:] MX, double [:,:] MY,
                               double[:,:] tk1,double[:,:] kt1, double [:,:] dtkn,
                               double[:,:] MXN, double[:,:] MYN, double dsubgridt,
                               double [:,:] wtnodes):
    
    ###start subgrid function body 
    # Interpolating parameters from markers to nodes
    cdef long int mm1;
    cdef long int xn,yn,dp;
    cdef double dx,dy,tkm,ktm,mwt,rhocpm,tdm,sdif,dtkm; 
    
    
    # Marker cycle
    for mm1 in range(1,marknum):
            
        # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

            #  yn    T[yn,xn]--------------------T[yn,xn+1]
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o Mrho[mm1,0]       ?
            #           ?                              ?
            #           ?                              ?
            #  yn+1  T[yn+1,xn]-------------------V[yn+1,xn+1]
            #
            #
            # Interpolating temperature changes from basic nodes
            #
            # Define indexes for upper left node in the cell where the marker is
            xn=int(MXN[mm1,0]);
            yn=int(MYN[mm1,0]);

            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
            dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

            # Compute marker weight koefficient from cell dimensions
            # Number of markers in a cell is in invert proportion to the cell volume
            mwt=1.0;#/xstp1(xn)/ystp1(yn);


            # Interpolate old nodal temperature for the marker
            tkm=0;
            tkm=tkm+(1.0-dx)*(1.0-dy)*tk1[yn,xn];
            tkm=tkm+(1.0-dx)*dy*tk1[yn+1,xn];
            tkm=tkm+dx*(1.0-dy)*tk1[yn,xn+1];
            tkm=tkm+dx*dy*tk1[yn+1,xn+1];
            # Calculate Nodal-Marker subgrid temperature difference
            dtkm=tkm-MTK[mm1,0];
            # Compute nodal k and RHO*Cp for the marker 
            # k
            ktm=0;
            ktm=ktm+(1.0-dx)*(1.0-dy)*kt1[yn,xn];
            ktm=ktm+(1.0-dx)*dy*kt1[yn+1,xn];
            ktm=ktm+dx*(1.0-dy)*kt1[yn,xn+1];
            ktm=ktm+dx*dy*kt1[yn+1,xn+1];
            # RHO*Cp
            rhocpm=0;
            rhocpm=rhocpm+(1.0-dx)*(1.0-dy)*rhocp1[yn,xn];
            rhocpm=rhocpm+(1.0-dx)*dy*rhocp1[yn+1,xn];
            rhocpm=rhocpm+dx*(1.0-dy)*rhocp1[yn,xn+1];
            rhocpm=rhocpm+dx*dy*rhocp1[yn+1,xn+1];

            # Compute local thermal diffusion timescale for the marker
            tdm=rhocpm/ktm/(2/xstp**2+2/ystp**2);

            # Computing subgrid diffusion
            sdif=-dsubgridt*timestep/tdm;
            if(sdif<-30): 
                sdif=-30;
            
            dtkm=dtkm*(1-np.exp(sdif));    

            # Correcting old temperature for the marker
            MTK[mm1,0]=MTK[mm1,0]+dtkm;

            # Interpolating subgrid temperature changes to 4 nodes
            dtkn[yn,xn]=dtkn[yn,xn]+(1.0-dx)*(1.0-dy)*dtkm*mwt;
            wtnodes[yn,xn]=wtnodes[yn,xn]+(1.0-dx)*(1.0-dy)*mwt;

            dtkn[yn+1,xn]=dtkn[yn+1,xn]+(1.0-dx)*dy*dtkm*mwt;
            wtnodes[yn+1,xn]=wtnodes[yn+1,xn]+(1.0-dx)*dy*mwt;

            dtkn[yn,xn+1]=dtkn[yn,xn+1]+dx*(1.0-dy)*dtkm*mwt;
            wtnodes[yn,xn+1]=wtnodes[yn,xn+1]+dx*(1.0-dy)*mwt;

            dtkn[yn+1,xn+1]=dtkn[yn+1,xn+1]+dx*dy*dtkm*mwt;
            wtnodes[yn+1,xn+1]=wtnodes[yn+1,xn+1]+dx*dy*mwt;


    # Computing subgrid diffusion for nodes
    for i in range(0,ynum):
        for j in range(0,xnum):
            # Density
            if (wtnodes[i,j]!=0):
                # Compute new value interpolated from markers
                dtkn[i,j]=dtkn[i,j]/wtnodes[i,j];
    
    
    return dtkn, MTK, wtnodes
    ###end of function body


###start of function body
@boundscheck(False)
@wraparound(False)
cdef Compute_subgrid_stress(int marknum,double[:] gridx,double[:] gridy,Py_ssize_t xnum,
                               Py_ssize_t ynum, double [:,:] MX, double [:,:] MY,
                               DTYPE_t[:,:] MI, double[:] MMU,double [:,:] META,
                               double dsubgrids,double timestep,double [:,:] MSXY, 
                               double [:,:] MSXX,double[:,:] wtetas,double[:,:] wtetan,
                               double[:,:] dsxyn):
    
   
    
    # Interpolating parameters from markers to nodes
    cdef long int mm1;
    cdef long int xn,yn,dp;
    cdef double sdm,sdif,sxym,dsxym,sxxm,dsxxm,sxyn,mwt,dx,dy; 
    
    ###start subgrid stress function body
    
    # Marker cycle
    for mm1 in range(1,marknum):

    # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

            # Compute local stress relaxation timescale (Maxwell time) for the marker
            sdm=META[mm1,0]/MMU[MI[mm1,0]];
            # Computing degree of subgrid stress relaxation
            sdif=-dsubgrids*timestep/sdm;
            if(sdif<-30): 
                sdif=-30;
            
            sdif=(1-np.exp(sdif));

            #  yn    sxy[yn,xn]--------------------sxy[yn,xn+1]
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o MSXY[mm1,0]      ?
            #           ?                              ?
            #           ?                              ?
            #  yn+1  sxy[yn+1,xn]-------------------sxy[yn+1,xn+1]
            #
            #
            # Interpolating old shear stress from Sxy nodes
            #
            # Define indexes for upper left node in the cell where the marker is
            xn=int(MXN[mm1,0]);
            yn=int(MYN[mm1,0]);

            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
            dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

            
            # Compute marker weight koefficient from cell dimensions
            # Number of markers in a cell is in invert proportion to the cell volume
            mwt=1;#/xstp1(xn)/ystp1(yn);


            # Interpolate old Sxy stress for the marker
            sxym=0;
            sxym=sxym+(1.0-dx)*(1.0-dy)*sxy1[yn,xn];
            sxym=sxym+(1.0-dx)*dy*sxy1[yn+1,xn];
            sxym=sxym+dx*(1.0-dy)*sxy1[yn,xn+1];
            sxym=sxym+dx*dy*sxy1[yn+1,xn+1];
            # Calculate Nodal-Marker subgrid Sxy stress difference
            dsxym=sxym-MSXY[mm1,0];
            # Relaxing Nodal-Marker subgrid Sxy stress difference
            dsxym=dsxym*sdif;    

            # Correcting old stress for the marker
            MSXY[mm1,0]=MSXY[mm1,0]+dsxym;

            # Interpolating subgrid Sxy stress changes to 4 nodes
            # only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 and dy<=0.5):
                dsxyn[yn,xn]=dsxyn[yn,xn]+(1.0-dx)*(1.0-dy)*dsxym*mwt;
                wtetas[yn,xn]=wtetas[yn,xn]+(1.0-dx)*(1.0-dy)*mwt;
            
            if(dx<=0.5 and dy>=0.5):
                dsxyn[yn+1,xn]=dsxyn[yn+1,xn]+(1.0-dx)*dy*dsxym*mwt;
                wtetas[yn+1,xn]=wtetas[yn+1,xn]+(1.0-dx)*dy*mwt;
            
            if(dx>=0.5 and dy<=0.5):
                dsxyn[yn,xn+1]=dsxyn[yn,xn+1]+dx*(1.0-dy)*dsxym*mwt;
                wtetas[yn,xn+1]=wtetas[yn,xn+1]+dx*(1.0-dy)*mwt;
            
            if(dx>=0.5 and dy>=0.5):
                dsxyn[yn+1,xn+1]=dsxyn[yn+1,xn+1]+dx*dy*dsxym*mwt;
                wtetas[yn+1,xn+1]=wtetas[yn+1,xn+1]+dx*dy*mwt;
            
            
            # Computing marker weight for the center of current
            # basic cell where Sxx stress is located
            mwt=mwt*(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy));
            #  yn    sxx[yn,xn]--------------------sxx[yn,xn+1]
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o MSXX[mm1,0]       ?
            #           ?                              ?
            #           ?                              ?
            #  yn+1  sxx[yn+1,xn]-------------------sxx[yn+1,xn+1]
            #
            #
            # Interpolating old normal stress from Sxx nodes
            #
            # Define, check indexes for upper left node in the Sxx cell where the marker is
            if (MX[mm1,0]<gridcx[xn+1]):
                xn=xn-1;
            
            if(xn<0):
                xn=0;
            
            if(xn>xnum-3):
                xn=xnum-3;
            
            if (MY[mm1,0]<gridcy[yn+1]):
                yn=yn-1;
            
            if(yn<0):
                yn=0;
            
            if(yn>ynum-3):
                yn=ynum-3;
            
            
            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridcx[xn+1])/xstpc1[xn+1];
            dy=(MY[mm1,0]-gridcy[yn+1])/ystpc1[yn+1];

            # Interpolate old Sxx stress for the marker
            sxxm=0;
            sxxm=sxxm+(1.0-dx)*(1.0-dy)*sxx1[yn,xn];
            sxxm=sxxm+(1.0-dx)*dy*sxx1[yn+1,xn];
            sxxm=sxxm+dx*(1.0-dy)*sxx1[yn,xn+1];
            sxxm=sxxm+dx*dy*sxx1[yn+1,xn+1];
            # Calculate Nodal-Marker subgrid Sxx stress difference
            dsxxm=sxxm-MSXX[mm1,0];
            # Relaxing Nodal-Marker subgrid Sxx stress difference
            dsxxm=dsxxm*sdif;    

            # Correcting old stress for the marker
            MSXX[mm1,0]=MSXX[mm1,0]+dsxxm;

            # Interpolating subgrid Sxx stress changes for the center of current basic cell
            xn=int(MXN[mm1,0]);
            yn=int(MYN[mm1,0]);
            dsxxn[yn,xn]=dsxxn[yn,xn]+dsxxm*mwt;
            wtetan[yn,xn]=wtetan[yn,xn]+mwt;
          
    # Computing subgrid stress changes for nodes
    for i in range(0,ynum):
        for j in range(0,xnum):
            # Density
            if (wtetas[i,j]!=0):
                # Compute new value interpolated from markers
                dsxyn[i,j]=dsxyn[i,j]/wtetas[i,j];
            
            if (j<xnum and i<ynum and wtetan[i,j]!=0):
                # Compute new value interpolated from markers
                dsxxn[i,j]=dsxxn[i,j]/wtetan[i,j];
    
    return np.asarray(wtetan,dtype=DTYPE),np.asarray(wtetas,dtype=DTYPE),np.asarray(dsxxn,dtype=DTYPE),\
                    np.asarray(dsxyn,dtype=DTYPE),np.asarray(MSXX,dtype=DTYPE),np.asarray(MSXY,dtype=DTYPE)
    ###end subgrid stress function body
    
    
    # Marker cycle
###end of function body



@boundscheck(False)
@wraparound(False)
cdef Move_markers_velocity_field(int marknum,double[:] gridx,double[:] gridy,Py_ssize_t xnum,
                               Py_ssize_t ynum, double [:,:] MX, double [:,:] MY,
                               double [:,:] vx1, double [:,:] vy1,
                               double [:,:] esp,
                               double[:,:] MXN, double[:,:] MYN, 
                               double timestep,double [:,:] MSXY, 
                               double [:,:] MSXX,
                               double[:,:] MGII,double[:,:] MBII, double[:,:] MEXX, double[:,:] MEXY, int markmove):
    
    ###start subgrid function body 
    # Interpolating parameters from markers to nodes
    cdef long int mm1;
    cdef long int xn,yn,dp,rk;
    cdef double dx,dy,xcur,ycur,tkm,ktm,mwt,rhocpm,tdm,sdif,dtkm; 
    
    cdef double[:] vxm=np.zeros([5],dtype=DTYPE); #(5,1) - spare index for row count 
    cdef double[:] vym=np.zeros([5],dtype=DTYPE);
    cdef double[:] espm=np.zeros([5],dtype=DTYPE);
    
    # Marker cycle
    for mm1 in range(1,marknum):

        # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

            # Save marker coordinates
            xcur=MX[mm1,0];
            ycur=MY[mm1,0];
            # Defining number of Runge-Kutta cycles
            for rk in range(1,markmove+1):

                #  xn    V(xn,yn)--------------------V(xn+1,yn)
                #           |           ^                  |
                #           |           ?                  |
                #           |          dy                  |
                #           |           ?                  |
                #           |           v                  |
                #           |<----dx--->o Mrho(xm,ym)      |
                #           |                              |
                #           |                              |
                #  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

                # Define indexes for upper left BASIC node in the cell where the marker is
                # using bisection
                # Load horizontal and vertical indexes
                if(rk==1):
                    xnmin=MXN[mm1,0];
                    ynmin=MYN[mm1,0];
                else:
                    # Find horizontal index
                    xnmin=0;
                    xnmax=xnum;
                    while ((xnmax-xnmin)>1):
                        # !!! SUBTRACT 0.5 since int16(0.5)=1
                        xn=int((xnmax+xnmin)/2);
                        if(gridx[xn]>xcur):
                            xnmax=xn;
                        else:
                            xnmin=xn;
                        
                    # Check horizontal index
                    if (xnmin<1):
                        xnmin=0;
                    
                    if (xnmin>xnum-1):
                        xnmin=xnum-1;
                    
                    # Find vertical index
                    ynmin=0;
                    ynmax=ynum;
                    while ((ynmax-ynmin)>1):
                        # !!! SUBTRACT 0.5 since int16(0.5)=1
                        yn=int((ynmax+ynmin)/2);
                        if(gridy[yn]>ycur):
                            ynmax=yn;
                        else:
                            ynmin=yn;

                    # Check vertical index
                    if (ynmin<1):
                        ynmin=0;

                    if (ynmin>ynum-1):
                        ynmin=ynum-1;

                # Define indexes for upper left node in the Vx-cell where the marker is
                # Horizontal Vx index
                xn=int(xnmin);
                # Vertical Vx index
                yn=int(ynmin);
                if(ycur>gridcy[yn+1]):
                    yn=yn+1;
                
                if (yn>ynum):
                    yn=ynum;

                # Define and check normalized distances from marker to the upper left VX-node;
                dx=(xcur-gridx[xn])/xstp1[xn];
                dy=(ycur-gridcy[yn])/ystpc1[yn];

                # Calculate Marker velocity from four surrounding Vx nodes
                vxm[rk]=0;
                vxm[rk]=vxm[rk]+(1.0-dx)*(1.0-dy)*vx1[yn,xn];
                vxm[rk]=vxm[rk]+(1.0-dx)*dy*vx1[yn+1,xn];
                vxm[rk]=vxm[rk]+dx*(1.0-dy)*vx1[yn,xn+1];
                vxm[rk]=vxm[rk]+dx*dy*vx1[yn+1,xn+1];

                # Define indexes for upper left node in the VY-cell where the marker is
                # Vertical Vy index
                yn=int(ynmin);
                # Horizontal Vy index
                xn=int(xnmin);
                if(xcur>gridcx[xn+1]):
                    xn=xn+1;
                
                if (xn>xnum):
                    xn=xnum;
                

                # Define and check normalized distances from marker to the upper left VX-node;
                dx=(xcur-gridcx[xn])/xstpc1[xn];
                dy=(ycur-gridy[yn])/ystp1[yn];

                # Calculate Marker velocity from four surrounding nodes
                vym[rk]=0;
                vym[rk]=vym[rk]+(1.0-dx)*(1.0-dy)*vy1[yn,xn];
                vym[rk]=vym[rk]+(1.0-dx)*dy*vy1[yn+1,xn];
                vym[rk]=vym[rk]+dx*(1.0-dy)*vy1[yn,xn+1];
                vym[rk]=vym[rk]+dx*dy*vy1[yn+1,xn+1];

                # Define indexes for upper left node in the Espin cell where the marker is
                xn=int(xnmin);
                yn=int(ynmin);

                # Define normalized distances from marker to the upper left node;
                dx=(xcur-gridx[xn])/xstp1[xn];
                dy=(ycur-gridy[yn])/ystp1[yn];

                # Interpolate old Sxy stress for the marker
                espm[rk]=0;
                espm[rk]=espm[rk]+(1.0-dx)*(1.0-dy)*esp[yn,xn];
                espm[rk]=espm[rk]+(1.0-dx)*dy*esp[yn+1,xn];
                espm[rk]=espm[rk]+dx*(1.0-dy)*esp[yn,xn+1];
                espm[rk]=espm[rk]+dx*dy*esp[yn+1,xn+1];


                # Update coordinates for the next cycle
                if(rk<4):
                    if (rk<3):
                        xcur=MX[mm1,0]+timestep/2*vxm[rk];
                        ycur=MY[mm1,0]+timestep/2*vym[rk];
                    else:
                        xcur=MX[mm1,0]+timestep*vxm[rk];
                        ycur=MY[mm1,0]+timestep*vym[rk];
                    
            # Recompute velocity and spin using 4-th order Runge_Kutta
            if (markmove==4):
                vxm[1]=(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])/6;
                vym[1]=(vym[1]+2*vym[2]+2*vym[3]+vym[4])/6;
                espm[1]=(espm[1]+2*espm[2]+2*espm[3]+espm[4])/6;
            

            # Displacing Marker according to its velocity
            MX[mm1,0]=MX[mm1,0]+timestep*vxm[1];
            MY[mm1,0]=MY[mm1,0]+timestep*vym[1];
            
            # Rotate stress on marker according to its spin
            # Compute amount of rotation from spin rate:
            # Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
            # (when x axis is directed rightward and y axis is directed downward) 
            espm[1]=espm[1]*timestep;
            # Save old stresses
            msxxold=MSXX[mm1,0];
            msxyold=MSXY[mm1,0];
            # SxyNEW=0.5(Sxx-Syy)*sin(2*Espin*dt)+Sxy*cos(2*Espin*dt)
            # where Sxx-Syy=2Sxx
            MSXY[mm1,0]=msxxold*np.sin(2*espm[1])+msxyold*np.cos(2*espm[1]);
            # SxxNEW=Sxx*(cos(Espin*dt))^2+Syy*(sin(Espin*dt))^2-Sxy*sin(2*Espin*dt)
            # where Sxx=-Syy
            MSXX[mm1,0]=msxxold*((np.cos(espm[1]))**2-(np.sin(espm[1]))**2)-msxyold*np.sin(2*espm[1]);

            # Adding marker plastic strain based on grid strain rates
            if (MGII[mm1,0]>0):
                MGII[mm1,0]=MGII[mm1,0]+timestep*(MEXX[mm1,0]**2+MEXY[mm1,0]**2)**0.5;
            
            # Adding marker bulk strain based on grid strain rates
            MBII[mm1,0]=MBII[mm1,0]+timestep*(MEXX[mm1,0]**2+MEXY[mm1,0]**2)**0.5;
            
    
    return  np.asarray(MX,dtype=DTYPE),np.asarray(MY,dtype=DTYPE),\
            np.asarray(MGII,dtype=DTYPE),np.asarray(MBII,dtype=DTYPE),\
            np.asarray(MEXX,dtype=DTYPE),np.asarray(MEXY,dtype=DTYPE),\
            np.asarray(MSXX,dtype=DTYPE),np.asarray(MSXY,dtype=DTYPE),\
            np.asarray(espm,dtype=DTYPE),np.asarray(vxm,dtype=DTYPE),np.asarray(vym,dtype=DTYPE)
    ###end of fuction body 


#exy = np.zeros([ynum,xnum]);
#exx = np.zeros([ynum-1,xnum-1]);
#esp = np.zeros([ynum,xnum]);
#eii = np.zeros([ynum-1,xnum-1]);

@boundscheck(False)
@wraparound(False)
cdef Compute_strain_rate_pressure_markers(int marknum, double[:] gridx,double[:] gridy,Py_ssize_t xnum,
                               Py_ssize_t ynum, double [:,:] MX, double [:,:] MY,
                               double[:,:] MXN, double[:,:] MYN, double[:,:] gridcx,double[:,:] gridcy,
                               double[:,:] xstpc1,double[:,:] ystp1,double timestep,
                               double [:,:] exy,double [:,:] exx, double [:,:] esp,double [:,:] eii,
                               double[:,:] MEXX,double[:,:] MPR, double[:,:] MEXY, double[:,:] META, 
                               double[:] MMU, double[:,:] MRAT):
    
    # Interpolating parameters from markers to nodes
    cdef long int mm1;
    cdef long int xn,yn;
    cdef double dx,dy,eiim,eiimg,exxm,sxxm,dsxxm; 
    
    cdef double[:] vxm=np.zeros([5],dtype=DTYPE); #(5,1) - spare index for row count 
    cdef double[:] vym=np.zeros([5],dtype=DTYPE);
    cdef double[:] espm=np.zeros([5],dtype=DTYPE);
    
    
    for mm1 in range(1,marknum):
             
        # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

            #  xn    V(xn,yn)--------------------V(xn+1,yn)
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o Mrho(xm,ym)       ?
            #           ?                              ?
            #           ?                              ?
            #  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

            # Define indexes for upper left BASIC node in the cell where the marker is
            xnmin=MXN[mm1,0];
            ynmin=MYN[mm1,0];

            # Calculating strain rate for marker
            #
            # Interpolating squares of EPS'xx=-EPS'yy from cell centers
            # EPS'xx-nodes are displaced rightward and downward for 1/2 of gridsteps
            # Horizontal EPS'xx index
            xn=int(xnmin);
            if(MX[mm1,0]<gridcx[xn+1,0]):
                xn=xn-1;
            
            if (xn<0):
                xn=0;
            
            if (xn>xnum-3):
                xn=xnum-3;
            
            # Vertical EPS'xx index
            yn=int(ynmin);
            if(MY[mm1,0]<gridcy[yn+1,0]):
                yn=yn-1;
            
            if (yn<0):
                yn=0;
            
            if (yn>ynum-3):
                yn=ynum-3;
            

            # Define and check normalized distances from marker to the upper left EPS'xx-node;
            dx=(MX[mm1,0]-gridcx[xn+1,0])/xstpc1[xn+1,0];
            dy=(MY[mm1,0]-gridcy[yn+1,0])/ystpc1[yn+1,0];

            # Calculate and save Marker EPS'xx from four surrounding nodes
            exxm=0;
            exxm=exxm+(1.0-dx)*(1.0-dy)*exx[yn,xn];
            exxm=exxm+(1.0-dx)*dy*exx[yn+1,xn];
            exxm=exxm+dx*(1.0-dy)*exx[yn,xn+1];
            exxm=exxm+dx*dy*exx[yn+1,xn+1];
            MEXX[mm1,0]=exxm;
            # Calculate Marker SIG'xx from four surrounding nodes
            sxxm=0;
            sxxm=sxxm+(1.0-dx)*(1.0-dy)*sxx2[yn,xn];
            sxxm=sxxm+(1.0-dx)*dy*sxx2[yn+1,xn];
            sxxm=sxxm+dx*(1.0-dy)*sxx2[yn,xn+1];
            sxxm=sxxm+dx*dy*sxx2[yn+1,xn+1];
            # Calculate Marker dSIG'xx from four surrounding nodes
            dsxxm=0;
            dsxxm=dsxxm+(1.0-dx)*(1.0-dy)*dsxx[yn,xn];
            dsxxm=dsxxm+(1.0-dx)*dy*dsxx[yn+1,xn];
            dsxxm=dsxxm+dx*(1.0-dy)*dsxx[yn,xn+1];
            dsxxm=dsxxm+dx*dy*dsxx[yn+1,xn+1];


            # Calculate and save Marker pressure from four surrounding nodes
            prm=0;
            prm=prm+(1.0-dx)*(1.0-dy)*pr1[yn,xn];
            prm=prm+(1.0-dx)*dy*pr1[yn+1,xn];
            prm=prm+dx*(1.0-dy)*pr1[yn,xn+1];
            prm=prm+dx*dy*pr1[yn+1,xn+1];
            MPR[mm1,0]=prm;


            # Interpolating EPSxy=EPSyx from basic nodes
            # Horizontal EPSxy index
            xn=int(xnmin);
            # Vertical EPSxy index
            yn=int(ynmin);

            # Define and check normalized distances from marker to the upper left VX-node;
            dx=(MX[mm1,0]-gridx[xn])/xstp1[xn,0];
            dy=(MY[mm1,0]-gridy[yn])/ystp1[yn,0];

            # Calculate and save Marker EPSxy from four surrounding nodes
            exym=0;
            exym=exym+(1.0-dx)*(1.0-dy)*exy[yn,xn];
            exym=exym+(1.0-dx)*dy*exy[yn+1,xn];
            exym=exym+dx*(1.0-dy)*exy[yn,xn+1];
            exym=exym+dx*dy*exy[yn+1,xn+1];
            MEXY[mm1,0]=exym;
            # Calculate Marker SIGxy from four surrounding nodes
            sxym=0;
            sxym=sxym+(1.0-dx)*(1.0-dy)*sxy2[yn,xn];
            sxym=sxym+(1.0-dx)*dy*sxy2[yn+1,xn];
            sxym=sxym+dx*(1.0-dy)*sxy2[yn,xn+1];
            sxym=sxym+dx*dy*sxy2[yn+1,xn+1];
            # Calculate Marker SIGxy from four surrounding nodes
            dsxym=0;
            dsxym=dsxym+(1.0-dx)*(1.0-dy)*dsxy[yn,xn];
            dsxym=dsxym+(1.0-dx)*dy*dsxy[yn+1,xn];
            dsxym=dsxym+dx*(1.0-dy)*dsxy[yn,xn+1];
            dsxym=dsxym+dx*dy*dsxy[yn+1,xn+1];

            # Computing second strain rate invariant 
            # for the marker using grid values
            eiimg=(MEXX[mm1,0]**2+MEXY[mm1,0]**2)**0.5;
            # Correcting strain rate for the marker using Maxwell model  
            if (eiimg>0):
                # Computing second strain rate invariant for the marker
                # from stresses using Maxwell model
                eiim=((sxxm/2/META[mm1,0]+dsxxm/2/timestep/MMU[MI[mm1,0]])**2+(sxym/2/META[mm1,0]+dsxym/2/timestep/MMU[MI[mm1,0]])**2)**0.5;
                # Computing EiiMarker/EiiGrid ratio
                MRAT[mm1,0]=(eiim/eiimg);
            else:
                MRAT[mm1,0]=1;
    
    return np.asarray(MX,dtype=DTYPE),np.asarray(MY,dtype=DTYPE),np.asarray(MXN,dtype=DTYPE),\
                            np.asarray(MYN,dtype=DTYPE),np.asarray(gridcx,dtype=DTYPE),\
                                np.asarray(gridcy,dtype=DTYPE),np.asarray(xstpc1,dtype=DTYPE),\
                                    np.asarray(ystp1,dtype=DTYPE),timestep,np.asarray(exy,dtype=DTYPE),\
                                        np.asarray(exx,dtype=DTYPE),np.asarray(esp,dtype=DTYPE),\
                                            np.asarray(eii,dtype=DTYPE),np.asarray(MEXX,dtype=DTYPE),\
                                                np.asarray(MPR,dtype=DTYPE),np.asarray(MEXY,dtype=DTYPE),\
                                                    np.asarray(META,dtype=DTYPE),np.asarray(MMU,dtype=DTYPE),np.asarray(MRAT,dtype=DTYPE)

#function for demonstration of markers
def show_markers(xsize,ysize,marknum,MX,MY,MI,MTK,gridx,gridy):
    # Visualizing marker type
    # Pixel grid resolution
    xresol=int(xsize/1000)+1;
    yresol=int(ysize/1000)+1;
    ngrid=2;
    sxstp=xsize/(xresol-1);
    systp=ysize/(yresol-1);
    # Process markers
    markcom=np.nan*np.ones([yresol,xresol]);
    markdis=1e+20*np.ones([yresol,xresol]);
    markgii=np.nan*np.ones([yresol,xresol]);
    markbii=np.nan*np.ones([yresol,xresol]);
    markmel=np.nan*np.ones([yresol,xresol]);
    marktk=np.nan*np.ones([yresol,xresol]);
    for mm1 in range(0,marknum):
        # Define pixel cell
        m1=int(np.int16((MX[mm1,0]-gridx[0])/sxstp-0.5))+1;
        m2=int(np.int16((MY[mm1,0]-gridy[0])/systp-0.5))+1;
        if (m1<0):
            m1=0;
    
        if (m1>xresol-2):
            m1=xresol-2;
    
        if (m2<0):
            m2=0;
        
        if (m2>yresol-2):
            m2=yresol-2;
        
        # Define indexes of surrounding pixels
        m10min=m1-ngrid;
        if (m10min<0):
            m10min=0;
        
        m10max=m1+1+ngrid;
        if (m10max>xresol-1):
            m10max=xresol-1;
        
        m20min=m2-ngrid;
        if (m20min<0):
            m20min=0;
        
        m20max=m2+1+ngrid;
        if (m20max>yresol-1):
            m20max=yresol-1;
        
        # Update pixels around the marker
        for m10 in range(m10min,m10max+1):
            for m20 in range(m20min,m20max+1): 
                # Check distance to current pixel
                dx=(MX[mm1,0]-gridx[0])-(m10-1)*sxstp;
                dy=(MY[mm1,0]-gridy[0])-(m20-1)*systp;
                dd=(dx*dx+dy*dy)**0.5;
                if(dd<markdis[m20,m10]):
                    markcom[m20,m10]=MI[mm1,0];
                    marktk[m20,m10]=MTK[mm1,0]; #check it
                
    # Draw composition
    markcom[0,0:9]=np.arange(0,9);

    # Save pic file
    fig3=plt.figure();
    
    print(get_now()+'Visualize test markers distribution')   
    #plt.imshow(pr1,alpha=1);
    #plt.imshow(res_tk1,alpha=1,cmap='jet');
    plt.imshow(marktk,alpha=1,cmap='jet');
    plt.colorbar();
    plt.imshow(markcom,alpha=0.75,cmap='gray');
    plt.title('marker distribution image');
    plt.savefig('out/pic_mark_{}.png'.format('test'),dpi=150);
    #plt.savefig('out/pic1_{}.svg'.format(ntimestep),dpi=150);
    plt.close();

    
def custom_print(message_to_print, log_file=start_params['log_file']):
    #print(message_to_print)
    if not os.path.isdir(os.path.dirname(log_file)):
        os.mkdir(os.path.dirname(log_file))
    with open(log_file, 'a') as of:
        of.write(message_to_print + '\n')

def sec_to_mdhs(sec): #service.sec_to_mdhs(1000*365.25*24*3600)
    years=sec//(365.25*24*3600);
    month=sec%365.25*24*3600//(29.5*24*3600);
    day=((sec%(365.25*24*3600))%(29.5*24*3600))//(24*3600);
    hours=sec%(365.25*24*3600)%(29.5*24*3600)%(24*3600)//3600;
    mins=((((sec%(365.25*24*3600)%(29.5*24*3600))%(24*3600))%3600)//60);
    secs=(((((sec%(365.25*24*3600))%(29.5*24*3600))%(24*3600))%3600)%60);
    res_string='Years: '+str(years)+' Month: '+str(month)+' Day: ' + str(day)+\
               ' Hours: '+str(hours)+' Minute: '+str(mins)+' Seconds: '+str(secs);       
    return res_string
###############end of functions block


#parse command line launch arguments
#cdef dict start_params={};
cdef str arg,key,val;

if len(sys.argv)>1:
    for arg in sys.argv[1:]:
        arg=arg.replace('-','');
        key,val=arg.split('=');
        start_params.update({key:val});

custom_print(get_now()+'Computations were started...');

cdef int start_time=time.time();
cdef str out_folder=start_params['outdirname'];
#cdef str log_file='log.out';
cdef str vis_on_off='off'; #show pictures or not

if not os.path.exists(out_folder): #if no folder
    os.mkdir(out_folder);


#if not os.path.isfile(log_file): #if log.out exists, delete it
#    os.remove(log_file);

custom_print(get_now()+'Assigning variables...');

#TODO load model with command line argument
#filename = 'test_magmatic2.cfg'
if not os.path.isfile(start_params['modelfile']):
    raise FileNotFoundError('Can not find JSON model file! Stopped')
#get parameters from model 
modelPhysics,bm = importModelJson(start_params['modelfile'])

# Acceleration of Gravity, m/s^2
cdef double gx=modelPhysics['gx'] #0;
cdef double gy=modelPhysics['gy'] #9.81;

# Temperature at the top, and bottom of the model, K
cdef int ttop=modelPhysics['ttop'] #273;
cdef int tbottom=modelPhysics['tbottom']# 1600;

# Gas constant J/mol/K
cdef double RGAS=8.314;

# Initial model size, m
cdef int xsize0,xsize,ysize0,ysize,ntimestep;
xsize0=modelPhysics['widthkm']*1000 ##100000;  #TODO convert from km of the visual editor
ysize0=modelPhysics['heightkm']*1000 ##120000;
xsize=xsize0;
ysize=ysize0;

# Initial water level
cdef int waterlev0,waterlev;
waterlev0=modelPhysics['seadepth']*1000 #8000;
waterlev=waterlev0;

# Defining grid resolution
cdef int xnum,ynum;
#xnum=201;
#ynum=61;
xnum=modelPhysics['sizepx'][1];
ynum=modelPhysics['sizepx'][0];


# Viscosity limits for rocks, Pa
cdef double etamin,etamax,etamelt,stressmin;
etamin=1e+14;   # Lower limit, Pa
etamax=1e+25;   # Upper limit, Pa
# Partally molten rock viscosity 
etamelt=1e+14; # Pa s
# Lower stress limit for power law, Pa
stressmin=1e+4;

# Viscoelastic timestep, s
cdef double timemax,markmax;
timemax=1e+2*365.25*24*3600; # 100 year
# Maximal marker displacement step, number of gridsteps
markmax=0.1;
# Moving Markers: 
# 0 = not moving at all
# 1 = simple advection
# 4 = 4-th order in space  Runge-Kutta
cdef int markmove=1;
# Velocity calculation
# 0 = by Solving momentum and continuity equations
# 1 = solid body rotation
cdef int movemod=0;
# Maximal temperature change, allowed for one timestep, K
cdef double tempmax=30;
# Amount of timesteps
cdef int stepmax=int(start_params['stepmax']);
# Saving mat file periodicity
cdef int savestep=int(start_params['savestep']);
cdef int savepicstep=int(start_params['savepicstep']);

cdef int dislocation_type=0; #0 - horizontal layers(no dislocation); 1 - monocline; 2 - folded layers
cdef int melt_type=1; #0- ultramafic (olivine); 1 - basaltic; 2-diorite; 3 - granitic
cdef str melt_word='basalt'

custom_print(get_now()+"Assigning Topography model...")

# Topography model
# Topography model size in horizontal direction
cdef int tsize;
tsize=xsize;
# Defining topography model resolution
cdef int tnum=301;
# Grid for topography profile
cdef cnp.ndarray gridt=np.zeros([6,tnum],dtype=DTYPE); #took 5 row matrix, cause it couldnt be done dynamically as in Matlab
cdef double tstp=tsize/(tnum-1); # topography grid step
cdef int i;

gridt[0,0]=0;
for i in range(1,tnum):
    gridt[0,i]=gridt[0,i-1]+tstp;

# Topography diffusion koefficient Ks, m^2/s
# dYt/dt=Ks*d2Yt/dx^2
# Define dYt/dt - erosion rate
cdef double dYtdt=1/(1000*365.25*24*3600); # 1 mm/yr
# Define d2Yt - max elevation
cdef int d2Yt=10*1000; # 10 km
# Define dx - transport lengthscale
cdef double dx=100*1000; # 100 km
cdef double Ks=dYtdt*dx**2/d2Yt
cdef cnp.ndarray topotime=np.zeros([stepmax,1],dtype=DTYPE);
cdef cnp.ndarray topohigh=np.zeros([stepmax,tnum],dtype=DTYPE);
cdef cnp.ndarray topowater=np.zeros([stepmax,1],dtype=DTYPE);

# Material properties
# MRHO = density (kg/m3): RHO*[1-ALP*(T-273)]*[1+BET*(P-1e+5)]
# MFLOW = power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
# MMU = shear modulus (Pa)
# MPL = Brittle/plastic strength (Pa): SIGMAyeild=C+sin(FI)*P
#       C=C0, FI=FI0 for strain<=GAM0
#       C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
#       C=C1, FI=FI1 for strain>=GAM0
# MCP = heat capacity (J/K/kg)
# MKT = thermal conductivity (W/m/K): k=k0+a/(T+77) 
# MHR = radiogenic heat production (W/m^3) 

#Materials
custom_print(get_now()+' Declare matrix for materials')
#Declare matrix for materials
cdef cnp.ndarray MRHO,MFLOW,MMU,MPL,MCP,MKT,MHR;

MRHO=np.zeros([10,5],dtype=DTYPE);
MFLOW=np.zeros([10,6],dtype=DTYPE);
MMU=np.zeros([10],dtype=DTYPE);
MPL=np.zeros([10,7],dtype=DTYPE);
MCP=np.zeros([10],dtype=DTYPE);
MKT=np.zeros([10,3],dtype=DTYPE);
MHR=np.zeros([10],dtype=DTYPE);

# 1 = Weak Layer ("sticky air/water")
MRHO[1,1]=1000;             # standard density, kg/m^3
MRHO[1,2]=0;                # thermal expansion, 1/K
MRHO[1,3]=0;                # compressibility, 1/Pa
MFLOW[1,1]=0;               # 0=constant viscosity
MFLOW[1,2]=1e+13;           # viscosity, Pa s
MMU[1]=1e+20;               # shear modulus, Pa
MPL[1,1]=0;                 # C0, Pa
MPL[1,2]=0;                 # C1, Pa
MPL[1,3]=0;                 # sin(FI0)
MPL[1,4]=0;                 # sin(FI1)
MPL[1,5]=0;                 # GAM0
MPL[1,6]=1;                 # GAM1
MCP[1]=3000;                # Cp, J/kg
MKT[1,1]=300;               # k0, W/m/K
MKT[1,2]=0;                 # a, W/m
MHR[1]=0;                   # radiogenic heat production, W/m^3
# 2 = Sediments
MRHO[2,1]=2700;             # standard density, kg/m^3
MRHO[2,2]=3e-5;             # thermal expansion, 1/K
MRHO[2,3]=1e-11;            # compressibility, 1/Pa
MRHO[2,4]=2400;             # melt density, kg/m^3
MFLOW[2,1]=1;               # 1=power law [wet quartzite: Ranalli, 1995]
MFLOW[2,2]=3.2e-4;          # AD, 1/s/MPa^n
MFLOW[2,3]=2.3;             # n
MFLOW[2,4]=154;             # Ea, kJ/mol
MFLOW[2,5]=0;               # Va, cm^3
MMU[2]=1e+10;               # shear modulus, Pa
MPL[2,1]=1e+6;              # C0, Pa
MPL[2,2]=1e+6;              # C1, Pa
MPL[2,3]=0.20;                 # sin(FI0)
MPL[2,4]=0.00;                 # sin(FI1)
MPL[2,5]=0;                 # GAM0
MPL[2,6]=0.1;                 # GAM1
MCP[2]=1000;                # Cp, J/kg
MKT[2,1]=0.64;              # k0, W/m/K
MKT[2,2]=807;               # a, W/m
MHR[2]=2.0e-6;              # radiogenic heat production, W/m^3
# 3 = Basalts
MRHO[3,1]=3000;             # standard density, kg/m^3
MRHO[3,2]=3e-5;             # thermal expansion, 1/K
MRHO[3,3]=1e-11;            # compressibility, 1/Pa
MRHO[3,4]=2400;             # melt density, kg/m^3
MFLOW[3,1]=1;               # 1=power law [wet quartzite: Ranalli, 1995]
MFLOW[3,2]=3.2e-4;          # AD, 1/s/MPa^n
MFLOW[3,3]=2.3;             # n
MFLOW[3,4]=154;             # Ea, kJ/mol
MFLOW[3,5]=0;               # Va, cm^3
MMU[3]=2.5e+10;             # shear modulus, Pa
MPL[3,1]=1e+6;              # C0, Pa
MPL[3,2]=1e+6;              # C1, Pa
MPL[3,3]=0.00;                 # sin(FI0)
MPL[3,4]=0.00;                 # sin(FI1)
MPL[3,5]=0;                 # GAM0
MPL[3,6]=0.1;                 # GAM1
MCP[3]=1000;                # Cp, J/kg
MKT[3,1]=1.18;              # k0, W/m/K
MKT[3,2]=474;               # a, W/m
MHR[3]=2.5e-7;              # radiogenic heat production, W/m^3
# 4 = Gabbro
MRHO[4,1]=3000;             # standard density, kg/m^3
MRHO[4,2]=3e-5;             # thermal expansion, 1/K
MRHO[4,3]=1e-11;            # compressibility, 1/Pa
MRHO[4,4]=2700;             # melt density, kg/m^3
MFLOW[4,1]=1;               # 1=power law [plagioclase An75: Ranalli, 1995]
MFLOW[4,2]=3.3e-4;          # AD, 1/s/MPa^n
MFLOW[4,3]=3.2;             # n
MFLOW[4,4]=238;             # Ea, kJ/mol
MFLOW[4,5]=0;               # Va, cm^3
MMU[4]=2.5e+10;             # shear modulus, Pa
MPL[4,1]=1e+6;              # C0, Pa
MPL[4,2]=1e+6;              # C1, Pa
MPL[4,3]=0.2;               # sin(FI0)
MPL[4,4]=0.00;               # sin(FI1)
MPL[4,5]=0;                 # GAM0
MPL[4,6]=0.1;                 # GAM1
MCP[4]=1000;                # Cp, J/kg
MKT[4,1]=1.18;              # k0, W/m/K
MKT[4,2]=474;               # a, W/m
MHR[4]=2.5e-7;              # radiogenic heat production, W/m^3
# 5 = Lithospheric mantle
MRHO[5,1]=3300;             # standard density, kg/m^3
MRHO[5,2]=3e-5;             # thermal expansion, 1/K
MRHO[5,3]=1e-11;            # compressibility, 1/Pa
MRHO[5,4]=2700;             # melt density, kg/m^3
MFLOW[5,1]=1;               # 1=power law [dry olivine: Ranalli, 1995]
MFLOW[5,2]=2.5e+4;          # AD, 1/s/MPa^n
MFLOW[5,3]=3.5;             # n
MFLOW[5,4]=532;             # Ea, kJ/mol
MFLOW[5,5]=10;               # Va, cm^3
MMU[5]=6.7e+10;             # shear modulus, Pa
MPL[5,1]=1e+6;              # C0, Pa
MPL[5,2]=1e+6;              # C1, Pa
MPL[5,3]=0.6;               # sin(FI0]
MPL[5,4]=0.00;               # sin(FI1]
MPL[5,5]=0;                 # GAM0
MPL[5,6]=0.1;                 # GAM1
MCP[5]=1000;                # Cp, J/kg
MKT[5,1]=0.73;              # k0, W/m/K
MKT[5,2]=1293;              # a, W/m
MHR[5]=2.2e-8;              # radiogenic heat production, W/m^3
# 6 = Asthenospheric mantle
MRHO[6,1]=3300;             # standard density, kg/m^3
MRHO[6,2]=3e-5;             # thermal expansion, 1/K
MRHO[6,3]=1e-11;            # compressibility, 1/Pa
MRHO[6,4]=2700;             # melt density, kg/m^3
MFLOW[6,1]=1;               # 1=power law (dry olivine: Ranalli, 1995]
MFLOW[6,2]=2.5e+4;          # AD, 1/s/MPa^n
MFLOW[6,3]=3.5;             # n
MFLOW[6,4]=532;             # Ea, kJ/mol
MFLOW[6,5]=10;               # Va, cm^3
MMU[6]=6.7e+10;             # shear modulus, Pa
MPL[6,1]=1e+6;              # C0, Pa
MPL[6,2]=1e+6;              # C1, Pa
MPL[6,3]=0.6;               # sin(FI0]
MPL[6,4]=0.00;               # sin(FI1)
MPL[6,5]=0;                 # GAM0
MPL[6,6]=0.1;                 # GAM1
MCP[6]=1000;                # Cp, J/kg
MKT[6,1]=0.73;              # k0, W/m/K
MKT[6,2]=1293;              # a, W/m
MHR[6]=2.2e-8;              # radiogenic heat production, W/m^3
# 7 = Hydrated mantle
MRHO[7,1]=3300;             # standard density, kg/m^3
MRHO[7,2]=3e-5;             # thermal expansion, 1/K
MRHO[7,3]=1e-11;            # compressibility, 1/Pa
MRHO[7,4]=2700;             # melt density, kg/m^3
MFLOW[7,1]=1;               # 1=power law (wet olivine: Ranalli, 1995]
MFLOW[7,2]=2.0e+3;          # AD, 1/s/MPa^n
MFLOW[7,3]=4.0;             # n
MFLOW[7,4]=471;             # Ea, kJ/mol
MFLOW[7,5]=0;               # Va, cm^3
MMU[7]=6.7e+10;             # shear modulus, Pa
MPL[7,1]=1e+6;              # C0, Pa
MPL[7,2]=1e+6;              # C1, Pa
MPL[7,3]=0.0;                 # sin(FI0)
MPL[7,4]=0.0;                 # sin(FI1)
MPL[7,5]=0;                 # GAM0
MPL[7,6]=0.1;                 # GAM1
MCP[7]=1000;                # Cp, J/kg
MKT[7,1]=0.73;              # k0, W/m/K
MKT[7,2]=1293;              # a, W/m
MHR[7]=2.2e-8;              # radiogenic heat production, W/m^3
# 8 = Upper continental crust (granodiorite]
MRHO[8,1]=2700;             # standard density, kg/m^3
MRHO[8,2]=3e-5;             # thermal expansion, 1/K
MRHO[8,3]=1e-11;            # compressibility, 1/Pa
MRHO[8,4]=2400;             # melt density, kg/m^3
MFLOW[8,1]=1;               # 1=power law [wet quartzite: Ranalli, 1995]
MFLOW[8,2]=3.2e-4;          # AD, 1/s/MPa^n
MFLOW[8,3]=2.3;             # n
MFLOW[8,4]=154;             # Ea, kJ/mol
MFLOW[8,5]=0;               # Va, cm^3
MMU[8]=1e+10;               # shear modulus, Pa
MPL[8,1]=1e+6;              # C0, Pa
MPL[8,2]=1e+6;              # C1, Pa
MPL[8,3]=0.2;               # sin(FI0)
MPL[8,4]=0.00;              # sin(FI1)
MPL[8,5]=0;                 # GAM0
MPL[8,6]=0.1;               # GAM1
MCP[8]=1000;                # Cp, J/kg
MKT[8,1]=0.64;              # k0, W/m/K
MKT[8,2]=807;               # a, W/m
MHR[8]=1.0e-6;              # radiogenic heat production, W/m^3
# 9 = Lower continental crust [diorite]
MRHO[9,1]=3000;             # standard density, kg/m^3
MRHO[9,2]=3e-5;             # thermal expansion, 1/K
MRHO[9,3]=1e-11;            # compressibility, 1/Pa
MRHO[9,4]=2700;             # melt density, kg/m^3
MFLOW[9,1]=1;               # 1=power law [plagioclase An75: Ranalli, 1995]
MFLOW[9,2]=3.3e-4;          # AD, 1/s/MPa^n
MFLOW[9,3]=3.2;             # n
MFLOW[9,4]=238;             # Ea, kJ/mol
MFLOW[9,5]=0;               # Va, cm^3
MMU[9]=2.5e+10;             # shear modulus, Pa
MPL[9,1]=1e+6;              # C0, Pa
MPL[9,2]=1e+6;              # C1, Pa
MPL[9,3]=0.2;               # sin(FI0)
MPL[9,4]=0.00;              # sin(FI1)
MPL[9,5]=0;                 # GAM0
MPL[9,6]=0.1;               # GAM1
MCP[9]=1000;                # Cp, J/kg
MKT[9,1]=1.18;              # k0, W/m/K
MKT[9,2]=474;               # a, W/m
MHR[9]=5.0e-7;              # radiogenic heat production, W/m^3

# Numerical Subgrid stress diffusion coefficient
cdef double dsubgrids=1;
# Numerical Subgrid temperature diffusion coefficient
cdef double dsubgridt=1;
# Shear heating on(1)/off(0)
cdef double frictyn=1;
# Adiabatic heating on(1)/off(0)
cdef double adiabyn=1;
# Weight for old visco-plastic viscosity
cdef double etawt=0;

# Pressure boundary conditions
# prfirst(1) = boundary condition mode:
# 0 - pressure in one cell definition
# 1 - pressure at the top and in the bottom of the channel
#cdef cnp.ndarray prfirst=np.zeros(2,dtype=DTYPE);
cdef double[:] prfirst=np.zeros(2,dtype=DTYPE);
prfirst[0]=0;
# prfirst(2) = boundary condition value
prfirst[1]=1e+5;

custom_print(get_now()+'Velocity Boundary condition specified by bleft,bright,btop,bbot ...')
# Velocity Boundary condition specified by bleft,bright,btop,bbot 
# are implemented from ghost nodes 
# directly into Stokes and continuity equations

# #These patterns could be adjusted for different geodynamic conditions
# # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
# # Upper boundary: Free slip
# # vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
# cdef cnp.ndarray stensil0=np.array([0,1,0,0],dtype=DTYPE);
# cdef cnp.ndarray btop=np.asarray(np.tile(stensil0,(xnum+1,1)),dtype=DTYPE);
# # Lower boundary: Free Slip  
# # vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
# cdef cnp.ndarray bbottom=np.asarray(np.tile(stensil0,(xnum+1,1)),dtype=DTYPE);

# custom_print(get_now()+'Assing marker arrays...')
# # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)

# cdef cnp.ndarray stensil1=np.array([0,0,0,1],dtype=DTYPE);
# # Left boundary: Free slip   
# # vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
# cdef cnp.ndarray bleft=np.asarray(np.tile(stensil1,(ynum+1,1)),dtype=DTYPE);
# # Right boundary: Free slip 
# # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
# cdef cnp.ndarray bright=np.asarray(np.tile(stensil1,(ynum+1,1)),dtype=DTYPE);

# Horizontal collision velocity, m/s
cdef double velcol=modelPhysics['horcolvel']; # 2 cm/yr
# Horizontal extension velocity, m/s
cdef double velext=modelPhysics['horextvel']; # 1 cm/yr
# Horizontal extension strain rate, dvx/dx, 1/s
cdef double epsext=10**modelPhysics['horextstrainpow'];

#TODO stationary
if modelPhysics['tmode'] == 'stationary':
    # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
    # Upper boundary: Free slip
    stensil00=np.array([0,1,0,0],dtype=DTYPE);
    stensil01=np.array([0,1,0,0],dtype=DTYPE);
    btop=np.asarray(np.tile(stensil00,(xnum+1,1)),dtype=DTYPE);
    # Lower boundary: Free Slip  
    bbottom=np.asarray(np.tile(stensil01,(xnum+1,1)),dtype=DTYPE);
    # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
    stensil10=np.array([0,0,0,1],dtype=DTYPE);
    stensil11=np.array([0,0,0,1],dtype=DTYPE);
    # Left boundary: Free slip   
    bleft=np.asarray(np.tile(stensil10,(ynum+1,1)),dtype=DTYPE);
    # Right boundary: Free slip 
    # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright=np.asarray(np.tile(stensil11,(ynum+1,1)),dtype=DTYPE);

#TODO extension
elif modelPhysics['tmode'] == 'extension':
   # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
    stensil00=np.array([0,1,0,0],dtype=DTYPE)
    btop=np.asarray(np.tile(stensil00,(xnum+1,1)),dtype=DTYPE);
    # Lower boundary: Free Slip  
    # vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    stensil01=np.array([0,1,-velext/xsize*ysize,0],dtype=DTYPE)
    bbottom=np.asarray(np.tile(stensil01,(xnum+1,1)),dtype=DTYPE);  

    # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
    stensil10=np.array([-velext/2,0,0,1],dtype=DTYPE);
    # Left boundary: Free slip   
    # vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft=np.asarray(np.tile(stensil10,(ynum+1,1)),dtype=DTYPE);
    # Right boundary: Free slip 
    stensil11=np.array([velext/2,0,0,1],dtype=DTYPE);
    # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright=np.asarray(np.tile(stensil11,(ynum+1,1)),dtype=DTYPE);

#TODO collision
elif modelPhysics['tmode'] == 'collision':
    
    # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
    stensil00=np.array([0,1,0,0],dtype=DTYPE)
    btop=np.asarray(np.tile(stensil00,(xnum+1,1)),dtype=DTYPE);
    # Lower boundary: Free Slip  
    # vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    stensil01=np.array([0,1,velcol/xsize*ysize,0],dtype=DTYPE)
    bbottom=np.asarray(np.tile(stensil01,(xnum+1,1)),dtype=DTYPE);  

    # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
    stensil10=np.array([0,0,0,1],dtype=DTYPE);
    # Left boundary: Free slip   
    # vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft=np.asarray(np.tile(stensil10,(ynum+1,1)),dtype=DTYPE);
    # Right boundary: Free slip 
    stensil11=np.array([-velcol,0,0,1],dtype=DTYPE);
    # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright=np.asarray(np.tile(stensil11,(ynum+1,1)),dtype=DTYPE); #TODO lateral pressure from left to right

#TODO
elif modelPhysics['tmode'] == 'subduction':
    
    # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
    stensil00=np.array([0,1,epsext*ysize/2,0],dtype=DTYPE)
    btop=np.asarray(np.tile(stensil00,(xnum+1,1)),dtype=DTYPE);
    # Lower boundary: Free Slip  
    # vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
    stensil01=np.array([0,1,-epsext*ysize/2,0],dtype=DTYPE)
    bbottom=np.asarray(np.tile(stensil01,(xnum+1,1)),dtype=DTYPE);  

    # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
    stensil10=np.array([-epsext*xsize/2,0,0,1],dtype=DTYPE);
    # Left boundary: Free slip   
    # vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
    bleft=np.asarray(np.tile(stensil10,(ynum+1,1)),dtype=DTYPE);
    # Right boundary: Free slip 
    stensil11=np.array([epsext*xsize/2,0,0,1],dtype=DTYPE);
    # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright=np.asarray(np.tile(stensil11,(ynum+1,1)),dtype=DTYPE);

else: #default regime - STATIC
    # Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
    # Upper boundary: Free slip
    stensil00=np.array([0,1,0,0],dtype=DTYPE);
    stensil01=np.array([0,1,0,0],dtype=DTYPE);
    btop=np.asarray(np.tile(stensil00,(xnum+1,1)),dtype=DTYPE);
    # Lower boundary: Free Slip  
    bbottom=np.asarray(np.tile(stensil01,(xnum+1,1)),dtype=DTYPE);
    # Left, Right boundaries: + Prescribed outward velocity (horizontal extension)
    stensil10=np.array([0,0,0,1],dtype=DTYPE);
    stensil11=np.array([0,0,0,1],dtype=DTYPE);
    # Left boundary: Free slip   
    bleft=np.asarray(np.tile(stensil10,(ynum+1,1)),dtype=DTYPE);
    # Right boundary: Free slip 
    # vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
    bright=np.asarray(np.tile(stensil11,(ynum+1,1)),dtype=DTYPE);




# Internal boundary condition: prescribed velocity of "mobile wall"
cdef cnp.ndarray bintern=np.zeros(8,dtype=DTYPE);
bintern[0]=-1;      # Horizontal position of vx nodes with prescrbed velocity (no susch condition if negative)
bintern[1]=0;       # Min vertical position
bintern[2]=0;       # Max vertical position
bintern[3]=-0;      # Prescribed shortening velocity, m/s 
bintern[4]=-1;      # Horizontal position of vy nodes with prescrbed velocity (no susch condition if negative)
bintern[5]=0;       # Min vertical position
bintern[6]=0;       # Max vertical position
bintern[7]=0;       # Prescribed vertical velocity, m/s

# Defining average initial gridsteps
cdef double xstp=xsize/(xnum-1);
cdef double ystp=ysize/(ynum-1);

# Defining gridline positions for irregular basic grid

# Horizontal grid
cdef cnp.ndarray gridx=np.zeros([xnum],dtype=DTYPE);
#for i=2:1:xnum:
#    gridx(i)=gridx(i-1)+xstp;

#gridx=np.arange(0,xnum*xstp,xstp).reshape(len(gridx),1)   #??????????
gridx=np.arange(0,xnum*xstp,xstp)

# Vertical grid
cdef int b,nn;
b=1000; # grid spacing in high resolution area
nn=31;  # number of nodes in high resolution area

cdef cnp.ndarray gridy=np.zeros([ynum],dtype=DTYPE);
# Define regular step in high resolution area
for i in range(0,nn):
    #gridy[i,0]=gridy[i-1,0]+b;
    gridy[i]=i*b;

# Define factor of grid spacing increase from the bottom
# of high resolution area
cdef float D=ysize-gridy[nn-1]; # distance to be covered by non-uniform grid
cdef int N=ynum-nn; # number of grid steps to be used in the grid
# Iterative search of F
cdef double F=1.1;
for i in range(0,100):
    F=(1+D/b*(1-1/F))**(1/N);

# Define position of nodal points
for i in range(nn,ynum):
    gridy[i-1]=gridy[i-2]+b*F**(i-nn);

gridy[ynum-1]=ysize;

# Thermal boundary conditions
# Upper, Lower boundaries: constant temperature

#TODO adjust boundary conditions regarding the file Collision (there is shortening) or Extension (extension )

cdef cnp.ndarray btopt=np.zeros([xnum,2],dtype=DTYPE);
cdef cnp.ndarray bbottomt=btopt.copy();
cdef int j;
for j in range(0,xnum):
    # Upper boundary
    # tk(1,j)=btopt(j,1)+tk(2,j)*btop(j,2)
    btopt[j,0]=ttop;
    btopt[j,1]=0;
    # Lower boundary
    # tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
    bbottomt[j,0]=tbottom;
    bbottomt[j,1]=0;

cdef cnp.ndarray bleftt=np.zeros([ynum,2],dtype=DTYPE)
cdef cnp.ndarray brightt=bleftt.copy();

# Left, Right boundaries: symmetry
for i in range(0,ynum):
    # Left boundary
    # tk(i,1)=bleftt(i,1)+bleftt(i,2)*tk(i,2);
    bleftt[i,0]=0;
    bleftt[i,1]=1;
    # Right boundary
    # tk(i,xnum)=brightt(i,1)+brightt(i,2)*tk(i,xnum-1);
    brightt[i,0]=0;
    brightt[i,1]=1;

#TODO initialize model workflow

# Defining number of markers and steps between them in the horizontal and vertical direction
cdef int mxnum,mynum,marknum;
cdef double mxstep,mystep;

mxnum=400; #number of markers in horizontal direction
mxstep=xsize/mxnum; #step between markers in horizontal direction
mynum=480;  #number of markers in vertical direction   
mystep=ysize/mynum; #step between markers in vertical direction
marknum=mxnum*mynum

# Creating markers arrays
cdef cnp.ndarray MX,MY,MTK,MXN,MYN,MSXX,MSXY,META,MEXX,MEXY,MPR,MBII,MGII,MRAT,MXM,MI 

cdef int layer_num=10;

MX=np.zeros([mynum*mxnum,1],dtype=DTYPE);   # X coordinate, m
MY=np.zeros([mynum*mxnum,1],dtype=DTYPE);   # Y coordinate, m
MTK=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # Temperature, K
MI=np.zeros([mynum*mxnum,1],dtype=DTYPE2);   # Type
MXN=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # Horizontal index
MYN=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # Vertical index
MSXX=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # SIGMAxx - deviatoric normal stress, Pa
MSXY=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # SIGMAyy - shear stress, Pa
META=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # viscosity, Pa s
MEXX=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # EPSILONxx - normal strain rate, 1/s
MEXY=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # EPSILONyy - shear strain rate, 1/s
MPR=np.zeros([mynum*mxnum,1],dtype=DTYPE);   # Pressure, Pa
MBII=np.zeros([mynum*mxnum,1],dtype=DTYPE);  # Accumulated bulk strain
MGII=-1e-20*np.ones([mynum*mxnum,1],dtype=DTYPE);  # Accumulated plastic strain
MRAT=np.ones([mynum*mxnum,1],dtype=DTYPE);   # EiiMarker/EiiGrid Ratio
MXM=np.zeros([mynum*mxnum,1],dtype=DTYPE);   # Melt fraction

cdef int is_reloaded=0 #default - no save file found
cdef str is_load='n'
if os.path.exists(os.path.join(start_params['outdirname'],'save.p')):
    #if 'is_load' not in [*key_val_dict]: #check command line input
    if 'is_load' not in [*start_params]: #check command line input
        while True:
          custom_print(get_now()+'Save file was found, want to load it? y/n...');  
          is_load=input('Save file was found, want to load it? y/n...').lower()
          if is_load in ['y','n']:
              break;
          else:
              print('You need to enter exactly y/n')
              continue;
    else:
       is_load='y' if (start_params['is_load']=='1') or (start_params['is_load']=='y') else 'n'
else:
    ntimestep = 1



#TODO declare more variables, loaded with pickle load

if is_load=='y':
    custom_print(get_now()+'Loading...')
    
    (timestep,timestept,timesum,ntimestep,MX,MY,MTK,MI,MXN,MYN,MSXX,MSXY,META,MEXX,MEXY,MPR,MGII,MRAT,\
    gridt,topotime,topohigh,topowater,prfirst,etas0,etan0,etas1,etan1,gridx,gridy,\
    mus0,mus1,mun0,mun1,sxy0,sxx0,sxy1,sxx1,rho0,rho1,tk0,tk2,rhocp0,rhocp1,\
    kt0,kt1,hr0,hr1,ha0,ha1)= pickle.load(open(os.path.join(start_params['outdirname'],'save.p'),"rb"))
    
    custom_print(get_now()+'Save file sucesfully loaded...');
    custom_print(get_now()+f'Reloaded ntimestep={ntimestep}');
    is_reloaded=1
else:
    custom_print(get_now()+'Save file was not foung or not loaded...');
    #TODO setup dislocation and melt type
    if 'dislocation_type' not in [*start_params]: #check command line input
        custom_print(get_now()+'Please enter preferences manually...');
        dislocation_type=int(input('Please, enter dislocation type (0 - horizontal layers(no dislocation); 1 - monocline; 2 - folded layers)'))
    else: 
        dislocation_type=int(start_params['dislocation_type'])
    
    if 'melt_type' not in [*start_params]: #check command line input
        custom_print(get_now()+'Please enter preferences manually...');
        melt_type=int(input('Please, enter melt type (0- ultramafic (olivine); 1 - basaltic; 2-diorite; 3 - granitic)'))
    else: 
        dislocation_type=int(start_params['melt_type'])
        
    #dislocation 
    #switch dislocation_type
    if dislocation_type==0:
        custom_print(get_now()+'Horizontal layering was choosen.')
    elif dislocation_type==1:
        custom_print(get_now()+'Monocline dislocation was choosen.')
    elif dislocation_type==2:
        custom_print(get_now()+'Folded dislocation was choosen.')
    else:
        custom_print(get_now()+'Number for dislocaiton not found.')
        custom_print(get_now()+'Set horizontal as default.')
        dislocation_type=0
    #melt 
    #switch melt_type
    if melt_type==0:
        custom_print(get_now()+'Olivine melt was choosen.')
        melt_word='olivine'
    elif melt_type==1:
        custom_print(get_now()+'Basaltic metl was choosen.')
        melt_word='basalt'
    elif melt_type==2:
        custom_print(get_now()+'Diorite melt was choosen.')
        melt_word='diorite'
    elif melt_type==3:
        custom_print(get_now()+'Granitic melt was choosen.')
        melt_word='granite'
    else:
        custom_print(get_now()+'Number for melt was not found.')
        custom_print(get_now()+'Set basaltic as default')
        melt_type=1   
        melt_word='basalt'
    # Initial time, s
    ntimestep=1; #start time from beginning
    timesum=0;
    timestepd=copy.copy(timemax);
    ntimestep=1;

cdef int dip_angle, h, n, folds, resolution,start_depth,end_depth,thickness,height,num_layers;
cdef double hv,length;
cdef cnp.ndarray sine;
cdef list pnts_poly,layers;

#TODO SEDIMENTS learn how to apply these structures to your polygons
if is_reloaded==0: #no save file found - need to redefine everything since the start
    custom_print(get_now()+'Defining intial position of markers...')
      
    
    
    # TODO Defining intial position of markers
    #markers are distributed uniformly. Need to check their position within poygons 
    # TODO convert pixels of marker coordinate into percent and find wchich marker they belong to
    # Defining lithological structure of the model
    
    # filename = 'test_magmatic.cfg'
    # modelPhysics,bm = importModelJson(filename)

    #types conversion is necessary
    mpxxstep:float=float(xnum)/float(mxnum) #step between markers in horizontal direction 
    mpxystep:float=float(ynum)/float(mynum)

    print('modelPhysics',modelPhysics)
    print('modelPhysics[sizepx][0]',modelPhysics['sizepx'][0])

    model = drawModel4(xnum,ynum,bm)

    print(np.max(np.max(model)))
    print(model)
    print("min model=",np.min(np.min(model)))

    if np.min(np.min(model))<=0:
        custom_print(get_now()+'replace 0 and negative values with default marker type (6 = Asthenosphere)')
        model[model<=0]=6 


    mm1=0;
    for xm in range(1,mxnum+1):
        for ym in range(1,mynum+1):
            #Coordinates with small random displacement
            MX[mm1,0]=xm*mxstep-mxstep/2+(np.random.rand()-0.5)*mxstep
            MY[mm1,0]=ym*mystep-mystep/2+(np.random.rand()-0.5)*mystep
            #default marker type (6 = Asthenosphere)
            MI[mm1,0]=6;

            # print('xm=',xm)
            # print('ym=',ym)
            
            #time.sleep(0.05)

            Mindex_x = int(xm*mpxxstep-mpxxstep/2+(np.random.rand()-0.5)*mpxxstep)
            Mindex_y = int(ym*mpxystep-mpxystep/2+(np.random.rand()-0.5)*mpxystep)

            # print('x=',Mindex_x)
            # print('y=',Mindex_y)
            # print('val=',model[Mindex_y,Mindex_x])

            MI[mm1,0] = model[Mindex_y,Mindex_x]
            
            #TODO Initial temperature structure
            # Adiabatic temperature gradient in the asthenosphere = 0.5 K/km
            dtdy=0.5/1000; # K/m
            MTK[mm1,0]=tbottom-dtdy*(ysize-MY[mm1,0]);
            # Sticky water
            if(MI[mm1,0]==1):
                MTK[mm1,0]=ttop;
        
            # TODO Linear sectionned geotherm  
            yast=modelPhysics['lithbottom']*1000; # Bottom of the lithosphere
            tast=tbottom-dtdy*(ysize-yast); # T (K) of asthenosphere at y=yast
            ymoho=modelPhysics['ymohokm']*1000; # Bottom of the crust
            tmoho=modelPhysics['tmoho']; # T (K) of the crust at y=ymoho
            # Mantle lithosphere
            if(MY[mm1,0]>ymoho and MY[mm1,0]<yast):
                MTK[mm1,0]=tmoho+(tast-tmoho)*(MY[mm1,0]-ymoho)/(yast-ymoho);
        
            # Crust
            if(MY[mm1,0]>7000 and MY[mm1,0]<ymoho):
                MTK[mm1,0]=ttop+(tmoho-ttop)*(MY[mm1,0]-7000)/(ymoho-7000);  #TODO apply sea depth here



            #Intrusion molten rocks 
            if MI[mm1,0] in [10,11,12]:
                if(rand()<25): #with less then 0.25 probability it is selected component, 0.75 it is hydrated mantle 
                    if labels[MI[mm1,0]-1] == 'Basalt(melt)': #TODO make intrusion postfix blocks for melts TODO GRANIteS!
                        MI[mm1,0] = 3
                    elif labels[MI[mm1,0]-1] == 'Diorite(melt)': #TODO make intrusion postfix blocks for melts TODO GRANIteS!
                        MI[mm1,0] = 4
                    elif labels[MI[mm1,0]-1] == 'Granodiorite(melt)': #TODO make intrusion postfix blocks for melts TODO GRANIteS!
                        MI[mm1,0] = 8
                else:
                    MI[mm1,0]=7  #hydrated mantle 75% of any melt
                MTK[mm1,0] = modelPhysics['tintrus']
            # else:
            # #if it is not and intrusion (geothermal gradient)
            #     MTK[mm1,0]=ttop+(tmoho-ttop)*(MY[mm1,0]-start_depth)/(ymoho-start_depth)

            #after assignment of type MI 
            mm1+=1;
            #print(mm1,'out of',mxnum*mynum)

    # Save Number of markers
    marknum=mm1
     
    
custom_print(get_now()+'Material data were implemented succesfully.');

# Initial elevation for topography profile
for i in range(0,tnum):
    # Above continental crust
    gridt[1,i]=7000; #(2,i)

# Save initial topography
topotime[0,0]=0;
topohigh[0,:]=gridt[1,:];
topowater[0,0]=waterlev;


# Density, viscosity, shear modulus, temperature, thermal conductivity, RHO*Cp arrays
etas1 = np.zeros([ynum,xnum],dtype=DTYPE);       # Viscosity for shear stress
etan1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);   # Viscosity for normal stress
mus1 = np.zeros([ynum,xnum],dtype=DTYPE);        # Shear modulus for shear stress
mun1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);    # Shear modulus for normal stress
sxy1 = np.zeros([ynum,xnum],dtype=DTYPE);        # Shear stress
sxx1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);    # Normal stress
rho1 = np.zeros([ynum,xnum],dtype=DTYPE);        # Density
tk1 = np.zeros([ynum,xnum],dtype=DTYPE);         # Old temperature
tk2=tk1.copy();                        # New temperature
rhocp1 = np.zeros([ynum,xnum],dtype=DTYPE);      # RHO*Cp (for temperature equation)
kt1 = np.zeros([ynum,xnum],dtype=DTYPE);         # Thermal conductivity
hr1 = np.zeros([ynum,xnum],dtype=DTYPE);         # Radiogenic heat production
ha1 = np.zeros([ynum,xnum],dtype=DTYPE);         # Adiabatic heat production/consuming

etas0 = np.zeros([ynum,xnum],dtype=DTYPE);       # Viscosity for shear stress
etan0 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);   # Viscosity for normal stress
mus0 = np.zeros([ynum,xnum],dtype=DTYPE);        # Shear modulus for shear stress
mun0 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);    # Shear modulus for normal stress
sxy0 = np.zeros([ynum,xnum],dtype=DTYPE);        # Shear stress
sxx0 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);    # Normal stress
rho0 = np.zeros([ynum,xnum],dtype=DTYPE);        # Density
tk0 = np.zeros([ynum,xnum],dtype=DTYPE);         # Old temperature
rhocp0 = np.zeros([ynum,xnum],dtype=DTYPE);      # RHO*Cp (for temperature equation)
kt0 = np.zeros([ynum,xnum],dtype=DTYPE);         # Thermal conductivity
hr0 = np.zeros([ynum,xnum],dtype=DTYPE);         # Radiogenic heat production
ha0 = np.zeros([ynum,xnum],dtype=DTYPE);         # Adiabatic heat production/consuming

custom_print(get_now()+'Creating matrix for Stokes and continuity equations.');
# Creating matrix for Stokes and continuity equations
L=lil_matrix(((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3),dtype=DTYPE); #scipy.sparse._lil.lil_matrix 
#L=np.zeros(((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3),dtype=np.float64);
R=np.zeros([(xnum-1)*(ynum-1)*3],dtype=DTYPE);
# Creating matrix for temperature equation
LT=lil_matrix((xnum*ynum,xnum*ynum),dtype=DTYPE);
#LT=np.zeros((xnum*ynum,xnum*ynum),dtype=np.float64);
RT=np.zeros([xnum*ynum,1],dtype=DTYPE);

# Initial time, s
# if not 'timesum' in locals() or not 'ntimestep' in locals():
#   custom_print(get_now()+'Resetting time variables...');  
#   timesum=0;
#   timestepd=copy.copy(timemax);
#   ntimestep=1;

#timesum=0;
#timestepd=copy.copy(timemax);
#ntimestep=1;
time_of_cycle=[]
# Main Time cycle
if is_reloaded == 1:
    ntimestep_start = ntimestep + 1 #start time iteration from loaded timestep
else:
    ntimestep_start = ntimestep
for ntimestep in range(ntimestep_start,stepmax+1):
    start_time_step=time.time()
    # Defining viscoelastic timestep
    timestep=copy.copy(timemax); # viscoelastic timestep
    # Plastic yeilding mark
    plastyn=0;

    # Backup transport properties arrays
    
    etas0 = np.asarray(etas1,dtype=DTYPE);
    etan0 =  np.asarray(etan1,dtype=DTYPE);
    mus0 =  np.asarray(mus1,dtype=DTYPE);
    mun0 =  np.asarray(mun1,dtype=DTYPE);
    sxy0 =  np.asarray(sxy1,dtype=DTYPE);
    sxx0 = np.asarray(sxx1,dtype=DTYPE);
    rho0 = np.asarray(rho1,dtype=DTYPE);
    tk0= np.asarray(tk2,dtype=DTYPE);
    rhocp0= np.asarray(rhocp1,dtype=DTYPE);
    kt0= np.asarray(kt1,dtype=DTYPE);
    hr0= np.asarray(hr1,dtype=DTYPE);
    ha0= np.asarray(ha1,dtype=DTYPE);
    
    # Clear transport properties arrays
    etas1 = np.zeros([ynum,xnum],dtype=DTYPE);
    etan1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);
    mus1 = np.zeros([ynum,xnum],dtype=DTYPE);
    mun1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);
    sxy1 = np.zeros([ynum,xnum],dtype=DTYPE);
    sxx1 = np.zeros([ynum-1,xnum-1],dtype=DTYPE);
    rho1 = np.zeros([ynum,xnum],dtype=DTYPE);
    tk1 = np.zeros([ynum,xnum],dtype=DTYPE);
    rhocp1 = np.zeros([ynum,xnum],dtype=DTYPE);
    kt1 = np.zeros([ynum,xnum],dtype=DTYPE);
    hr1 = np.zeros([ynum,xnum],dtype=DTYPE);
    ha1 = np.zeros([ynum,xnum],dtype=DTYPE);
    # Clear wights for basic nodes
    wtnodes=np.zeros([ynum,xnum],dtype=DTYPE);
    # Clear wights for etas
    wtetas=np.zeros([ynum,xnum],dtype=DTYPE);
    # Clear wights for etan
    wtetan=np.zeros([ynum-1,xnum-1],dtype=DTYPE);

    # Computing grid steps for basic nodes
    #TODO refactor index iteration
    
    xstp1=np.zeros([xnum-1,1],dtype=DTYPE);
    ystp1=np.zeros([ynum-1,1],dtype=DTYPE);
    
    
    for i in range(0,xnum-1):
        xstp1[i]=gridx[i+1]-gridx[i];
    
    for i in range(0,ynum-1):
        ystp1[i]=gridy[i+1]-gridy[i];
    
    #xstp1=gridx[1:]-gridx[:-1];
    #ystp1=gridy[1:]-gridy[:-1];
    
    custom_print(get_now()+'Computing grids and grid steps for Vx, Vy nodes.');    
    # Computing grids and grid steps for Vx, Vy nodes
    # Horizontal (for Vy)
    gridcx=np.zeros([xnum+1,1],dtype=DTYPE);
    xstpc1=np.zeros([xnum,1],dtype=DTYPE);
    # Vertical (for Vx)
    gridcy=np.zeros([ynum+1,1],dtype=DTYPE);
    ystpc1=np.zeros([ynum,1],dtype=DTYPE);
    # First and last nodes and steps (for external nodes)
    # Horizontal (for Vy)
    gridcx[0]=gridx[0]-xstp1[0]/2;
    xstpc1[0]=xstp1[0];
    gridcx[xnum]=gridx[xnum-1]+xstp1[xnum-2]/2;
    xstpc1[xnum-1]=xstp1[xnum-2];
    # Vertical (for Vx)
    gridcy[0]=gridy[0]-ystp1[0]/2;
    ystpc1[0]=ystp1[0];
    gridcy[ynum]=gridy[ynum-1]+ystp1[ynum-2]/2;
    ystpc1[ynum-1]=ystp1[ynum-2];
    
    
    # Internal nodes
    for i in range(1,xnum):
        gridcx[i]=(gridx[i]+gridx[i-1])/2;
    
    for i in range(1,ynum):
        gridcy[i]=(gridy[i]+gridy[i-1])/2;

    # Internal grid steps
    for i in range(1,xnum-1):
        xstpc1[i]=(gridx[i+1]-gridx[i-1])/2;
    
    for i in range(1,ynum-1):
        ystpc1[i]=(gridy[i+1]-gridy[i-1])/2;
    
    
    
    custom_print(get_now()+'Interpolating parameters from markers to nodes.'); 
    #Interpolating parameters from markers to nodes
    rho1,tk1,kt1,rhocp1,hr1, ha1, wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum=Interpolate_markers_nodes(marknum,gridx,gridy,xnum,ynum,MX,MY,
                                  MI, MRAT, MGII,MXN,MYN,MFLOW,MPL,MXM,
                                  tstp,tnum,gridt,waterlev,stressmin,ttop,etamelt,
                                  rho1,tk1,kt1,rhocp1,hr1,ha1,wtnodes,
                                  etas1,mus1,sxy1,wtetas,
                                  etan1,mun1,sxx1,wtetan,timesum)
    
    
    #end of outer funciton
    custom_print(get_now()+'Pefrorm rendering checkups...'); 

    #check phisical matrices for broken (nans ans zeros values)
    
    arr_list = {'rho1':rho1,'tk1':tk1,'kt1':kt1,'rhocp1':rhocp1,'hr1':hr1, 'ha1':ha1, 'wtnodes':wtnodes,
                'etas1':etas1,'mus1':mus1,'sxy1':sxy1,'wtetas':wtetas,'etan1':etan1,'mun1':mun1,
                'sxx1':sxx1,'wtetan':wtetan}
    for key in arr_list:
        #mask = ((arr_list[key] == 0) | (np.isnan(arr_list[key])))
        mask = np.isnan(arr_list[key])
        if mask.any():
            custom_print(get_now()+f'Interpolated matrix {key} has incostintency.'); 
            custom_print(get_now()+'Replace inf, 0 and nan values with mean rowwise'); 
            #means = np.ma.array(arr_list[key], mask = mask).mean(0) #means columnwise
            #np.asarray(arr_list[key])[mask] = means #assign mean values to broken pixels
            #arr_list[key]=np.asarray(arr_list[key])
            #arr_list[key] = np.where(~mask,arr_list[key],means)
            arr_list[key]=np.asarray(arr_list[key])
            masked = np.ma.masked_where(np.isinf(arr_list[key]) | (arr_list[key] == 0) | np.isnan(arr_list[key]), arr_list[key])
            means = np.ma.array(arr_list[key], mask = masked.mask).mean(1)
            arr_list[key]=np.where(masked.mask,masked.mask*means.data[:,None],arr_list[key])
            exec(f'{key}=arr_list["{key}"]') #replace script level variable with fixed
    del arr_list        
    

    #time.sleep(100);
    
    custom_print(get_now()+'computing  Viscosity, density, rock type for nodal points:');        
    # Computing  Viscosity, density, rock type for nodal points
    for i in range(0,ynum):
        for j in range(0,xnum):
            # Density
            if (wtnodes[i,j]!=0):
                # Compute new value interpolated from markers
                rho1[i,j]=rho1[i,j]/wtnodes[i,j];
                tk1[i,j]=tk1[i,j]/wtnodes[i,j];
                kt1[i,j]=kt1[i,j]/wtnodes[i,j];
                rhocp1[i,j]=rhocp1[i,j]/wtnodes[i,j];
                hr1[i,j]=hr1[i,j]/wtnodes[i,j];
                ha1[i,j]=ha1[i,j]/wtnodes[i,j];
            else:
                # If no new value is interpolated from markers old value is used
                rho1[i,j]=rho0[i,j];
                tk1[i,j]=tk0[i,j];
                kt1[i,j]=kt0[i,j];
                rhocp1[i,j]=rhocp0[i,j];
                hr1[i,j]=hr0[i,j];
                ha1[i,j]=ha0[i,j];
            # Shear viscosity
            if (wtetas[i,j]!=0):
                # Compute new value interpolated from markers
                etas1[i,j]=etas1[i,j]/wtetas[i,j];
                mus1[i,j]=1/(mus1[i,j]/wtetas[i,j]);
                sxy1[i,j]=sxy1[i,j]/wtetas[i,j];
            else:
                # If no new value is interpolated from markers old value is used
                etas1[i,j]=etas0[i,j];
                mus1[i,j]=mus0[i,j];
                sxy1[i,j]=sxy0[i,j];
            
            # Flatten density distribution for "air/water" boundary
            # in order to avoid perturbations in the weak layer
            if (rho1[i,j]<=1000 and gridy[i]<waterlev and hr1[i,j]==0):
                rho1[i,j]=1;
            
            if (rho1[i,j]<1000 and gridy[i]>=waterlev and hr1[i,j]==0):
                rho1[i,j]=1000;

            # Normal viscosity  
            if (i<ynum-1 and j<xnum-1):   
                if (wtetan[i,j]!=0):
                    # Compute new value interpolated from markers
                    etan1[i,j]=etan1[i,j]/wtetan[i,j];
                    mun1[i,j]=1/(mun1[i,j]/wtetan[i,j]);
                    sxx1[i,j]=sxx1[i,j]/wtetan[i,j];
                else:
                    # If no new value is interpolated from markers old value is used
                    etan1[i,j]=etan0[i,j];
                    mun1[i,j]=mun0[i,j];
                    sxx1[i,j]=sxx0[i,j];
    #TODO - Vectorization
    # Computing  Viscosity, density, rock type for nodal points

    #End of vectorization for Viscosity, density, rock type for nodal points                
    
    
    custom_print(get_now()+'Applying thermal boundary conditions for interpolated temperature')    
    # Applying thermal boundary conditions for interpolated temperature
    # Upper, Lower boundaries
    for j in range(1,xnum-1):
        # Upper boundary
        tk1[0,j]=btopt[j,0]+btopt[j,1]*tk1[1,j];
        # Lower boundary
        tk1[ynum-1,j]=bbottomt[j,0]+bbottomt[j,1]*tk1[ynum-2,j];
        
    # Left, Right boundaries: constant temperature
    for i in range(0,ynum):
        # Left boundary
        tk1[i,0]=bleftt[i,0]+bleftt[i,1]*tk1[i,1];
        # Right boundary
        tk1[i,xnum-1]=brightt[i,0]+brightt[i,1]*tk1[i,xnum-2];
    
    
    if np.fix(ntimestep/savepicstep)*savepicstep==ntimestep or ntimestep==1:
        custom_print(get_now()+' Visualise physical properties plots')
        fig0=plt.figure();
        
        x1=plt.subplot(131);
        ax = plt.gca()
        #print('gridx',gridx)
        #print('gridy',gridy)
        im=plt.pcolor(gridx[:]/1000,gridy[:]/1000,np.asarray(tk1)-273,cmap='jet',shading='auto');
        #plt.colorbar();           # showing a colorbar for the map
        plt.title('T C '+modelPhysics['modelname']);
        plt.gca().invert_yaxis()
        x1.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        
        
        x2=plt.subplot(132);
        ax = plt.gca()
        im=plt.pcolor(gridx[:]/1000,gridy[:]/1000,np.log10(etas1),cmap='jet',shading='auto');
        plt.title('log viscosity');
        plt.gca().invert_yaxis()
        x2.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        
        
        x3=plt.subplot(133);
        ax = plt.gca()
        im=plt.pcolor(gridx[:]/1000,gridy[:]/1000,rho1,shading='auto');
        #plt.colorbar();           # showing a colorbar for the map
        plt.title('density');
        plt.gca().invert_yaxis()
        x3.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.savefig(os.path.join(start_params['outdirname'],f'pic0_{ntimestep}.png'),dpi=150);
        
        plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.8, 
                        hspace=0.4)
        if vis_on_off=='off':
            plt.close();
        else:
            plt.show();

        #save csv files if needed
        if int(start_params['csv'])==1:
            custom_print(get_now()+'CSV output tk1 etas1 rho1')
            savetk1fn = os.path.join(start_params['outdirname'],f'tk1_{ntimestep}.csv')
            saveetas1fn = os.path.join(start_params['outdirname'],f'etas1_{ntimestep}.csv')
            saverho1fn = os.path.join(start_params['outdirname'],f'rho1_{ntimestep}.csv')
            ui.mat2csv(tk1,savetk1fn,heading='temperature, K')
            ui.mat2csv(etas1,saveetas1fn,heading='viscosity, Pa')
            ui.mat2csv(rho1,saverho1fn,heading='density, kg/m**3')
    
    # Interpolating initial nodal temperatures back to markers
    # to avoid initial discrepancies between markers and nodes 
    if (timesum==0):
        # Marker cycle
        for mm1 in range(1,marknum):

            # Check markers inside the grid
            if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

                #  yn    T[yn,xn]--------------------T[yn,xn+1]
                #           ?           ^                  ?
                #           ?           ?                  ?
                #           ?          dy                  ?
                #           ?           ?                  ?
                #           ?           v                  ?
                #           ?<----dx--->o Mrho[mm1,0]       ?
                #           ?                              ?
                #           ?                              ?
                #  yn+1  T[yn+1,xn]-------------------V[yn+1,xn+1]
                #
                #
                # Interpolating temperature changes from basic nodes
                #
                # Define indexes for upper left node in the cell where the marker is
                xn=int(MXN[mm1,0]);
                yn=int(MYN[mm1,0]);

                # Define normalized distances from marker to the upper left node;
                dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
                dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

                # Interpolate nodal temperature for the marker
                tkm=0;
                tkm=tkm+(1.0-dx)*(1.0-dy)*tk1[yn,xn];
                tkm=tkm+(1.0-dx)*dy*tk1[yn+1,xn];
                tkm=tkm+dx*(1.0-dy)*tk1[yn,xn+1];
                tkm=tkm+dx*dy*tk1[yn+1,xn+1];
                # Reset marker temperature
                MTK[mm1,0]=tkm;
    
    # Computing viscoelastic (numerical) viscosity and stress
    # Shear stress
    for i in range(0,ynum):
        for j in range(0,xnum):
            #Viscoelasticity factor
            try:
                xelvis=etas1[i,j]/(etas1[i,j]+timestep*mus1[i,j]);
                #print('xelvis=',xelvis)
            except ZeroDivisionError:
                xelvis = 0.0001
            # Viscoelastic viscosity = (1-xelvis)*ETA
            etas0[i,j]=etas1[i,j]*(1-xelvis);
            # Vsicoelastic stress = xelvis*Sxy
            sxy0[i,j]=sxy1[i,j]*xelvis;
        
    # Normal stress
    for i in range(0,ynum-2):
        for j in range(0,xnum-2):
            #Viscoelasticity factor
            try:
                xelvis=etan1[i,j]/(etan1[i,j]+timestep*mun1[i,j]);
                #print('xelvis=',xelvis)
            except ZeroDivisionError:
                xelvis = 0.0001
            # Viscoelastic viscosity = (1-xelvis)*ETA
            etan0[i,j]=etan1[i,j]*(1-xelvis);
            # Vsicoelastic stress = xelvis*Sxx
            sxx0[i,j]=sxx1[i,j]*xelvis;
    custom_print(get_now()+' Computing right part of mechanical viscoelastic equations')    
    # Computing right part of mechanical viscoelastic equations
    # x-Stokes
    RX1=np.zeros([ynum+1,xnum],dtype=DTYPE);
    # y-Stokes
    RY1=np.zeros([ynum,xnum+1],dtype=DTYPE);
    # continuity
    RC1=np.zeros([ynum-1,xnum-1],dtype=DTYPE);
    # Grid points cycle
    for i in range(1,ynum-1):
        for j in range(1,xnum-1):
            # Right part of x-Stokes Equation
            if(j<xnum):
                RX1[i,j]=-gx*(rho1[i,j]+rho1[i-1,j])/2;
                # Adding xelvis*dSxx0/dx
                RX1[i,j]=RX1[i,j]-(sxx0[i-1,j]-sxx0[i-1,j-1])/xstpc1[j];
                # Adding xelvis*dSxy0/dy
                RX1[i,j]=RX1[i,j]-(sxy0[i,j]-sxy0[i-1,j])/ystp1[i-1];
            
            # Right part of y-Stokes Equation
            if(i<ynum):
                RY1[i,j]=-gy*(rho1[i,j]+rho1[i,j-1])/2;
                # Adding xelvis*dSyy0/dy using that Syy0=-Sxx0 (deviatoric stress)
                RY1[i,j]=RY1[i,j]+(sxx0[i,j-1]-sxx0[i-1,j-1])/ystpc1[i];
                # Adding xelvis*dSyx0/dx using that Syx0=Sxy0 
                RY1[i,j]=RY1[i,j]-(sxy0[i,j]-sxy0[i,j-1])/xstp1[j-1];
            
    # Computing velocity field
    if (movemod==0):
        # Solving of Stokes and Continuity equations on nodes
        # and computing residuals
        # by calling function Stokes_Continuity_solver_grid() 
        # with viscoelastic numerical viscosity
        # and modified right parts
        #print('bintern shape',np.shape(bintern))  #
        
        L,R,vx1,resx1,vy1,resy1,pr1,resc1=Stokes_Continuity_solver_sandbox(L,R,prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern);
    
    custom_print(get_now()+'Solid body rotation')
    # Solid body rotation
    if (movemod==1):
        
        # Velocity at the right wall for solid body rotation
        vyright=1e-9;
        
        for i in range(0,ynum):
            for j in range(0,xnum):
                # Vx
                if(j<xnum):
                    # Relative distance of vx node from the model center
                    dx=((j-1)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1.5)*ystp-ysize/2)/(xsize/2);
                    dr=(dx**2+dy**2)**0.5;
                    # Set vx
                    vx1[i,j]=-vyright*dy;
                
                # Vy
                if(i<ynum):
                    # Relative distance of vy node from the model center
                    dx=((j-1.5)*xstp-xsize/2)/(xsize/2);
                    dy=((i-1)*ystp-ysize/2)/(xsize/2);
                    dr=(dx**2+dy**2)**0.5;
                    # Set vy
                    vy1[i,j]=vyright*dx;
                

    # Computing EPS'xx=-EPS'yy, EPSxy=EPSyx deviatoric strain rate tensor components from vx, vy
    # Computing spin tensor Espin
    exy = np.zeros([ynum,xnum]);
    exx = np.zeros([ynum-1,xnum-1]);
    esp = np.zeros([ynum,xnum]);
    eii = np.zeros([ynum-1,xnum-1]);
    
    # Grid points cycle
    for i in range(0,ynum):
        for j in range(0,xnum):
            # EPS'xx=-EPS'yy=1/2(dvx/dx-dvy/dy)
            if(i<ynum-1 and j<xnum-1):
                exx[i,j]=0.5*((vx1[i+1,j+1]-vx1[i+1,j])/xstp1[j]-(vy1[i+1,j+1]-vy1[i,j+1])/ystp1[i]);
            # EPSxy=EPSyx=1/2(dvx/dy+dvy/dx)
            exy[i,j]=0.5*((vx1[i+1,j]-vx1[i,j])/ystpc1[i]+(vy1[i,j+1]-vy1[i,j])/xstpc1[j]);
            # Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
            # (when x axis is directed rightward and y axis is directed downward) 
            esp[i,j]=0.5*((vy1[i,j+1]-vy1[i,j])/xstpc1[j]-(vx1[i+1,j]-vx1[i,j])/ystpc1[i]);
            # EPSii=(EPS'xx^2+EPSxy^2)^0.5
            if(i>0 and j>0):
                eii[i-1,j-1]=(exx[i-1,j-1]**2+(exy[i-1,j-1]**2+exy[i,j-1]**2+exy[i-1,j]**2+exy[i,j]**2)/4)**0.5;
            
    # Check maximal velocity
    vxmax=max(abs(max(vx1.flatten())),abs(min(vx1.flatten())));
    vymax=max(abs(max(vy1.flatten())),abs(min(vy1.flatten())));
    # Check marker displacement step
    if (vxmax>0):
        if (timestep>markmax*xstp/vxmax):
            timestep=markmax*xstp/vxmax;
        
    if (vymax>0):
        if (timestep>markmax*ystp/vymax):
            timestep=markmax*ystp/vymax;
        
    # Defining displacement timestep
    timestep=timestep # final displacement step
    timestepd=timestep;
    
    # Computing new stresses and stress change using the displacement timestep
    sxy2 = np.zeros([ynum,xnum]);
    sxx2 = np.zeros([ynum-1,xnum-1]);
    # Shear stress
    for i in range(0,ynum):
        for j in range(0,xnum):
            #Viscoelasticity factor
            try:
                xelvis=etas1[i,j]/(etas1[i,j]+timestep*mus1[i,j]);
            except ZeroDivisionError:
                xelvis = 0.0001
            # New viscoelastic stress = (1-xelvis)*2*ETA*EPSxy + xelvis*Sxy0
            sxy2[i,j]=(1-xelvis)*2*etas1[i,j]*exy[i,j]+xelvis*sxy1[i,j];
        
    # Normal stress
    for i in range(0,ynum-1):
        for j in range(0,xnum-1):
            #Viscoelasticity factor
            try:
                xelvis=etan1[i,j]/(etan1[i,j]+timestep*mun1[i,j]);
            except ZeroDivisionError:
                xelvis = 0.0001
            # New viscoelastic stress = (1-xelvis)*2*ETA*EPSxx + xelvis*Sxx0
            sxx2[i,j]=(1-xelvis)*2*etan1[i,j]*exx[i,j]+xelvis*sxx1[i,j];
        
    # Stress change
    dsxy = sxy2-sxy1;
    dsxx = sxx2-sxx1;
    
    custom_print(get_now()+'Computing strain rate and pressure for markers')
    # Computing strain rate and pressure for markers
    
    MX,MY,MXN,MYN,gridcx,gridcy,xstpc1,ystp1,timestep,exy,exx,esp,eii,\
    MEXX,MPR,MEXY,META,MMU,MRAT = Compute_strain_rate_pressure_markers( marknum, gridx,gridy,xnum,\
                                    ynum, MX, MY, MXN, MYN, gridcx, gridcy,xstpc1,ystp1,timestep,\
                                   exy,exx,  esp, eii,MEXX,MPR, MEXY, META, MMU,  MRAT)
    
    
                
            
    custom_print(get_now()+'Computing subgrid stress changes for markers')        
    # Computing subgrid stress changes for markers
    if (dsubgrids>0):
        # Clear subgrid stress changes for nodes
        dsxyn=np.zeros([ynum,xnum],dtype=DTYPE);
        dsxxn=np.zeros([ynum-1,xnum-1],dtype=DTYPE);
        # Clear wights for Sxy
        wtetas=np.zeros([ynum,xnum],dtype=DTYPE);
        # Clear wights for Sxx
        wtetan=np.zeros([ynum,xnum],dtype=DTYPE);
        
        
        ###start subgrid stress function body
        wtetan,wtetas,dsxxn,dsxyn,MSXX,MSXY = Compute_subgrid_stress(marknum,gridx,gridy,xnum,ynum, MX, MY,
                                       MI, MMU,META, dsubgrids,timestep,MSXY, 
                                       MSXX, wtetas, wtetan,dsxyn)
        
        
        
        
        # Subtracting subgrid stress change part from nodal stress changes
        dsxy=dsxy-dsxyn;
        dsxx=dsxx-dsxxn;
    
    custom_print(get_now()+'Updating stress for markers')
    ###start of funciton body
    # Updating stress for markers
    for mm1 in range(1,marknum):

        # Check markers inside the grid
        if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

            #  yn    sxy[yn,xn]--------------------sxy[yn,xn+1]
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o MSXY[mm1,0]       ?
            #           ?                              ?
            #           ?                              ?
            #  yn+1  sxy[yn+1,xn]-------------------sxy[yn+1,xn+1]
            #
            #
            # Interpolating old shear stress changes from Sxy nodes
            #
            # Define indexes for upper left node in the cell where the marker is
            xn=int(MXN[mm1,0]);
            yn=int(MYN[mm1,0]);

            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
            dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

            # Interpolate old Sxy stress change for the marker
            dsxym=0;
            dsxym=dsxym+(1.0-dx)*(1.0-dy)*dsxy[yn,xn];
            dsxym=dsxym+(1.0-dx)*dy*dsxy[yn+1,xn];
            dsxym=dsxym+dx*(1.0-dy)*dsxy[yn,xn+1];
            dsxym=dsxym+dx*dy*dsxy[yn+1,xn+1];

            # Update stress for the marker
            MSXY[mm1,0]=MSXY[mm1,0]+dsxym;

            #  yn    sxx[yn,xn]--------------------sxx[yn,xn+1]
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o MSXX[mm1,0]       ?
            #           ?                              ?
            #           ?                              ?
            #  yn+1  sxx[yn+1,xn]-------------------sxx[yn+1,xn+1]
            #
            #
            # Interpolating old normal stress changes from Sxx nodes
            #
            # Define, check indexes for upper left node in the Sxx cell where the marker is
            if (MX[mm1,0]<gridcx[xn+1]):
                xn=xn-1;
            
            if(xn<0):
                xn=0;
            
            if(xn>xnum-3):
                xn=xnum-3;
            
            if (MY[mm1,0]<gridcy[yn+1]):
                yn=yn-1;
            
            if(yn<0):
                yn=0;
            
            if(yn>ynum-3):
                yn=ynum-3;
            

            # Define normalized distances from marker to the upper left node;
            dx=(MX[mm1,0]-gridcx[xn+1])/xstpc1[xn+1];
            dy=(MY[mm1,0]-gridcy[yn+1])/ystpc1[yn+1];

            # Interpolate old Sxx stress for the marker
            dsxxm=0;
            dsxxm=dsxxm+(1.0-dx)*(1.0-dy)*dsxx[yn,xn];
            dsxxm=dsxxm+(1.0-dx)*dy*dsxx[yn+1,xn];
            dsxxm=dsxxm+dx*(1.0-dy)*dsxx[yn,xn+1];
            dsxxm=dsxxm+dx*dy*dsxx[yn+1,xn+1];

            # Correcting old stress for the marker
            MSXX[mm1,0]=MSXX[mm1,0]+dsxxm;
    
    ###end of funciton body
    
    custom_print(get_now()+'Solving Temperature equation')
    # Solving Temperature equation
    if (timestep>0 and tempmax>0):
        
        # Computing right part of temperature equation
        RT1=hr1.copy();
        # Grid points cycle
        for i in range(1,ynum-1):
            for j in range(1,xnum-1):
                # Adiabatic heating on(1)/off(0)
                if(adiabyn==1):
                    # Adding alp*T*DP/dt where DP/dt ~ vx*gx*rho+vy*gy*rho
                    RT1[i,j]=RT1[i,j]+ha1[i,j]*rho1[i,j]*(gx*(vx1[i,j]+vx1[i+1,j])+gy*(vy1[i,j]+vy1[i,j+1]))/2;
                
                # Computing viscoelastic shear heating for Temperature nodes
                # Hs=2*Sxx*Sxx/2/etan+2*Sxy*Sxy/2/etas
                # Shear heating on(1)/off(0)
                if(frictyn==1):
                    # Adding 2*Sxy*Sxy/2/etas
                    RT1[i,j]=RT1[i,j]+sxy2[i,j]**2/etas1[i,j];
                    # Computing and adding 2*Sxx*Sxx/2/etan
                    RT1[i,j]=RT1[i,j]+(sxx2[i-1,j-1]**2/etan1[i-1,j-1]+\
                                       sxx2[i,j-1]**2/etan1[i,j-1]+sxx2[i-1,j]**2/etan1[i-1,j]+\
                                           sxx2[i,j]**2/etan1[i,j])/4;

        # Solving temperature equation making (if needed) several thermal
        # timesteps for one displacement timestep
        # Set current thermal timestep 
        timestept=timestep
        # Set total thermal timestep
        timesteps=0;
        
        
        
        # Set old Temperature
        tk0=tk1.copy();
        #custom_print('tk0 was backuped');
        #custom_print('tk0=');
        #custom_print(tk0);
        while (timesteps<timestep):
            # Solving Temperature equation with thermal timestep
            LT,RT,tk2,rest=Temperature_solver_grid(LT,RT,timestept,xnum,ynum,gridx,gridy,kt1,rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt);
            #custom_print('cheking temperature solver performance...');
            #custom_print('tk2=');
            #custom_print(tk2);
            #custom_print('tk0=');
            #custom_print(tk0);
            #if np.isnan(tk2).any():
            #    raise ValueError('Uvaga! tk2 contains NaN values. Do somethig!');
            # Computing temperature changes
            dtk1=tk2-tk0;
            # Checking temperature changes
            dtkmax=max(abs(dtk1.flatten()));
            # Repeating temperature solution if temperature changes are too big
            if(dtkmax>tempmax):
                # Computing reduced timestep
                timestept=timestept*tempmax/dtkmax;
                # Solving Temperature equation with reduced timestep
                LT,RT,tk2,rest=Temperature_solver_grid(LT,RT,timestept,xnum,ynum,gridx,gridy,kt1,rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt);
                # Computing temperature changes
            
            # Add total thermal timestep
            timesteps=timesteps+timestept;
            # Compute current thermal timestep
            if (timestept>timestep-timesteps):
                timestept=timestep-timesteps
            else:
                timestept=timestept
        
            # Update old temperature
            tk0=tk2.copy();
        
        #check tk2 and tk1 for consistensy 
        arr_list = {'tk2':tk2,'tk1':tk1}
        for key in arr_list:
            mask = ((arr_list[key] == 0) | (np.isnan(arr_list[key])))
            if mask.any():
                custom_print(get_now()+f'Temperature matrix {key} has incostintency.'); 
                custom_print(get_now()+'Replace inf, 0 and nan values with mean rowwise'); 
                #means = np.ma.array(arr_list[key], mask = mask).mean(0) #means columnwise
                #np.asarray(arr_list[key])[mask] = means #assign mean values to broken pixels
                arr_list[key]=np.asarray(arr_list[key])
                masked = np.ma.masked_where(np.isinf(arr_list[key]) | (arr_list[key] == 0) | np.isnan(arr_list[key]), arr_list[key])
                means = np.ma.array(arr_list[key], mask = masked.mask).mean(1)
                arr_list[key]=np.where(masked.mask,masked.mask*means.data[:,None],arr_list[key])
                #arr_list[key] = np.where(~mask,arr_list[key],means)
                exec(f'{key}=arr_list["{key}"]') #replace script level variable with fixed
        del arr_list  
        
        
        # Compute temperature changes
        dtk1=tk2-tk1;

        custom_print(get_now()+'Computing subgrid diffusion for markers')
        # Computing subgrid diffusion for markers
        if (dsubgridt>0):
            # Clear subgrid temperature changes for nodes
            dtkn=np.zeros([ynum,xnum],dtype=np.float64);
            # Clear wights for basic nodes
            wtnodes=np.zeros([ynum,xnum]);
            
            #считалось на чистом Python 16 секунд, после преобразования в процедуру за 1 сек
            
            dtkn, MTK, wtnodes = Compute_subgrid_diffusion(marknum,gridx,gridy,xnum,
                                           ynum, MX, MY,tk1,kt1, dtkn, MXN, MYN, dsubgridt,
                                           wtnodes)
            
            
            
            # Subtracting subgrid diffusion part from nodal temperature changes
            dtk1=dtk1-dtkn;
        custom_print(get_now()+'Updating temperature for markers')
        # Updating temperature for markers
        for mm1 in range(1,marknum):

            # Check markers inside the grid
            if (MX[mm1,0]>=gridx[0] and MX[mm1,0]<=gridx[xnum-1] and MY[mm1,0]>=gridy[0] and MY[mm1,0]<=gridy[ynum-1]): 

                #  yn    T[yn,xn]--------------------T[yn,xn+1]
                #           ?           ^                  ?
                #           ?           ?                  ?
                #           ?          dy                  ?
                #           ?           ?                  ?
                #           ?           v                  ?
                #           ?<----dx--->o Mrho[mm1,0]       ?
                #           ?                              ?
                #           ?                              ?
                #  yn+1  T[yn+1,xn]-------------------V[yn+1,xn+1]
                #
                #
                # Interpolating temperature changes from basic nodes
                #
                # Define indexes for upper left node in the cell where the marker is
                xn=int(MXN[mm1,0]);
                yn=int(MYN[mm1,0]);

                # Define normalized distances from marker to the upper left node;
                dx=(MX[mm1,0]-gridx[xn])/xstp1[xn];
                dy=(MY[mm1,0]-gridy[yn])/ystp1[yn];

                # Calculate Marker temperature change from four surrounding nodes
                dtkm=0;
                dtkm=dtkm+(1.0-dx)*(1.0-dy)*dtk1[yn,xn];
                dtkm=dtkm+(1.0-dx)*dy*dtk1[yn+1,xn];
                dtkm=dtkm+dx*(1.0-dy)*dtk1[yn,xn+1];
                dtkm=dtkm+dx*dy*dtk1[yn+1,xn+1];
                #
                #Computing new temperature for the marker
                MTK[mm1,0]=MTK[mm1,0]+dtkm;
                
    custom_print(get_now()+'Moving Markers by velocity field')            
    # Moving Markers by velocity field
    if(markmove>0):
        # Create arrays for velocity and spin of markers
        #vxm=np.zeros([5,1],dtype=DTYPE); #(5,1) - spare index for row count 
        #vym=np.zeros([5,1],dtype=DTYPE);
        #espm=np.zeros([5,1],dtype=DTYPE);
        
        
        MX,MY,MGII,MBII,MEXX,MEXY,MSXX,MSXY,espm,vxm,vym=Move_markers_velocity_field(marknum,gridx,gridy,xnum,ynum,MX,MY,vx1,vy1,esp,
                                       MXN, MYN, timestep,MSXY, MSXX,MGII, MBII, MEXX,MEXY,markmove)
        
    
    custom_print(get_now()+' Recomputing topography surface') 
    # Recomputing topography surface
    # Set velocity initially to zero
    gridt[3:7,:]=0;
    # Defining advection velocity at topography points
    for i in range(0,tnum):
        # Check topography nodes inside the grid
        if (gridt[0,i]>=gridx[0] and gridt[0,i]<=gridx[xnum-1] and gridt[1,i]>=gridy[0] and gridt[1,i]<=gridy[ynum-1]): 
            #  xn    V(xn,yn)--------------------V(xn+1,yn)
            #           ?           ^                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o Mrho(xm,ym)       ?
            #           ?                              ?
            #           ?                              ?
            #  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

            # Define indexes for upper left BASIC node in the cell where the topograhy node is
            # using bisection
            xcur=gridt[0,i];
            ycur=gridt[1,i];
            # Find horizontal index
            xnmin=0;
            xnmax=xnum-1;
            while ((xnmax-xnmin)>1):
                # !!! SUBTRACT 0.5 since int16(0.5)=1
                xn=int((xnmax+xnmin)/2);
                if(gridx[xn]>xcur):
                    xnmax=xn;
                else:
                    xnmin=xn;
                
            # Check horizontal index
            if (xnmin<0):
                xnmin=0;
            
            if (xnmin>xnum-2):
                xnmin=xnum-2;
            
            # Find vertical index
            ynmin=0;
            ynmax=ynum-1;
            while ((ynmax-ynmin)>1):
                # !!! SUBTRACT 0.5 since int16(0.5)=1
                yn=int((ynmax+ynmin)/2);
                if(gridy[yn]>ycur):
                    ynmax=yn;
                else:
                    ynmin=yn;

            # Check vertical index
            if (ynmin<0):
                ynmin=0;
                
            if (ynmin>ynum-2):
                ynmin=ynum-2;
            
            # Define indexes for upper left node in the Vx-cell where topography node is
            # Horizontal Vx index
            xn=int(xnmin);
            # Vertical Vx index
            yn=int(ynmin);
            if(ycur>gridcy[yn+1]):
                yn=yn+1;
            
            if (yn>ynum):
                yn=ynum;
            
            # Define and check normalized distances from topography node to the upper left VX-node;
            dx=(xcur-gridx[xn])/xstp1[xn];
            dy=(ycur-gridcy[yn])/ystpc1[yn];

            # Calculate topography point velocity from four surrounding Vx nodes
            gridt[3,i]=gridt[3,i]+(1.0-dx)*(1.0-dy)*vx1[yn,xn];
            gridt[3,i]=gridt[3,i]+(1.0-dx)*dy*vx1[yn+1,xn];
            gridt[3,i]=gridt[3,i]+dx*(1.0-dy)*vx1[yn,xn+1];
            gridt[3,i]=gridt[3,i]+dx*dy*vx1[yn+1,xn+1];

            # Define indexes for upper left node in the VY-cell where the topography node is
            # Vertical Vy index
            yn=int(ynmin);
            # Horizontal Vy index
            xn=int(xnmin);
            if(xcur>gridcx[xn+1]):
                xn=xn+1;
            
            if (xn>xnum):
                xn=xnum;

            # Define and check normalized distances from topography node to the upper left VX-node;
            dx=(xcur-gridcx[xn])/xstpc1[xn];
            dy=(ycur-gridy[yn])/ystp1[yn];

            # Calculate topography node velocity from four surrounding nodes
            gridt[3,i]=gridt[3,i]+(1.0-dx)*(1.0-dy)*vy1[yn,xn];
            gridt[3,i]=gridt[3,i]+(1.0-dx)*dy*vy1[yn+1,xn];
            gridt[3,i]=gridt[3,i]+dx*(1.0-dy)*vy1[yn,xn+1];
            gridt[3,i]=gridt[3,i]+dx*dy*vy1[yn+1,xn+1];
    custom_print(get_now()+' Advect topography vertically') 
    # Advect topography vertically
    # Diffuse Topography (downhill diffusion)
    # Build topography diffusion matrix
    LA=lil_matrix((tnum,tnum));
    RA=np.zeros([tnum,1]);
    # First point: symmetry
    LA[0,0]=1;
    LA[0,1]=-1;
    RA[0,0]=0;
    # Intermediate points: dYt/dt=d(Ks*dYt/dx)/dx
    # Yt(i-1)-----Yt(i)-----Yt(i+1)
    # FD representation Yt(i)-Ks*dt*(Yt(i-1)-2*Yt(i)+Yt(i+1))/tstp^2
    for i in range(1,tnum-1):
        # Internal points
        if (gridt[0,i]>=gridx[0] and gridt[0,i]<=gridx[xnum-1]):
            # Left part 
            LA[i,i-1]=-Ks*timestep/tstp**2;
            LA[i,i]=1+2*Ks*timestep/tstp**2;
            LA[i,i+1]=-Ks*timestep/tstp**2;
            # Right part, advect topography vertically
            RA[i,0]=gridt[1,i]+gridt[4,i]*timestep;
        # External points
        else:
            # To the left from the left boundary: symmetry
            if (gridt[0,i]<gridx[0]):
                LA[i,i]=1;
                LA[i,i+1]=-1;
                RA[i,0]=0;
            # To the right from the right boundary: symmetry
            else:
                LA[i,i]=1;
                LA[i,i-1]=-1;
                RA[i,0]=0;

    # Last point: symmetry
    LA[tnum-1,tnum-1]=1;
    LA[tnum-1,tnum-2]=-1;
    RA[tnum-1,0]=0;
    # Solve Matrix
    SA=spsolve(LA,RA);
    #SA=np.linalg.solve(LA,RA);
    
    if np.isnan(SA).any():
        #raise ValueError('SA has NaN inside. Execution should be stopped!');
        custom_print(get_now()+'SA has NaN inside. Implement nan_to_mean convertion...')
        SA = np.where(np.isnan(SA), np.ma.array(SA, 
               mask = np.isnan(SA)).mean(axis = 0), SA)
        #SA = np.nan_to_num(SA)
        custom_print(get_now()+'nan_to_mean convertion is done')
    
    # Reload solutions
    for i in range(0,tnum):
        gridt[1,i]=SA[i];
    
    # Advect topography horizontally
    # Define maximal horizontal velocity at topography nodes
    vxmax=max(abs(gridt[3,:]));
    # Defining topography advection timestep
    dt=timestep;
    nt=1;
    if(vxmax>0):
        dt=tstp/vxmax;
        if (dt<timestep):
            nt=np.float64(np.int16(timestep/dt-0.5))+1;
            dt=timestep/nt;
        else:
            dt=timestep;
        
    # Defining FCT parameter MU
    mu=1/8;
    # Advect topography with FCT
    for t in range(0,nt):
        # Step 0: Set new profile
        gridt[2,:]=gridt[2,:];
        # Step 1: Transport+numerical diffusion stage
        for i in range(1,tnum-2):
            # Defining FCT parameters EPS and NU
            eps=gridt[3,i]*dt/tstp;
            nu=1/8+(eps**2)/2;
            # Change topography
            gridt[2,i]=gridt[1,i]-eps/2*(gridt[1,i+1]-gridt[1,i-1])+nu*(gridt[1,i+1]-2*gridt[1,i]+gridt[1,i-1]);
        
        # Step 2: Antidiffusion stage
        # Antidiffusion flow for the first cell
        gridt[4,0]=0;
        for i in range(1,tnum-3):
            # Corrected antidiffusion flow for current cell
            delt0=gridt[2,i]-gridt[2,i-1];
            delt1=gridt[2,i+1]-gridt[2,i];
            delt2=gridt[2,i+2]-gridt[2,i+1];
            s=np.sign(delt1);
            gridt[5,i]=s*max(0,min(min(s*delt2,s*delt0),mu*abs(delt1)));
            gridt[1,i]=gridt[2,i]-gridt[5,i]+gridt[5,i-1];
    

    # Interpolating vx, vy for the basic grid
    vxb=np.zeros([ynum,xnum]);
    vyb=np.zeros([ynum,xnum]);
    for j in range(0,xnum):
        for i in range(0,ynum):
            vxb[i,j]=(vx1[i,j]+vx1[i+1,j])/2;
            vyb[i,j]=(vy1[i,j]+vy1[i,j+1])/2;
    
            
    
    
    if np.fix(ntimestep/savepicstep)*savepicstep==ntimestep or ntimestep==1:
        custom_print(get_now()+'Visualize stress and strain')        
        fig2=plt.figure();
        
        x1=plt.subplot(221);
        ax = plt.gca()
        #plt.pcolor(gridcx[0,xnum-1]/1000,gridcy[1,ynum-1]/1000,sxx2,cmap='jet',shading='auto');
        im=plt.pcolor(gridcx[:-1,0]/1000,gridcy[:-1,0]/1000,sxx2,cmap='jet',shading='auto');
        plt.title('Sxx'+modelPhysics['modelname']);
        plt.gca().invert_yaxis()
        x1.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
            
        x2=plt.subplot(222);
        ax = plt.gca()
        plt.pcolor(gridx[1:-1]/1000,gridy[1:-1]/1000,sxy2[1:-1,1:-1],cmap='jet',shading='auto');
        plt.title('Sxy');
        plt.gca().invert_yaxis()
        x2.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        
        x3=plt.subplot(223);
        ax = plt.gca()
        plt.pcolor(gridcx[:-1,0]/1000,gridcy[:-1,0]/1000,pr1,cmap='jet',shading='auto');
        plt.title('Pr');
        plt.gca().invert_yaxis()
        x3.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        #plt.colorbar();           # showing a colorbar for the map
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        
        x4=plt.subplot(224);
        ax = plt.gca()
        plt.pcolor(gridcx[:-1,0]/1000,gridcy[:-1,0]/1000,np.log10(eii),cmap='jet',shading='auto');
        #plt.colorbar();           # showing a colorbar for the map
        plt.title('Eii');
        plt.gca().invert_yaxis()
        x4.set_aspect('equal', 'box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.savefig(os.path.join(start_params['outdirname'],f'pic1_{ntimestep}.png'),dpi=150);
        if vis_on_off=='off':
            plt.close();
        else:
            plt.show();
        del fig2

        fig3=plt.figure();
        
        x1=plt.subplot(121);
        ax = plt.gca()
        plt.pcolor(gridx[:]/1000,gridy[:]/1000,np.log10(etas0),cmap='jet',shading='auto');
        # xx       =   np.arange(0,xsize+xstp,xstp); # Horizontal
        # yy       =   np.arange(0,ysize+ystp,ystp); # Vertical
        xx = np.linspace(0, xsize+xstp, num=xnum) # Horizontal
        yy = np.linspace(0, ysize+ystp, num=ynum) # Vertical

        
        try:
            plt.streamplot(xx/1000,yy/1000,vxb,vyb);
        except:
            print('ValueError: y must be strictly increasing')
        plt.title('Step='+str(ntimestep)+' Myr='+str(timesum*1e-6/(365.25*24*3600))+' '+modelPhysics['modelname']);
        plt.gca().invert_yaxis();
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        x1.set_aspect('equal', 'box')    
        
        x2=plt.subplot(122);
        ax = plt.gca()
        plt.pcolor(gridx[:]/1000,gridy[:]/1000,np.log10(etas1),cmap='jet',shading='auto');
        #plt.colorbar();           # showing a colorbar for the map
        plt.title('log10(eta1)');
        plt.gca().invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        x2.set_aspect('equal', 'box')  
        plt.savefig(os.path.join(start_params['outdirname'],f'pic2_{ntimestep}.png'),dpi=150);
        if vis_on_off=='off':
            plt.close();
        else:
            plt.show();
        del fig3 

        #print('start_params',start_params)

        if int(start_params['csv'])==1:
            custom_print(get_now()+'CSV output Sxx Sxy Eii')
            saveSxxfn = os.path.join(start_params['outdirname'],f'sxx_{ntimestep}.csv')
            saveSxyfn = os.path.join(start_params['outdirname'],f'sxy_{ntimestep}.csv')
            saveeiifn = os.path.join(start_params['outdirname'],f'eii_{ntimestep}.csv')
            ui.mat2csv(sxx2,saveSxxfn,heading='normal stress Sxx, Pa')
            ui.mat2csv(sxy2,saveSxyfn,heading='shear stress Sxy, Pa')
            ui.mat2csv(eii,saveeiifn,heading='strain rate invariant from power law, AD*exp(-Ea/RT)')

        # TODO Visualizing marker type
        # Pixel grid resolution
        xresol=int(np.int16(xsize/500))+1;
        yresol=int(np.int16(ysize/500))+1;
        ngrid=2;
        sxstp=xsize/(xresol-1);
        systp=ysize/(yresol-1);
        mx=np.arange(0,(xsize/1000),sxstp/1000);
        my=np.arange(0,(ysize/1000),systp/1000);
        # Process markers
        markcom=np.nan*np.ones([yresol,xresol]);
        markdis=1e+20*np.ones([yresol,xresol]);
        markgii=np.nan*np.ones([yresol,xresol]);
        markbii=np.nan*np.ones([yresol,xresol]);
        markmel=np.nan*np.ones([yresol,xresol]);
        for mm1 in range(1,marknum):
            # Define pixel cell
            m1=int(np.int16((MX[mm1,0]-gridx[0])/sxstp-0.5))+1;
            m2=int(np.int16((MY[mm1,0]-gridy[0])/systp-0.5))+1;
            if (m1<0):
                m1=0;
    
            if (m1>xresol-2):
                m1=xresol-2;
    
            if (m2<0):
                m2=0;
            
            if (m2>yresol-2):
                m2=yresol-2;
            
            # Define indexes of surrounding pixels
            m10min=m1-ngrid;
            if (m10min<0):
                m10min=0;
            
            m10max=m1+1+ngrid;
            if (m10max>xresol-1):
                m10max=xresol-1;
            
            m20min=m2-ngrid;
            if (m20min<0):
                m20min=0;
            
            m20max=m2+1+ngrid;
            if (m20max>yresol-1):
                m20max=yresol-1;
            
            # Update pixels around the marker
            for m10 in range(m10min,m10max+1):
                for m20 in range(m20min,m20max+1): 
                    # Check distance to current pixel
                    dx=(MX[mm1,0]-gridx[0])-(m10-1)*sxstp;
                    dy=(MY[mm1,0]-gridy[0])-(m20-1)*systp;
                    dd=(dx*dx+dy*dy)**0.5;
                    if(dd<markdis[m20,m10]):
                        markcom[m20,m10]=MI[mm1,0];
                        markgii[m20,m10]=max(MGII[mm1,0],1e-5);
                        markbii[m20,m10]=MBII[mm1,0];
                        markmel[m20,m10]=MXM[mm1,0];
                        markdis[m20,m10]=dd;
                    
        
        # C_map=np.array([
        #     [0.50196, 0.50196, 0.50196],
        #     [0.75294,0.75294,0.75294],
        #     [1, 0, 0],
        #     [0.486, 0.388, 0.898],
        #     [0, 0, 0.71765],       
        #     [0, 0.84314, 0],    
        #     [1, 1, 1],
        #     [0.68235, 0.34118, 0],
        #     [0.50588, 0.99608, 0.78824]  
        #     ]);

        # Draw composition
        markcom[0,0:9]=np.arange(0,9);
        res_pr1=resize(pr1,[markcom.shape[0],markcom.shape[1]]); #rescale for overlay
        
        
        fig3=plt.figure();
        
        custom_print(get_now()+'Visualize markers distribution')   
        #fig,ax = plt.subplots()
        #fig=Figure() #for garbage collection
        #canvas = FigureCanvas(fig)
        ax = fig3.add_subplot(111)
        im = ax.imshow(res_pr1,alpha=1,cmap='jet');
        fig3.colorbar(im,orientation='vertical');
        ax.imshow(markcom,alpha=0.75,cmap='gray');
        #fig.canvas.draw()
        xlabels = [int(item) for item in ax.get_xticks()]
        ylabels = [int(item) for item in ax.get_yticks()]
        ax.xaxis.set_major_locator(mticker.FixedLocator(xlabels)) #fix old places
        ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels)) #fix old places
        xlabels = [int(item)//2 for item in ax.get_xticks()]    #set new values
        ylabels = [int(item)//2 for item in ax.get_yticks()]
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        plt.title('Step='+str(ntimestep)+' Myr='+str(timesum*1e-6/(365.25*24*3600))+' '+modelPhysics['modelname']);
        plt.savefig(os.path.join(start_params['outdirname'],f'pic3_{ntimestep}.png'),dpi=150);
        if vis_on_off=='off':
            plt.close(fig3);
        else:
            plt.show(fig3);

        #clears memory
        del fig3,ax

        custom_print(get_now()+'Draw strain of markers') 
        # Draw strain
        fig4=plt.figure();
        x41=plt.subplot(131);
        ax = plt.gca()
        plt.gca().invert_yaxis()
        #plt.pcolor(mx,my,-np.log10(abs(markgii)),cmap='jet',shading='auto');
        plt.pcolor(-np.log10(abs(markgii)),cmap='jet',shading='auto');
        plt.title('markgii');
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        x41.set_aspect('equal', 'box')  
        
        
        # Draw strain
        x42=plt.subplot(132);
        ax = plt.gca()
        #plt.pcolor(mx,my,-np.log10(abs(markbii)),cmap='jet',shading='auto');
        plt.pcolor(-np.log10(abs(markbii)),cmap='jet',shading='auto');
        plt.gca().invert_yaxis()
        plt.title('markbii');
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        x42.set_aspect('equal', 'box')  
        custom_print(get_now()+'Draw melt fraction')
        # Draw melt fraction
        x43=plt.subplot(133);
        ax = plt.gca()
        plt.pcolor(markmel,cmap='jet',shading='auto');
        plt.gca().invert_yaxis()
        plt.title('melt fraction');
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        x43.set_aspect('equal', 'box') 
        plt.savefig(os.path.join(start_params['outdirname'],f'pic4_{ntimestep}.png'),dpi=150);
        if vis_on_off=='off':
            plt.close();
        else:
            plt.show();
        
        del x41,x42,x43,ax,fig4

    custom_print(get_now()+'timesum={}'.format(timesum))
    custom_print(get_now()+'timestep={}'.format(timestep))
    
    # Advance in time
    timesum=timesum+timestep

    custom_print(get_now()+'ntimestep={}'.format(ntimestep))
    #ntimestep=ntimestep
    
    # Save topography
    topotime[ntimestep-1,0]=timesum;
    topohigh[ntimestep-1,:]=gridt[1,:];
    topowater[ntimestep-1,0]=waterlev;
    
    # Save mat file
    if(np.fix(ntimestep/savestep)*savestep==ntimestep):
        custom_print(get_now()+'Saving into save{}.p file...'.format(ntimestep))
        #namemat    =  ['intrusion1_'+str(ntimestep)];
        #save(namemat);
        #pickle.dump( namemat, open( "save.p", "wb" ) )
        varslist = [timestep,timestept,timesum,ntimestep,MX,MY,MTK,MI,MXN,MYN,MSXX,MSXY,META,MEXX,MEXY,MPR,MGII,MRAT,\
        gridt,topotime,topohigh,topowater,prfirst,etas0,etan0,etas1,etan1,gridx,gridy,\
        mus0,mus1,mun0,mun1,sxy0,sxx0,sxy1,sxx1,rho0,rho1,tk0,tk2,rhocp0,rhocp1,\
        kt0,kt1,hr0,hr1,ha0,ha1]
        
        var2save = tuple([np.asarray(el) for el in varslist])
            
        os.path.join(start_params['outdirname'],f'save{ntimestep}.p')

        with open(os.path.join(start_params['outdirname'],f'save{ntimestep}.p'), 'wb') as f:
            pickle.dump((var2save), f)
        with open(os.path.join(start_params['outdirname'],'save.p'), 'wb') as f:  #last save
            pickle.dump((var2save), f)
        custom_print(get_now()+'Saved.')
       
    time_of_cycle.append(time.time()-start_time_step);
    mean_time_cycle=sum(time_of_cycle)/len(time_of_cycle);
    ETA=mean_time_cycle*(stepmax-ntimestep);
    custom_print(get_now()+f'Timesum: {timesum}');
    custom_print(get_now()+f'Timestep: {ntimestep}');
    custom_print(get_now()+f'Estimated time left: {sec_to_mdhs(ETA)}');
