#cython: language_level=2
import sys
import numpy as np
import json

#parse command line arguments
def argparse():
    key_val_dict = {}
    for arg in sys.argv[1:]:
        arg=arg.replace('--','');
        if arg[0] == '-':
            arg = arg[1:]
        key,val=arg.split('=');
        try:
            key_val_dict.update({key:eval(val)});
        except:
            key_val_dict.update({key:val});
    return key_val_dict

#save command line args into startup file
def args2post(outfile):
    pargv = []
    for arg in sys.argv[1:]:
        arg=arg.replace('--','');
        if arg[0] == '-':
             pargv.append(arg[1:])
        else:
            pargv.append(arg)
    out = '?'+'&'.join(pargv)
    with open(outfile,'w') as f:
        f.write(out)

def parse_startup(infile):
    params={}
    with open(infile,'r') as f:
        txt = f.read()[1:]
    keyvals = txt.split('&')
    for keyval in keyvals:
        key,val=keyval.split('=');
        try:
            params.update({key:eval(val)});
        except:
            params.update({key:val});
    return params

def generate_polygons(dislocation_type:int,xsize:int):
    """
    how to apply returned result:
    
    #dislocation_type==1 or 2 folds or monocline
    #         else:
    #             if((MY[mm1,0]>=7000) and (MY[mm1,0]<32000)):  #only count sedimentary cover
    #                 #use that processor 
    #                 for pnts in pnts_poly:
    #                     poly = Polygon(pnts);
    #                     pnt=Point(float(MX[mm1,0]),float(MY[mm1,0]));
    #                     if (poly.contains(pnt)):
    #                         MI[mm1,0]=2;

    """
    if dislocation_type!=0: #build polygons for sedimentary dislocations
        #custom_print(get_now()+'Not horizontal layer structure was selected, building polygons...');
        if dislocation_type==1: #monocline
            dip_angle=35; 
            h=5000;      #monocline layer thickness
            hv=h/np.sin(np.deg2rad(dip_angle));
            dx=(32e3-7e3)/np.tan(np.deg2rad(dip_angle));
      
            pnts_poly=[]
            for n in range(1,int(xsize/h*2)):
                if not n%2==0:
                    hv_r=int(np.random.randint(5,10)*0.1*hv) #randomized hv thickness
                    #h_r=h
                    pnts_poly.append([(-8*h+n*h,7e3),(-8*h+hv_r+n*h,7e3),(-8*h+dx+hv_r+n*h,32e3),(-8*h+dx+n*h,32e3)])

        if dislocation_type==2: #folds
            #dislocation description
            dip_angle=35
            
            folds = 3 # how many folds cycles
            resolution = 50 # how many datapoints to generate
            start_depth=7000   #interval start
            end_depth= 32000 #interval end
            thickness=2000 #layer thickness
            height=1000 #vertical amplitude of folds
            num_layers=(end_depth-start_depth)//thickness
            
            length = np.pi * 2 * folds
            sine=np.sin(np.arange(0, length, length / resolution))*height  #wave
            
            layers=[] #two lines, lower and upper margin
            
            for n in range(num_layers):
                if n%2 !=0:
                    my_wave = start_depth+n*thickness + sine
                    my_wave2=my_wave+thickness*np.random.randint(8,10)*0.1
                    x_points=np.linspace(0, xsize,len(my_wave))
                    #layers.append([x_points+x_points[::-1],my_wave+my_wave2[::-1]])
                    layers.append((np.concatenate((x_points,np.flip(x_points))),\
                                   np.concatenate((my_wave,np.flip(my_wave2)))))
    
            #create polygons for layers polylines
            pnts_poly=[]
            for i in range(len(layers)):
                poly=[]
                for x,y in zip(layers[i][0],layers[i][1]):
                    poly.append((x,y))
                pnts_poly.append(poly)
    return pnts_poly

def csv2mat(fname:str,sep=';',eol='\n')->tuple[str,np.array]:
    with open(fname,'r') as f:
        lines = f.readlines()
        out = []
        for line in lines[1:]:
            out.append([float(i) for i in line.replace(eol, '').split(sep)])
        return (lines[0],np.array(out))
            
def mat2csv(mat:np.array,fname:str,heading='matrix of physical discribution parameter',sep=';',eol='\n')->None:
    rows,cols = np.shape(mat)
    with open(fname,'w') as f:
        f.write(heading+eol)
        #f.write('>>>'+eol) 
        for i in range(rows):
            #f.write('\n')
            mat_line = [str(j) for j in mat[i,:]]
            f.write(sep.join(mat_line)+eol)

#TODO mat2json
def mat2json(mat:np.array,fname:str='parameter.json', name:str='matrix input',description:str='description',timestep=-1,myr=-1):
    matlist = []
    #check if value is float conversible
    isint = isinstance(mat[0,0],int)
    for i in range(mat.shape[0]):
        row = []
        for j in range(mat.shape[1]):
            if isint:
                row.append(int(mat[i,j]))
            else:
                row.append(float(mat[i,j]))
        matlist.append(row)
    out = {
        'fname':fname,              #output filename
        'name':name,                #name of the matrix to be saved
        'description':description,  #description of the matrix
        'timestep':timestep,        #number of the timesep
        'myr':myr,                  #age of the model
        'resx':mat.shape[1],      #resolution in pixels per columns
        'resy':mat.shape[0],      #resoluiton in pixels per rows
        'matrix':matlist            #matrix output
    }
    #запись в json
    try:
        with open(fname,mode='w',encoding='utf-8') as f:
            f.write(json.dumps(out,ensure_ascii=False,indent=4))
            print(f'saved to {fname}')
    except FileNotFoundError:
        print('file not found')


def json2mat(filename):
        #TODO convert json file to matrix
        #читаем из json
        info_read = []
        with open(filename,encoding='utf-8') as f:
            info_read = json.loads(f.read())


if __name__ == '__main__':
    #python uilib.py -modelfile=subduction.cfg -outdirname=out_model -csv=0 -is_load=n dislocation_type=0 melt_type=0   
    #saves cl args into a file
    if len(sys.argv)>1:
        args2post('startup.p')
        #loads arguments
        sa = parse_startup('startup.p')
        print(sa)