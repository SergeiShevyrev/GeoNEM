
import sys, os
import uilib as ui


file_dir = os.path.dirname(__file__) 
sys.path.append(file_dir)

#TODO launch launcher with command-line variables, create txt file startup.p in a folder with parameters to be read by the
# launcher.py That startup.p can be read from main_b.so and taken into work 

#python launcher.py -modelfile=subduction.cfg -outdirname=out_model log_file=model_log.out -csv=0 


#create file with startup parameters (if any)
if len(sys.argv)>1:
    #dump args into startup.p file
    ui.args2post('startup.p')
    
    #load startup parameters
    #startup=ui.parse_startup('startup.p')


#launch Cython compiled file
import main_b


