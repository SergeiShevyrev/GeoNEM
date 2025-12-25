GeoNEM (Geodynamic Numeric Environmental Modeling)

That folder hosts the marker in cell (MIC) computational model initially developed in Matlab by Taras Gerya and then ported into Python/Cython by Sergei L Shevyrev (fegi.ru).

Picture legeng:
Sxx  second invariant of future marker stress for given viscoelastic timestep
Eii  compute strain rate invariant from power law
Sxy  second invariant of marker stress
markgii accumulated plastic strain
markbii accumulated bulk strain

1) Compile mylib.c using: 
gcc -Wall -pedantic -shared -fPIC -o mylib.so mylib.c

2) compile main_b.pyx and uilib.pyx using: 
python compile.py

3) Use command line args to load model file. Number of steps can be passed through the command line (stepmax):

python launcher.py -modelfile=closed_basin_fault_magma_rising.cfg --outdirname=model_out --csv=0 --stepmax=15000 -savestep=1000 --log_file=model_out/model.out

You may track progress with log file and image output
