import numpy as np

DTYPE2 = np.int64
DTYPE = np.double 

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
# 8 = Upper continental crust (granodiorite)
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