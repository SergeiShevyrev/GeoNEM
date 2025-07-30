// 
// main.c
//compile it with:
//gcc -Wall -pedantic -shared -fPIC -o mylib.so mylib.c

//array as funciton parameter:
//https://stackoverflow.com/questions/5862915/passing-numpy-arrays-to-a-c-function-for-input-and-output
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//function below is for indexing testing only
//numpy pointer to double c array convertion (otherwise can not convert retun types back to numpy )
int get_value(double (*arr_m)[3],int (*arr_i)[3],int mm1, int mat_num){
    double res;
    printf("%s","MRHO[1][1]=");
    printf("%f",arr_m[1][1]);
    
    printf("material type:");
    printf("\n");
    printf("%i",arr_i[mm1][0]);
    printf("property value:");
    printf("\n");
    printf("%f",arr_m[arr_i[mm1][0]][mat_num]);
    printf("\n");
    printf("\n");
    res = arr_m[arr_i[mm1][0]][mat_num];
return res;
}

//2d array numpy in c
//https://stackoverflow.com/questions/22425921/pass-a-2d-numpy-array-to-c-using-ctypes
//passing c array of variable size
//https://stackoverflow.com/questions/14548753/passing-a-multidimensional-variable-length-array-to-a-function

double** get_value2(const int col_num, double arr_m[][col_num],size_t rows, size_t cols){
    
    
    printf("\n im already tired...");
    printf("%s","arr_m[1][1]=");
    printf("%f",arr_m[1][1]);
    printf("\n\n\n\n");
    size_t i,j;
        
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%f" ,arr_m[i][j]); //round 
            printf("|");
            }
            printf("\n");
            }
    //cooking something to return 
    double **res = malloc(rows * sizeof(double));

    for (i = 0; i < rows; i++) {
        res[i] = malloc(cols * sizeof(double));
        
        for (j = 0; j < cols; j++) {
            res[i][j] = arr_m[i][j]; //round 
            }
            }        
            
            return res;

}

double** ptr2double_arr(double (*arr)[3],size_t m, size_t n){
    double **res = malloc(m * sizeof(double*));
    int i,j; 
    for (i = 0; i < m; i++) {
        res[i] = malloc(n * sizeof(double));
        for (j = 0; j < n; j++) {
            res[i][j] = arr[i][j]; //round 
            }
            }
return res;
}

double** ptr2double_arr_new(const int cols, double arr[][cols],size_t rows){
    double **res = malloc(rows * sizeof(double*));
    int i,j; 
    for (i = 0; i < rows; i++) {
        res[i] = malloc(cols * sizeof(double));
        for (j = 0; j < cols; j++) {
            res[i][j] = arr[i][j]; //round 
            }
            }
return res;
}



struct result_melt{
    double xmelt;
    double hlat;
};

struct result_melt Melt_fraction(double ppa, double mtk, int rock){
    /*
    printf("ppa=");
    printf("\n");
    printf("%f",ppa);
    printf("\n");
    printf("mtk=");
    printf("\n");
    printf("%f",mtk);
    printf("\n");
    printf("rock=");
    printf("\n");
    printf("%i",rock);
    printf("\n");
    */
    
    
    struct result_melt res;
    // Calculate melt fraction using marker type
    double P=ppa*1e-6; // Pressure, MPa
    double tl=0; // Liquidus temperature
    double xmelt,hlat,ts;
    long int HL;
    
    //switch rock
    
    // 1 = Sticky air/water, no melting
    if (rock==1){
        tl=0;
    }
        // 2 = Sediments: latent heat 300 kJ/kg
    else if (rock==2){
        // Solidus Temperature
        if (P<1200){ 
            ts=889+17900/(P+54)+20200/pow((P+54),2); 
            }
        else{
            ts=831+0.06*P;
        }
        // Liquidus temperature
        tl=1262+0.09*P;
        // Latent heat
        HL=300000;
    }
        // 3 = Basalt: latent heat 380 kJ/kg
    else if (rock==3){
        // Solidus Temperature
        if (P<1600){ 
            ts=973-70400/(P+354)+77800000/pow((P+354),2); 
            }
        else{
            ts=935+0.0035*P+0.0000062*pow(P,2);
        }
        // Liquidus temperature
        tl=1423+0.105*P;
        // Latent heat
        HL=380000;
    }
    // 4 = Gabbro: latent heat 380 kJ/kg
    else if (rock==4){
        // Solidus Temperature
        if (P<1600){ 
            ts=973-70400/(P+354)+77800000/pow((P+354),2); 
            }
        else{
            ts=935+0.0035*P+0.0000062*pow(P,2);
        }
        // Liquidus temperature
        tl=1423+0.105*P;
        // Latent heat
        HL=380000;
    }
    // 5 = Lithospheric mantle (dry): latent heat 400 kJ/kg
    else if (rock==5){
        // Solidus Temperature
        if (P<10000){ 
            ts=1394+0.132899*P-0.000005104*pow(P,2); 
            }
        else{
            ts=2212+0.030819*(P-10000);
        }
        // Liquidus temperature
        tl=2073+0.114*P;
        // Latent heat
        HL=400000;
    }
    // 6 = Asthenospheric mantle (dry): latent heat 400 kJ/kg
    else if (rock==6){
        // Solidus Temperature
        if (P<10000){ 
            ts=1394+0.132899*P-0.000005104*pow(P,2); 
            }
        else{
            ts=2212+0.030819*(P-10000);
        }
        
        // Liquidus temperature
        tl=2073+0.114*P;
        // Latent heat
        HL=400000;
    }
    // 7 = Hydrated mantle (wet): latent heat 400 kJ/kg
    else if (rock==7){
        // Solidus Temperature
        if (P<2400){ 
            ts=1240+49800/(P+323); 
            }
        else{
            ts=1266-0.0118*P+0.0000035*pow(P,2);
            }
        // Liquidus temperature
        tl=2073+0.114*P;
        // Latent heat
        HL=400000;
    }
    // 8 = Upper continental crust: latent heat 300 kJ/kg
    else if (rock==8){
    // Solidus Temperature
        if (P<1200){ 
            ts=889+17900/(P+54)+20200/pow((P+54),2); 
            }
        else{
            ts=831+0.06*P;
        }
        // Liquidus temperature
        tl=1262+0.09*P;
        // Latent heat
        HL=300000;
    }
    
    // 9 = Lower continental crust: latent heat 380 kJ/kg
    else if (rock==9){
    // Solidus Temperature
        if (P<1600){
            ts=973-70400/(P+354)+77800000/pow((P+354),2); 
            }
        else{
            ts=935+0.0035*P+0.0000062*pow(P,2);
            }
        // Liquidus temperature
        tl=1423+0.105*P;
        // Latent heat
        HL=380000;
        }
    else{
        printf("error - unknown rock type\n");
        printf("%i",rock);
        }
    
    
    //exit(0);
        
    // Melt fraction and latent heat calc, check
    xmelt=0;
    hlat=0;
    if (tl>0){
    	// Solidus and liquidus must not entersect
        // in the extrapolation region
        if (ts>tl-100){ 
            ts=tl-100;
            }
        // Melt fraction
        xmelt=(mtk-ts)/(tl-ts);
        if (xmelt<0){ 
            xmelt=0;
            }
        if (xmelt>1){
            xmelt=1;
            }
        	// Latent heat calc 
        hlat=HL*xmelt;
 	}
 	//return structure
 	
 	res.xmelt=xmelt;
 	res.hlat=hlat;
 	
    return res;
}




struct result_t {
  double** rho1;
  double** tk1;
  double** kt1;
  double** rhocp1;
  double** hr1;
  double** ha1;
  double** wtnodes;
  double** etas1;
  double** mus1;
  double** sxy1;
  double** wtetas;
  double** etan1;
  double** mun1;
  double** sxx1;
  double** wtetan;
  double timesum;
  int plastyn;
 
  };

struct result_t interpolate_markers_nodes(  //TODO REDO INPUTS IN Python code
                const int xnum, //same as xnum NEW input parameter
                const int ynum, //same as xnum NEW input parameter
                //const int mark_num, //same as mark_num NEW input parameter
                const int layer_num, //same as layer_num NEW input parameter
                const int marknum,
                double gridx[xnum],
                double gridy[ynum],
                //size_t xnum,      //mark for removal
                //size_t ynum,
                double MX[marknum][1], // MX[][xnum]
                double MY[marknum][1],
                long int MI[marknum][1],//int (*MI)[3],     //right shape for vector column
                double MRAT[marknum][1], 
                double MGII[marknum][1],
                double MXN[marknum][1],
                double MYN[marknum][1],
                double MFLOW[][6],
                double MPL[][7],
                double MXM[marknum][1],
                double tstp,
                const int tnum,
                double gridt[][tnum], //
                int waterlev,
                double stressmin,
                int ttop,   
                double etamelt,
                double rho1[][xnum],
                double tk1[][xnum],
                double kt1[][xnum],
                double rhocp1[][xnum],
                double hr1[][xnum],
                double ha1[][xnum],
                
                double wtnodes[][xnum],
                double etas1[][xnum],
                double mus1[][xnum],
                double sxy1[][xnum],
                double wtetas[][xnum],
                double etan1[][xnum-1],
                double mun1[][xnum-1],
                double sxx1[][xnum-1],
                double wtetan[][xnum-1],
                double timesum, //new parameters are added from here
                double xstp1[xnum-1][1],
                double ystp1[xnum-1][1],
                double MRHO[][5],
                double MKT[][3],
                double MTK[marknum][1],
                double MPR[marknum][1],
                double MCP[layer_num],
                double MHR[layer_num],
                double MSXX[marknum][1],
                double MSXY[marknum][1],
                double MEXX[marknum][1],
                double MEXY[marknum][1],
                double META[marknum][1],
                double MMU[layer_num],
                double timestep,
                double etamin,
                double etamax,
                int ntimestep,
                double etawt,
                int plastyn
                )

                    
{   
//Interpolating parameters from markers to nodes
long int mm1;
long int xn,yn,dp;
double dx, dy,mwt;
long int xnmin,xnmax,ynmin,ynmax;
long int plawiter,grid_ind;



struct result_t res; //return value structure

struct result_melt res_m,res_m0,res_m1; //structure for local result melt calls

//!!
//TODO do memory allocation for every matrix created inside of function
//!!


double MRHOCUR,MRHOCPCUR,MKTCUR,MHACUR,METACUR,METACURP,sii0,
            plawexp,xelvis, sxxnew,sxynew,
            siinewe,sii1,eii0,eii1,siicur,siinew,
            eiicur,MCOHES,MFRICT,
            xmelt,hlat,hlat1,hlat0,eta0,sxxnewe,sxynewe,siiyeld,MMUCUR;

double RGAS = 8.314; //Gas constant J/mol/K

//printf("Mark #1");


//TODO НАДО ПЕРЕД ПЕРЕДАЧЕЙ В C КОНВЕРТИРОВАТЬ NUMPY МАССИВЫ В ФОРМАТ, КОТОРЫЙ МОЖЕТ ИНДЕКСИРВОАТЬ C
// printf("%s","MRHO[1][1]=");
// printf("%f",MRHO[1][1]);
// printf("\n");
// 
// printf("%s","MI[100][0]=");
// printf("%i",MI[100][0]);
// printf("\n");
//MI[mm1][0]

//exit(0);
/*
for (mm1 = 0; mm1 < marknum; mm1++) {
    printf("MX[mm1][0]=%f\n",MX[mm1][0]);
    printf("MY[mm1][0]=%f\n",MY[mm1][0]);
    printf("MI[mm1][0]=%d\n",MI[mm1][0]);

}

exit(0);
*/
for (mm1 = 0; mm1 < marknum; mm1++) {
        // Check markers inside the grid
        
        
        //printf("mark num=%i",mm1);
        //printf(" out of %i \n",marknum);
        
        /*
        printf("rock type=%i",MI[mm1][0]); //каждый второй маркер 0 НУЛЕВОЙ?! ПОЧЕМУ!
        printf("\n");

        continue;
        */
        
        
        
        if (MX[mm1][0]>=gridx[0] && MX[mm1][0]<=gridx[xnum-1] && 
            MY[mm1][0]>=gridy[0] && MY[mm1][0]<=gridy[ynum-1]){ 
            //Erosion-sedimentation
            //Find topography node for the marker
            xn=floor((MX[mm1][0]-gridt[0][0])/tstp-0.5)+1;
            
            if (xn<0){
                xn=0;
            }
            if (xn>tnum-2){
                xn=tnum-2;
            }
            
            // Compute relative distance to topography node
            dx=(MX[mm1][0]-gridt[0][xn])/tstp;
            // Compute topograhy elevation above the marker
            dy=gridt[1][xn]*(1-dx)+gridt[1][xn+1]*dx;
            // water/air to sediments transformation
            if (MI[mm1][0]==1 && MY[mm1][0]>dy){
                MI[mm1][0]=2; // Change marker type
                MRAT[mm1][0]=1; // Reset strain rate ratio                  
                MGII[mm1][0]=0; // Reset strain
            }
            // Rocks to water/air transformation
            if (MI[mm1][0]>1 && MY[mm1][0]<dy){
                MI[mm1][0]=1; // Change marker type
                MRAT[mm1][0]=1; // Reset strain rate ratio
                MGII[mm1][0]=0; // Reset strain
           }
            //  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            //           ?           ^                  ?
            //           ?           ?                  ?
            //           ?          dy                  ?
            //           ?           ?                  ?
            //           ?           v                  ?
            //           ?<----dx--->o Mrho(xm,ym)       ?
            //           ?                              ?
            //           ?                              ?
            //  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            //
            // Define indexes for upper left node in the cell where the marker is
            // using bisection
            // Find horizontal index
            xnmin=0;
            
            //continue;
            
            xnmax=xnum;
            while ((xnmax-xnmin)>1){            //TODO THIS CAUSES HANGING!!!!!! WHY!?!
                // !!! SUBTRACT 0.5 since int16(0.5)=1
                //xn=rintf((xnmax+xnmin)/2-0.5);
                //xn=floor((xnmax+xnmin)/2);
                xn=floor((xnmax+xnmin)/2);
                //grid_ind = (int)xn;     //??????
                //grid_ind = (int)xn;
                if(gridx[xn]>MX[mm1][0]){
                    xnmax=xn;
                    }
                else{
                    xnmin=xn;
                    }
                }
                                               
            
            //continue;               
            
            xn=xnmin;
            // Check horizontal index
            if (xn<0){
                xn=0;
            }
            if (xn>xnum-2){
                xn=xnum-2;
            }
            
            
            // Save horizontal index
            MXN[mm1][0]=xn;
            // Find vertical index
            ynmin=0;
            ynmax=ynum;
            while ((ynmax-ynmin)>1){
                // !!! SUBTRACT 0.5 since int16(0.5)=1
                //yn=rintf((ynmax+ynmin)/2-0.5);
                yn=floor((ynmax+ynmin)/2);
                if(gridy[yn]>MY[mm1][0]){
                    ynmax=yn;
                }
                else{
                    ynmin=yn;
                    }
            }

            yn=ynmin;
            // Check vertical index
            
            
            if (yn<0){
                yn=0;
                }
            if (yn>ynum-2){
                yn=ynum-2;
                }
            // Save Vertical index
            MYN[mm1][0]=yn;

            // Define normalized distances from marker to the upper left node;
            dx=(MX[mm1][0]-gridx[xn])/xstp1[xn][0];
            dy=(MY[mm1][0]-gridy[yn])/ystp1[yn][0];

            //Compute marker weight koefficient from cell dimensions
            //Number of markers in a cell is in invert proportion to the cell volume
            mwt=1; // /xstp1(xn)/ystp1(yn);

            // Compute density from marker temperature
            MRHOCUR = MRHO[MI[mm1][0]][1]*(1-MRHO[MI[mm1][0]][2]*(MTK[mm1][0]-273))*(1+MRHO[MI[mm1][0]][3]*(MPR[mm1][0]-1e+5));

            // Compute rho*Cp for marker 
            MRHOCPCUR= MRHOCUR * MCP[MI[mm1][0]];
            
            //custom_print('Change density for \"air\"')
            // Change density for "air"
            if (MI[mm1][0]==1 && MY[mm1][0]<waterlev){
                MRHOCUR=1;
            }
            
            // Compute thermal conductivity from marker temperature
            // Rock thermal conductivity (W/m/K): k=k0+a/(T+77)
            MKTCUR=MKT[MI[mm1][0]][1]+MKT[MI[mm1][0]][2]/(MTK[mm1][0]+77);
            
            // Compute adiabatic heating term (alp*T*DP/Dt)
            MHACUR=MRHO[MI[mm1][0]][2]*MTK[mm1][0];
            
            // Computing Marker Viscosity
            //printf("Mark Computing Marker Viscosity"); 
            
            if(MFLOW[MI[mm1][0]][1]==0){
                // Constant viscosity
                METACUR=MFLOW[MI[mm1][0]][2];
                }
            else{
                // Power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
                // Iterate for viscosity
                // First viscosity value
                // Compute and check old marker stress invariant in Pa
                sii0=pow((pow(MSXX[mm1][0],2)+pow(MSXY[mm1][0],2)),0.5);
                // Check old marker stress invariant (should be allways positive to be used in power law)
                if(sii0<stressmin){
                    sii0=stressmin;
                    }
                }
                // Check marker temperature
                plawexp=MTK[mm1][0];
                if(plawexp<ttop){
                    plawexp=ttop;
                }
                // Compute exponential term: 
                // Ea is in J/mol(=1000*kJ/mol)
                // Va is in J/Pa (=1e-6*cm^3) 
                // Cut if too big (at cold temperature);
                plawexp=(MFLOW[MI[mm1][0]][4]*1000+MFLOW[MI[mm1][0]][5]*1e-6*MPR[mm1][0])/RGAS/plawexp;
                if(plawexp>150){
                    plawexp=150;
                }
                // Compute AD*exp[-Ea/RT)
                plawexp=MFLOW[MI[mm1][0]][2]*exp(-plawexp);
                // Compute strain rate invariant from power law
                eii0=plawexp*pow((1e-6*sii0),MFLOW[MI[mm1][0]][3]);
                // Compute effective viscosity
                eta0=sii0/2/eii0;
                // Forcasting second invariant of future marker stress for given viscoelastic timestep
                xelvis=eta0/(MMU[MI[mm1][0]]*timestep+eta0);
                sxxnew=MSXX[mm1][0]*xelvis+2*eta0*MEXX[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                sxynew=MSXY[mm1][0]*xelvis+2*eta0*MEXY[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                sii1=pow((pow(sxxnew,2)+pow(sxynew,2)),0.5);
                // Check new marker stress invariant (should be allways positive to be used in power law)
                if(sii1<stressmin){
                    sii1=stressmin;
                }
                // Compute strain rate invariant from power law
                eii1=plawexp*pow((1e-6*sii1),MFLOW[MI[mm1][0]][3]);
                // Compute effective viscosity
                METACUR=sii1/2/eii1;
                // Iterate for viscosity which corresponds to future stress invariant using bisection
                // Iteration counter
                plawiter=0;
                while(plawiter<20 && fabs(sii1-sii0)>1){
                    // Add iteration counter
                    plawiter=plawiter+1;
                    // Compute middle stress
                    siicur=(sii0+sii1)/2;
                    // Compute strain rate invariant from power law
                    eiicur=plawexp*pow((1e-6*siicur),MFLOW[MI[mm1][0]][3]);
                    // Compute effective viscosity
                    METACUR=siicur/2/eiicur;
                    // Forcasting second invariant of future marker stress for given viscoelastic timestep
                    xelvis=METACUR/(MMU[MI[mm1][0]]*timestep+METACUR);
                    sxxnew=MSXX[mm1][0]*xelvis+2*METACUR*MEXX[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                    sxynew=MSXY[mm1][0]*xelvis+2*METACUR*MEXY[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                    siinew=pow((pow(sxxnew,2)+pow(sxynew,2)),0.5);
                    // Changing bisection limits
                    if((sii0<sii1 && siicur<siinew) || (sii0>sii1 && siicur>siinew)){
                        sii0=siicur;
                        }
                    else{
                        sii1=siicur;
                   }
                }
                // Limiting viscosity for the power law
                if (METACUR<etamin){
                    METACUR=etamin;
                }
                if (METACUR>etamax){ 
                    METACUR=etamax;
                    }
            
            
            }
            //printf("Mark Check if any plastic yeiding condition is present ");        
            
            // Check if any plastic yeiding condition is present 
            if (ntimestep>1 && (MPL[MI[mm1][0]][1]>0 || MPL[MI[mm1][0]][3]>0)){
                // Checking for plastic yeilding
                // Compute second invariant for a purely elastic stress build-up
                sxxnewe=MSXX[mm1][0]+2*MMU[MI[mm1][0]]*timestep*MEXX[mm1][0]*MRAT[mm1][0];
                sxynewe=MSXY[mm1][0]+2*MMU[MI[mm1][0]]*timestep*MEXY[mm1][0]*MRAT[mm1][0];
                siinewe=pow((pow(sxxnewe,2)+pow(sxynewe,2)),0.5);
                // Checking yeilding criterion for strain weakening/hardening
                // C=C0, FI=FI0 for strain<=GAM0
                // C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
                // C=C1, FI=FI1 for strain>=GAM0
                MCOHES=MPL[MI[mm1][0]][1];
                MFRICT=MPL[MI[mm1][0]][3];
                if (MGII[mm1][0]>=MPL[MI[mm1][0]][6]){
                    MCOHES=MPL[MI[mm1][0]][2];
                    MFRICT=MPL[MI[mm1][0]][4];
                }
                if (MGII[mm1][0]>MPL[MI[mm1][0]][5] && MGII[mm1][0]<MPL[MI[mm1][0]][6]){
                    MCOHES=MPL[MI[mm1][0]][1]+(MPL[MI[mm1][0]][2]-MPL[MI[mm1][0]][1])/(MPL[MI[mm1][0]][6]-MPL[MI[mm1][0]][5])*(MGII[mm1][0]-MPL[MI[mm1][0]][5]);
                    MFRICT=MPL[MI[mm1][0]][3]+(MPL[MI[mm1][0]][4]-MPL[MI[mm1][0]][3])/(MPL[MI[mm1][0]][6]-MPL[MI[mm1][0]][5])*(MGII[mm1][0]-MPL[MI[mm1][0]][5]);
                }
                // Computing yelding stress for the marker
                siiyeld=MCOHES+MFRICT*MPR[mm1][0];
                if (siiyeld<0){ 
                    siiyeld=0;
                    }
                // Correcting rock properties for yeilding 
                if (siiyeld<siinewe){
                    // Bringing marker viscosity to yeilding stress
                    METACURP=MMU[MI[mm1][0]]*timestep*siiyeld/(siinewe-siiyeld);
                    METACURP=pow(METACURP,(1-etawt))*pow(META[mm1][0],etawt);
                    if(METACURP<METACUR){
                        METACUR=METACURP;
                        // Limiting viscosity for the yeilding
                        if (METACUR<etamin){ 
                            METACUR=etamin;
                        }
                        if (METACUR>etamax){
                            METACUR=etamax;
                        }
                        // Mark that plastic yeildind occur
                        plastyn=1;
                        // Mark that plastic strain needs to be accumulated
                        MGII[mm1][0]=fabs(MGII[mm1][0]);
                        }
                    else{
                        // Reset plastic strain if no yelding
                        MGII[mm1][0]=-1e-20;
                        }
                    }
            }
            // Compute 1/MU values (MU is shear modulus) 
            MMUCUR=1/MMU[MI[mm1][0]];
            
            //custom_print(get_now()+'Molten rocks.'); 
            //Molten rocks
            
            //printf("Mark #3"); 
            /*                         double      double     long int
            
            printf("mm1=%d",mm1);
            printf("\n");
            printf("MPR[mm1][0]=");
            printf("\n");
            printf("%f",MPR[mm1][0]);
            printf("\n");
            printf("MTK[mm1][0]=");
            printf("\n");
            printf("%f",MTK[mm1][0]);
            printf("\n");
            printf("MI[mm1][0]=");
            printf("\n");
            printf("%i",MI[mm1][0]);
            printf("\n");
            */
            //continue;
            //exit(0);
            /////////////////////
            
            res_m=Melt_fraction(MPR[mm1][0],MTK[mm1][0],MI[mm1][0]);
            xmelt = res_m.xmelt;
            hlat = res_m.hlat;
            // Save marker melting
            if(timesum>0){
                //custom_print('Save marker melt ratio');
                MXM[mm1][0]=xmelt;
                }
            //custom_print('xmelt={}'.format(xmelt));
            if(xmelt>0 && timesum>0){
                // Reset creep parameters for molten rocks                   
                MRAT[mm1][0]=1; // Reset strain rate ratio
                MGII[mm1][0]=0; // Reset strain
                // Viscosity of partially molten rocks
                if(xmelt>0.1){
                    METACUR=etamelt;
                }
                // Density
                MRHOCUR=MRHOCUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // RHO*CP
                MRHOCPCUR=MRHOCPCUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // Compute adiabatic heating term (alp*T*DP/Dt)
                MHACUR=MHACUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // Latent heating: effective adiabatic term, RHOCP
                if(xmelt<1){
                    // Melting adiabatic term: alpham=-rho*(dHlat/dP)/T
                    // Numerical differentiation
                    dp=1000; // Pressure increment, Pa
                    
                    res_m0=Melt_fraction(MPR[mm1][0]-dp,MTK[mm1][0],MI[mm1][0]);
                    res_m1=Melt_fraction(MPR[mm1][0]+dp,MTK[mm1][0],MI[mm1][0]);
                    
                    hlat0 = res_m0.hlat;
                    xmelt = res_m1.xmelt;
                    hlat1 = res_m1.hlat;
                    
                    MHACUR=MHACUR-(hlat1-hlat0)/(2.0*dp);
                    // Melting heat capacity term: cpm=dHlat/dT 
                    // Numerical differentiation 
                    double dt=1.0; // Temperature increment, Pa
                    res_m0=Melt_fraction(MPR[mm1][0],MTK[mm1][0]-dt,MI[mm1][0]);
                    res_m1=Melt_fraction(MPR[mm1][0],MTK[mm1][0]+dt,MI[mm1][0]);
                    
                    hlat0 = res_m0.hlat;
                    xmelt = res_m1.xmelt;
                    hlat1 = res_m1.hlat;
                    
                    //printf("hlat0=%f\n",hlat0);
                    //printf("xmelt=%f\n",xmelt);
                    
                    MRHOCPCUR=MRHOCPCUR+MRHOCUR*(hlat1-hlat0)/(2.0*dt);
                }
                }
            // Save marker viscosity
            META[mm1][0]=METACUR;

            // Add properties to 4 surrounding nodes
            rho1[yn][xn]=rho1[yn][xn]+(1.0-dx)*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn][xn]=tk1[yn][xn]+(1.0-dx)*(1.0-dy)*MTK[mm1][0]*mwt;
            kt1[yn][xn]=kt1[yn][xn]+(1.0-dx)*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn][xn]=rhocp1[yn][xn]+(1.0-dx)*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn][xn]=hr1[yn][xn]+(1.0-dx)*(1.0-dy)*MHR[MI[mm1][0]]*mwt;
            ha1[yn][xn]=ha1[yn][xn]+(1.0-dx)*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn][xn]=wtnodes[yn][xn]+(1.0-dx)*(1.0-dy)*mwt;

            rho1[yn+1][xn]=rho1[yn+1][xn]+(1.0-dx)*dy*MRHOCUR*mwt;
            tk1[yn+1][xn]=tk1[yn+1][xn]+(1.0-dx)*dy*MTK[mm1][0]*mwt;
            kt1[yn+1][xn]=kt1[yn+1][xn]+(1.0-dx)*dy*MKTCUR*mwt;
            rhocp1[yn+1][xn]=rhocp1[yn+1][xn]+(1.0-dx)*dy*MRHOCPCUR*mwt;
            hr1[yn+1][xn]=hr1[yn+1][xn]+(1.0-dx)*dy*MHR[MI[mm1][0]]*mwt;
            ha1[yn+1][xn]=ha1[yn+1][xn]+(1.0-dx)*dy*MHACUR*mwt;
            wtnodes[yn+1][xn]=wtnodes[yn+1][xn]+(1.0-dx)*dy*mwt;

            rho1[yn][xn+1]=rho1[yn][xn+1]+dx*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn][xn+1]=tk1[yn][xn+1]+dx*(1.0-dy)*MTK[mm1][0]*mwt;
            kt1[yn][xn+1]=kt1[yn][xn+1]+dx*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn][xn+1]=rhocp1[yn][xn+1]+dx*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn][xn+1]=hr1[yn][xn+1]+dx*(1.0-dy)*MHR[MI[mm1][0]]*mwt;
            ha1[yn][xn+1]=ha1[yn][xn+1]+dx*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn][xn+1]=wtnodes[yn][xn+1]+dx*(1.0-dy)*mwt;

            rho1[yn+1][xn+1]=rho1[yn+1][xn+1]+dx*dy*MRHOCUR*mwt;
            tk1[yn+1][xn+1]=tk1[yn+1][xn+1]+dx*dy*MTK[mm1][0]*mwt;
            kt1[yn+1][xn+1]=kt1[yn+1][xn+1]+dx*dy*MKTCUR*mwt;
            rhocp1[yn+1][xn+1]=rhocp1[yn+1][xn+1]+dx*dy*MRHOCPCUR*mwt;
            hr1[yn+1][xn+1]=hr1[yn+1][xn+1]+dx*dy*MHR[MI[mm1][0]]*mwt;
            ha1[yn+1][xn+1]=ha1[yn+1][xn+1]+dx*dy*MHACUR*mwt;
            wtnodes[yn+1][xn+1]=wtnodes[yn+1][xn+1]+dx*dy*mwt;

            //Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
            // only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 && dy<=0.5){
                etas1[yn][xn]=etas1[yn][xn]+(1.0-dx)*(1.0-dy)*METACUR*mwt;
                mus1[yn][xn]=mus1[yn][xn]+(1.0-dx)*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn][xn]=sxy1[yn][xn]+(1.0-dx)*(1.0-dy)*MSXY[mm1][0]*mwt;
                wtetas[yn][xn]=wtetas[yn][xn]+(1.0-dx)*(1.0-dy)*mwt;
            }
            if(dx<=0.5 && dy>=0.5){
                etas1[yn+1][xn]=etas1[yn+1][xn]+(1.0-dx)*dy*METACUR*mwt;
                mus1[yn+1][xn]=mus1[yn+1][xn]+(1.0-dx)*dy*MMUCUR*mwt;
                sxy1[yn+1][xn]=sxy1[yn+1][xn]+(1.0-dx)*dy*MSXY[mm1][0]*mwt;
                wtetas[yn+1][xn]=wtetas[yn+1][xn]+(1.0-dx)*dy*mwt;
                }
            
            if(dx>=0.5 && dy<=0.5){
                etas1[yn][xn+1]=etas1[yn][xn+1]+dx*(1.0-dy)*METACUR*mwt;
                mus1[yn][xn+1]=mus1[yn][xn+1]+dx*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn][xn+1]=sxy1[yn][xn+1]+dx*(1.0-dy)*MSXY[mm1][0]*mwt;
                wtetas[yn][xn+1]=wtetas[yn][xn+1]+dx*(1.0-dy)*mwt;
            }
            if(dx>=0.5 && dy>=0.5){
                etas1[yn+1][xn+1]=etas1[yn+1][xn+1]+dx*dy*METACUR*mwt;
                mus1[yn+1][xn+1]=mus1[yn+1][xn+1]+dx*dy*MMUCUR*mwt;
                sxy1[yn+1][xn+1]=sxy1[yn+1][xn+1]+dx*dy*MSXY[mm1][0]*mwt;
                wtetas[yn+1][xn+1]=wtetas[yn+1][xn+1]+dx*dy*mwt;
            }

            // Add viscosity etan(), normal stress sxx() and shear modulus mun() to the center of current cell
            etan1[yn][xn]=etan1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*METACUR*mwt;
            mun1[yn][xn]=mun1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*MMUCUR*mwt;
            sxx1[yn][xn]=sxx1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*MSXX[mm1][0]*mwt;
            wtetan[yn][xn]=wtetan[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*mwt;
    }
    
    //printf("Mark #4");        
    //todo have to return structure
    //end of interpolating parameters from markers to nodes
    //return rho1,tk1,kt1,rhocp1,hr1, ha1, wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum
    /*
    res.rho1=ptr2double_arr(rho1,ynum,xnum);
    res.tk1=ptr2double_arr(tk1,ynum,xnum);
    res.kt1=ptr2double_arr(kt1,ynum,xnum);
    res.rhocp1=ptr2double_arr(rhocp1,ynum,xnum);
    res.hr1=ptr2double_arr(hr1,ynum,xnum);
    res.ha1=ptr2double_arr(ha1,ynum,xnum);
    res.wtnodes=ptr2double_arr(wtnodes,ynum,xnum);
    res.etas1=ptr2double_arr(etas1,ynum,xnum);
    res.mus1=ptr2double_arr(mus1,ynum,xnum);
    res.sxy1=ptr2double_arr(sxy1,ynum,xnum);
    res.wtetas=ptr2double_arr(wtetas,ynum,xnum);
    res.etan1=ptr2double_arr(etan1,ynum-1,xnum-1);
    res.mun1=ptr2double_arr(mun1,ynum-1,xnum-1);
    res.sxx1=ptr2double_arr(sxx1,ynum-1,xnum-1);
    res.wtetan=ptr2double_arr(wtetan,ynum-1,xnum-1);
    */
    res.rho1=ptr2double_arr_new(xnum,rho1,ynum);
    res.tk1=ptr2double_arr_new(xnum,tk1,ynum);
    res.kt1=ptr2double_arr_new(xnum,kt1,ynum);
    res.rhocp1=ptr2double_arr_new(xnum,rhocp1,ynum);
    res.hr1=ptr2double_arr_new(xnum,hr1,ynum);
    res.ha1=ptr2double_arr_new(xnum,ha1,ynum);
    res.wtnodes=ptr2double_arr_new(xnum,wtnodes,ynum);
    res.etas1=ptr2double_arr_new(xnum,etas1,ynum);
    res.mus1=ptr2double_arr_new(xnum,mus1,ynum);
    res.sxy1=ptr2double_arr_new(xnum,sxy1,ynum);
    res.wtetas=ptr2double_arr_new(xnum,wtetas,ynum);
    res.etan1=ptr2double_arr_new(xnum-1,etan1,ynum-1);
    res.mun1=ptr2double_arr_new(xnum-1,mun1,ynum-1);
    res.sxx1=ptr2double_arr_new(xnum-1,sxx1,ynum-1);
    res.wtetan=ptr2double_arr_new(xnum-1,wtetan,ynum-1);
    res.timesum=timesum;
    res.plastyn = plastyn;
    
    return res;
    }    
    
//below is old interpolation function for attention and comparison

struct result_t interpolate_markers_nodes_old(
                int marknum,
                double *gridx,
                double *gridy,
                size_t xnum,
                size_t ynum,
                double (*MX)[2],
                double (*MY)[2],
                int (*MI)[3],//int (*MI)[3],     //right shape for vector column
                double (*MRAT)[2], 
                double (*MGII)[2],
                double (*MXN)[2],
                double (*MYN)[2],
                double (*MFLOW)[2],
                double (*MPL)[2],
                double (*MXM)[2],
                double tstp,
                long int tnum,
                double (*gridt)[2],
                int waterlev,
                double stressmin,
                int ttop,
                double etamelt,
                double (*rho1)[3],
                double (*tk1)[3],
                double (*kt1)[3],
                double (*rhocp1)[3],
                double (*hr1)[3],
                double (*ha1)[3],
                double (*wtnodes)[3],
                double (*etas1)[3],
                double (*mus1)[3],
                double (*sxy1)[3],
                double (*wtetas)[3],
                double (*etan1)[3],
                double (*mun1)[3],
                double (*sxx1)[3],
                double (*wtetan)[3],
                double timesum, //new parameters are added from here
                long int *xstp1,
                long int *ystp1,
                double (*MRHO)[3],
                double (*MKT)[3],
                double (*MTK)[3],
                double (*MPR)[2],
                double *MCP,
                double *MHR,
                double (*MSXX)[2],
                double (*MSXY)[2],
                double (*MEXX)[2],
                double (*MEXY)[2],
                double (*META)[2],
                double *MMU,
                double timestep,
                double etamin,
                double etamax,
                int ntimestep,
                double etawt,
                int plastyn
                )
                
                    //double (*v)[3], size_t m, size_t n)
{   
//Interpolating parameters from markers to nodes
long int mm1;
long int xn,yn,dp;
double dx, dy,mwt;
long int xnmin,xnmax,ynmin,ynmax;
long int plawiter,grid_ind;



struct result_t res; //return value structure

struct result_melt res_m,res_m0,res_m1; //structure for local result melt calls

//!!
//TODO do memory allocation for every matrix created inside of function
//!!


double MRHOCUR,MRHOCPCUR,MKTCUR,MHACUR,METACUR,METACURP,sii0,
            plawexp,xelvis, sxxnew,sxynew,
            siinewe,sii1,eii0,eii1,siicur,siinew,
            eiicur,MCOHES,MFRICT,
            xmelt,hlat,hlat1,hlat0,eta0,sxxnewe,sxynewe,siiyeld,MMUCUR;

double RGAS = 8.314; //Gas constant J/mol/K

//printf("Mark #1");


//TODO НАДО ПЕРЕД ПЕРЕДАЧЕЙ В C КОНВЕРТИРОВАТЬ NUMPY МАССИВЫ В ФОРМАТ, КОТОРЫЙ МОЖЕТ ИНДЕКСИРВОАТЬ C
/*
printf("%s","MRHO[1][1]=");
printf("%f",MRHO[1][1]);
printf("\n");

printf("%s","MI[100][0]=");
printf("%i",MI[100][0]);
printf("\n");

exit(0);

*/
//MI[mm1][0]


for (mm1 = 0; mm1 < marknum; mm1++) {
        // Check markers inside the grid
        
        
        //printf("mark num=%i",mm1);
        //printf(" out of %i",marknum);
        /*
        printf("rock type=%i",MI[mm1][0]); //каждый второй маркер 0 НУЛЕВОЙ?! ПОЧЕМУ!
        printf("\n");

        continue;
        */
 
        if (MX[mm1][0]>=gridx[0] && MX[mm1][0]<=gridx[xnum-1] && 
            MY[mm1][0]>=gridy[0] && MY[mm1][0]<=gridy[ynum-1]){ 
            //Erosion-sedimentation
            //Find topography node for the marker
            xn=floor((MX[mm1][0]-gridt[0][0])/tstp-0.5)+1;
            
            if (xn<0){
                xn=0;
            }
            if (xn>tnum-2){
                xn=tnum-2;
            }
            
            // Compute relative distance to topography node
            dx=(MX[mm1][0]-gridt[0][xn])/tstp;
            // Compute topograhy elevation above the marker
            dy=gridt[1][xn]*(1-dx)+gridt[1][xn+1]*dx;
            // water/air to sediments transformation
            if (MI[mm1][0]==1 && MY[mm1][0]>dy){
                MI[mm1][0]=2; // Change marker type
                MRAT[mm1][0]=1; // Reset strain rate ratio                  
                MGII[mm1][0]=0; // Reset strain
            }
            // Rocks to water/air transformation
            if (MI[mm1][0]>1 && MY[mm1][0]<dy){
                MI[mm1][0]=1; // Change marker type
                MRAT[mm1][0]=1; // Reset strain rate ratio
                MGII[mm1][0]=0; // Reset strain
           }
            //  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            //           ?           ^                  ?
            //           ?           ?                  ?
            //           ?          dy                  ?
            //           ?           ?                  ?
            //           ?           v                  ?
            //           ?<----dx--->o Mrho(xm,ym)       ?
            //           ?                              ?
            //           ?                              ?
            //  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            //
            // Define indexes for upper left node in the cell where the marker is
            // using bisection
            // Find horizontal index
            xnmin=0;
            
            
            xnmax=xnum;
            while ((xnmax-xnmin)>1){
                // !!! SUBTRACT 0.5 since int16(0.5)=1
                xn=floor((xnmax+xnmin)/2);
                //grid_ind = (int)xn;
                if(gridx[xn]>MX[mm1][0]){
                    xnmax=xn;
                    }
                else{
                    xnmin=xn;
                    }
                }
            }
            xn=xnmin;
            // Check horizontal index
            if (xn<0){
                xn=0;
            }
            if (xn>xnum-2){
                xn=xnum-2;
            }
            
            
            // Save horizontal index
            MXN[mm1][0]=xn;
            // Find vertical index
            ynmin=0;
            ynmax=ynum;
            while ((ynmax-ynmin)>1){
                // !!! SUBTRACT 0.5 since int16(0.5)=1
                yn=floor((ynmax+ynmin)/2);
                if(gridy[yn]>MY[mm1][0]){
                    ynmax=yn;
                }
                else{
                    ynmin=yn;
                    }
            }

            yn=ynmin;
            // Check vertical index
            
            
            
            if (yn<0){
                yn=0;
                }
            if (yn>ynum-2){
                yn=ynum-2;
                }
            // Save Vertical index
            MYN[mm1][0]=yn;

            // Define normalized distances from marker to the upper left node;
            dx=(MX[mm1][0]-gridx[xn])/xstp1[xn];
            dy=(MY[mm1][0]-gridy[yn])/ystp1[yn];

            //Compute marker weight koefficient from cell dimensions
            //Number of markers in a cell is in invert proportion to the cell volume
            mwt=1; // /xstp1(xn)/ystp1(yn);

            // Compute density from marker temperature
            MRHOCUR = MRHO[MI[mm1][0]][1]*(1-MRHO[MI[mm1][0]][2]*(MTK[mm1][0]-273))*(1+MRHO[MI[mm1][0]][3]*(MPR[mm1][0]-1e+5));

            // Compute rho*Cp for marker 
            MRHOCPCUR= MRHOCUR * MCP[MI[mm1][0]];
            
            //custom_print('Change density for \"air\"')
            // Change density for "air"
            if (MI[mm1][0]==1 && MY[mm1][0]<waterlev){
                MRHOCUR=1;
            }
            
            // Compute thermal conductivity from marker temperature
            // Rock thermal conductivity (W/m/K): k=k0+a/(T+77)
            MKTCUR=MKT[MI[mm1][0]][1]+MKT[MI[mm1][0]][2]/(MTK[mm1][0]+77);
            
            // Compute adiabatic heating term (alp*T*DP/Dt)
            MHACUR=MRHO[MI[mm1][0]][2]*MTK[mm1][0];
            
            // Computing Marker Viscosity
            //printf("Mark Computing Marker Viscosity"); 
            
            if(MFLOW[MI[mm1][0]][1]==0){
                // Constant viscosity
                METACUR=MFLOW[MI[mm1][0]][2];
                }
            else{
                // Power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
                // Iterate for viscosity
                // First viscosity value
                // Compute and check old marker stress invariant in Pa
                sii0=pow((pow(MSXX[mm1][0],2)+pow(MSXY[mm1][0],2)),0.5);
                // Check old marker stress invariant (should be allways positive to be used in power law)
                if(sii0<stressmin){
                    sii0=stressmin;
                    }
                
                // Check marker temperature
                plawexp=MTK[mm1][0];
                if(plawexp<ttop){
                    plawexp=ttop;
                }
                // Compute exponential term: 
                // Ea is in J/mol(=1000*kJ/mol)
                // Va is in J/Pa (=1e-6*cm^3) 
                // Cut if too big (at cold temperature);
                plawexp=(MFLOW[MI[mm1][0]][4]*1000+MFLOW[MI[mm1][0]][5]*1e-6*MPR[mm1][0])/RGAS/plawexp;
                if(plawexp>150){
                    plawexp=150;
                }
                // Compute AD*exp[-Ea/RT)
                plawexp=MFLOW[MI[mm1][0]][2]*exp(-plawexp);
                // Compute strain rate invariant from power law
                eii0=pow(plawexp*(1e-6*sii0),MFLOW[MI[mm1][0]][3]);
                // Compute effective viscosity
                eta0=sii0/2/eii0;
                // Forcasting second invariant of future marker stress for given viscoelastic timestep
                xelvis=eta0/(MMU[MI[mm1][0]]*timestep+eta0);
                sxxnew=MSXX[mm1][0]*xelvis+2*eta0*MEXX[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                sxynew=MSXY[mm1][0]*xelvis+2*eta0*MEXY[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                sii1=pow((pow(sxxnew,2)+pow(sxynew,2)),0.5);
                // Check new marker stress invariant (should be allways positive to be used in power law)
                if(sii1<stressmin){
                    sii1=stressmin;
                }
                // Compute strain rate invariant from power law
                eii1=plawexp*pow((1e-6*sii1),MFLOW[MI[mm1][0]][3]);
                // Compute effective viscosity
                METACUR=sii1/2/eii1;
                // Iterate for viscosity which corresponds to future stress invariant using bisection
                // Iteration counter
                plawiter=0;
                while(plawiter<20 && fabs(sii1-sii0)>1){
                    // Add iteration counter
                    plawiter=plawiter+1;
                    // Compute middle stress
                    siicur=(sii0+sii1)/2;
                    // Compute strain rate invariant from power law
                    eiicur=plawexp*pow((1e-6*siicur),MFLOW[MI[mm1][0]][3]);
                    // Compute effective viscosity
                    METACUR=siicur/2/eiicur;
                    // Forcasting second invariant of future marker stress for given viscoelastic timestep
                    xelvis=METACUR/(MMU[MI[mm1][0]]*timestep+METACUR);
                    sxxnew=MSXX[mm1][0]*xelvis+2*METACUR*MEXX[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                    sxynew=MSXY[mm1][0]*xelvis+2*METACUR*MEXY[mm1][0]*MRAT[mm1][0]*(1-xelvis);
                    siinew=pow((pow(sxxnew,2)+pow(sxynew,2)),0.5);
                    // Changing bisection limits
                    if((sii0<sii1 && siicur<siinew) || (sii0>sii1 && siicur>siinew)){
                        sii0=siicur;
                        }
                    else{
                        sii1=siicur;
                   }
                }
                // Limiting viscosity for the power law
                if (METACUR<etamin){
                    METACUR=etamin;
                }
                if (METACUR>etamax){ 
                    METACUR=etamax;
                    }
            
             }
            
            //printf("Mark Check if any plastic yeiding condition is present ");        
            
            // Check if any plastic yeiding condition is present 
            if (ntimestep>1 && (MPL[MI[mm1][0]][1]>0 || MPL[MI[mm1][0]][3]>0)){
                // Checking for plastic yeilding
                // Compute second invariant for a purely elastic stress build-up
                sxxnewe=MSXX[mm1][0]+2*MMU[MI[mm1][0]]*timestep*MEXX[mm1][0]*MRAT[mm1][0];
                sxynewe=MSXY[mm1][0]+2*MMU[MI[mm1][0]]*timestep*MEXY[mm1][0]*MRAT[mm1][0];
                siinewe=pow((pow(sxxnewe,2)+pow(sxynewe,2)),0.5);
                // Checking yeilding criterion for strain weakening/hardening
                // C=C0, FI=FI0 for strain<=GAM0
                // C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
                // C=C1, FI=FI1 for strain>=GAM0
                MCOHES=MPL[MI[mm1][0]][1];
                MFRICT=MPL[MI[mm1][0]][3];
                if (MGII[mm1][0]>=MPL[MI[mm1][0]][6]){
                    MCOHES=MPL[MI[mm1][0]][2];
                    MFRICT=MPL[MI[mm1][0]][4];
                }
                if (MGII[mm1][0]>MPL[MI[mm1][0]][5] && MGII[mm1][0]<MPL[MI[mm1][0]][6]){
                    MCOHES=MPL[MI[mm1][0]][1]+(MPL[MI[mm1][0]][2]-MPL[MI[mm1][0]][1])/(MPL[MI[mm1][0]][6]-MPL[MI[mm1][0]][5])*(MGII[mm1][0]-MPL[MI[mm1][0]][5]);
                    MFRICT=MPL[MI[mm1][0]][3]+(MPL[MI[mm1][0]][4]-MPL[MI[mm1][0]][3])/(MPL[MI[mm1][0]][6]-MPL[MI[mm1][0]][5])*(MGII[mm1][0]-MPL[MI[mm1][0]][5]);
                }
                // Computing yelding stress for the marker
                siiyeld=MCOHES+MFRICT*MPR[mm1][0];
                if (siiyeld<0){ 
                    siiyeld=0;
                    }
                // Correcting rock properties for yeilding 
                if (siiyeld<siinewe){
                    // Bringing marker viscosity to yeilding stress
                    METACURP=MMU[MI[mm1][0]]*timestep*siiyeld/(siinewe-siiyeld);
                    METACURP=pow(METACURP,(1-etawt))*pow(META[mm1][0],etawt);
                    
                    if(METACURP<METACUR){
                        METACUR=METACURP;
                        // Limiting viscosity for the yeilding
                        if (METACUR<etamin){ 
                            METACUR=etamin;
                        }
                        if (METACUR>etamax){
                            METACUR=etamax;
                        }
                        // Mark that plastic yeildind occur
                        plastyn=1;
                        // Mark that plastic strain needs to be accumulated
                        MGII[mm1][0]=fabs(MGII[mm1][0]);
                        }
                    else{
                        // Reset plastic strain if no yelding
                        MGII[mm1][0]=-1e-20;
                        }
                    }
            }
            // Compute 1/MU values (MU is shear modulus) 
            MMUCUR=1/MMU[MI[mm1][0]];
            
            //custom_print(get_now()+'Molten rocks.'); 
            //Molten rocks
            
            //printf("Mark #3"); 
            //                         double      double     long int
            
            /*
            printf("MPR[mm1][0]=");
            printf("\n");
            printf("%f",MPR[mm1][0]);
            printf("\n");
            printf("MTK[mm1][0]=");
            printf("\n");
            printf("%f",MTK[mm1][0]);
            printf("\n");
            printf("MI[mm1][0]=");
            printf("\n");
            printf("%i",MI[mm1][0]);
            printf("\n");
            
            exit(0);
            */
            res_m=Melt_fraction(MPR[mm1][0],MTK[mm1][0],MI[mm1][0]);
            xmelt = res_m.xmelt;
            hlat = res_m.hlat;
            // Save marker melting
            if(timesum>0){
                //custom_print('Save marker melt ratio');
                MXM[mm1][0]=xmelt;
                }
            //custom_print('xmelt={}'.format(xmelt));
            if(xmelt>0 && timesum>0){
                // Reset creep parameters for molten rocks                   
                MRAT[mm1][0]=1; // Reset strain rate ratio
                MGII[mm1][0]=0; // Reset strain
                // Viscosity of partially molten rocks
                if(xmelt>0.1){
                    METACUR=etamelt;
                }
                // Density
                MRHOCUR=MRHOCUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // RHO*CP
                MRHOCPCUR=MRHOCPCUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // Compute adiabatic heating term (alp*T*DP/Dt)
                MHACUR=MHACUR*((1-xmelt)+MRHO[MI[mm1][0]][4]/MRHO[MI[mm1][0]][1]*xmelt);
                // Latent heating: effective adiabatic term, RHOCP
                if(xmelt<1){
                    // Melting adiabatic term: alpham=-rho*(dHlat/dP)/T
                    // Numerical differentiation
                    dp=1000; // Pressure increment, Pa
                    
                    res_m0=Melt_fraction(MPR[mm1][0]-dp,MTK[mm1][0],MI[mm1][0]);
                    res_m1=Melt_fraction(MPR[mm1][0]+dp,MTK[mm1][0],MI[mm1][0]);
                    
                    hlat0 = res_m0.hlat;
                    xmelt = res_m1.xmelt;
                    hlat1 = res_m1.hlat;
                    
                    MHACUR=MHACUR-(hlat1-hlat0)/(2.0*dp);
                    // Melting heat capacity term: cpm=dHlat/dT 
                    // Numerical differentiation 
                    double dt=1.0; // Temperature increment, Pa
                    res_m0=Melt_fraction(MPR[mm1][0],MTK[mm1][0]-dt,MI[mm1][0]);
                    res_m1=Melt_fraction(MPR[mm1][0],MTK[mm1][0]+dt,MI[mm1][0]);
                    
                    hlat0 = res_m0.hlat;
                    xmelt = res_m1.xmelt;
                    hlat1 = res_m1.hlat;
                    
                    MRHOCPCUR=MRHOCPCUR+MRHOCUR*(hlat1-hlat0)/(2.0*dt);
                }
                }
            // Save marker viscosity
            META[mm1][0]=METACUR;

            // Add properties to 4 surrounding nodes
            rho1[yn][xn]=rho1[yn][xn]+(1.0-dx)*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn][xn]=tk1[yn][xn]+(1.0-dx)*(1.0-dy)*MTK[mm1][0]*mwt;
            kt1[yn][xn]=kt1[yn][xn]+(1.0-dx)*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn][xn]=rhocp1[yn][xn]+(1.0-dx)*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn][xn]=hr1[yn][xn]+(1.0-dx)*(1.0-dy)*MHR[MI[mm1][0]]*mwt;
            ha1[yn][xn]=ha1[yn][xn]+(1.0-dx)*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn][xn]=wtnodes[yn][xn]+(1.0-dx)*(1.0-dy)*mwt;

            rho1[yn+1][xn]=rho1[yn+1][xn]+(1.0-dx)*dy*MRHOCUR*mwt;
            tk1[yn+1][xn]=tk1[yn+1][xn]+(1.0-dx)*dy*MTK[mm1][0]*mwt;
            kt1[yn+1][xn]=kt1[yn+1][xn]+(1.0-dx)*dy*MKTCUR*mwt;
            rhocp1[yn+1][xn]=rhocp1[yn+1][xn]+(1.0-dx)*dy*MRHOCPCUR*mwt;
            hr1[yn+1][xn]=hr1[yn+1][xn]+(1.0-dx)*dy*MHR[MI[mm1][0]]*mwt;
            ha1[yn+1][xn]=ha1[yn+1][xn]+(1.0-dx)*dy*MHACUR*mwt;
            wtnodes[yn+1][xn]=wtnodes[yn+1][xn]+(1.0-dx)*dy*mwt;

            rho1[yn][xn+1]=rho1[yn][xn+1]+dx*(1.0-dy)*MRHOCUR*mwt;
            tk1[yn][xn+1]=tk1[yn][xn+1]+dx*(1.0-dy)*MTK[mm1][0]*mwt;
            kt1[yn][xn+1]=kt1[yn][xn+1]+dx*(1.0-dy)*MKTCUR*mwt;
            rhocp1[yn][xn+1]=rhocp1[yn][xn+1]+dx*(1.0-dy)*MRHOCPCUR*mwt;
            hr1[yn][xn+1]=hr1[yn][xn+1]+dx*(1.0-dy)*MHR[MI[mm1][0]]*mwt;
            ha1[yn][xn+1]=ha1[yn][xn+1]+dx*(1.0-dy)*MHACUR*mwt;
            wtnodes[yn][xn+1]=wtnodes[yn][xn+1]+dx*(1.0-dy)*mwt;

            rho1[yn+1][xn+1]=rho1[yn+1][xn+1]+dx*dy*MRHOCUR*mwt;
            tk1[yn+1][xn+1]=tk1[yn+1][xn+1]+dx*dy*MTK[mm1][0]*mwt;
            kt1[yn+1][xn+1]=kt1[yn+1][xn+1]+dx*dy*MKTCUR*mwt;
            rhocp1[yn+1][xn+1]=rhocp1[yn+1][xn+1]+dx*dy*MRHOCPCUR*mwt;
            hr1[yn+1][xn+1]=hr1[yn+1][xn+1]+dx*dy*MHR[MI[mm1][0]]*mwt;
            ha1[yn+1][xn+1]=ha1[yn+1][xn+1]+dx*dy*MHACUR*mwt;
            wtnodes[yn+1][xn+1]=wtnodes[yn+1][xn+1]+dx*dy*mwt;

            //Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
            // only using markers located at <=0.5 gridstep distances from nodes
            if(dx<=0.5 && dy<=0.5){
                etas1[yn][xn]=etas1[yn][xn]+(1.0-dx)*(1.0-dy)*METACUR*mwt;
                mus1[yn][xn]=mus1[yn][xn]+(1.0-dx)*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn][xn]=sxy1[yn][xn]+(1.0-dx)*(1.0-dy)*MSXY[mm1][0]*mwt;
                wtetas[yn][xn]=wtetas[yn][xn]+(1.0-dx)*(1.0-dy)*mwt;
            }
            if(dx<=0.5 && dy>=0.5){
                etas1[yn+1][xn]=etas1[yn+1][xn]+(1.0-dx)*dy*METACUR*mwt;
                mus1[yn+1][xn]=mus1[yn+1][xn]+(1.0-dx)*dy*MMUCUR*mwt;
                sxy1[yn+1][xn]=sxy1[yn+1][xn]+(1.0-dx)*dy*MSXY[mm1][0]*mwt;
                wtetas[yn+1][xn]=wtetas[yn+1][xn]+(1.0-dx)*dy*mwt;
                }
            
            if(dx>=0.5 && dy<=0.5){
                etas1[yn][xn+1]=etas1[yn][xn+1]+dx*(1.0-dy)*METACUR*mwt;
                mus1[yn][xn+1]=mus1[yn][xn+1]+dx*(1.0-dy)*MMUCUR*mwt;
                sxy1[yn][xn+1]=sxy1[yn][xn+1]+dx*(1.0-dy)*MSXY[mm1][0]*mwt;
                wtetas[yn][xn+1]=wtetas[yn][xn+1]+dx*(1.0-dy)*mwt;
            }
            if(dx>=0.5 && dy>=0.5){
                etas1[yn+1][xn+1]=etas1[yn+1][xn+1]+dx*dy*METACUR*mwt;
                mus1[yn+1][xn+1]=mus1[yn+1][xn+1]+dx*dy*MMUCUR*mwt;
                sxy1[yn+1][xn+1]=sxy1[yn+1][xn+1]+dx*dy*MSXY[mm1][0]*mwt;
                wtetas[yn+1][xn+1]=wtetas[yn+1][xn+1]+dx*dy*mwt;
            }

            // Add viscosity etan(), normal stress sxx() and shear modulus mun() to the center of current cell
            etan1[yn][xn]=etan1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*METACUR*mwt;
            mun1[yn][xn]=mun1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*MMUCUR*mwt;
            sxx1[yn][xn]=sxx1[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*MSXX[mm1][0]*mwt;
            wtetan[yn][xn]=wtetan[yn][xn]+(1.0-fabs(0.5-dx))*(1.0-fabs(0.5-dy))*mwt;
    }
    
    //printf("Mark #4");        
    //todo have to return structure
    //end of interpolating parameters from markers to nodes
    //return rho1,tk1,kt1,rhocp1,hr1, ha1, wtnodes,etas1,mus1,sxy1,wtetas,etan1,mun1,sxx1,wtetan,timesum
    res.rho1=ptr2double_arr(rho1,ynum,xnum);
    res.tk1=ptr2double_arr(tk1,ynum,xnum);
    res.kt1=ptr2double_arr(kt1,ynum,xnum);
    res.rhocp1=ptr2double_arr(rhocp1,ynum,xnum);
    res.hr1=ptr2double_arr(hr1,ynum,xnum);
    res.ha1=ptr2double_arr(ha1,ynum,xnum);
    res.wtnodes=ptr2double_arr(wtnodes,ynum,xnum);
    res.etas1=ptr2double_arr(etas1,ynum,xnum);
    res.mus1=ptr2double_arr(mus1,ynum,xnum);
    res.sxy1=ptr2double_arr(sxy1,ynum,xnum);
    res.wtetas=ptr2double_arr(wtetas,ynum,xnum);
    res.etan1=ptr2double_arr(etan1,ynum-1,xnum-1);
    res.mun1=ptr2double_arr(mun1,ynum-1,xnum-1);
    res.sxx1=ptr2double_arr(sxx1,ynum-1,xnum-1);
    res.wtetan=ptr2double_arr(wtetan,ynum-1,xnum-1);
    res.timesum=timesum;
    res.plastyn = plastyn;
    
    return res;
    }    

//function for testing only

struct result_test {
  double** Mat;
   };

   struct  result_test matrix_test(  //TODO REDO INPUTS IN Python code
                size_t xnum, //same as xnum NEW input parameter
                size_t ynum, //same as xnum NEW input parameter
                double matrix[][xnum])
                {
    //Do something with the matrix
   struct result_test res;
    
    //initialize empty array
    double arr[ynum][xnum];

    //double arr[][xnum] = malloc( ynum*xnum*sizeof(double) );
    
   //cooking something to return 
   //res.Mat = ptr2double_arr_new(xnum,matrix,ynum);

//    double **arr = malloc(ynum * sizeof(double));
 for (size_t i = 0; i < ynum; i++) {
     for (size_t j = 0; j < xnum; j++) {
         printf("%f ",matrix[i][j]);
         //res.Mat[i][j]=matrix[i][j];
         arr[i][j] = matrix[i][j];  //copy matrix to new array
     }
     printf("\n");
 }

 //create result matrix
 res.Mat = ptr2double_arr_new(xnum,arr,ynum);

//    for (int i = 0; i < ynum; i++) {
//        double arr[i] = malloc(xnum * sizeof(double));
       
//         for (int j = 0; j < xnum; j++) {
//             arr[i][j] = matrix[i][j]; //round 
//             }
//             }           
  

    //res.Mat = ptr2double_arr_new(xnum,matrix,ynum);
    //res.Mat = arr; 

    return res;     
}


int main()
{
printf("hello");
}
