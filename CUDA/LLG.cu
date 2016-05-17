#include <stdio.h>      
#include <math.h>    
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include </u/local/cuda/5.0/include/cuda.h>
#include </u/local/cuda/5.0/include/cuda_runtime.h>
#include </u/local/cuda/5.0/include/curand_kernel.h>
using namespace std;
//#define VARIATION
#define shaodi_pi 3.1415926
#define CUDA_CALL(x) do { if( (x) ! =  cudaSuccess ){\
	printf("Error at %s:%d\n",__FILE__,__LINE__ );\
	exit(1);} } while(0) 
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) {\
	printf("Error at %s:%d\n",__FILE__,__LINE__);\
	exit(1);}} while(0)


__global__ void LLG(double* g_v_pulse, double* finalR, int initial_state, double t_step, int trials_p_thread,int n_sim, bool enableVCMA, double ori_length, double sigma_l, double ori_width, double sigma_w, double ori_tfl, double sigma_tfl, double ori_Nx, double ori_Ny, double ori_Nz, double* g_lin_dep_factor  ){
/* -------------------------------------------
 Input Parameters From User
 -------------------------------------------*/
int this_id = (blockIdx.x * blockDim.x + threadIdx.x) ;
finalR[this_id] = 0 ;
//initiate state for following random generation
curandState_t localState;
curand_init(this_id, this_id, 0, &localState);

double Nx = ori_Nx;//origin:                                                // x Demagnetization factor
double Ny = ori_Ny;                                                // y Demagnetization factor
double Nz = ori_Nz;                                                // z Demagnetization factor
double length = ori_length; //nominal length of MTJ
double width = ori_width; //nominal width of MTJ
double tfl = ori_tfl; //nominal thickness of free layer
double areamtj = shaodi_pi*ori_length*ori_width/4;            // nominal area without variation
double R0 = 1e3;  				//Nominal MTJ resistance
double dMgO_a = 1.54e-3, dMgO_b = 1.1537278e10;  
double dMgO = (log(R0 * areamtj * 10e12) - log(dMgO_a)) / dMgO_b;     // Nominal MgO thickness [m]
double Temperature = 27+273;                       // Temperature
double TMR = 1;                                    // TMR at zero bias voltage

/******************Monte-Carlo simulation*****************/
for( int i_trial = 0; i_trial < trials_p_thread; i_trial++){ // Simulation for each trials
#ifdef VARIATION
//MTJ Dimention variation
length = ori_length + sigma_l*curand_normal_double(&localState); 
width = ori_width + sigma_w*curand_normal_double(&localState);
tfl = ori_tfl + sigma_tfl*curand_normal_double(&localState);
double temp_Nx = ori_Nx + g_lin_dep_factor[0]*(length-ori_length) + g_lin_dep_factor[1] * ( width - ori_width) + g_lin_dep_factor[2] * (tfl - ori_tfl) ;
double temp_Ny = ori_Ny + g_lin_dep_factor[3]*(length-ori_length) + g_lin_dep_factor[4] * ( width - ori_width) + g_lin_dep_factor[5] * (tfl - ori_tfl) ;
double temp_Nz = ori_Nz + g_lin_dep_factor[6]*(length-ori_length) + g_lin_dep_factor[7] * ( width - ori_width) + g_lin_dep_factor[8] * (tfl - ori_tfl) ;
Nx = temp_Nx / (temp_Nx + temp_Ny + temp_Nz);
Ny = temp_Ny / (temp_Nx + temp_Ny + temp_Nz);
Nz = temp_Nz / (temp_Nx + temp_Ny + temp_Nz); 
#endif

areamtj = shaodi_pi*length*width/4;            // Area with variation
double Rp = exp(dMgO * dMgO_b)*dMgO_a / (areamtj * 10e12);         // Parallel resistance [Ohms] with variation
double Rap = (1+TMR)*Rp;                                           // Anti-parallel resistance [Ohms]
double B1 = 0;                                   // Field-like torque linear parameter [unitless]
double B2 = 0;                                     // Field-like torque quadratic parameter [1/A]
//int initial_state is input;                          // Inital state [0 = parallel, 1 = anti-parallel]
double P [3] = {0, 0, -1};                             // Direction of polarization

/* -------------------------------------------
 Constants
 -------------------------------------------*/
double hbar = 1.05457173e-34;                                      // Reduced Planck constant, [J*s]
double k = 1.3806488e-23;                                          // Boltzmann constant, [J/K]
double u0 = 4e-7*shaodi_pi;                                               // Vacuum permeability, [V�s/(A�m)]
double q = 1.60217657e-19;                                         // Electron charge, [C]
double alphac = 0.02;                                              // LLGE damping factor
double gammap = (221276/(1+pow(alphac,2)));                             // Gyromagnetic ratio [m/(A x s)]
double T0 = 1120;
double Ms0 = 1393128.323;
double Ms = Ms0 * ( 1 - pow(Temperature/T0,1.5));                  // Saturation magnetization [A/m] - 1e6 A/m = 1000 emu/cc 
double dstray = 20e-9, tstray = 1.164656e-9;
double Ext[3]	 = {-Ms*length*width/4/shaodi_pi*((dstray+tstray)/(pow(length/2,2)*sqrt(pow(length/2,2)+pow(dstray+tstray,2)))-(dstray-tstray)/(pow(length/2,2)*sqrt(pow(length/2,2)+pow(dstray-tstray,2)))),0,0}; // External magnetic field [A/m] - 1 oersted [Oe] = 79.5774715459424 ampere/meter [A/m]
double Ki0 =1.479036e-3;
double Ki = Ki0 * pow(Ms/Ms0, 2.18);                          // Anisotropy field constant [J/m^2]
double Xi0 =53.39247e-15; 				// Energy barrier change slope for VCMA effect
double Xi = Xi0* pow(Ms/Ms0, 2.83);                                // VCMA field constant [J/(V x m)]
double Gt = 1/(Rp*(1+(TMR/(TMR+2))));                              // Direct elastic tunneling conductance [S]
double KiPF = (2*Ki)/(tfl*u0*Ms);                                  // Prefactor for interface anisotropy effective field
double VCMAPF = (2*Xi)/(u0*Ms*dMgO*tfl);                           // Prefactor for VCMA effective field
if (!enableVCMA){ // if VCMA is dsable
	VCMAPF = 0;
}
double Gsi	= 0;                                                    // Conductance due to imperfections in Mgo [S]

double P_tunnel = 0.2;                                  // the polarization of the tunnel currentdouble

double Pol = 0.6;                                                  // Polarization for Spin Torque

double volume = areamtj*tfl;                                       // MTJ volume [m^3]
double Hth = sqrt((2*k*Temperature*alphac)/(u0*gammap*Ms*volume*t_step));    // Amplitude of Thermal Field

/* -------------------------------------------
 Internal Variables
 -------------------------------------------*/

double costheta = 0;                                       // the angle between the magnization  of free and reference layers
double g_sv = 0;                                        // the polarization efficiency in spin valve
double g_tunnel = 0;                                    // the polarization efficiency in tunnel current


//double m_old [3] = {0, 0, 0};                              // Normalized previous magnetization
double Heff_old [3] = {0, 0, 0};                           // Previous Heff components [A/m]
double m_int [3] = {0, 0, 0};                              // Intermediate normalized magnetization
double dm_int [3] = {0, 0, 0};                             // Intermediate derivative of magnetization
double M_int [3] = {0, 0, 0};                              // Intermediate denormalized magnetization                              
//double Heff_int [3] = {0, 0, 0};                           // Intermediate Heff components [A/m]
double dm [3] = {0, 0, 0};                                 // Time derivative of magnetization [1/s]
double M [3] = {0, 0, 0};                                  // Denormalized magnetization
double mcrossp_int [3] = {0, 0, 0};                        // Intermediate cross product components (m x p)
double mcrossHeff_int [3] = {0, 0, 0};                     // Intermediate cross product components (m x Heff)
double mcrossHth_int [3] = {0, 0, 0};                      // Intermediate cross product components (m x Hth)
double mcrossmcrossp_int [3] = {0, 0, 0};                  // Intermediate double cross product components (m x m x p)
double mcrossmcrossHeff_int [3] = {0, 0, 0};               // Intermediate double cross product components (m x m x Heff)
double mcrossp [3] = {0, 0, 0};                            // Cross product components (m x p)
double mcrossHeff [3] = {0, 0, 0};                         // Cross product components (m x Heff)
double mcrossHth [3] = {0, 0, 0};                          // Cross product components (m x Hth)
double mcrossmcrossp [3] = {0, 0, 0};                      // Cross product components (m x m x p)
double mcrossmcrossHeff [3] = {0, 0, 0};                   // Cross product components (m x m x Heff)
double randomHth [3] = {0, 0, 0};                          // Vector of random variables
double STT  = 0;                                        // Strenght of STT term
double FLT  = 0;                                        // Strenght of FLT term

// -------------------------------------------
// Initialize Variables
// -------------------------------------------
double m [3] = {0, 0, 1};                             // Normalized mangetization
double R  = Rap;                						// MTJ resistance [Ohms]
if(initial_state != 1){
                       
    R  = Rp;                                    // MTJ resistance [Ohms]
    m[2]  = -1;                             // Normalized mangetization
}

double J = 0;                                          // Current density [A/m^2]
double V = 0;                                          // MTJ Voltage [V]



for(int i=1;i<=n_sim;i++){

    // Update values
//    double t_current = i*t_step;
    double m_old [3] = {m[0], m[1], m[2]};
    
    
    // Update voltage/current density
//    if ( (t_current>=t_delay) && (t_current<=(t_delay+t_pulse))){
    V = g_v_pulse[i-1] + 1e-10 * (this_id * trials_p_thread+i_trial);

    // Update effective magnetic field Heff_old
    Heff_old[0] = Ext[0]-Ms*Nx*m_old[0];
    Heff_old[1] = Ext[1]-Ms*Ny*m_old[1];
    Heff_old[2] = Ext[2]-Ms*Nz*m_old[2]+(KiPF*m_old[2]-VCMAPF*m_old[2]*V);
    
    //Calculate STT factor
    J = V/(R*areamtj);
    costheta = m_old[0]*P[0] + m_old[1]*P[1] + m_old[2]*P[2];
    g_tunnel = 1/2 * P_tunnel / ( 1 + pow(P_tunnel,2)*costheta);
    g_sv = 1 / ( -4 + pow(( 1 / sqrt(Pol) + sqrt(Pol) ), 3) * (3 + costheta) / 4);
    STT = gammap*J* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0); // STT model, polarization factor is (g_tunnel+g_sv)
    FLT = STT*B1+STT*B2*areamtj*J;

    // Calculate m x Hth
    mcrossHth_int[0]=m_old[1]*randomHth[2]-m_old[2]*randomHth[1];
    mcrossHth_int[1]=m_old[2]*randomHth[0]-m_old[0]*randomHth[2];
    mcrossHth_int[2]=m_old[0]*randomHth[1]-m_old[1]*randomHth[0];

    // Calculate m x p and m x m x p
    mcrossp_int[0]=m_old[1]*P[2]-m_old[2]*P[1];
    mcrossp_int[1]=m_old[2]*P[0]-m_old[0]*P[2];
    mcrossp_int[2]=m_old[0]*P[1]-m_old[1]*P[0];
    mcrossmcrossp_int[0]=m_old[1]*mcrossp_int[2]-m_old[2]*mcrossp_int[1];
    mcrossmcrossp_int[1]=m_old[2]*mcrossp_int[0]-m_old[0]*mcrossp_int[2];
    mcrossmcrossp_int[2]=m_old[0]*mcrossp_int[1]-m_old[1]*mcrossp_int[0];

    // Calculate m x Heff and m x m x Heff
    mcrossHeff_int[0]=m_old[1]*Heff_old[2]-m_old[2]*Heff_old[1];
    mcrossHeff_int[1]=m_old[2]*Heff_old[0]-m_old[0]*Heff_old[2];
    mcrossHeff_int[2]=m_old[0]*Heff_old[1]-m_old[1]*Heff_old[0];
    mcrossmcrossHeff_int[0]=m_old[1]*mcrossHeff_int[2]-m_old[2]*mcrossHeff_int[1];
    mcrossmcrossHeff_int[1]=m_old[2]*mcrossHeff_int[0]-m_old[0]*mcrossHeff_int[2];
    mcrossmcrossHeff_int[2]=m_old[0]*mcrossHeff_int[1]-m_old[1]*mcrossHeff_int[0];
    // Use the LLG equation w/ Heun's Method to update the magnetization
    dm_int[0] = -gammap*(mcrossHeff_int[0]+mcrossHth_int[0]) - gammap*alphac*mcrossmcrossHeff_int[0] + STT*mcrossmcrossp_int[0] + FLT*mcrossp_int[0];
    dm_int[1] = -gammap*(mcrossHeff_int[1]+mcrossHth_int[1]) - gammap*alphac*mcrossmcrossHeff_int[1] + STT*mcrossmcrossp_int[1] + FLT*mcrossp_int[1];
    dm_int[2] = -gammap*(mcrossHeff_int[2]+mcrossHth_int[2]) - gammap*alphac*mcrossmcrossHeff_int[2] + STT*mcrossmcrossp_int[2] + FLT*mcrossp_int[2];
    M_int[0] = m_old[0] + (dm_int[0]*t_step);
    M_int[1] = m_old[1] + (dm_int[1]*t_step);
    M_int[2] = m_old[2] + (dm_int[2]*t_step);
    m_int[0] = M_int[0]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[1] = M_int[1]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[2] = M_int[2]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);

    // Update the thermal field and current values (time evolves)
    
	double2 gen_x12;
	double gen_x3;
	gen_x12 = curand_normal2_double(&localState);
      randomHth[0] = Hth*gen_x12.x;
      randomHth[1] = Hth*gen_x12.y;
	gen_x3 = curand_normal_double(&localState);
      randomHth[2] = Hth*gen_x3;
    
//STT calculation
    costheta = m_int[0]*P[0] + m_int[1]*P[1] + m_int[2]*P[2];
    g_tunnel = 1/2 * P_tunnel / ( 1 + pow(P_tunnel,2)*costheta);
    g_sv = 1 / ( -4 + pow(( 1 / sqrt(Pol) + sqrt(Pol) ), 3) * (3 + costheta) / 4); 
    STT = gammap*J* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0);                  // STT model, polarization factor is (g_tunnel+g_sv)
    //STT = gammap*J/Jc0;
    FLT = STT*B1+STT*B2*areamtj*J;

    // Update intermediate effective magnetic field Heff
    double Heff_int [3] = {Ext[0]-Ms*Nx*m_int[0], Ext[1]-Ms*Ny*m_int[1], Ext[2]-Ms*Nz*m_int[2]+(KiPF*m_int[2]-VCMAPF*m_int[2]*V)};

    // Calculate m x Hth
    mcrossHth[0]=m_int[1]*randomHth[2]-m_int[2]*randomHth[1];
    mcrossHth[1]=m_int[2]*randomHth[0]-m_int[0]*randomHth[2];
    mcrossHth[2]=m_int[0]*randomHth[1]-m_int[1]*randomHth[0];
    // Calculate m x p and m x m x p
    mcrossp[0]=m_int[1]*P[2]-m_int[2]*P[1];
    mcrossp[1]=m_int[2]*P[0]-m_int[0]*P[2];
    mcrossp[2]=m_int[0]*P[1]-m_int[1]*P[0];
    mcrossmcrossp[0]=m_int[1]*mcrossp[2]-m_int[2]*mcrossp[1];
    mcrossmcrossp[1]=m_int[2]*mcrossp[0]-m_int[0]*mcrossp[2];
    mcrossmcrossp[2]=m_int[0]*mcrossp[1]-m_int[1]*mcrossp[0];

    // Calculate m x Heff and m x m x Heff
    mcrossHeff[0]=m_int[1]*Heff_int[2]-m_int[2]*Heff_int[1];
    mcrossHeff[1]=m_int[2]*Heff_int[0]-m_int[0]*Heff_int[2];
    mcrossHeff[2]=m_int[0]*Heff_int[1]-m_int[1]*Heff_int[0];
    mcrossmcrossHeff[0]=m_int[1]*mcrossHeff[2]-m_int[2]*mcrossHeff[1];
    mcrossmcrossHeff[1]=m_int[2]*mcrossHeff[0]-m_int[0]*mcrossHeff[2];
    mcrossmcrossHeff[2]=m_int[0]*mcrossHeff[1]-m_int[1]*mcrossHeff[0];

    // Now use intermediate value in final value computation 
    dm[0] = -gammap*(mcrossHeff[0]+mcrossHth[0]) - gammap*alphac*mcrossmcrossHeff[0] + STT*mcrossmcrossp[0] + FLT*mcrossp[0];
    dm[1] = -gammap*(mcrossHeff[1]+mcrossHth[1]) - gammap*alphac*mcrossmcrossHeff[1] + STT*mcrossmcrossp[1] + FLT*mcrossp[1];
    dm[2] = -gammap*(mcrossHeff[2]+mcrossHth[2]) - gammap*alphac*mcrossmcrossHeff[2] + STT*mcrossmcrossp[2] + FLT*mcrossp[2];
    M[0] = m_old[0] + (t_step/2)*(dm[0] + dm_int[0]);
    M[1] = m_old[1] + (t_step/2)*(dm[1] + dm_int[1]);
    M[2] = m_old[2] + (t_step/2)*(dm[2] + dm_int[2]);
    m[0] = M[0]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[1] = M[1]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[2] = M[2]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    // Update final values for next step
    R = 1/(Gt*(1+(TMR/(TMR+2))*(m[0]*P[0]+m[1]*P[1]+m[2]*P[2]))+Gsi);

}
   	if( initial_state ==0){
   		if( R >= Rp*(1+TMR/2)){
   			finalR[this_id]++;
   		}
   	}
   	else {
   		if( R <= Rp*(1+TMR/2)){
   			finalR[this_id]++;
   		}
   	}
}
                  
}


