#include <iostream>
#include <stdio.h>      
#include <math.h>    
#include <stdlib.h>
#include "./Mythread.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include "./Demagnetization_factors.cpp"
#include <fstream>

using namespace std;
//extern double t_sim, t_step;
//#define shaodi_pi 3.1415926
//#define OUTCURVE
//#define DEBUG

#ifdef OUTCURVE
std::fstream fout ("curve.log", std::fstream::out);
#endif
template<class T>
T gen_normal(T generator,
            double res[][3], int row)
{
  for(int i=0; i<row; ++i){
	  for (int j=0; j<3;j++ ){
		  res[i][j]=generator();
	  }
  }
  // Note the generator is returned back
  return  generator;
}

void Mythread::Run(){
/* -------------------------------------------
 Input Parameters From User
 -------------------------------------------*/

boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    	generator(boost::mt19937( static_cast<unsigned int>(this_id)),
    			boost::normal_distribution<>());
 
double length = 60e-9;                            // MTJ length [m]
double width = 60e-9;                              // MTJ width [m]
//double dMgO = 1.3e-9;                              // MgO thickness [m]
double tfl = 1.22e-9;                               // Free layer thickness [m]
double Rp = 100e3;                                 // Parallel resistance [Ohms]
double TMR = 1;                                    // TMR at zero bias voltage
double T = 300;                                    // Temperature [K]
double Ext [3] = {-11000, 0, 0};                        // External magnetic field [A/m] - 1 oersted [Oe] = 79.5774715459424 ampere/meter [A/m]
double B1 = 0.2;                                   // Field-like torque linear parameter [unitless]
double B2 = 0;                                     // Field-like torque quadratic parameter [1/A]
int initial_state = initialstate;                          // Inital state [0 = parallel, 1 = anti-parallel]
double P [3] = {0, 0, -1};                             // Direction of polarization

double t_delay = 2e-9;                             // Time to initiate pulse application [s]

/* -------------------------------------------
 Constants
 -------------------------------------------*/

double hbar = 1.05457173e-34;                                      // Reduced Planck constant, [J*s]
double k = 1.3806488e-23;                                          // Boltzmann constant, [J/K]
double u0 = 4e-7*shaodi_pi;                                               // Vacuum permeability, [V�s/(A�m)]
double q = 1.60217657e-19;                                         // Electron charge, [C]
double alphac = 0.02;                                              // LLGE damping factor
double gammap = (221276/(1+pow(alphac,2)));                             // Gyromagnetic ratio [m/(A x s)]
double Ms0 = 1393128.323;                                                 // Saturation magnetization at 1120K [A/m] - 1e6 A/m = 1000 emu/cc  
double Ms = Ms0*(1- pow(T/1120,1.5));
double Xi = 130.39247e-15* pow(Ms/Ms0, 2.83); //origin:53.39247e-15* pow(Ms/Ms0, 2.83)
double Ki = 1.479036e-3*pow(Ms/Ms0,2.18);//1.0056364e-3;                          // Anisotropy field constant [J/m^2]
double Pol = 0.6;                                                  // Polarization for Spin Torque
double P_tunnel = 0.2;                                  // the polarization of the tunnel currentdouble
double Jc0 = (2*Ms*tfl*q*u0)/(hbar*Pol);                           // Normalization Constant for Current Density
double Nx = 0;                                                // x Demagnetization factor
double Ny = 0;                                                // y Demagnetization factor
double Nz = 0;                                                // z Demagnetization factor
Demagnetization_factors(length,width,tfl,Nx,Ny,Nz);               // Calculate demagnetization_factor
double areamtj = shaodi_pi*length*width/4;                                // MTJ area [m^2]
double volume = areamtj*tfl;                                       // MTJ volume [m^3]
double dMgO_a = 1.54e-3, dMgO_b = 1.1537278e10;//origin:9.24e9;
double dMgO = (log(Rp * areamtj * 10e12) - log(dMgO_a)) / dMgO_b;     // MgO thickness [m]
double Rap = (1+TMR)*Rp;                                           // Anti-parallel resistance [Ohms]
double Gt = 1/(Rp*(1+(TMR/(TMR+2))));                              // Direct elastic tunneling conductance [S]
double KiPF = (2*Ki)/(tfl*u0*Ms);                                  // Prefactor for interface anisotropy effective field
double VCMAPF = (2*Xi)/(u0*Ms*dMgO*tfl);                           // Prefactor for VCMA effective field
//std::cout<<"VCMAPF: "<<VCMAPF/KiPF<<std::endl;
double Hth = sqrt((2*k*T*alphac)/(u0*gammap*Ms*volume*t_step));    // Amplitude of Thermal Field
//std::cout<<"Hth: "<<Hth<<std::endl;
double Gsi	= 0;                                                    // Conductance due to imperfections in Mgo [S]
//Thermal calculation
#ifdef DEBUG
double Energy = 0;
double mean_Jc0 = 0;
double zeron = 0;
double Ms1 = Ms0*(1 - pow((T+50)/1120,1.5));
double Ki1 = 1.479036e-3*pow(Ms1/Ms0,2.18);
double Hkeff = 0.99612 * 2*Ki/(tfl*Ms)/(u0) - Ms*(Nz - max(Nx,Ny));
double Hkeff1 = 0.99612 * 2*Ki1/(tfl*Ms1)/(u0) - Ms1*(Nz - max(Nx,Ny));
double thermS = Ms*Hkeff*(length*width*shaodi_pi/4*tfl)/(2*k*T)*u0;
double thermS1 = Ms1*Hkeff1*(length*width*shaodi_pi/4*tfl)/(2*k*(T+50))*u0;
std::cout<<"Thermal stability: "<< thermS<<" "<< thermS1 <<" retention time:"<<1e-9*exp(thermS)/3600/24<<" "<<1e-9*exp(thermS1)/3600/24<<std::endl;
std::cout<<"VCMA voltage: "<<(0.99612*KiPF - (Nz-max(Nx,Ny))*Ms)/VCMAPF << std::endl;
double maxJc0=0, maxTheta=0;
#endif

for (int i_trial =0; i_trial<n_trials; i_trial++){
double Vmtj = base_Vmtj + (n_trials*this_id + i_trial)*1e-10;
double randomArray [int(t_sim/t_step)][3];
generator=gen_normal(generator, randomArray, int(t_sim/t_step) );
/* -------------------------------------------
 Internal Variables
 -------------------------------------------*/

double m_old [3] = {0, 0, 0};                              // Normalized previous magnetization
double Heff_old [3] = {0, 0, 0};                           // Previous Heff components [A/m]
double m_int [3] = {0, 0, 0};                              // Intermediate normalized magnetization
double dm_int [3] = {0, 0, 0};                             // Intermediate derivative of magnetization
double M_int [3] = {0, 0, 0};                              // Intermediate denormalized magnetization                              
double Heff_int [3] = {0, 0, 0};                           // Intermediate Heff components [A/m]
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
double costheta = 0;                                       // the angle between the magnization  of free and reference layers
double g_sv = 0;                                        // the polarization efficiency in spin valve
double g_tunnel = 0;                                    // the polarization efficiency in tunnel current

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
double _current = 0;                                  // Current simulation time [s]
//srand((unsigned)time(NULL));                               // Random number generator

// -------------------------------------------
// Uncomment this block if T = 0 or if thermal noise is disabled
// -------------------------------------------

// if (T == 0)
//     delta = 0.001;
//     if(m(3) == 1)
//           m = [delta; delta; 1-2*delta];
//     end
// 
//     if(m(3) == -1)
//           m = [-delta; -delta; -1+2*delta];
//     end
// end

// -------------------------------------------
// Cycle for the simulation time
// -------------------------------------------
 

for(int i=1;i<=int(t_sim/t_step);i++){

    // Update values
    double t_current = i*t_step;
    double m_old [3] = {m[0], m[1], m[2]};
    
    
    // Update voltage/current density
    if ( (t_current>=t_delay) && (t_current<=(t_delay+t_pulse))){
        V = Vmtj;
    }
    else{
        V = 0;
    }

    // Update effective magnetic field Heff_old
    Heff_old[0] = Ext[0]-Ms*Nx*m_old[0];
    Heff_old[1] = Ext[1]-Ms*Ny*m_old[1];
    Heff_old[2] = Ext[2]-Ms*Nz*m_old[2]+(KiPF*m_old[2]-VCMAPF*m_old[2]*V);
#ifdef OUTCURVE
    fout<<i*t_step<<"\t"<<V<<"\t"<<Heff_old[0] << "\t" << m_old[2]<< "\t" << V/R<<std::endl;
#endif
//    std::cout<< Heff_old[2] << " " <<-Ms*Nz*m_old[2]<<" "<<KiPF*m_old[2]<<" "<<-VCMAPF*m_old[2]*V<<" "<<-Ms*Nz*m_old[2]+(KiPF*m_old[2]-VCMAPF*m_old[2]*V)<<std::endl;
//    cout<< Heff_old[2]<<endl;
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
    //cout<<" m_old[0] "<<m_old[0]<<" m_old[1] "<<m_old[1]<<" m_old[2] "<<m_old[2]<<endl;
    // Use the LLG equation w/ Heun's Method to update the magnetization
    dm_int[0] = -gammap*(mcrossHeff_int[0]+mcrossHth_int[0]) - gammap*alphac*mcrossmcrossHeff_int[0] + STT*mcrossmcrossp_int[0] + FLT*mcrossp_int[0];
    dm_int[1] = -gammap*(mcrossHeff_int[1]+mcrossHth_int[1]) - gammap*alphac*mcrossmcrossHeff_int[1] + STT*mcrossmcrossp_int[1] + FLT*mcrossp_int[1];
    dm_int[2] = -gammap*(mcrossHeff_int[2]+mcrossHth_int[2]) - gammap*alphac*mcrossmcrossHeff_int[2] + STT*mcrossmcrossp_int[2] + FLT*mcrossp_int[2];
    //std::cout<<V<<" "<<dm_int[0] << " "<<dm_int[1] << " " <<dm_int[2] << std::endl;
    M_int[0] = m_old[0] + (dm_int[0]*t_step);
    M_int[1] = m_old[1] + (dm_int[1]*t_step);
    M_int[2] = m_old[2] + (dm_int[2]*t_step);
    //cout<<"dm_int[1] "<<dm_int[1]<<" t_step "<<t_step<<endl;
    //cout<<"M_int "<<M_int[0]<<" "<<M_int[1]<<" "<<M_int[2]<<endl;;
    m_int[0] = M_int[0]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[1] = M_int[1]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[2] = M_int[2]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);

    // Update the thermal field and current values (time evolves)
    
    
    randomHth[0] = Hth*randomArray[i-1][0];
    randomHth[1] = Hth*randomArray[i-1][1];
    randomHth[2] = Hth*randomArray[i-1][2];
//    randomHth[0] = Hth*0.5;
    //randomHth[1] = Hth*(-0.7);
    //randomHth[2] = Hth*0.3;
    //if(i>=0)	cout<<" randomHth[0] "<<randomHth[0]<<" randomHth[1] "<<randomHth[1]<<" randomHth[2] "<<randomHth[2]<<endl;
    //Calculate theta, the angle between the magnization of the free and reference layer
    costheta = m_int[0]*P[0] + m_int[1]*P[1] + m_int[2]*P[2];
    g_tunnel = 1/2 * P_tunnel / ( 1 + pow(P_tunnel,2)*costheta);
    g_sv = 1 / ( -4 + pow(( 1 / sqrt(Pol) + sqrt(Pol) ), 3) * (3 + costheta) / 4);
    J = V/(R*areamtj);
   // cout<<t_current<<" "<<V<<" "<<J<<endl;
//    STT = gammap*J/Jc0;
    //STT = gammap*J/Jc0;
    STT = gammap*J* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0);
#ifdef DEBUG
    Energy += V/R * V *t_step;
    double absJc00 =abs((gammap*(mcrossHeff_int[0]+mcrossHth_int[0]) + gammap*alphac*mcrossmcrossHeff_int[0])/  (mcrossmcrossp_int[0] * gammap* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0)) );
    double absJc01 =abs((gammap*(mcrossHeff_int[1]+mcrossHth_int[1]) + gammap*alphac*mcrossmcrossHeff_int[1])/  (mcrossmcrossp_int[1] * gammap* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0)) );
    double absJc02 =abs((gammap*(mcrossHeff_int[2]+mcrossHth_int[2]) + gammap*alphac*mcrossmcrossHeff_int[2])/  (mcrossmcrossp_int[2] * gammap* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0)) );
    double absJc0 = min(absJc01, absJc00);
    absJc0 = min(absJc02, absJc0);
    double absJc1 = abs(Heff_old[2]/ ((g_tunnel+g_sv) * hbar /(2*q*Ms*tfl*u0)));  
    if(absJc1 > maxJc0 && !isnan(absJc0)){
        maxJc0 = absJc1;
    //    std::cout<<"absJc02:" <<absJc02<<std::endl;
        maxTheta = costheta;
        mean_Jc0 += absJc0;
//        cout<<"m_int[0]: "<<m_int[0]<<" m_int[1]: "<<m_int[1]<<" m_int[2]: "<<m_int[2]<<std::endl;
//        cout<<"Heff_int[0]: "<<Heff_old[0]<<" Heff_int[1]: "<<Heff_old[1]<<" Heff_int[2]: "<<Ext[2]-Ms*Nz*m_old[2]<<std::endl;
//        cout<<"mHth_int[0]: "<<mcrossHth_int[0]<<" mHth_int[1]: "<<mcrossHth_int[1]<<" mHth_int[2]: "<<mcrossHth_int[2]<<std::endl;
      //  std::cout<<"mcrossHeff_int[2]: "<<mcrossHeff_int[2]<<" mcrossHth_int[2]: "<<mcrossHth_int[2]<<" alphac*mcrossmcrossHeff_int[2]: "<<alphac*mcrossmcrossHeff_int[2]<<" mcrossmcrossp_int[2]: "<<mcrossmcrossp_int[2]<<" hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0): "<< hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0)<<" Heff_old[2]: "<<Heff_old[2]<<std::endl;

    }
    //if( abs(costheta) < 0.001)
    if( isnan(absJc0))
        zeron = mcrossmcrossp_int[2];
#endif
    FLT = STT*B1+STT*B2*areamtj*J;

    // Update intermediate effective magnetic field Heff
    double Heff_int [3] = {Ext[0]-Ms*Nx*m_int[0], Ext[1]-Ms*Ny*m_int[1], Ext[2]-Ms*Nz*m_int[2]+(KiPF*m_int[2]-VCMAPF*m_int[2]*V)};

    // Calculate m x Hth
    mcrossHth[0]=m_int[1]*randomHth[2]-m_int[2]*randomHth[1];
    mcrossHth[1]=m_int[2]*randomHth[0]-m_int[0]*randomHth[2];
    mcrossHth[2]=m_int[0]*randomHth[1]-m_int[1]*randomHth[0];
    //cout<<" m_int[1] "<<m_int[1]<<" randomHth[2] "<<randomHth[2]<<" m_int[2] "<< m_int[2] << " randomHth[1] "<<randomHth[1]<<endl;
    // Calculate m x p and m x m x p
    mcrossp[0]=m_int[1]*P[2]-m_int[2]*P[1];
    //cout<<"m_int[1] "<<m_int[1]<<" m_int[2] "<<m_int[2]<<" P[2] "<<P[2]<<" P[1] "<<P[1]<<endl;
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
    //cout<<" gammap "<<gammap<<" mcrossHeff[0] "<<mcrossHeff[0]<<" mcrossHth[0] "<<mcrossHth[0]<<" mcrossmcrossHeff[0] "<<mcrossmcrossHeff[0]<<" mcrossmcrossp[0] "<<mcrossmcrossp[0]<<" mcrossp[0] "<<mcrossp[0]<<" STT "<<STT<<" FLT "<<FLT<<endl;
    //cout<< -gammap*(mcrossHeff[0]+mcrossHth[0]) << " "<< gammap*alphac*mcrossmcrossHeff[0]<<" "<<STT*mcrossmcrossp[0]<<" "<<FLT*mcrossp[0]<<endl;
    dm[0] = -gammap*(mcrossHeff[0]+mcrossHth[0]) - gammap*alphac*mcrossmcrossHeff[0] + STT*mcrossmcrossp[0] + FLT*mcrossp[0];
    dm[1] = -gammap*(mcrossHeff[1]+mcrossHth[1]) - gammap*alphac*mcrossmcrossHeff[1] + STT*mcrossmcrossp[1] + FLT*mcrossp[1];
    dm[2] = -gammap*(mcrossHeff[2]+mcrossHth[2]) - gammap*alphac*mcrossmcrossHeff[2] + STT*mcrossmcrossp[2] + FLT*mcrossp[2];
    M[0] = m_old[0] + (t_step/2)*(dm[0] + dm_int[0]);
    M[1] = m_old[1] + (t_step/2)*(dm[1] + dm_int[1]);
    M[2] = m_old[2] + (t_step/2)*(dm[2] + dm_int[2]);
    //cout<<"m_old[0] "<<m_old[0]<<" dm[0] "<<dm[0]<<" dm_int[0] "<<dm_int[0]<<endl;
    m[0] = M[0]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[1] = M[1]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[2] = M[2]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    //cout<<"M[0] "<<M[0]<<" M[1] "<<M[1]<<" M[2] "<<M[2]<<endl;
    // Update final values for next step
    R = 1/(Gt*(1+(TMR/(TMR+2))*(m[0]*P[0]+m[1]*P[1]+m[2]*P[2]))+Gsi);

    //cout<<"Gt "<<Gt<<" m[0] "<<m[0]<<" P[0] "<<P[0]<<" m[1] "<<m[1]<<" P[1] "<<P[1]<<" m[2] "<<m[2]<<" P[2] "<<P[2]<<" Gsi "<< Gsi<<endl;
    //cout<<R<<endl;
}
    
    if( initial_state ==0){

		if( R >= Rp*(1+TMR/2)){
			sum++;
        }
	}
	else {
		 if( R <= Rp*(1+TMR/2)){
		   	sum++;
		 }
    }

    }
#ifdef DEBUG
    std::cout<<"Switching Energy: "<<Energy<<std::endl;
    std::cout<<"Max Jc0: "<<maxJc0*length*width<<" maxTheta: "<<maxTheta<<" mean Jc0: "<<mean_Jc0/t_sim*t_step<<" zero para: "<<zeron<<std::endl;
#endif
#ifdef OUTCURVE
    fout.close();
#endif
}


