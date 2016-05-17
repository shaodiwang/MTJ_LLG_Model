
#define shaodi_pi 3.1415926
void ellipke(const double kf, double& K, double&E){
    K = 0;
    E = 0;
    const double step = 0.00000025;
    for (double t = 0; t< 1; t+=step){
        K += step / sqrt((1- t*t)*(1-kf*kf*t*t));
        E += sqrt( (1-kf*kf*t*t)/(1-t*t) ) * step;
    }
//    std::cout<<K<<" "<<E<<std::endl;
}


void Demagnetization_factors( const double length, const double width, const double thickness, double& Nx, double &Ny, double& Nz){
    double    a = length/2;
    double    b = width/2;
    double    t = thickness;
    double    d = t/2;
    double    Beta = b/a;
    double    Epsilon = sqrt(1-Beta*Beta);
    double    Tao = d/a;
    double    Zita = d/b;
    double    KTao = 1/sqrt(1+Tao*Tao);
    double    KZita = 1/sqrt(1+Zita*Zita);
    double K = 0, E = 0;
    ellipke(KZita*KZita,K,E);

    // Calculate Nz

    double u0 = 1+4*(1-((1-Zita*Zita)*E+Zita*Zita*K)/KZita)/(3*shaodi_pi*Zita);
    double u1 = 1*(1-((1+2*Zita*Zita)*E-2*Zita*Zita*K)/KZita)/(3*shaodi_pi*Zita);
    double u2 = 3*(1-(E-KZita*KZita*Zita*Zita*K)/KZita)/(16*shaodi_pi*Zita);
    double u3 = 5*(5-KZita*((5+4*Zita*Zita)*E-4*Zita*Zita*K))/(192*shaodi_pi*Zita);
    double u4 = 35*(35-pow(KZita,3)*((5+3*Zita*Zita)*(7+8*Zita*Zita)*E-Zita*Zita*(25+24*Zita*Zita)*K))/(12288*shaodi_pi*Zita);
    double u5 = 21*(315-pow(KZita,5)*((315+857*Zita*Zita+742*pow(Zita,4)+196*pow(Zita,6))*E-2*Zita*Zita*(105+203*Zita*Zita+96*pow(Zita,4))*K))/(81920*shaodi_pi*Zita);
//    std::cout<< u0 << " "<<u1 <<" "<<u2 <<" "<<u3 <<" "<< u4 << " "<< u5 <<" "<< Epsilon<<std::endl;
    double Eps0_u = u0; //pow(Epsilon,0);
    double Eps1_u = u1*pow(Epsilon,2);
    double Eps2_u = u2*pow(Epsilon,4);
    double Eps3_u = u3*pow(Epsilon,6);
    double Eps4_u = u4*pow(Epsilon,8);
    double Eps5_u = u5*pow(Epsilon,10);

    Nz = Eps0_u + Eps1_u + Eps2_u + Eps3_u + Eps4_u + Eps5_u;

// Calculate Nx

    double v0 = -2*(1-((1-Zita*Zita)*E+Zita*Zita*K)/KZita)/(3*shaodi_pi*Zita);
    double v1 = -1*(1-((1+8*Zita*Zita)*E-8*Zita*Zita*K)/KZita)/(12*shaodi_pi*Zita);
    double v2 = -1*(1-(E-5*KZita*KZita*Zita*Zita*K)/KZita)/(32*shaodi_pi*Zita);
    double v3 = -5*(5-KZita*((5-2*Zita*Zita)*E-22*Zita*Zita*K))/(1536*shaodi_pi*Zita);
    double v4 = -7*(35-pow(KZita,3)*((35-3*Zita*Zita-56*pow(Zita,4))*E-Zita*Zita*(145+136*Zita*Zita)*K))/(24576*shaodi_pi*Zita);
    double v5 = -7*(315-pow(KZita,5)*((315+147*Zita*Zita-928*pow(Zita,4)-848*pow(Zita,6))*E-4*Zita*Zita*(315+594*Zita*Zita+268*pow(Zita,4))*K))/(327680*shaodi_pi*Zita);

    double Eps0_v = v0;//Epsilon^0;
    double Eps1_v = v1*pow(Epsilon,2);
    double Eps2_v = v2*pow(Epsilon,4);
    double Eps3_v = v3*pow(Epsilon,6);
    double Eps4_v = v4*pow(Epsilon,8);
    double Eps5_v = v5*pow(Epsilon,10);

    Nx = Eps0_v + Eps1_v + Eps2_v + Eps3_v + Eps4_v + Eps5_v;

// Calculate Ny

    Ny = 1 - Nx - Nz;
//    std::cout<<Nx<<" "<<Ny<<" "<<Nz<<std::endl;
}


