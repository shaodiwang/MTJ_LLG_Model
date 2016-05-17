// -------------------------------------------
// Updated May 6 - 2016
// Vector coordinates: [x; y; z]
// -------------------------------------------
#include <iostream>
#include <math.h>
#include "LLG.cu"
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;
// -------------------------------------------
// Calculate Wall time
// -------------------------------------------
#include <sys/time.h>
#include "./Demagnetization_factors.cu"
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}



// -------------------------------------------
// WER Calculation
// -------------------------------------------
static cudaError_t crc;

void g_allocate_1D(double **g_f, int nsize, int *irc) {
/* allocate global double memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(double)*nsize);
   if (crc) {
      printf("cudaMalloc double Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_f = (double *)gptr;
   return;
}

void g_allocate_2D(double ***g_f, size_t * pitch, int row, int column, int *irc) {
/* allocate global double memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMallocPitch( &gptr, pitch, sizeof(double)*column, row);
   if (crc) {
      printf("cudaMalloc double Error=%d:%s,row=%d,column=%d\n",crc,
              cudaGetErrorString(crc),row,column);
      *irc = 1;
   }
   *g_f = (double **)gptr;
   return;
}


void copyin_gmemptr_1D (double *f, double *g_f, int nsize) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(double)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice1D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyin_gmemptr_2D (double **f, double **g_f, size_t& pitch, int row, int column) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy2D(g_f,pitch,f,sizeof(double)*column,sizeof(double)*column, row, cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice2D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyout_gmemptr_1D(double *f, double *g_f, int nsize) {
/* copy double array from global GPU memory to main memory */
   crc = cudaMemcpy(f,g_f,sizeof(double)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost1D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyout_gmemptr_2D(double **f, double **g_f, size_t &pitch, int row, int column) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy2D(f,sizeof(double)*column,g_f,pitch,sizeof(double)*column, row, cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost2D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

int main(int    argc, char* argv[])
{
	if(argc < 7) {
		cout<<" arguments: [total write trials] [number of blocks (32n) in GPU] [number of threads in block (32m) in GPU]  [MTJ initial state: initial state 0 = parallel, 1 = anti-parallel] [input pulse shape file, e.g., example_pulse.txt, every line contains two parameters: time, voltage on MTJ ] [Enable VCMA effect?, 1:enable, 0: disable]"<<endl;
	return 1;
	}
	//dimention 
        double length = 40e-9;                            // MTJ length [m]
	double width = 40e-9;                            // MTJ width [m]
	double tfl = 1.2e-9;                              // Free layer thickness [m]
	double sigma_l = 2e-9;				 //standard deviation of length in process variation
	double sigma_w = 2e-9;				 //std of w in process variation
	double sigma_tfl = 1e-10;			// std of free layer thickness in process variation
	//Demagnetization calculation
	double start_t = get_wall_time(); // Recording computation time
	std::cout<<"Demagnitization calculation"<<std::endl;
	double Nx = 0, Ny = 0, Nz = 0;			//Demagnetization factors
	double* lin_dep_factor = new double[9];		//The dependence of demagnetization factors on MTJ dimension
	Extract_linear_dependent(length, sigma_l,  width, sigma_w, tfl, sigma_tfl, Nx, Ny, Nz,lin_dep_factor);	
	
	double t_step = 5e-12;                             // Time step [s]
	double t_sim = 3e-9;                              // preset simulation time [s]
	
	int trials = atoi(argv[1]);			//input trials
	int GridSize = atoi(argv[2]);
	int GridLength = sqrt(GridSize);
	int BlockSize = atoi(argv[3]);
	string input_pulse_file(argv[5]);
	bool enableVCMA=true; // is it precessional switching or STT
	if ( atoi(argv[6]) == 0){
		 enableVCMA = false;
	}
	int BlockLength = sqrt(BlockSize);
	int trials_p_thread = ceil(double(trials)/ double(GridSize*BlockSize));
	int real_trials = trials_p_thread * GridSize*BlockSize; //Real trials in gpu simulations
	dim3 dimBlock(BlockSize, 1);
	dim3 dimGrid(GridSize,1);
	double* global_finalR = new double[GridSize * BlockSize];
	int initial_state = atoi(argv[4]);			//initial_state of MTJ
	//Read pulse shape and open output file
	fstream fs,fi;
	fi.open(argv[5],std::fstream::in);
	vector<double> v_time, v_voltage;
	double readBuffer,voltage=0, pulse_start_t=0, pulse_end_t=0;
	
	double Demag_t = get_wall_time(); // Recording computation time
	cout<<"Read in pulse shape file: "<<argv[5]<<endl;
	while(fi >> readBuffer){
		v_time.push_back(readBuffer);
		if( readBuffer > 1e-7){
			cout<<"Error: simulation time > 100ns"<<endl;
			return 1;
		}
		if (readBuffer > t_sim){ // Pulse time is longer than pre-set simulation time
			t_sim = readBuffer;
		}
		fi >> readBuffer;
		v_voltage.push_back(readBuffer);
		if(readBuffer > voltage) voltage = readBuffer;
	}
	fi.close();
	if(v_time.size() != v_voltage.size()){
		cout<<"Error: number of times > number of voltages"<<endl;
		return 1;
	}
	int n_sim = ceil(t_sim/t_step); // total simulation steps in a simulation
	//Move pulse shape to 1D array to copy to GPU
	double *v_pulse;
	v_pulse = new double[n_sim];
	int p_v=-1;
	double current_time=0;
	//Characterize pulse votlageinto every time steps
	for (int i_v = 0; i_v < n_sim; i_v ++ ){
		current_time = (i_v+1) * t_step;
		while( (v_time[p_v+1] < current_time) && (p_v < (int)(v_time.size()-1)) ) p_v++;
		if( p_v < (int) (v_time.size() -1)  && p_v >= 0 ){
			v_pulse[i_v] = v_voltage[p_v] + (v_voltage[p_v+1] - v_voltage[p_v])/(v_time[p_v+1]-v_time[p_v]) * ( current_time - v_time[p_v]) ;
		}
		else{
			v_pulse[i_v] = 0;
		}
		if ( v_pulse[i_v] > voltage*0.9 && pulse_start_t == 0)
			pulse_start_t = current_time;
		if ( v_pulse[i_v] < voltage*0.9 && pulse_start_t > 0 && pulse_end_t == 0)
			pulse_end_t = current_time;
	}
	double pulse_width = pulse_end_t - pulse_start_t;
	v_time.clear();
	v_voltage.clear();
	vector<double>().swap(v_time);
	vector<double>().swap(v_voltage);	
	if( initial_state == 1){
		fs.open("ap2p.txt",std::fstream::out | std::fstream::app);
	}
	else{
		fs.open("p2ap.txt", std::fstream::out | std::fstream::app);
	}
	
	double readin_t = get_wall_time();
	
	cout<<"Copy data from host memory to GPU..."<<endl;
	double * g_lin_dep_factor;
	int irc = 0;
	g_allocate_1D(&g_lin_dep_factor, 9 ,&irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D(  lin_dep_factor , g_lin_dep_factor, 9);

	double * g_finalR, *g_v_pulse;
	g_allocate_1D(&g_finalR,GridSize * BlockSize ,&irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( global_finalR ,g_finalR, GridSize * BlockSize);
	
	g_allocate_1D(&g_v_pulse, n_sim, &irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( v_pulse ,g_v_pulse, n_sim);

	cout<<"GPU calculation starts...";
	double copyin_t = get_wall_time();
   	LLG<<<dimGrid,dimBlock>>>(g_v_pulse, g_finalR, initial_state, t_step, trials_p_thread, n_sim, enableVCMA, length, sigma_l,  width, sigma_w, tfl, sigma_tfl, Nx, Ny, Nz, g_lin_dep_factor );
	cudaDeviceSynchronize();
	double calculate_t = get_wall_time();
	cout<<"ends"<<endl;
	
	cout<<"Copy GPU to host memory"<<endl;
	double copyout_t = get_wall_time();
	copyout_gmemptr_1D(global_finalR, g_finalR, GridSize * BlockSize);


	int sum=0; // total switched MTJ numbers
   	for( int i =0 ; i<GridSize * BlockSize; i++){
   		sum+= global_finalR[i];
   	}
	double couting_t = get_wall_time();

	cout<<"Complete calculation\n"<<endl;
   	cout << "**********************\nMonte-Carlo Simulation Results:\n  switching from "<< ((initial_state==0)? "P" : "AP")<<" to "<< ((initial_state==1)? "P" : "AP")<<"\n"<<"  Pulse voltage: "<<voltage<<" V, Pulse width: "<<pulse_width<<"s\n  Simulation time: "<<t_sim<<"s, time step: "<<t_step<<"s\n"
   	<<"  total number of trials: "<<real_trials<<", "<<sum<<" trials success,"<<" switching rate: "<<std::setprecision(9)<<(double(sum) / double(real_trials)) << endl;
	cout<<"Runtime summary:\n"<<"  Demagetization calculation: "<< Demag_t - start_t<<"s\n  Readin pulse shape: "<< readin_t - Demag_t << " s\n  copy memory to GPU: "<<copyin_t - readin_t<<" s\n  GPU calculation: "<<calculate_t - copyin_t<<" s\n  copy memory out to CPU: "<<copyout_t - calculate_t <<" s\n***********************\n"<<endl;
	fs << voltage <<" "<<pulse_width<<" "<<std::setprecision(9)<<(double(sum) / double(real_trials))<<" "<<couting_t-start_t<<" s"<<endl;
	fs.close();
	return 0;
}





