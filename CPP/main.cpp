// -------------------------------------------
// Precessional Switching Old Parameters
// Updated November 21 - 2013
// Vector coordinates: [x; y; z]
// -------------------------------------------
#include <iostream>
//#include <math.h>
#include <vector>
#include "Mythread.h"
//#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;

#include <sys/time.h>
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

//std::vector<int> CountRfinal;
                            // Total simulation time [s]

double t_step = 7.5e-13;                             // Time step [s]
double t_sim = 25e-9;

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

double * global_finalR;

int main(int    argc, char* argv[])
{
  
	if(argc < 6) {
		cout<<" arguments: [number of trials] [number of threads] [pulse voltage e.g. 0.9V] [pulse width in ns eg 1] [initial_state 1:AP or 0:P]"<<endl;
		return 1;
	}
    double start_t = get_wall_time();
	int trials = atoi(argv[1]);
	int n_threads = atoi(argv[2]);
	double pulse_t = atof(argv[4]) * 1e-9;
	double voltage = atof(argv[3]);
	int initial_state = atoi(argv[5]);

//Setup random variable

	

	double Rp = 100e3;
	double TMR = 1;

	vector<Mythread*> threads;
	threads.resize(n_threads);
	int sum=0;
//	double randomArray [int(t_sim/t_step)][3];
	
		for (int inner_l=0; inner_l<n_threads; inner_l ++){
//			generator=gen_normal(generator, randomArray, int(t_sim/t_step) );
			Mythread* t_mythread = new Mythread( voltage, pulse_t, initial_state, int(trials/n_threads),inner_l, 1024*1024, true);
			threads[inner_l] =  t_mythread;
			t_mythread->Start();
		}
		for (int inner_l=0; inner_l<n_threads; inner_l ++){
			threads[inner_l]->Join();
		}
		for( int i =0 ; i<n_threads; i++){

			sum += threads[i]->getsum();
		}
		for(int inner_l = 0; inner_l<n_threads; inner_l++){
			delete threads[inner_l];
		}
		
	  	
    double end_t = get_wall_time();
   	cout << "Monte-Carlo Results for "<<initial_state<<" to "<< (initial_state == 0)<<" :\n"<<"Pulse voltage: "
   	<<voltage<<" V, Pulse Width: "<<pulse_t<<"\nnumber of threads: "<<n_threads<<" , total number of trials: "<<trials
   	<<"\n"<<sum<<" success of "<<trials<<" trials, percentage:"<<(double(sum) / double(int(trials/n_threads)*n_threads)) <<"\n"<<"Run time: "<<end_t - start_t << endl;
   	return 0;
}





