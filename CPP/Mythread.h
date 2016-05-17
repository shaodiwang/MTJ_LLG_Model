#include "stdlib.h"
#include <pthread.h>
#include "./pthread_library/thread.h"

using namespace Threads;
extern double t_sim, t_step;
#ifndef _Mythread
#define _Mythread
class Mythread : public Thread { // Derived from Thread class
    public:
    Mythread( double _Vmtj, double _t_pulse, int _initialstate, int _n_trials, int _this_id, int stack_size, bool joinable ) : 
    			Thread(Thread::Options().set_stack_size(stack_size).set_joinable(joinable)), 
    			base_Vmtj(_Vmtj), t_pulse(_t_pulse), initialstate(_initialstate), n_trials(_n_trials), this_id(_this_id) {
/*    	randomArray = new double* [int(t_sim/t_step)];
    	for (int i=0; i< int(t_sim/t_step); i++){
    		randomArray[i] = new double[3];
    		randomArray[i][0] = _randomArray[i][0];
    		randomArray[i][1] = _randomArray[i][1];
    		randomArray[i][2] = _randomArray[i][2];
    	}*/
    }     
    ~Mythread()  {
/*    	for (int i=0; i< int(t_step/t_sim); i++){
    		delete randomArray[i];
    	}
    	delete[] randomArray;*/
    };
    void Run(); // Redefine Run function in Mythread.cpp
    int getsum(){
    	return sum;
    }
    protected:
    double base_Vmtj, t_pulse;
    int initialstate,n_trials,this_id,sum;
};


#endif
