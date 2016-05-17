
To compile
	download boost and extract under the same directory
	type make to compile


To run WERsim:
	./main [number of simulations] [number of threads] [pulse voltage e.g. 0.9V] [pulse width in ns ] 
	[initial state 0 = parallel, 1 = anti-parallel]

	e.g. ./main 1000 10 0.9 1 
	1000 is the #. of Monte-Carlo simulations, 10 is the #. of threads
	0.9 is the pulse width (0.9ns), 1 means initial state is anti-parallel (high resistance),
	
	For STT effect:	p to ap requires positive STT current (positive voltage), ap to p
	requires negative STT current. 
	For VCMA effect: positive voltage reduces coersivity (thermal 
	stability), and negative voltage increases coersivity.

Output:
	see terminal

To modify MTJ dimension, external file, MTJ resistance
	Mythread.cpp: lines 43-51
To modify MTJ's STT polarization efficiency:
	Mythread.cpp: lines 73,73
To modify MTJ's VCMA factor:
	Mythread.cpp: line 70
	
