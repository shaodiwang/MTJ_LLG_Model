
To compile
	change CCFLAG to include your own cuda directory
	type make to compile


To run WERsim:
	./WER  [total write trials] [number of blocks (32n) in GPU] [number of threads in block (32m) in GPU]  
	[MTJ initial state: initial state 0 = parallel, 1 = anti-parallel] 
	[input pulse shape file, e.g., example_pulse.txt, every line contains two parameters: time, voltage on MTJ ] 
	[Enable VCMA effect?, 1:enable, 0: disable]
	
	e.g. ./WER 1048576 128 256 0 ./example_pulse.txt 1  
	1048576 is the number of Monte-Carlo simulations, 128 is the #. of utilized blocks 
	in GPU, 256 is the #. of threads per block, 0 means initial state is parallel (low resistance),
	and example_pulse.txt is the write pulse in the formats [ time (second)  voltage (V)] for each
	line. The last 1 means VCMA effect is enabled
	
	For STT effect:	p to ap requires positive STT current (positive voltage), ap to p
	requires negative STT current. 
	For VCMA effect: positive voltage reduces coersivity (thermal 
	stability), and negative voltage increases coersivity.

Output files:
	ap2p.txt or p2ap.txt

To modify MTJ dimension and variation:
	WERSim.cu: lines 113-118
To modify MTJ's external field:
	LLG.cu: line 78
To modify MTJ's STT polarization efficiency:
	LLG.cu: lines 166,214
To modify MTJ's VCMA factor:
	LLG.cu: line 82, 85
To modify MTJ resistance;
	LLG.cu: line 36
	
