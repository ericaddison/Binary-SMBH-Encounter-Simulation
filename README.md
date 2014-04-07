Binary-SMBH-Encounter-Simulation
================================

This repository contains a gravitational N-body code based on few-body regularization (Mikkola & Tanikawa 1999) and the CHAIN structure (Mikkola & Aarseth 1993). The integrator consists of Burlish-Stoer extrapolation on top of a leap-frog integrator, and evolves the Newtonian equations of motion for three bodies. Standard 4th order Runge-Kutta with adaptive time stepping is also available.

The code is currently configured to evolve a binary system with center of mass in orbit around a third mass. By default the code will run until either the binary is disrupted or it passes the third mass and reaches a new stable orbital configuration. Default outputs include a file with position information for each timestep, gravitational wave output for each timestep, and a results file with the intial and final parameters.

A sample command line call to the program looks like:

	./BEMRI -ang -bs2 -so -ips 200 -N 5 -fname ../data/out.dat
	
This command runs the simulation with randomized binary orientations using the Burlish-Stoer integrator, suppressing screen output, sampling over 200 evenly spaced initial binary phases per orientation, running 5 total orientations (no position or GW output), and outputting to the file ../data/out.dat.

	
The following command line options are available:

	//	-c					run to completion
	//	-s [value]			use value for seed
	//	-t [value]			run to tmax = value*Ph
	//	-nd					do not stop sim when BEMRI is disrupted
	//	-N [value]			run simulation N times, do not output position or quadrupole data
	//	-so					suppress any screen output
  	//  -PbN [value]		set Pb = value*tnode
  	//	-th0 [value]		set theta0 = value
  	//	-bs2				use BS integrator instead of RK4
  	//	-fname [string]		use string as file name instead of BEMRI.dat
  	//	-fpars [string] 	use string as filename to find input parameters
  	//	-beta [value]		use beta = value
	//	-gam [value]		use gamma = value
  	//	-ips [value]		run [value] runs per set of parameters, sampling initial phases
  	//	-eh [value]			set BH orbit eccentricity
	//	-Htest				Test Heggie values
	//	-Hrat [value]		Use [value] for rp/a ratio in Heggie test
	//	-geo				Use geometricized units -- G = c = 1
	//	-ang				Randomize binary orientation angles
