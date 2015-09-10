function [rho] = st_cross_coef(spktr1, spktr2,time_window, T);


%{


	dt = 0.1			% time resolution [ms]
	t = 0:dt:1000*T+dt
    boxkernel = ones(1,ceil(length(time_win/dt)))
    

	t = np.arange(0,1000*T+dt,dt)	
	# time window is implemented as a boxkernel of length time_window/dt
	boxkernel = np.ones(np.ceil(time_window/dt))
	# create arrays that contain 1s for spiketimes, and 0s otherwise 	
	hist1 = np.zeros(len(t))
	hist1[(spktr1*1000/dt).astype(int)]=1
	hist2 = np.zeros(len(t))
	hist2[(spktr2*1000/dt).astype(int)]=1

	# convolution with the boxkernel is an efficient way to count the events in the 
	# interval time_window 
	n1 = np.convolve(hist1,boxkernel,'same') %% conv
	n2 = np.convolve(hist2,boxkernel,'same')

	# use a built-in numpy function to compute the correlation coefficient between the 
	# count sequences:
	rho = np.corrcoef(n1,n2)[0,1]
	return rho

%}


dt = 0.1;

