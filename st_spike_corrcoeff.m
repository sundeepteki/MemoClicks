function [rho,p] = st_spike_corrcoeff(spktr1,spktr2,time_window,tON,T)
%{
  Function to compute cross-correlation of spike trains.
  Trial-by-trial across simultaneously recorded pairs of neurons.

  usage: [r,p] = st_corrcoeff(spktr1,spktr2,time_window,dur)
  spktr1: contains spike times in seconds
  time_window: 50ms; can try different time scales
  dur: vary according to pre-stim, stim, rep1/2/3 durations (in sec)

  Based on Martin's tutorial at CAMP NCBS 2015.

  Sundeep Teki.
  Created:           09.07.15
  Last modified:     09.07.15
  Backup:            st_spike_corrcoeff_copy
  Last backup saved: 09.07.15

%}

%%

dt = 0.0001; % sec    
t  = tON:dt:T; 

% time window is implemented as a boxkernel of length time_window/dt
boxkernel = ones(1,ceil(length(time_window/dt))); 

% create arrays that contain 1s for spiketimes, and 0s otherwise
hist1 = zeros(1,length(t)); 
tmp1  =  ceil(spktr1/dt);
hist1(tmp1(tmp1>0))=1;
% tmp1 = tmp1(tmp1<T/dt);
% hist1(tmp1>0) = 1;

hist2 = zeros(1,length(t));
tmp2  =  ceil(spktr2/dt);
hist2(tmp2(tmp2>0))=1;
% tmp2 = tmp2(tmp2<T/dt);
% hist2(tmp2>0) = 1;

% convolution with the boxkernel is an efficient way to count the events in the interval time_window
n1 = conv(hist1,boxkernel,'same'); 
n2 = conv(hist2,boxkernel,'same');

% use a built-in numpy function to compute the correlation coefficient between the count sequences:
[rho,p] = corrcoef(n1,n2);
rho = rho(1,2);
p   = p(1,2);

