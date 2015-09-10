function [x,y,z] = st_getMMCclicktimes(varargin)

%{
  Function to get correct clicktimes by using seed
  
Sundeep Teki
Created            07.07.15
Last modified:     07.07.15
Backup:            st_getMMCclicktimes_copy
Last backup saved: 0x.07.15

%}

stim     = {};
events   = {};
trialnum = 100; %input('Enter no. of trials: ');
o        = MemoClicks();
o.maxgap = MMCsorted.MMCinfo.MaxGap;
o.Key    = MMCsorted.MMCinfo.RefRCSeed;

for i = 1:trialnum
  [w, ev,o] = waveform (o,1,[],[],i);
  stim{i}   = w;
  events{i} = ev;
end

x=get(o,'stimulus');
xC=(find(x==0));
xRC=(find(x==1));
xRefRC=(find(x==2));

y=get(o,'Seeds');

z = get(o,'ClickTimes');
