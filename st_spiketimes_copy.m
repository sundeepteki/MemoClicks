function [spiketimes,trials] = st_spiketimes(varargin)
clear all; close all; clc;

%{
  Function to get spiketimes (in sec) for each trial sorted by condition
  usage: [spiketimes] = st_spiketimes(MMC)
  
  Sundeep Teki.
  Created:           09.07.15
  Last modified:     09.07.15
  Backup:            st_spiketimes_copy
  Last backup saved: 09.07.15
%}

%%

datapath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/Results/';
savepath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/SpikeTimes/';
addpath(datapath);
addpath(savepath);

%%

tmpdir = dir(datapath);
spikefs     = 31250;

for ij = 1:(length(tmpdir)-3)
    
    data{ij} =  [tmpdir(ij+3).name];
    load(data{ij});
    fname = data{ij}; fname = fname(1:end-2);
    disp(['Loading file - ' fname]);
    
    for i = 1:100
        spiketimes.all{i} = MMCsorted.spiketimes(MMCsorted.trials == i)/spikefs;
    end
    
    spiketimes.C     = spiketimes.all(MMCsorted.Ctrials);
    spiketimes.RC    = spiketimes.all(MMCsorted.RCtrials);
    spiketimes.RefRC = spiketimes.all(MMCsorted.RefRCtrials);
    
    trials.C         = MMCsorted.Ctrials;
    trials.RC        = MMCsorted.RCtrials;
    trials.RefRC     = MMCsorted.RefRCtrials;
    
end

save([savepath MMCsorted.filename(10:end) '_' num2str(MMCsorted.maxgap*1000,'%0.2dmsMaxGap.mat')],'spiketimes','trials');
