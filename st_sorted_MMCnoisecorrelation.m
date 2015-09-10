function [r,p] = st_sorted_MMCnoisecorrelation(varargin)
% clear all; close all; clc;

%{
Sundeep Teki.
Created:           03.07.15
Last modified:     11.07.15
Backup:            st_sorted_MMCnoisecorrelation_copy
Last backup saved: 11.07.15

- script to run noise correlation analysis
- on pairs_within of simultaneously recorded neurons
- use different time windows (pre-stim, stim of different windows)
- use no. of spikes on each trial
- compare across conditions (and time windows)

- shuffle trials not spike times


####
1. DEFINE WITHIN POPULATION pairs_within
2. DEFINE DIFFERENT POPULATIONS
high input rates control correlations

TO DO:
done. assign pair IDs
perform for different time windows
perform for within and between populations
mean(cell2mat(r.RefRC))

change savepath according to time window of correlation
skipped pair 1 (maroiller21) as it has different click timings
[r,fnames] = .... for group analysis

%}

%% Paths

datapath    = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/SpikeTimes/';
savepath    = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/Stim/';
savepathbtw = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/Stim_Btw/';
addpath(datapath);
addpath(savepath);
addpath(savepathbtw);
tmpdir      = dir(datapath);

%% Choose analyses

corrwithinpairs = 1;
analysecorrwithinpairs = 0;
corrbtwpairs    = 0;
analysecorrbtwpairs = 0;

%% Define pairs

    
    pairs_within = {{16,17}, {18,19}, {20,21}, {22,23}, {24,25}, {26,27}, {28,29}, ...
        {30,31}, {32,33}, {34,35}, {36,37}, {40,41}, {42,43}, {44,45},...
        {46,47}, {48,49}, {50,51}, {52,53}, {54,55}, {56,57}, {58,59},...
        {59,60}, {58,60}, {61,62}, {62,63}, {61,63}, {64,65}, {64,66},...
        {64,67}, {65,66}, {65,67}, {66,67}, {68,69}, {68,70}, {68,71},...
        {69,70}, {69,71}, {70,71}, {80,81}, {80,82}, {80,83}, {80,84},...
        {80,85}, {81,82}, {81,83}, {81,84}, {81,85}, {82,83}, {82,84},...
        {82,85}, {83,84}, {83,85}, {84,85}, {86,87}, {86,88}, {86,89},...
        {86,90}, {86,91}, {87,88}, {87,89}, {87,90}, {87,91}, {88,89},...
        {88,90}, {88,91}, {89,90}, {89,91}, {90,91}, {92,93}, {92,94},...
        {92,95}, {92,96}, {92,97}, {93,94}, {93,95}, {93,96}, {93,97},...
        {94,95}, {94,96}, {94,97}, {95,96}, {95,97}, {96,97}, {98,99},...
        {98,100}, {98,101}, {98,102}, {98,103}, {99,100}, {99,101},...
        {99,102}, {99,103}, {100,101}, {100,102}, {100,103}, {101,102},...
        {101,103}, {102,103}, {104,105}, {104,106}, {105,106}, {108,109},...
        {110,111}, {112,113}, {114,115}, {116,117}, {118,119}, {120,121},...
        {122,123}, {124,125}, {126,127}}; % pair ids from 1-127
    
    
    %% Analyse corr within each pair

if(corrwithinpairs)
    
    for j = 1:length(pairs_within) % for j = 111; set + 2 instead of +4
        
        % load data for both pairs_within
        data1{j} =  [tmpdir(pairs_within{j}{1}+4).name];
        data2{j} =  [tmpdir(pairs_within{j}{2}+4).name];
        fname{1} =  [tmpdir(pairs_within{j}{1}+4).name];
        fname{2} =  [tmpdir(pairs_within{j}{2}+4).name];
        
        load(data1{j});
        spktr1C              = spiketimes.C;
        spktr1RC             = spiketimes.RC;
        spktr1RefRC          = spiketimes.RefRC;
        shuffledspktr1C      = tshuffle(spiketimes.C);
        shuffledspktr1RC     = tshuffle(spiketimes.RC);
        shuffledspktr1RefRC  = tshuffle(spiketimes.RefRC);


        clear spiketimes;
        
        load(data2{j});
        spktr2C              = spiketimes.C;
        spktr2RC             = spiketimes.RC;
        spktr2RefRC          = spiketimes.RefRC;
        shuffledspktr2C      = tshuffle(spiketimes.C);
        shuffledspktr2RC     = tshuffle(spiketimes.RC);
        shuffledspktr2RefRC  = tshuffle(spiketimes.RefRC);        
        clear spiketimes;
        
        disp(['Loading file - ' tmpdir(pairs_within{j}{1}+2).name(1:end-4)])
        disp(['Loading file - ' tmpdir(pairs_within{j}{2}+2).name(1:end-4)])
        
        %% within pair analysis
        
        time_window = 0.050; % 50ms, try different time scales
        % time_window = 0.025;
        % time_window = 0.1;
        
        % dur_prestim    = [0 0.6];  tON = 0;
        % dur_stim     = [0.6 2.1]; tON = 0.6;
        % dur_rep1     = [0.6 1.1]; tON = 0.6;
        % dur_rep2     = [1.1 1.6]; tON = 1.1;
        % dur_rep3     = [1.6 2.1]; tON = 1.6;
        % dur_poststim = [2.1 2.45]; tON = 2.1;
        % dur_torc     = [2.45 2.95]; tON = 2.45;
        dur_all      = [0 3.75]; tON = 0;
        
        
        %% compute pairwise correlations for C
        
        for i = 1:length(trials.C)
            
            if length(spktr1C{i}) < length(spktr2C{i})
                spktr1C{i}(length(spktr1C{i})+1 : length(spktr2C{i})) = 0;
            elseif length(spktr1C{i}) > length(spktr2C{i})
                spktr2C{i}(length(spktr2C{i})+1 : length(spktr1C{i})) = 0;
            end
            
            if length(shuffledspktr1C{i}) < length(shuffledspktr2C{i})
                shuffledspktr1C{i}(length(shuffledspktr1C{i})+1 : length(shuffledspktr2C{i})) = 0;
            elseif length(shuffledspktr1C{i}) > length(shuffledspktr2C{i})
                shuffledspktr2C{i}(length(shuffledspktr2C{i})+1 : length(shuffledspktr1C{i})) = 0;
            end

            [r.C{i},p.C{i}]    = st_spike_corrcoeff(spktr1C{i},spktr2C{i},time_window,tON,tON + diff(dur_all));
            [r.shC{i},p.shC{i}] = st_spike_corrcoeff(shuffledspktr1C{i},shuffledspktr2C{i},time_window,tON,tON + diff(dur_all));

        end
        

        %% compute pairwise correlations for RC
        
        for i = 1:length(trials.RC)
                       
            if length(spktr1RC{i}) < length(spktr2RC{i})
                spktr1RC{i}(length(spktr1RC{i})+1 : length(spktr2RC{i})) = 0;
            elseif length(spktr1RC{i}) > length(spktr2RC{i})
                spktr2RC{i}(length(spktr2RC{i})+1 : length(spktr1RC{i})) = 0;
            end
            
            if length(shuffledspktr1RC{i}) < length(shuffledspktr2RC{i})
                shuffledspktr1RC{i}(length(shuffledspktr1RC{i})+1 : length(shuffledspktr2RC{i})) = 0;
            elseif length(shuffledspktr1RC{i}) > length(shuffledspktr2RC{i})
                shuffledspktr2RC{i}(length(shuffledspktr2RC{i})+1 : length(shuffledspktr1RC{i})) = 0;
            end
            
            [r.RC{i},p.RC{i}] = st_spike_corrcoeff(spktr1RC{i},spktr2RC{i},time_window,tON,tON + diff(dur_all));
            [r.shRC{i},p.shRC{i}] = st_spike_corrcoeff(shuffledspktr1RC{i},shuffledspktr2RC{i},time_window,tON,tON + diff(dur_all));

        end
                
        
        %% compute pairwise correlations for RefRC
        
        for i = 1:length(trials.RefRC)
            
            if length(spktr1RefRC{i}) < length(spktr2RefRC{i})
                spktr1RefRC{i}(length(spktr1RefRC{i})+1 : length(spktr2RefRC{i})) = 0;
            elseif length(spktr1RefRC{i}) > length(spktr2RefRC{i})
                spktr2RefRC{i}(length(spktr2RefRC{i})+1 : length(spktr1RefRC{i})) = 0;
            end
            
            if length(shuffledspktr1RefRC{i}) < length(shuffledspktr2RefRC{i})
                shuffledspktr1RefRC{i}(length(shuffledspktr1RefRC{i})+1 : length(shuffledspktr2RefRC{i})) = 0;
            elseif length(shuffledspktr1RefRC{i}) > length(shuffledspktr2RefRC{i})
                shuffledspktr2RefRC{i}(length(shuffledspktr2RefRC{i})+1 : length(shuffledspktr1RefRC{i})) = 0;
            end
            
            [r.RefRC{i},p.RefRC{i}] = st_spike_corrcoeff(spktr1RefRC{i},spktr2RefRC{i},time_window,tON,tON + diff(dur_all));
            [r.shRefRC{i},p.shRefRC{i}] = st_spike_corrcoeff(shuffledspktr1RefRC{i},shuffledspktr2RefRC{i},time_window,tON,tON + diff(dur_all));

        end
        
       
        save
        save([savepath 'PairShuffled' num2str(j) '_corr_stim.mat'],'r','p','fname');
        
    end
    
end

%% Group analysis of within pair correlation

if(analysecorrwithinpairs)
    
    ravg.RefRC   = [];
    ravg.shRefRC = [];
    ravg.RC      = [];
    ravg.shRC    = [];
    ravg.C       = [];
    ravg.shC     = [];
    
    
    for j = 1:length(pairs_within) % for j = 111; set + 2 instead of +4
        
        % load data for both pairs_within
        load(['Pair' num2str(j) '_corr_stim']);
        
        ravg.shRefRC = [ravg.shRefRC nanmean(cell2mat(r.shRefRC))];
        ravg.RefRC   = [ravg.RefRC nanmean(cell2mat(r.RefRC))];
        
        ravg.shRC = [ravg.shRC nanmean(cell2mat(r.shRC))];
        ravg.RC   = [ravg.RC nanmean(cell2mat(r.RC))];
        
        ravg.shC = [ravg.shC nanmean(cell2mat(r.shC))];
        ravg.C   = [ravg.C nanmean(cell2mat(r.C))];
        
        
    end
    
end


%% Between pair analysis

    pairs_btw = {{1:5}, {6:10}, {11:13},{14:17},{18,19},{20:37}, {38:100},{101}, {102:105},{106:110}};

if(corrbtwpairs)
        
    for j = 1:length(pairs_btw)
        
        rpop.C = []; rpop.RC = []; rpop.RefRC = []; fnames = {};
        rpop.shC = []; rpop.shRC = []; rpop.shRefRC = [];
        
        for k=1:length(pairs_btw{j}{1})
            load([savepath 'Pair' num2str(pairs_btw{j}{1}(k)) '_corr_stim']);
            
            rpop.C = [rpop.C  nanmean(cell2mat(r.C))];
            rpop.RC = [rpop.RC  nanmean(cell2mat(r.RC))];
            rpop.RefRC = [rpop.RefRC  nanmean(cell2mat(r.RefRC))];
            
            rpop.shC = [rpop.shC  nanmean(cell2mat(r.shC))];
            rpop.shRC = [rpop.shRC  nanmean(cell2mat(r.shRC))];
            rpop.shRefRC = [rpop.shRefRC  nanmean(cell2mat(r.shRefRC))];
            
            clear r p;
            fnames{k} = {fname{1}(1:14)};
            
            % pool correlations across groups - don't take mean
            
        end
        
        save([savepathbtw 'BtwPair' num2str(j) '_corr_stim.mat'],'rpop','fname');
    end
end

%% Group analysis of between pair correlation

if(analysecorrbtwpairs)
    
    ravgbtw.RefRC   = [];
    ravgbtw.shRefRC = [];
    ravgbtw.RC      = [];
    ravgbtw.shRC    = [];
    ravgbtw.C       = [];
    ravgbtw.shC     = [];
    
    
    for j = 1:length(pairs_btw) 
        
        load(['BtwPair' num2str(j) '_corr_stim']);
        
        ravgbtw.shRefRC = [ravgbtw.shRefRC nanmean(rpop.shRefRC)];
        ravgbtw.RefRC   = [ravgbtw.RefRC nanmean(rpop.RefRC)];
        
        ravgbtw.shRC = [ravgbtw.shRC nanmean(rpop.shRC)];
        ravgbtw.RC   = [ravgbtw.RC nanmean(rpop.RC)];
        
        ravgbtw.shC = [ravgbtw.shC nanmean(rpop.shC)];
        ravgbtw.C   = [ravgbtw.C nanmean(rpop.C)];
        
        
    end
    
end

