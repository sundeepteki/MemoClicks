function st_sorted_getspiketimes(varargin)
clear all; close all; clc;

% function to load spike times and trial ids from sorted data
% ST. v1:            28.06.15
% Last modified:     01.07.15
% Backup:            st_sorted_getspiketimes_copy
% Last backup saved: 01.07.15

    
%% separate single units from multi-units and noise

% session -> channels -> units
% check all dir sizes for dir(datafolder)
% su/mu
% sort trial ids acc. to conditions
% save maxgap and electrode number
% get clicktimes as well
% save data in same folder
% set m = 1 in for loop if only one unit per session

%% Initialize
datapath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/';
conditions = {'MMC','MTC','MMCexp','MTCexp'};
disp(conditions);
cond = input('Enter condition to analyze - 1 to 4: ');

MMCinfo.datapath    = [datapath '#MMC'];
MMCinfo.savepath    = [MMCinfo.datapath filesep '#analyzed'];
MTCinfo.datapath    = [datapath '#MTC'];
MTCinfo.savepath    = [MTCinfo.datapath filesep '#analyzed'];
MMCexpinfo.datapath = [datapath '#MMCexp'];
MMCexpinfo.savepath = [MMCexpinfo.datapath filesep '#analyzed'];
MTCexpinfo.datapath = [datapath '#MTCexp'];
MTCexpinfo.savepath = [MTCexpinfo.datapath filesep '#analyzed'];
addpath(MMCinfo.datapath);
addpath(MTCinfo.datapath);
addpath(MMCexpinfo.datapath);
addpath(MTCexpinfo.datapath);

%% MMC

if cond == 1
    
    tmpdir = dir(MMCinfo.datapath);
    
    for i = 1:(length(tmpdir)-4)
        data{i} =  [tmpdir(i+4).name];
        load(data{i});
        fname = data{i}; fname = fname(1:end-4);
        MMCinfo.filename = fname;
        disp(['Loading file - ' fname]);
        
        chan = [];

        for j = 1:length(sortextras) % electrode number
            chan = [chan ~isempty(sortextras{j})];
            elec = find(chan==1);
            
            for k = 1:length(elec) % number of units
                num_units{j} = length(sortinfo{elec(k)}{1});
                
                for m = 1:num_units{j} % spike times and trial id
                    MMC.trials     = sortinfo{elec(k)}{1}(m).unitSpikes(1,:);
                    MMC.spiketimes = sortinfo{elec(k)}{1}(m).unitSpikes(2,:);
                    MMCinfo.elec = elec(k);
                    MMCinfo.unit = m;
                    save([MMCinfo.savepath filesep 'analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m) '.mat'],'MMC','MMCinfo');
                    savename = ['analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m)];
                    disp(['Saving file - ' savename]);
                end
                
                
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%% MTC %%%%%%%%%%%%%%
    
elseif cond == 2
    
    tmpdir = dir(MTCinfo.datapath);
    
    for i = 1:(length(tmpdir)-3)
        data{i} =  [tmpdir(i+3).name];
        load(data{i});
        fname = data{i}; fname = fname(1:end-4);
        MTCinfo.filename = fname;
        disp(['Loading file - ' fname]);
        
        chan = [];

        for j = 1:length(sortextras) % electrode number
            chan = [chan ~isempty(sortextras{j})];
            elec = find(chan==1);
            
            for k = 1:length(elec) % number of units
                num_units{j} = length(sortinfo{elec(k)}{1});
                
                for m = 1:num_units{j} % spike times and trial id
                    MTC.trials     = sortinfo{elec(k)}{1}(m).unitSpikes(1,:);
                    MTC.spiketimes = sortinfo{elec(k)}{1}(m).unitSpikes(2,:);
                    MTCinfo.elec = elec(k);
                    MTCinfo.unit = m;
                    save([MTCinfo.savepath filesep 'analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m) '.mat'],'MTC','MTCinfo');
                    savename = ['analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m)];
                    disp(['Saving file - ' savename]);
                end
                                   
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%% MMCexp %%%%%%%%%%%%%%%%%%%%
    
elseif cond == 3
    
    tmpdir = dir(MMCexpinfo.datapath);
    
    for i = 1:(length(tmpdir)-3)
        data{i} =  [tmpdir(i+3).name];
        load(data{i});
        fname = data{i}; fname = fname(1:end-4);
        MMCexpinfo.filename = fname;
        disp(['Loading file - ' fname]);
        
        chan = [];
        
        for j = 1:length(sortextras) % electrode number
            chan = [chan ~isempty(sortextras{j})];
            elec = find(chan==1);
            
            for k = 1:length(elec) % number of units
                num_units{j} = length(sortinfo{elec(k)}{1});
                
                for m = 1:num_units{j} % spike times and trial id
                    MMCexp.trials     = sortinfo{elec(k)}{1}(m).unitSpikes(1,:);
                    MMCexp.spiketimes = sortinfo{elec(k)}{1}(m).unitSpikes(2,:);
                    MMCexpinfo.elec = elec(k);
                    MMCexpinfo.unit = m;
                    save([MMCexpinfo.savepath filesep 'analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m) '.mat'],'MMCexp','MMCexpinfo');
                    savename = ['analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m)];
                    disp(['Saving file - ' savename]);
                end
                
                
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%% MTCexp %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif cond == 4
    
    tmpdir = dir(MTCexpinfo.datapath);
        
    for i = 1:(length(tmpdir)-4)
        data{i} =  [tmpdir(i+4).name];
        load(data{i});
        fname = data{i}; fname = fname(1:end-4);
        MTCexpinfo.filename = fname;
        disp(['Loading file - ' fname]);
        
        chan = [];
        
        for j = 1:length(sortextras) % electrode number
            chan = [chan ~isempty(sortextras{j})];
            elec = find(chan==1);
            
            for k = 1:length(elec) % number of units
                num_units{j} = length(sortinfo{elec(k)}{1});
                
                for m = 1:num_units{j} % spike times and trial id
                    MTCexp.trials     = sortinfo{elec(k)}{1}(m).unitSpikes(1,:);
                    MTCexp.spiketimes = sortinfo{elec(k)}{1}(m).unitSpikes(2,:);
                    MTCexpinfo.elec = elec(k);
                    MTCexpinfo.unit = m;
                    save([MTCexpinfo.savepath filesep 'analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m) '.mat'],'MTCexp','MTCexpinfo');
                    savename = ['analyzed_' fname '_elec' num2str(elec(k)) '_unit' num2str(m)];
                    disp(['Saving file - ' savename]);
                end
                
                
            end
            
        end
    end
    
end

