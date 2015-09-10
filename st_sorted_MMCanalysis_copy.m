function [MMCsorted,MMCsortedall,MMCsortedavg,MMCgroupstats] = st_sorted_MMCanalysis(varargin)
clear all; close all; clc;

%{
function to analyze MMC sorted data

ST. v1:            01.07.15
Last modified:     09.07.15
Backup:            st_analyze_MMCsorted_copy
Last backup saved: 08.07.15


% analyses
raster,
psth,
firing rate as a function of: trials, click id, click interval
compare baseline to stimulus
compare: repeats 1 - 3
compare: first 5 vs. last 5 trials

% spikes to analyse
54, 55: 40e04
77: 45b04e3
83,84,85: 46d03e1

% click analysis left: 39

%}

%% Set paths

MMCsorted.datapath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/';
MMCsorted.savepath = [MMCsorted.datapath '#analyzed/'];
addpath(MMCsorted.datapath);
addpath(MMCsorted.savepath);

%% analysis flags

loaddata                    = 1;
loadMMCunsorteddata         = 1;

getallsorteddata            = 0;
getclicktimes               = 0;
getspiketimes               = 0;

plotrasterstim              = 0;
plotspikeratestim           = 0;
plotspikeraterepeat         = 0;
plotpsth                    = 0;
plotclickpsth               = 0;
plotclickspikeratepertrial  = 0;
plotclickspikerateperclick  = 0;
plotclickintervalpsth       = 0;

dogroupanalysis             = 1;
dogroupstats                = 1;
savedataflag                = 0;
savegroupdataflag           = 0;


%% load sorted data

tmpdir = dir(MMCsorted.datapath);
q = [];

for ij = [1:53 56:76 78:82 86:133]%1:(length(tmpdir)-4)
    
    close all;
    
    data{ij} =  [tmpdir(ij+4).name];
    load(data{ij});
    fname = data{ij}; fname = fname(1:end-4);    
    disp(['Loading file - ' fname]);
    
    MMCsorted.filename   = fname;
    MMCsorted.trials     = MMC.trials;
    MMCsorted.spiketimes = MMC.spiketimes;
    MMCsorted.elec       = MMCinfo.elec;
    MMCsorted.unit       = MMCinfo.unit;    
    
    
    %% load MMC unsorted analyzed data
    
    if(loadMMCunsorteddata)
        
        MMCanalyzed.path   = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#Data/#Analysis/MMC/';
        tmpdir1 = dir(MMCanalyzed.path);
        
        for i1 = 1:(length(tmpdir1)-9)
            data1{i1} =  [tmpdir1(i1+9).name];
            fname1 = data1{i1};
            longfiles{i1} = fname1;
            fname2 = fname1(1:14);
            shortfiles{i1} = fname2;
        end
        
        % load corresponding data
        match = 0;
        matchid = [];
        for i2 = 1:length(tmpdir1)-9
            if(strcmpi(MMCsorted.filename(10:23),shortfiles{i2}))
                match = 1;
                matchid = [matchid i2];
                break;
            end
        end
        
        %% Save relevant parameters - trial id, clicktimes, MMCsortedgap
        
        load([MMCanalyzed.path longfiles{matchid}]);
        disp(['Loading file - ' longfiles{matchid}]);
        MMCsorted.MMCinfo       = MMC;
        MMCsorted.Cunsorted     = C;
        MMCsorted.RCunsorted    = RC;
        MMCsorted.RefRCunsorted = RefRC;
        MMCsorted.clickunsorted = click;
        MMCsorted.noiseunsorted = noise;
        MMCsorted.Ctrials     = MMCsorted.Cunsorted.trial;
        MMCsorted.RCtrials    = MMCsorted.RCunsorted.trial;
        MMCsorted.RefRCtrials = MMCsorted.RefRCunsorted.trial;
        MMCsorted.maxgap      = MMCsorted.MMCinfo.MaxGap;
        MMCsorted.seed        = MMCsorted.MMCinfo.RefRCSeed;
        
    end
    
    %% Load saved data for subsequent analysis
    
    if(loaddata)
        load([MMCsorted.savepath 'Results' filesep MMCsorted.filename(10:end) '_' num2str(MMCsorted.maxgap*1000,'%0.2dmsMaxGap.mat')])
    end
    
    %% collect data for group analysis
        
    
    if(getallsorteddata)
        
%         q = [q MMCsorted.MMCinfo.MaxGap];
%         MMCsortedall.maxgaps = q;
        % baseline, stimulus, post-baseline
        MMCsortedall.Cavgspikeratebaseline{ij}     = MMCsorted.Cavgspikeratebaseline;
        MMCsortedall.RCavgspikeratebaseline{ij}    = MMCsorted.RCavgspikeratebaseline;
        MMCsortedall.RefRCavgspikeratebaseline{ij} = MMCsorted.RefRCavgspikeratebaseline;
        
        MMCsortedall.Cavgspikerateclick{ij}     = MMCsorted.Cavgspikerateclick;
        MMCsortedall.RCavgspikerateclick{ij}    = MMCsorted.RCavgspikerateclick;
        MMCsortedall.RefRCavgspikerateclick{ij} = MMCsorted.RefRCavgspikerateclick;
        
        MMCsortedall.Cavgspikeratepostclick{ij}     = MMCsorted.Cavgspikeratepostclick;
        MMCsortedall.RCavgspikeratepostclick{ij}    = MMCsorted.RCavgspikeratepostclick;
        MMCsortedall.RefRCavgspikeratepostclick{ij} = MMCsorted.RefRCavgspikeratepostclick;
        
        % repeat 1,2,3
        MMCsortedall.Cavgspikerateclickrep1{ij}     = MMCsorted.Cavgspikerateclickrep1;
        MMCsortedall.RCavgspikerateclickrep1{ij}    = MMCsorted.RCavgspikerateclickrep1;
        MMCsortedall.RefRCavgspikerateclickrep1{ij} = MMCsorted.RefRCavgspikerateclickrep1;
        
        MMCsortedall.Cavgspikerateclickrep2{ij}     = MMCsorted.Cavgspikerateclickrep2;
        MMCsortedall.RCavgspikerateclickrep2{ij}    = MMCsorted.RCavgspikerateclickrep2;
        MMCsortedall.RefRCavgspikerateclickrep2{ij} = MMCsorted.RefRCavgspikerateclickrep2;
        
        MMCsortedall.Cavgspikerateclickrep3{ij}     = MMCsorted.Cavgspikerateclickrep3;
        MMCsortedall.RCavgspikerateclickrep3{ij}    = MMCsorted.RCavgspikerateclickrep3;
        MMCsortedall.RefRCavgspikerateclickrep3{ij} = MMCsorted.RefRCavgspikerateclickrep3;
        
        % click-evoked spikerate
        MMCsortedall.Cclickavgspikerate{ij}     = MMCsorted.Cclickavgspikerate;
        MMCsortedall.RCclickavgspikerate{ij}    = MMCsorted.RCclickavgspikerate;
        MMCsortedall.RefRCclickavgspikerate{ij} = MMCsorted.RefRCclickavgspikerate;
        
        MMCsortedall.Cclickavgspikerateperclickfirst5{ij}     = MMCsorted.Cclickavgspikerateperclickfirst5;
        MMCsortedall.RCclickavgspikerateperclickfirst5{ij}    = MMCsorted.RCclickavgspikerateperclickfirst5;
        MMCsortedall.RefRCclickavgspikerateperclickfirst5{ij} = MMCsorted.RefRCclickavgspikerateperclickfirst5;
        
        MMCsortedall.Cclickavgspikerateperclickfirst{ij}     = MMCsorted.Cclickavgspikerateperclickfirst;
        MMCsortedall.RCclickavgspikerateperclickfirst{ij}    = MMCsorted.RCclickavgspikerateperclickfirst;
        MMCsortedall.RefRCclickavgspikerateperclickfirst{ij} = MMCsorted.RefRCclickavgspikerateperclickfirst;
        
        MMCsortedall.Cclickavgspikerateperclicklast5{ij}     = MMCsorted.Cclickavgspikerateperclicklast5;
        MMCsortedall.RCclickavgspikerateperclicklast5{ij}    = MMCsorted.RCclickavgspikerateperclicklast5;
        MMCsortedall.RefRCclickavgspikerateperclicklast5{ij} = MMCsorted.RefRCclickavgspikerateperclicklast5;
        
    end
    
    %% get click times
    
    if(getclicktimes)
        
        stim     = {};
        events   = {};
        trialnum = 100;
        o        = MemoClicks();
        o        = set(o,'maxgap',MMCsorted.MMCinfo.MaxGap);
        o        = set(o,'Key',   MMCsorted.MMCinfo.RefRCSeed);
        
        for i = 1:trialnum
            [w, ev,o] = waveform (o,1,[],[],i);
            stim{i}   = w;
            events{i} = ev;
        end
        
        MMCsorted.stimulus   = get(o,'stimulus');
        MMCsorted.seedall    = get(o,'Seeds');
        MMCsorted.clicktimes = get(o,'ClickTimes');
        
    end
    
    %% Sorted data - spiketimes, spiketrialids
    
    xC = [];
    xRC = [];
    xRefRC = [];
    
    for i = 1:length(MMCsorted.trials)
        for iC = 1:length(MMCsorted.Ctrials)
            if MMCsorted.trials(i)==MMCsorted.Ctrials(iC)
                xC = [xC i];
            end
        end
        
        for iRC = 1:length(MMCsorted.RCtrials)
            if MMCsorted.trials(i)==MMCsorted.RCtrials(iRC)
                xRC = [xRC i];
            end
        end
        
        for iRefRC = 1:length(MMCsorted.RefRCtrials)
            if MMCsorted.trials(i)==MMCsorted.RefRCtrials(iRefRC)
                xRefRC = [xRefRC i];
            end
        end
    end
    
    MMCsorted.Cspiketrialid = MMCsorted.trials(xC);
    MMCsorted.RCspiketrialid = MMCsorted.trials(xRC);
    MMCsorted.RefRCspiketrialid = MMCsorted.trials(xRefRC);
    
    MMCsorted.Cspiketime = MMCsorted.spiketimes(xC)/MMCsorted.MMCinfo.spikefs;
    MMCsorted.RCspiketime = MMCsorted.spiketimes(xRC)/MMCsorted.MMCinfo.spikefs;
    MMCsorted.RefRCspiketime = MMCsorted.spiketimes(xRefRC)/MMCsorted.MMCinfo.spikefs;
    
    
    
    %% Get spike count vectors for variance analysis, cf. Churchland et al. 2010 Nat Neurosci
    
    if(getspiketimes)
        
        binstep = 1;     % 1ms
        binsize = 3750;  % 3.5s of stimulus including post-noise baseline
        spikes = zeros(MMCsorted.MMCinfo.NumTrials,binsize);
        spikenum = zeros(1,MMCsorted.MMCinfo.NumTrials);
        Cspikes = zeros(MMCsorted.MMCinfo.NumTrials/2,binsize);
        RCspikes = zeros(MMCsorted.MMCinfo.NumTrials/4,binsize);
        RefRCspikes = zeros(MMCsorted.MMCinfo.NumTrials/4,binsize);
        
        for i = 1:1:MMCsorted.MMCinfo.NumTrials
            
            for j = 1:binstep:binsize
                hit = ceil((MMCsorted.spiketimes(MMCsorted.trials==i)*1000/31250));
                spikes(i,hit) = 1;
            end
            
            spikenum(i) = length(find(spikes(i,:)==1));
            
        end
        
        clickonset  = MMCsorted.clickunsorted;
        noiseonset  = MMCsorted.noiseunsorted;
        Cspikes     = spikes(MMCsorted.Ctrials,:);
        RCspikes    = spikes(MMCsorted.RCtrials,:);
        RefRCspikes = spikes(MMCsorted.RefRCtrials,:);
        trialid.C     = MMCsorted.Ctrials;
        trialid.RC    = MMCsorted.RCtrials;
        trialid.RefRC = MMCsorted.RefRCtrials;
        
        save([MMCsorted.savepath 'SpikeTimes' filesep MMCsorted.filename(10:end) '_spikes.mat'],...
            'spikes','spikenum','clickonset','noiseonset','Cspikes','RCspikes','RefRCspikes','trialid');
        
    end
    
    
    %% Plot rasters
    
    if(plotrasterstim)
        
        %%%% Raster - C %%%%%
        
        fig1 = figure;
        subplot(311); hold on;
        set(gca,'Position',[0.13 0.63 0.7750 0.35])
        
        for ii = 1:length(MMCsorted.Ctrials)
            
            tmp = [find(MMCsorted.Cspiketrialid == MMCsorted.Ctrials(ii))];
            
            for j = 1:length(tmp)
                plot(MMCsorted.Cspiketime(tmp),[ii*ones(1,length(tmp))],'k.','MarkerSize',1);
            end
        end
        
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.5:length(MMCsorted.Ctrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.Ctrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + 2*MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.Ctrials)],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.5:length(MMCsorted.Ctrials)],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.5:length(MMCsorted.Ctrials)],'b.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.5:length(MMCsorted.Ctrials)],'b.');
        ylim([0 length(MMCsorted.Ctrials)]);
        ylabel('C','FontSize',12,'FontWeight','bold','Color',[0 0 0]);
        xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color',[0 0 0]);
        
        %%%%% Raster - RC %%%%%%
        
        subplot(312); hold on;
        set(gca,'Position',[0.13 0.375 0.7750 0.175])
        
        for ii = 1:length(MMCsorted.RCtrials)
            
            tmp = [find(MMCsorted.RCspiketrialid == MMCsorted.RCtrials(ii))];
            
            for j = 1:length(tmp)
                plot(MMCsorted.RCspiketime(tmp),[ii*ones(1,length(tmp))],'k.','MarkerSize',1);
            end
        end
        
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.5:length(MMCsorted.RCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.RCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + 2*MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.RCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.5:length(MMCsorted.RCtrials)],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.5:length(MMCsorted.RCtrials)],'b.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.5:length(MMCsorted.RCtrials)],'b.');
        ylim([0 length(MMCsorted.RCtrials)]);
        ylabel('RC','FontSize',12,'FontWeight','bold','Color',[0 0 0]);
        
        %%%%%% Raster - RefRC %%%%%%%
        
        subplot(313); hold on;
        set(gca,'Position',[0.13 0.1 0.7750 0.175])
        
        for ii = 1:length(MMCsorted.RefRCtrials)
            
            tmp = [find(MMCsorted.RefRCspiketrialid == MMCsorted.RefRCtrials(ii))];
            
            for j = 1:length(tmp)
                plot(MMCsorted.RefRCspiketime(tmp),[ii*ones(1,length(tmp))],'k.','MarkerSize',1);
            end
        end
        
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.5:length(MMCsorted.RefRCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.RefRCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.start + 2*MMCsorted.MMCinfo.ClickDur,2),[0:2:length(MMCsorted.RefRCtrials)],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.5:length(MMCsorted.RefRCtrials)],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.5:length(MMCsorted.RefRCtrials)],'b.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.5:length(MMCsorted.RefRCtrials)],'b.');
        ylim([0 length(MMCsorted.RefRCtrials)]);
        ylabel('RefRC','FontSize',12,'FontWeight','bold','Color',[0 0 0]);
        xlabel(MMCsorted.filename,'FontSize',12,'FontWeight','bold','Color',[0 0 0]);
        
        %%%% save raster %%%%%
        %         saveas(fig1,[MMCsorted.savepath filesep 'Raster/Sorted_1_Raster_' MMCsorted.filename(10:end) '.tiff']);
        
    end
    
    %% Plot PSTH
    
    if(plotpsth)
        
        bin      = 0.025; %25ms bin
        baseline = MMCsorted.clickunsorted.start;
        edges    = [0:bin:MMCsorted.noiseunsorted.post];
        
        % PSTH - C
        MMCsorted.Cpsth           = zeros(length(edges),1);
        MMCsorted.Cspikerate      = zeros(length(edges),1);
        MMCsorted.Cnormspikerate  = zeros(length(edges),1);
        MMCsorted.Chistspiketime  = zeros(length(edges),1);
        
        for j = 1:length(MMCsorted.Ctrials)
            tmpC1{j}         = find(MMCsorted.Cspiketrialid == MMCsorted.Ctrials(j));
            MMCsorted.Chistspiketime = histc( (MMCsorted.Cspiketime(tmpC1{j}) ) ,edges);
            
            if size(MMCsorted.Cpsth,1)==size(MMCsorted.Chistspiketime,1)
                MMCsorted.Cpsth          = MMCsorted.Cpsth + MMCsorted.Chistspiketime;
            else
                MMCsorted.Cpsth          = MMCsorted.Cpsth + MMCsorted.Chistspiketime';
            end
            
            MMCsorted.Cspikerate     = MMCsorted.Cpsth/(bin*length(MMCsorted.Ctrials));
            tempC           = MMCsorted.Cspikerate - MMCsorted.Cspikerate./mean(MMCsorted.Cspikerate(1:baseline/bin));
            MMCsorted.Cnormspikerate = tempC/max(tempC);
        end
        
        fig2 = figure;
        ax1 = subplot(311); hold on;
        bar(edges,MMCsorted.Cspikerate);
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start + MMCsorted.MMCinfo.ClickDur,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start + 2*MMCsorted.MMCinfo.ClickDur,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'c.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.1:max(MMCsorted.Cspikerate)/2],'c.');
        title('C','FontSize',12,'FontWeight','bold');
        
        % PSTH - RC
        MMCsorted.RCpsth           = zeros(length(edges),1);
        MMCsorted.RCspikerate      = zeros(length(edges),1);
        MMCsorted.RCnormspikerate  = zeros(length(edges),1);
        MMCsorted.RChistspiketime  = zeros(length(edges),1);
        
        for j = 1:length(MMCsorted.RCtrials)
            tmpRC1{j}         = find(MMCsorted.RCspiketrialid == MMCsorted.RCtrials(j));
            MMCsorted.RChistspiketime = histc( (MMCsorted.RCspiketime(tmpRC1{j}) ) ,edges);
            
            if size(MMCsorted.RCpsth,1)==size(MMCsorted.RChistspiketime,1)
                MMCsorted.RCpsth          = MMCsorted.RCpsth + MMCsorted.RChistspiketime;
            else
                MMCsorted.RCpsth          = MMCsorted.RCpsth + MMCsorted.RChistspiketime';
            end
            
            MMCsorted.RCspikerate     = MMCsorted.RCpsth/(bin*length(MMCsorted.RCtrials));
            tempRC           = MMCsorted.RCspikerate - MMCsorted.RCspikerate./mean(MMCsorted.RCspikerate(1:baseline/bin));
            MMCsorted.RCnormspikerate = tempRC/max(tempRC);
        end
        
        ax2 = subplot(312); hold on;
        bar(edges,MMCsorted.RCspikerate);
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start+MMCsorted.MMCinfo.ClickDur,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start+MMCsorted.MMCinfo.ClickDur*2,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'c.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.1:max(MMCsorted.RCspikerate)/2],'c.');
        ylabel('Spike rate (Hz)','FontSize',12,'FontWeight','bold')
        title('RC','FontSize',12,'FontWeight','bold');
        
        % PSTH - RefRC
        MMCsorted.RefRCpsth           = zeros(length(edges),1);
        MMCsorted.RefRCspikerate      = zeros(length(edges),1);
        MMCsorted.RefRCnormspikerate  = zeros(length(edges),1);
        MMCsorted.RefRChistspiketime  = zeros(length(edges),1);
        
        for j = 1:length(MMCsorted.RefRCtrials)
            tmpRefRC1{j}         = find(MMCsorted.RefRCspiketrialid == MMCsorted.RefRCtrials(j));
            MMCsorted.RefRChistspiketime = histc( (MMCsorted.RefRCspiketime(tmpRefRC1{j}) ) ,edges);
            
            if size(MMCsorted.RefRCpsth,1)==size(MMCsorted.RefRChistspiketime,1)
                MMCsorted.RefRCpsth          = MMCsorted.RefRCpsth + MMCsorted.RefRChistspiketime;
            else
                MMCsorted.RefRCpsth          = MMCsorted.RefRCpsth + MMCsorted.RefRChistspiketime';
            end
            
            MMCsorted.RefRCspikerate     = MMCsorted.RefRCpsth/(bin*length(MMCsorted.RefRCtrials));
            tempRefRC           = MMCsorted.RefRCspikerate - MMCsorted.RefRCspikerate./mean(MMCsorted.RefRCspikerate(1:baseline/bin));
            MMCsorted.RefRCnormspikerate = tempRefRC/max(tempRefRC);
        end
        
        ax3 = subplot(313); hold on;
        bar(edges,MMCsorted.RefRCspikerate);
        plot(round(MMCsorted.clickunsorted.start,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start+MMCsorted.MMCinfo.ClickDur,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.start+MMCsorted.MMCinfo.ClickDur*2,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'g.');
        plot(round(MMCsorted.clickunsorted.end,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'g.');
        plot(round(MMCsorted.noiseunsorted.start,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'c.');
        plot(round(MMCsorted.noiseunsorted.end,2),[0:0.1:max(MMCsorted.RefRCspikerate)/2],'c.');
        xlabel('Time (s)','FontSize',12,'FontWeight','bold')
        title('RefRC','FontSize',12,'FontWeight','bold');
        
        
        %%%% save psth %%%%%
        %         saveas(fig2,[MMCsorted.savepath filesep 'PSTH/Sorted_2_PSTH_' MMCsorted.filename(10:end) '.tiff']);
        
    end
    
    %% Average spike rates across baseline and stimulus
    
    %  bins based on RefRC
    bininv          = 1/0.025;
    bin.preclick    = 1 : round(MMCsorted.clickunsorted.start*bininv);
    bin.click       = bin.preclick(end)+1 : round(MMCsorted.clickunsorted.end*bininv);
    bin.postclick   = bin.click(end)+1 : round(MMCsorted.clickunsorted.post*bininv);
    
    if(plotspikeratestim)
        
        % Average spikerate during pre-click baseline
        MMCsorted.Cavgspikeratebaseline     = mean(MMCsorted.Cspikerate(bin.preclick));
        MMCsorted.RCavgspikeratebaseline    = mean(MMCsorted.RCspikerate(bin.preclick));
        MMCsorted.RefRCavgspikeratebaseline = mean(MMCsorted.RefRCspikerate(bin.preclick));
        
        % Average spikerate during clicktrain
        MMCsorted.Cavgspikerateclick     = mean(MMCsorted.Cspikerate(bin.click));
        MMCsorted.RCavgspikerateclick    = mean(MMCsorted.RCspikerate(bin.click));
        MMCsorted.RefRCavgspikerateclick = mean(MMCsorted.RefRCspikerate(bin.click));
        
        % Average spikerate during post-click
        MMCsorted.Cavgspikeratepostclick     = mean(MMCsorted.Cspikerate(bin.postclick));
        MMCsorted.RCavgspikeratepostclick    = mean(MMCsorted.RCspikerate(bin.postclick));
        MMCsorted.RefRCavgspikeratepostclick = mean(MMCsorted.RefRCspikerate(bin.postclick));
        
        fig3 = figure; hold on;
        xticks = [0.5 1 1.5];
        plot(xticks, [MMCsorted.Cavgspikeratebaseline MMCsorted.Cavgspikerateclick MMCsorted.Cavgspikeratepostclick],'ko-',...
            xticks, [MMCsorted.RCavgspikeratebaseline MMCsorted.RCavgspikerateclick MMCsorted.RCavgspikeratepostclick],'bo-',...
            xticks, [MMCsorted.RefRCavgspikeratebaseline MMCsorted.RefRCavgspikerateclick MMCsorted.RefRCavgspikeratepostclick],'ro-',...
            'LineWidth',2,'MarkerSize',8);
        set(gca,'XTick', xticks);
        set(gca,'XTickLabel',[{'Pre-stim'} {'Stim'} {'Post-stim'}],'FontSize',12,'FontWeight','bold');
        
        xlim([0.4 1.6]);
        ylabel('Average Spike rate (Hz)','FontSize',12,'FontWeight','bold');
        legend('C','RC','RefRC','Location','Best');
        %         title(['SpikeRates: ' MMCdata.infoname{1}(1:end-4) 'Elec' num2str(MMC.elecID) 'MaxGap' num2str(MMC.MaxGap*1000,'%0.2dms')])
        saveas(fig3,[MMCsorted.savepath filesep 'SpikeRates/Sorted_3_SpikerateStim_' MMCsorted.filename(10:end) '.tiff']);
        
    end
    
    
    %% Average spike rates across repeats
    
    if(plotspikeraterepeat)
        
        bin.rep1 = bin.click(1)  : bin.click(1)+round(click.dur/3,1)*bininv;
        bin.rep2 = bin.rep1(end)+1  : bin.rep1(end)+1+round(click.dur/3,1)*bininv;
        bin.rep3 = bin.rep2(end)+1  : bin.rep2(end)+1+round(click.dur/3,1)*bininv;
        
        MMCsorted.Cavgspikerateclickrep1 = mean(MMCsorted.Cspikerate(bin.rep1));
        MMCsorted.Cavgspikerateclickrep2 = mean(MMCsorted.Cspikerate(bin.rep2));
        MMCsorted.Cavgspikerateclickrep3 = mean(MMCsorted.Cspikerate(bin.rep3));
        
        MMCsorted.RCavgspikerateclickrep1 = mean(MMCsorted.RCspikerate(bin.rep1));
        MMCsorted.RCavgspikerateclickrep2 = mean(MMCsorted.RCspikerate(bin.rep2));
        MMCsorted.RCavgspikerateclickrep3 = mean(MMCsorted.RCspikerate(bin.rep3));
        
        MMCsorted.RefRCavgspikerateclickrep1 = mean(MMCsorted.RefRCspikerate(bin.rep1));
        MMCsorted.RefRCavgspikerateclickrep2 = mean(MMCsorted.RefRCspikerate(bin.rep2));
        MMCsorted.RefRCavgspikerateclickrep3 = mean(MMCsorted.RefRCspikerate(bin.rep3));
        
        fig4 = figure; hold on;
        xticks = [0.5 1 1.5];
        plot(xticks, [MMCsorted.Cavgspikerateclickrep1 MMCsorted.Cavgspikerateclickrep2 MMCsorted.Cavgspikerateclickrep3],'ko-',...
            xticks, [MMCsorted.RCavgspikerateclickrep1 MMCsorted.RCavgspikerateclickrep2 MMCsorted.RCavgspikerateclickrep3],'bo-',...
            xticks, [MMCsorted.RefRCavgspikerateclickrep1 MMCsorted.RefRCavgspikerateclickrep2 MMCsorted.RefRCavgspikerateclickrep3],'ro-',...
            'LineWidth',2,'MarkerSize',8);
        
        set(gca,'XTick', xticks);
        set(gca,'XTickLabel',[{'Repeat1'} {'Repeat2'} {'Repeat3'}],'FontSize',12,'FontWeight','bold');
        xlim([0.4 1.6]);
        ylabel('Average Spike rate (Hz)','FontSize',12,'FontWeight','bold');
        legend('C','RC','RefRC','Location','Best');
        saveas(fig4,[MMCsorted.savepath filesep 'SpikeRates/Sorted_4_SpikerateRep_' MMCsorted.filename(10:end) '.tiff']);
    end
    
    
    %% Click-evoked PSTH %%%
    
    if(plotclickpsth)
        
        % ClickPSTH - C
        
        bin      = 0.001;
        binsize  = 0.020; %20ms
        tmp1C    = {};
        MMCsorted.Cnum_psthtrials = length(MMCsorted.Ctrials); % all, 5, 10 etc.
        
        for i = 1:MMCsorted.Cnum_psthtrials
            
            tmp2C = MMCsorted.Ctrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2C}) % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.Ctrials(i)}(j) - binsize;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.Ctrials(i)}(j) + binsize;
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.Cpsth           = zeros(length(edges),1);
                MMCsorted.Cspikerate      = zeros(length(edges),1);
                MMCsorted.Cnormspikerate  = zeros(length(edges),1);
                MMCsorted.Chistspiketime  = zeros(length(edges),1);
                
                tmp1C{i}{j}         = find(MMCsorted.Cspiketrialid == MMCsorted.Ctrials(i));
                MMCsorted.Chistspiketime     = histc((MMCsorted.Cspiketime(tmp1C{i}{j})),edges);
                
                if size(MMCsorted.Cpsth,1)==size(MMCsorted.Chistspiketime,1)
                    MMCsorted.Cpsth          = MMCsorted.Cpsth + MMCsorted.Chistspiketime;
                else
                    MMCsorted.Cpsth          = MMCsorted.Cpsth + MMCsorted.Chistspiketime';
                end
                
                MMCsorted.Cspikerate              = 1000/length(MMCsorted.Cpsth)*sum(MMCsorted.Cpsth);
                MMCsorted.Cclickpsth{i}{j}        = MMCsorted.Cpsth;
                MMCsorted.Cclickspikerate{i}{j}   = MMCsorted.Cspikerate;
                
            end
        end; clear i j
        
        
        % ClickPSTH - RC
        tmp1RC = {};
        MMCsorted.RCnum_psthtrials = length(MMCsorted.RCtrials);
        
        for i = 1:MMCsorted.RCnum_psthtrials
            
            tmp2RC = MMCsorted.RCtrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2RC}) % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RCtrials(i)}(j) - binsize;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RCtrials(i)}(j) + binsize;
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.RCpsth           = zeros(length(edges),1);
                MMCsorted.RCspikerate      = zeros(length(edges),1);
                MMCsorted.RCnormspikerate  = zeros(length(edges),1);
                MMCsorted.RChistspiketime  = zeros(length(edges),1);
                
                tmp1RC{i}{j}         = find(MMCsorted.RCspiketrialid == MMCsorted.RCtrials(i));
                MMCsorted.RChistspiketime     = histc( (MMCsorted.RCspiketime(tmp1RC{i}{j}) ) ,edges);
                
                if size(MMCsorted.RCpsth,1)==size(MMCsorted.RChistspiketime,1)
                    MMCsorted.RCpsth          = MMCsorted.RCpsth + MMCsorted.RChistspiketime;
                else
                    MMCsorted.RCpsth          = MMCsorted.RCpsth + MMCsorted.RChistspiketime';
                end
                
                MMCsorted.RCspikerate              = 1000/length(MMCsorted.RCpsth)*sum(MMCsorted.RCpsth);
                MMCsorted.RCclickpsth{i}{j}        = MMCsorted.RCpsth;
                MMCsorted.RCclickspikerate{i}{j}   = MMCsorted.RCspikerate;
                
            end
        end
        
        
        % ClickPSTH - RefRC
        tmp1RefRC = {};
        MMCsorted.RefRCnum_psthtrials = length(MMCsorted.RefRCtrials);
        
        for i = 1:MMCsorted.RefRCnum_psthtrials
            
            tmp2RefRC = MMCsorted.RefRCtrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2RefRC}) % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RefRCtrials(i)}(j) - binsize;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RefRCtrials(i)}(j) + binsize;
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.RefRCpsth           = zeros(length(edges),1);
                MMCsorted.RefRCspikerate      = zeros(length(edges),1);
                MMCsorted.RefRCnormspikerate  = zeros(length(edges),1);
                MMCsorted.RefRChistspiketime  = zeros(length(edges),1);
                
                tmp1RefRC{i}{j}         = find(MMCsorted.RefRCspiketrialid == MMCsorted.RefRCtrials(i));
                MMCsorted.RefRChistspiketime     = histc( (MMCsorted.RefRCspiketime(tmp1RefRC{i}{j}) ) ,edges);
                
                if size(MMCsorted.RefRCpsth,1)==size(MMCsorted.RefRChistspiketime,1)
                    MMCsorted.RefRCpsth          = MMCsorted.RefRCpsth + MMCsorted.RefRChistspiketime;
                else
                    MMCsorted.RefRCpsth          = MMCsorted.RefRCpsth + MMCsorted.RefRChistspiketime';
                end
                
                MMCsorted.RefRCspikerate              = 1000/length(MMCsorted.RefRCpsth)*sum(MMCsorted.RefRCpsth);
                MMCsorted.RefRCclickpsth{i}{j}        = MMCsorted.RefRCpsth;
                MMCsorted.RefRCclickspikerate{i}{j}   = MMCsorted.RefRCspikerate;
                
            end
        end
        
        
        %%%% avg click PSTH analysis %%%%%
        
        % C
        MMCsorted.Cclickavgpsth = [];
        
        for i=1:MMCsorted.Cnum_psthtrials
            for j=1:length(MMCsorted.Cclickpsth{i}) % no. of clicks per trial
                
                tmpCavg = [];
                for k=1:length(MMCsorted.Cclickpsth{i}{j}) % no. of bins
                    tmpCavg = [tmpCavg MMCsorted.Cclickpsth{i}{j}(k)];
                end
                MMCsorted.Cclickavgpsth = [MMCsorted.Cclickavgpsth; tmpCavg];
            end
        end
        MMCsorted.Cclickavgpsthall = mean(MMCsorted.Cclickavgpsth);
        
        % RC
        MMCsorted.RCclickavgpsth = [];
        
        for i=1:MMCsorted.RCnum_psthtrials
            for j=1:length(MMCsorted.RCclickpsth{i}) % no. of clicks per trial
                
                tmpRCavg = [];
                for k=1:length(MMCsorted.RCclickpsth{i}{j}) % no. of bins
                    tmpRCavg = [tmpRCavg MMCsorted.RCclickpsth{i}{j}(k)];
                end
                MMCsorted.RCclickavgpsth = [MMCsorted.RCclickavgpsth; tmpRCavg];
            end
        end
        MMCsorted.RCclickavgpsthall = mean(MMCsorted.RCclickavgpsth);
        
        
        % RefRC
        
        MMCsorted.RefRCclickavgpsth = [];
        
        for i=1:MMCsorted.RefRCnum_psthtrials
            for j=1:length(MMCsorted.RefRCclickpsth{i}) % no. of clicks per trial
                
                tmpRefRCavg1 = [];
                for k=1:length(MMCsorted.RefRCclickpsth{i}{j})
                    tmpRefRCavg1 = [tmpRefRCavg1 MMCsorted.RefRCclickpsth{i}{j}(k)];
                end
                MMCsorted.RefRCclickavgpsth = [MMCsorted.RefRCclickavgpsth; tmpRefRCavg1];
            end
        end
        
        MMCsorted.RefRCclickavgpsthall = mean(MMCsorted.RefRCclickavgpsth);
        
        
        % plot
        
        xbins = [-binsize:0.001:binsize];
        fig5  = figure; hold on;
        errorbar(xbins,MMCsorted.Cclickavgpsthall,std(MMCsorted.Cclickavgpsth)/sqrt(size(MMCsorted.Cclickavgpsth,1)),'ko-','LineWidth',2,'MarkerSize',4);
        errorbar(xbins,MMCsorted.RCclickavgpsthall,std(MMCsorted.RCclickavgpsth)/sqrt(size(MMCsorted.RCclickavgpsth,1)),'bo-','LineWidth',2,'MarkerSize',4);
        errorbar(xbins, MMCsorted.RefRCclickavgpsthall,std(MMCsorted.RefRCclickavgpsth)/sqrt(size(MMCsorted.RefRCclickavgpsth,1)),'ro-','LineWidth',2,'MarkerSize',4);
        xlabel('Peri-stimulus Time (s)','FontSize',12,'FontWeight','bold')
        ylabel('Click-evoked mean spike count','FontSize',12,'FontWeight','Bold')
        legend('C','RC','RefRC','Location','Best')
        set(gca,'XTick', [-0.02 -0.01 0 0.01 0.02]);
        saveas(fig5,[MMCsorted.savepath filesep 'ClickEvoked/PSTH/Sorted_5_ClickPSTH_' MMCsorted.filename(10:end) '.tiff']);
    end
    
    
    %%% Average click-evoked spike rate across trials %%%
    
    if(plotclickspikeratepertrial)
        
        %C
        MMCsorted.Cclickavgspikerate = [];
        MMCsorted.Cclickavgspikeratesem = [];
        
        for i=1:length(MMCsorted.Ctrials)
            tmpCavg = [];
            for j=1:length(MMCsorted.Cclickspikerate{i}) % no. of clicks per trial
                tmpCavg = [tmpCavg MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerate = [MMCsorted.Cclickavgspikerate mean(tmpCavg)];
            MMCsorted.Cclickavgspikeratesem = [MMCsorted.Cclickavgspikeratesem std(tmpCavg)/sqrt(length(MMCsorted.Cclickspikerate{i}))];
            
        end
        
        %RC
        MMCsorted.RCclickavgspikerate = [];
        MMCsorted.RCclickavgspikeratesem = [];
        
        for i=1:length(MMCsorted.RCtrials)
            tmpRCavg = [];
            for j=1:length(MMCsorted.RCclickspikerate{i}) % no. of clicks per trial
                tmpRCavg = [tmpRCavg MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerate = [MMCsorted.RCclickavgspikerate mean(tmpRCavg)];
            MMCsorted.RCclickavgspikeratesem = [MMCsorted.RCclickavgspikeratesem std(tmpRCavg)/sqrt(length(MMCsorted.RCclickspikerate{i}))];
            
        end
        
        % RefRC
        
        MMCsorted.RefRCclickavgspikerate = [];
        MMCsorted.RefRCclickavgspikeratesem = [];
        
        for i=1:length(MMCsorted.RefRCtrials) %
            tmpCavg = [];
            for j=1:length(MMCsorted.RefRCclickspikerate{i}) % no. of clicks per trial
                tmpCavg = [tmpCavg MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerate    = [MMCsorted.RefRCclickavgspikerate mean(tmpCavg)];
            MMCsorted.RefRCclickavgspikeratesem = [MMCsorted.RefRCclickavgspikeratesem std(tmpCavg)/sqrt(length(MMCsorted.RefRCclickspikerate{i}))];
            
        end
        
        % plot
        Cbins = [1:length(MMCsorted.Ctrials)];
        RCbins = [1:length(MMCsorted.RCtrials)];
        RefRCbins = [1:length(MMCsorted.RefRCtrials)];
        
        fig6 = figure; hold on;
        errorbar(MMCsorted.Cclickavgspikerate,MMCsorted.Cclickavgspikeratesem,'ko-','LineWidth',2,'MarkerSize',4);
        errorbar(MMCsorted.RCclickavgspikerate,MMCsorted.RCclickavgspikeratesem,'bo-','LineWidth',2,'MarkerSize',4);
        errorbar(MMCsorted.RefRCclickavgspikerate,MMCsorted.RefRCclickavgspikeratesem,'ro-','LineWidth',2,'MarkerSize',4);
        legend('C','RC','RefRC','Location','Best');
        xlabel('Trials','FontSize',12,'FontWeight','bold');
        xlim([0 51]);
        ylabel('Click-evoked spike rate (Hz)', 'FontSize',12,'FontWeight','bold');
        %         title(['ClickSpikeRateTrial: ' MMCdata.infoname{1}(1:end-4) 'Elec' num2str(MMC.elecID) 'MaxGap' num2str(MMC.MaxGap*1000,'%0.2dms')])
        saveas(fig6,[MMCsorted.savepath filesep 'ClickEvoked/SpikeRateTrial/Sorted_6_ClickRateTrial_' MMCsorted.filename(10:end) '.tiff']);
        
    end
    
    
    %%% average click-evoked spike rate per click across all trials %%%
    
    if(plotclickspikerateperclick)
        
        % min. clicks
        tmp = [];
        Ctmp = MMCsorted.clicktimes(MMCsorted.Ctrials);
        
        for i=1:length(MMCsorted.Ctrials)
            tmp = [tmp length(Ctmp{i})];
        end
        MMCsorted.Cnum_clicks = min(tmp); % to average
        
        tmp = []; RCtmp = MMCsorted.clicktimes(MMCsorted.RCtrials);
        for i=1:length(MMCsorted.RCtrials)
            tmp = [tmp length(RCtmp{i})];
        end
        MMCsorted.RCnum_clicks = min(tmp);
        
        x=MMCsorted.clicktimes(MMCsorted.RefRCtrials);
        MMCsorted.RefRCnum_clicks = length(x{1});
        MMCsorted.num_common_clicks = min([MMCsorted.Cnum_clicks MMCsorted.RCnum_clicks MMCsorted.RefRCnum_clicks]);
        
        %C
        MMCsorted.Cclickavgspikerateperclick = [];
        MMCsorted.Cclickavgspikerateperclicksem = [];
        MMCsorted.Cclickavgspikerateperclickfirst5 = [];
        MMCsorted.Cclickavgspikerateperclickfirst5sem = [];
        MMCsorted.Cclickavgspikerateperclickfirst = [];
        MMCsorted.Cclickavgspikerateperclickfirstsem = [];
        MMCsorted.Cclickavgspikerateperclicklast5 = [];
        MMCsorted.Cclickavgspikerateperclicklast5sem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            
            % average across all trials
            tmpCavg = [];
            
            for i=1:length(MMCsorted.Ctrials)
                tmpCavg    = [tmpCavg MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclick    = [MMCsorted.Cclickavgspikerateperclick; mean(tmpCavg)];
            MMCsorted.Cclickavgspikerateperclicksem = [MMCsorted.Cclickavgspikerateperclicksem; std(tmpCavg)/sqrt(length(tmpCavg))];
            
            % average across first 5 trials
            tmpCavgfirst5    = [];
            
            for i=1:5
                tmpCavgfirst5    = [tmpCavgfirst5 MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclickfirst5    = [MMCsorted.Cclickavgspikerateperclickfirst5; mean(tmpCavgfirst5)];
            MMCsorted.Cclickavgspikerateperclickfirst5sem = [MMCsorted.Cclickavgspikerateperclickfirst5sem; std(tmpCavgfirst5)/sqrt(length(tmpCavgfirst5))];
            
            % average across first trial
            tmpCavgfirst    = [];
            
            for i=1
                tmpCavgfirst    = [tmpCavgfirst MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclickfirst    = [MMCsorted.Cclickavgspikerateperclickfirst; mean(tmpCavgfirst)];
            MMCsorted.Cclickavgspikerateperclickfirstsem = [MMCsorted.Cclickavgspikerateperclickfirstsem; std(tmpCavgfirst)/sqrt(length(tmpCavgfirst))];
            
            % average across last 5 trials
            tmpCavglast5    = [];
            
            for i=length(MMCsorted.Ctrials)-4:length(MMCsorted.Ctrials)
                tmpCavglast5    = [tmpCavglast5 MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclicklast5    = [MMCsorted.Cclickavgspikerateperclicklast5; mean(tmpCavglast5)];
            MMCsorted.Cclickavgspikerateperclicklast5sem = [MMCsorted.Cclickavgspikerateperclicklast5sem; std(tmpCavglast5)/sqrt(length(tmpCavglast5))];
            
        end
        
        
        % RC
        MMCsorted.RCclickavgspikerateperclick = [];
        MMCsorted.RCclickavgspikerateperclicksem = [];
        MMCsorted.RCclickavgspikerateperclickfirst5 = [];
        MMCsorted.RCclickavgspikerateperclickfirst5sem = [];
        MMCsorted.RCclickavgspikerateperclickfirst = [];
        MMCsorted.RCclickavgspikerateperclickfirstsem = [];
        MMCsorted.RCclickavgspikerateperclicklast5 = [];
        MMCsorted.RCclickavgspikerateperclicklast5sem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            
            % average across all trials
            tmpRCavg1 = [];
            
            for i=1:length(MMCsorted.RCtrials)
                tmpRCavg1 = [tmpRCavg1 MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclick    = [MMCsorted.RCclickavgspikerateperclick; mean(tmpRCavg1)];
            MMCsorted.RCclickavgspikerateperclicksem = [MMCsorted.RCclickavgspikerateperclicksem; std(tmpRCavg1)/sqrt(length(tmpRCavg1))];
            
            % average across first 5 trials
            tmpRCavgfirst5    = [];
            
            for i=1:5
                tmpRCavgfirst5    = [tmpRCavgfirst5 MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclickfirst5    = [MMCsorted.RCclickavgspikerateperclickfirst5; mean(tmpRCavgfirst5)];
            MMCsorted.RCclickavgspikerateperclickfirst5sem = [MMCsorted.RCclickavgspikerateperclickfirst5sem; std(tmpRCavgfirst5)/sqrt(length(tmpRCavgfirst5))];
            
            % average across first trial
            tmpRCavgfirst    = [];
            
            for i=1
                tmpRCavgfirst    = [tmpRCavgfirst MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclickfirst    = [MMCsorted.RCclickavgspikerateperclickfirst; mean(tmpRCavgfirst)];
            MMCsorted.RCclickavgspikerateperclickfirstsem = [MMCsorted.RCclickavgspikerateperclickfirstsem; std(tmpRCavgfirst)/sqrt(length(tmpRCavgfirst))];
            
            % average across last 5 trials
            tmpRCavglast5    = [];
            
            for i=length(MMCsorted.RCtrials)-4:length(MMCsorted.RCtrials)
                tmpRCavglast5    = [tmpRCavglast5 MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclicklast5    = [MMCsorted.RCclickavgspikerateperclicklast5; mean(tmpRCavglast5)];
            MMCsorted.RCclickavgspikerateperclicklast5sem = [MMCsorted.RCclickavgspikerateperclicklast5sem; std(tmpRCavglast5)/sqrt(length(tmpRCavglast5))];
            
        end
        
        
        % RefRC
        
        MMCsorted.RefRCclickavgspikerateperclick = [];
        MMCsorted.RefRCclickavgspikerateperclicksem = [];
        MMCsorted.RefRCclickavgspikerateperclickfirst5 = [];
        MMCsorted.RefRCclickavgspikerateperclickfirst5sem = [];
        MMCsorted.RefRCclickavgspikerateperclickfirst = [];
        MMCsorted.RefRCclickavgspikerateperclickfirstsem = [];
        MMCsorted.RefRCclickavgspikerateperclicklast5 = [];
        MMCsorted.RefRCclickavgspikerateperclicklast5sem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            
            % average across all trials
            tmpRefRCavg = [];
            
            for i=1:length(MMCsorted.RefRCtrials)
                tmpRefRCavg = [tmpRefRCavg MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclick    = [MMCsorted.RefRCclickavgspikerateperclick; mean(tmpRefRCavg)];
            MMCsorted.RefRCclickavgspikerateperclicksem = [MMCsorted.RefRCclickavgspikerateperclicksem; std(tmpRefRCavg)/sqrt(length(tmpRefRCavg))];
            
            % average across first 5 trials
            tmpRefRCavgfirst5    = [];
            
            for i=1:5
                tmpRefRCavgfirst5    = [tmpRefRCavgfirst5 MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclickfirst5    = [MMCsorted.RefRCclickavgspikerateperclickfirst5; mean(tmpRefRCavgfirst5)];
            MMCsorted.RefRCclickavgspikerateperclickfirst5sem = [MMCsorted.RefRCclickavgspikerateperclickfirst5sem; std(tmpRefRCavgfirst5)/sqrt(length(tmpRefRCavgfirst5))];
            
            % average across first trial
            tmpRefRCavgfirst    = [];
            
            for i=1
                tmpRefRCavgfirst    = [tmpRefRCavgfirst MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclickfirst    = [MMCsorted.RefRCclickavgspikerateperclickfirst; mean(tmpRefRCavgfirst)];
            MMCsorted.RefRCclickavgspikerateperclickfirstsem = [MMCsorted.RefRCclickavgspikerateperclickfirstsem; std(tmpRefRCavgfirst)/sqrt(length(tmpRefRCavgfirst))];
            
            % average across last 5 trials
            tmpRefRCavglast5    = [];
            
            for i=length(MMCsorted.RefRCtrials)-4:length(MMCsorted.RefRCtrials)
                tmpRefRCavglast5    = [tmpRefRCavglast5 MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclicklast5    = [MMCsorted.RefRCclickavgspikerateperclicklast5; mean(tmpRefRCavglast5)];
            MMCsorted.RefRCclickavgspikerateperclicklast5sem = [MMCsorted.RefRCclickavgspikerateperclicklast5sem; std(tmpRefRCavglast5)/sqrt(length(tmpRefRCavglast5))];
            
        end
        
        % plot
        fig7 = figure;
        
        ax1 = subplot(411); hold on; % first trial
        errorbar(MMCsorted.Cclickavgspikerateperclickfirst,MMCsorted.Cclickavgspikerateperclickfirstsem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclickfirst,MMCsorted.RCclickavgspikerateperclickfirstsem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclickfirst,MMCsorted.RefRCclickavgspikerateperclickfirstsem,'ro-','LineWidth',2, 'MarkerSize',8);
        xlim([0 MMCsorted.num_common_clicks+1]);
        legend('C','RC','RefRC','Location','Best');
        
        ax2 = subplot(412); hold on; % first 5 trials
        errorbar(MMCsorted.Cclickavgspikerateperclickfirst5,MMCsorted.Cclickavgspikerateperclickfirst5sem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclickfirst5,MMCsorted.RCclickavgspikerateperclickfirst5sem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclickfirst5,MMCsorted.RefRCclickavgspikerateperclickfirst5sem,'ro-','LineWidth',2, 'MarkerSize',8);
        xlim([0 MMCsorted.num_common_clicks+1]);
        
        ax3 = subplot(413); hold on; % last 5 trials
        errorbar(MMCsorted.Cclickavgspikerateperclicklast5,MMCsorted.Cclickavgspikerateperclicklast5sem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclicklast5,MMCsorted.RCclickavgspikerateperclicklast5sem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclicklast5,MMCsorted.RefRCclickavgspikerateperclicklast5sem,'ro-','LineWidth',2, 'MarkerSize',8);
        ylabel('Click-evoked mean spike rate per click (Hz)','FontSize',12,'FontWeight','bold');
        xlim([0 MMCsorted.num_common_clicks+1]);
        
        ax4 = subplot(414); hold on; % plot 4 - all trials
        errorbar(MMCsorted.Cclickavgspikerateperclick,MMCsorted.Cclickavgspikerateperclicksem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclick,MMCsorted.RCclickavgspikerateperclicksem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclick,MMCsorted.RefRCclickavgspikerateperclicksem,'ro-','LineWidth',2, 'MarkerSize',8);
        xlim([0 MMCsorted.num_common_clicks+1]);
        xlabel('Click number','FontWeight','bold');
        
        tmpval = [MMCsorted.Cclickavgspikerateperclick MMCsorted.Cclickavgspikerateperclickfirst5 MMCsorted.Cclickavgspikerateperclickfirst,MMCsorted.Cclickavgspikerateperclicklast5,...
            MMCsorted.RCclickavgspikerateperclick MMCsorted.RCclickavgspikerateperclickfirst5 MMCsorted.RCclickavgspikerateperclickfirst,MMCsorted.RCclickavgspikerateperclicklast5,...
            MMCsorted.RefRCclickavgspikerateperclick MMCsorted.RefRCclickavgspikerateperclickfirst5 MMCsorted.RefRCclickavgspikerateperclickfirst,MMCsorted.RefRCclickavgspikerateperclicklast5];
        linkaxes([ax1,ax2,ax3,ax4],'y');
        ax4.YLim = ([0 max(max(tmpval))]);
        saveas(fig7,[MMCsorted.savepath filesep 'ClickEvoked/SpikeRateClick/Sorted_7_ClickRateTrial_' MMCsorted.filename(10:end) '.tiff']);
    end
    
    %%% average click-evoked spike rate per click histogram as a function of click interval %%%
    
    if(plotclickintervalpsth)
        
        % ClickPSTH - C
        
        bin      = 0.001;
        tmp1C    = {};
        MMCsorted.Cnum_psthtrials = length(MMCsorted.Ctrials);
        
        for i = 1:MMCsorted.Cnum_psthtrials
            
            tmp2C = MMCsorted.Ctrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2C})-1 % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.Ctrials(i)}(j) - bin;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.Ctrials(i)}(j+1);
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.Cintpsth           = zeros(length(edges),1);
                MMCsorted.Cintspikerate      = zeros(length(edges),1);
                MMCsorted.Cintnormspikerate  = zeros(length(edges),1);
                MMCsorted.Cinthistspiketime  = zeros(length(edges),1);
                
                tmp1C{i}{j}         = find(MMCsorted.Cspiketrialid == MMCsorted.Ctrials(i));
                MMCsorted.Cinthistspiketime     = histc((MMCsorted.Cspiketime(tmp1C{i}{j})),edges);
                
                if size(MMCsorted.Cintpsth,1)==size(MMCsorted.Cinthistspiketime,1)
                    MMCsorted.Cintpsth          = MMCsorted.Cintpsth + MMCsorted.Cinthistspiketime;
                else
                    MMCsorted.Cintpsth          = MMCsorted.Cintpsth + MMCsorted.Cinthistspiketime';
                end
                
                MMCsorted.Cintspikerate              = 1000/length(MMCsorted.Cintpsth)*sum(MMCsorted.Cintpsth);
                MMCsorted.Cintclickpsth{i}{j}        = MMCsorted.Cintpsth;
                MMCsorted.Cintclickspikerate{i}{j}   = MMCsorted.Cintspikerate;
                
            end
        end; clear i j
        
        
        % ClickPSTH - RC
        tmp1RC = {};
        MMCsorted.RCnum_psthtrials = length(MMCsorted.RCtrials);
        
        for i = 1:MMCsorted.RCnum_psthtrials
            
            tmp2RC = MMCsorted.RCtrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2RC})-1 % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RCtrials(i)}(j) - bin;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RCtrials(i)}(j+1);
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.RCintpsth           = zeros(length(edges),1);
                MMCsorted.RCintspikerate      = zeros(length(edges),1);
                MMCsorted.RCintnormspikerate  = zeros(length(edges),1);
                MMCsorted.RCinthistspiketime  = zeros(length(edges),1);
                
                tmp1RC{i}{j}         = find(MMCsorted.RCspiketrialid == MMCsorted.RCtrials(i));
                MMCsorted.RCinthistspiketime     = histc( (MMCsorted.RCspiketime(tmp1RC{i}{j}) ) ,edges);
                
                if size(MMCsorted.RCintpsth,1)==size(MMCsorted.RCinthistspiketime,1)
                    MMCsorted.RCintpsth          = MMCsorted.RCintpsth + MMCsorted.RCinthistspiketime;
                else
                    MMCsorted.RCintpsth          = MMCsorted.RCintpsth + MMCsorted.RCinthistspiketime';
                end
                
                MMCsorted.RCintspikerate              = 1000/length(MMCsorted.RCintpsth)*sum(MMCsorted.RCintpsth);
                MMCsorted.RCintclickpsth{i}{j}        = MMCsorted.RCintpsth;
                MMCsorted.RCintclickspikerate{i}{j}   = MMCsorted.RCintspikerate;
                
            end
        end
        
        
        % ClickPSTH - RefRC
        tmp1RefRC = {};
        MMCsorted.RefRCnum_psthtrials = length(MMCsorted.RefRCtrials);
        
        for i = 1:MMCsorted.RefRCnum_psthtrials
            
            tmp2RefRC = MMCsorted.RefRCtrials(i);
            
            for j = 1:length(MMCsorted.clicktimes{tmp2RefRC})-1 % no. of C clicks
                
                baseline_pre  = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RefRCtrials(i)}(j) - bin;
                baseline_post = MMCsorted.clickunsorted.start + MMCsorted.clicktimes{MMCsorted.RefRCtrials(i)}(j+1);
                edges         = [baseline_pre : bin: baseline_post];
                
                MMCsorted.RefRCintpsth           = zeros(length(edges),1);
                MMCsorted.RefRCintspikerate      = zeros(length(edges),1);
                MMCsorted.RefRCintnormspikerate  = zeros(length(edges),1);
                MMCsorted.RefRCinthistspiketime  = zeros(length(edges),1);
                
                tmp1RefRC{i}{j}         = find(MMCsorted.RefRCspiketrialid == MMCsorted.RefRCtrials(i));
                MMCsorted.RefRCinthistspiketime     = histc( (MMCsorted.RefRCspiketime(tmp1RefRC{i}{j}) ) ,edges);
                
                if size(MMCsorted.RefRCintpsth,1)==size(MMCsorted.RefRCinthistspiketime,1)
                    MMCsorted.RefRCintpsth          = MMCsorted.RefRCintpsth + MMCsorted.RefRCinthistspiketime;
                else
                    MMCsorted.RefRCintpsth          = MMCsorted.RefRCintpsth + MMCsorted.RefRCinthistspiketime';
                end
                
                MMCsorted.RefRCintspikerate              = 1000/length(MMCsorted.RefRCintpsth)*sum(MMCsorted.RefRCintpsth);
                MMCsorted.RefRCintclickpsth{i}{j}        = MMCsorted.RefRCintpsth;
                MMCsorted.RefRCintclickspikerate{i}{j}   = MMCsorted.RefRCintspikerate;
                
            end
        end
        
        
        % collect interval sizes and spike rates for corresponding interval
        
        %C
        for i1 = 1:length(MMCsorted.Ctrials)
            tmpintC = MMCsorted.clicktimes(MMCsorted.Ctrials);
            MMCsorted.Cint{i1} = diff(tmpintC{i1})*1000;
            MMCsorted.Cintspikeratehist{i1} = [];
            
            for j1 = 1:length(MMCsorted.Cint{i1})
                MMCsorted.Cintspikeratehist{i1}  = [MMCsorted.Cintspikeratehist{i1} MMCsorted.Cintclickspikerate{i1}{j1}];
            end
        end
        
        %RC
        for i1 = 1:length(MMCsorted.RCtrials)
            tmpintRC = MMCsorted.clicktimes(MMCsorted.RCtrials);
            MMCsorted.RCint{i1} = diff(tmpintRC{i1})*1000;
            MMCsorted.RCintspikeratehist{i1} = [];
            
            for j1 = 1:length(MMCsorted.RCint{i1})
                MMCsorted.RCintspikeratehist{i1}  = [MMCsorted.RCintspikeratehist{i1} MMCsorted.RCintclickspikerate{i1}{j1}];
            end
        end
        
        %RefRC
        for i1 = 1:length(MMCsorted.RefRCtrials)
            tmpintRefRC = MMCsorted.clicktimes(MMCsorted.RefRCtrials);
            MMCsorted.RefRCint{i1} = diff(tmpintRefRC{i1})*1000;
            MMCsorted.RefRCintspikeratehist{i1} = [];
            
            for j1 = 1:length(MMCsorted.RefRCint{i1})
                MMCsorted.RefRCintspikeratehist{i1}  = [MMCsorted.RefRCintspikeratehist{i1} MMCsorted.RefRCintclickspikerate{i1}{j1}];
            end
        end
        
        % collect all and plot scatter plot and compute correlations
        
        % C
        MMCsorted.Cintall = [];
        MMCsorted.Cintspikeratehistall = [];
        MMCsorted.rhoCall = [];
        MMCsorted.pCall   = [];
        
        for i2 = 1:length(MMCsorted.Ctrials)
            MMCsorted.Cintall = [MMCsorted.Cintall MMCsorted.Cint{i2}];
            MMCsorted.Cintspikeratehistall = [MMCsorted.Cintspikeratehistall MMCsorted.Cintspikeratehist{i2}];
            [MMCsorted.rhoC{i2},MMCsorted.pC{i2}]=corrcoef(MMCsorted.Cint{i2},MMCsorted.Cintspikeratehist{i2});
            MMCsorted.rhoCall         = [MMCsorted.rhoCall MMCsorted.rhoC{i2}(1,2)];
            MMCsorted.pCall         = [MMCsorted.pCall MMCsorted.pC{i2}(1,2)];
            
        end
        
        % RC
        MMCsorted.RCintall = [];
        MMCsorted.RCintspikeratehistall = [];
        MMCsorted.rhoRCall = [];
        MMCsorted.pRCall   = [];
        
        for i2 = 1:length(MMCsorted.RCtrials)
            MMCsorted.RCintall = [MMCsorted.RCintall MMCsorted.RCint{i2}];
            MMCsorted.RCintspikeratehistall = [MMCsorted.RCintspikeratehistall MMCsorted.RCintspikeratehist{i2}];
            [MMCsorted.rhoRC{i2},MMCsorted.pRC{i2}]=corrcoef(MMCsorted.RCint{i2},MMCsorted.RCintspikeratehist{i2});
            MMCsorted.rhoRCall         = [MMCsorted.rhoRCall MMCsorted.rhoRC{i2}(1,2)];
            MMCsorted.pRCall         = [MMCsorted.pRCall MMCsorted.pRC{i2}(1,2)];
        end
        
        % RefRC
        MMCsorted.RefRCintall = [];
        MMCsorted.RefRCintspikeratehistall = [];
        MMCsorted.rhoRefRCall = [];
        MMCsorted.pRefRCall   = [];
        
        for i2 = 1:length(MMCsorted.RefRCtrials)
            MMCsorted.RefRCintall = [MMCsorted.RefRCintall MMCsorted.RefRCint{i2}];
            MMCsorted.RefRCintspikeratehistall = [MMCsorted.RefRCintspikeratehistall MMCsorted.RefRCintspikeratehist{i2}];
            [MMCsorted.rhoRefRC{i2},MMCsorted.pRefRC{i2}]=corrcoef(MMCsorted.RefRCint{i2},MMCsorted.RefRCintspikeratehist{i2});
            MMCsorted.rhoRefRCall         = [MMCsorted.rhoRefRCall MMCsorted.rhoRefRC{i2}(1,2)];
            MMCsorted.pRefRCall         = [MMCsorted.pRefRCall MMCsorted.pRefRC{i2}(1,2)];
        end
        
        
        [MMCsorted.rhoRC,MMCsorted.pRC]=corrcoef(MMCsorted.RCintall,MMCsorted.RCintspikeratehistall);
        [MMCsorted.rhoRefRC,MMCsorted.pRefRC]=corrcoef(MMCsorted.RefRCintall,MMCsorted.RefRCintspikeratehistall);
        
    end
    
    
    %% Save individual data
    
    if(savedataflag)
        save([MMCsorted.savepath 'Results' filesep MMCsorted.filename(10:end) '_' num2str(MMCsorted.maxgap*1000,'%0.2dmsMaxGap.mat')],'MMCsorted','MMCsortedall');
    end
    
end


%% Do group analysis

if(dogroupanalysis)
    
    % 1 avgspikeratebaseline
    % 2 avgspikerateclick
    % 3 avgspikeratepostclick
    % 4 avgspikerateclickrep1
    % 5 avgspikerateclickrep2
    % 6 avgspikerateclickrep3
    % 7 clickavgspikerate
    % 8 clickavgspikerateperclickfirst5
    % 9 clickavgspikerateperclickfirst
    % 10 clickavgspikerateperclicklast5
    
    num_cells_avg = 127;
    
    MMCsortedavg.Cavgspikeratebaseline  = [];
    MMCsortedavg.Cavgspikerateclick     = [];
    MMCsortedavg.Cavgspikeratepostclick = [];
    MMCsortedavg.Cavgspikerateclickrep1 = [];
    MMCsortedavg.Cavgspikerateclickrep2 = [];
    MMCsortedavg.Cavgspikerateclickrep3 = [];
    MMCsortedavg.Cclickavgspikerate     = [];
    MMCsortedavg.Cclickavgspikerateperclickfirst5 = [];
    MMCsortedavg.Cclickavgspikerateperclickfirst  = [];
    
    MMCsortedavg.RCavgspikeratebaseline  = [];
    MMCsortedavg.RCavgspikerateclick     = [];
    MMCsortedavg.RCavgspikeratepostclick = [];
    MMCsortedavg.RCavgspikerateclickrep1 = [];
    MMCsortedavg.RCavgspikerateclickrep2 = [];
    MMCsortedavg.RCavgspikerateclickrep3 = [];
    MMCsortedavg.RCclickavgspikerate     = [];
    MMCsortedavg.RCclickavgspikerateperclickfirst5 = [];
    MMCsortedavg.RCclickavgspikerateperclickfirst  = [];
    
    MMCsortedavg.RefRCavgspikeratebaseline  = [];
    MMCsortedavg.RefRCavgspikerateclick     = [];
    MMCsortedavg.RefRCavgspikeratepostclick = [];
    MMCsortedavg.RefRCavgspikerateclickrep1 = [];
    MMCsortedavg.RefRCavgspikerateclickrep2 = [];
    MMCsortedavg.RefRCavgspikerateclickrep3 = [];
    MMCsortedavg.RefRCclickavgspikerate     = [];
    MMCsortedavg.RefRCclickavgspikerateperclickfirst5 = [];
    MMCsortedavg.RefRCclickavgspikerateperclickfirst  = [];
    
%     load pp; for jj=pp % pp = find(rate > threshold)
                for jj=1:num_cells_avg
        
        MMCsortedavg.Cavgspikeratebaseline      = [MMCsortedavg.Cavgspikeratebaseline nanmean(MMCsortedall.Cavgspikeratebaseline{jj})];
        MMCsortedavg.RCavgspikeratebaseline     = [MMCsortedavg.RCavgspikeratebaseline nanmean(MMCsortedall.RCavgspikeratebaseline{jj})];
        MMCsortedavg.RefRCavgspikeratebaseline  = [MMCsortedavg.RefRCavgspikeratebaseline nanmean(MMCsortedall.RefRCavgspikeratebaseline{jj})];
        
        MMCsortedavg.Cavgspikerateclick         = [MMCsortedavg.Cavgspikerateclick nanmean(MMCsortedall.Cavgspikerateclick{jj})];
        MMCsortedavg.RCavgspikerateclick        = [MMCsortedavg.RCavgspikerateclick nanmean(MMCsortedall.RCavgspikerateclick{jj})];
        MMCsortedavg.RefRCavgspikerateclick     = [MMCsortedavg.RefRCavgspikerateclick nanmean(MMCsortedall.RefRCavgspikerateclick{jj})];
        
        MMCsortedavg.Cavgspikeratepostclick     = [MMCsortedavg.Cavgspikeratepostclick nanmean(MMCsortedall.Cavgspikeratepostclick{jj})];
        MMCsortedavg.RCavgspikeratepostclick    = [MMCsortedavg.RCavgspikeratepostclick nanmean(MMCsortedall.RCavgspikeratepostclick{jj})];
        MMCsortedavg.RefRCavgspikeratepostclick = [MMCsortedavg.RefRCavgspikeratepostclick nanmean(MMCsortedall.RefRCavgspikeratepostclick{jj})];
        
        MMCsortedavg.Cavgspikerateclickrep1     = [MMCsortedavg.Cavgspikerateclickrep1 nanmean(MMCsortedall.Cavgspikerateclickrep1{jj})];
        MMCsortedavg.RCavgspikerateclickrep1    = [MMCsortedavg.RCavgspikerateclickrep1 nanmean(MMCsortedall.RCavgspikerateclickrep1{jj})];
        MMCsortedavg.RefRCavgspikerateclickrep1 = [MMCsortedavg.RefRCavgspikerateclickrep1 nanmean(MMCsortedall.RefRCavgspikerateclickrep1{jj})];
        
        MMCsortedavg.Cavgspikerateclickrep2     = [MMCsortedavg.Cavgspikerateclickrep2 nanmean(MMCsortedall.Cavgspikerateclickrep2{jj})];
        MMCsortedavg.RCavgspikerateclickrep2    = [MMCsortedavg.RCavgspikerateclickrep2 nanmean(MMCsortedall.RCavgspikerateclickrep2{jj})];
        MMCsortedavg.RefRCavgspikerateclickrep2 = [MMCsortedavg.RefRCavgspikerateclickrep2 nanmean(MMCsortedall.RefRCavgspikerateclickrep2{jj})];
        
        MMCsortedavg.Cavgspikerateclickrep3     = [MMCsortedavg.Cavgspikerateclickrep3 nanmean(MMCsortedall.Cavgspikerateclickrep3{jj})];
        MMCsortedavg.RCavgspikerateclickrep3    = [MMCsortedavg.RCavgspikerateclickrep3 nanmean(MMCsortedall.RCavgspikerateclickrep3{jj})];
        MMCsortedavg.RefRCavgspikerateclickrep3 = [MMCsortedavg.RefRCavgspikerateclickrep3 nanmean(MMCsortedall.RefRCavgspikerateclickrep3{jj})];
        
        MMCsortedavg.Cclickavgspikerate     = [MMCsortedavg.Cclickavgspikerate nanmean(MMCsortedall.Cclickavgspikerate{jj})];
        MMCsortedavg.RCclickavgspikerate    = [MMCsortedavg.RCclickavgspikerate nanmean(MMCsortedall.RCclickavgspikerate{jj})];
        MMCsortedavg.RefRCclickavgspikerate = [MMCsortedavg.RefRCclickavgspikerate nanmean(MMCsortedall.RefRCclickavgspikerate{jj})];
        
        MMCsortedavg.Cclickavgspikerateperclickfirst5     = [MMCsortedavg.Cclickavgspikerateperclickfirst5 nanmean(MMCsortedall.Cclickavgspikerateperclickfirst5{jj})];
        MMCsortedavg.RCclickavgspikerateperclickfirst5    = [MMCsortedavg.RCclickavgspikerateperclickfirst5 nanmean(MMCsortedall.RCclickavgspikerateperclickfirst5{jj})];
        MMCsortedavg.RefRCclickavgspikerateperclickfirst5 = [MMCsortedavg.RefRCclickavgspikerateperclickfirst5 nanmean(MMCsortedall.RefRCclickavgspikerateperclickfirst5{jj})];
        
        MMCsortedavg.Cclickavgspikerateperclickfirst     = [MMCsortedavg.Cclickavgspikerateperclickfirst nanmean(MMCsortedall.Cclickavgspikerateperclickfirst{jj})];
        MMCsortedavg.RCclickavgspikerateperclickfirst    = [MMCsortedavg.RCclickavgspikerateperclickfirst nanmean(MMCsortedall.RCclickavgspikerateperclickfirst{jj})];
        MMCsortedavg.RefRCclickavgspikerateperclickfirst = [MMCsortedavg.RefRCclickavgspikerateperclickfirst nanmean(MMCsortedall.RefRCclickavgspikerateperclickfirst{jj})];
        
        %            MMCsortedavg.Cclickavgspikerateperclicklast5     = [MMCsortedavg.Cclickavgspikerateperclicklast5 nanmean(MMCsortedall.Cclickavgspikerateperclicklast5{jj})];
        %            MMCsortedavg.RCclickavgspikerateperclicklast5    = [MMCsortedavg.RCclickavgspikerateperclicklast5 nanmean(MMCsortedall.RCclickavgspikerateperclicklast5{jj})];
        %            MMCsortedavg.RefRCclickavgspikerateperclicklast5 = [MMCsortedavg.RefRCclickavgspikerateperclicklast5 nanmean(MMCsortedall.RefRCclickavgspikerateperclicklast5{jj})];
        
    end
    
end
    
    %% do group stats
    
    if(dogroupstats)
        
        %% anova
        [MMCgroupstats.pavgspikeratebaseline,...
            MMCgroupstats.anovaavgspikeratebaseline,...
            MMCgroupstats.statsavgspikeratebaseline] = st_anova1(MMCsortedavg.Cavgspikeratebaseline,...
            MMCsortedavg.RCavgspikeratebaseline,...
            MMCsortedavg.RefRCavgspikeratebaseline);
        
        [MMCgroupstats.pavgspikerateclick,...
            MMCgroupstats.anovaavgspikerateclick,...
            MMCgroupstats.statsavgspikerateclick] = st_anova1(MMCsortedavg.Cavgspikerateclick,...
            MMCsortedavg.RCavgspikerateclick,...
            MMCsortedavg.RefRCavgspikerateclick);
        [MMCgroupstats.pavgspikeratepostclick,...
            MMCgroupstats.anovaavgspikeratepostclick,...
            MMCgroupstats.statsavgspikeratepostclick] = st_anova1(MMCsortedavg.Cavgspikeratepostclick,...
            MMCsortedavg.RCavgspikeratepostclick,...
            MMCsortedavg.RefRCavgspikeratepostclick);
        
        [MMCgroupstats.pavgspikerateclickrep1,...
            MMCgroupstats.anovaavgspikerateclickrep1,...
            MMCgroupstats.statsavgspikerateclickrep1] = st_anova1(MMCsortedavg.Cavgspikerateclickrep1,...
            MMCsortedavg.RCavgspikerateclickrep1,...
            MMCsortedavg.RefRCavgspikerateclickrep1);
        
        [MMCgroupstats.pavgspikerateclickrep2,...
            MMCgroupstats.anovaavgspikerateclickrep2,...
            MMCgroupstats.statsavgspikerateclickrep2] = st_anova1(MMCsortedavg.Cavgspikerateclickrep2,...
            MMCsortedavg.RCavgspikerateclickrep2,...
            MMCsortedavg.RefRCavgspikerateclickrep2);
        
        [MMCgroupstats.pavgspikerateclickrep3,...
            MMCgroupstats.anovaavgspikerateclickrep3,...
            MMCgroupstats.statsavgspikerateclickrep3] = st_anova1(MMCsortedavg.Cavgspikerateclickrep3,...
            MMCsortedavg.RCavgspikerateclickrep3,...
            MMCsortedavg.RefRCavgspikerateclickrep3);
        
        [MMCgroupstats.pclickavgspikerate,...
            MMCgroupstats.anovaclickavgspikerate,...
            MMCgroupstats.statsclickavgspikerate] = st_anova1(MMCsortedavg.Cclickavgspikerate,...
            MMCsortedavg.RCclickavgspikerate,...
            MMCsortedavg.RefRCclickavgspikerate);
        
        [MMCgroupstats.pclickavgspikerateperclickfirst5,...
            MMCgroupstats.anovaclickavgspikerateperclickfirst5,...
            MMCgroupstats.statsclickavgspikerateperclickfirst5] = st_anova1(MMCsortedavg.Cclickavgspikerateperclickfirst5,...
            MMCsortedavg.RCclickavgspikerateperclickfirst5,...
            MMCsortedavg.RefRCclickavgspikerateperclickfirst5);
        
        [MMCgroupstats.pclickavgspikerateperclickfirst,...
            MMCgroupstats.anovaclickavgspikerateperclickfirst,...
            MMCgroupstats.statsclickavgspikerateperclickfirst] = st_anova1(MMCsortedavg.Cclickavgspikerateperclickfirst,...
            MMCsortedavg.RCclickavgspikerateperclickfirst,...
            MMCsortedavg.RefRCclickavgspikerateperclickfirst);
        
        
        %% t-test
        
        [MMCgroupstats.ttestHavgspikeratebaseline,...
            MMCgroupstats.ttestpanovaavgspikeratebaseline,...
            MMCgroupstats.ttestCIavgspikeratebaseline,...
            MMCgroupstats.ttestStatsavgspikeratebaseline] = ttest2(MMCsortedavg.RCavgspikeratebaseline,...
            MMCsortedavg.RefRCavgspikeratebaseline);
        
        [MMCgroupstats.ttestHavgspikerateclick,...
            MMCgroupstats.ttestpanovaavgspikerateclick,...
            MMCgroupstats.ttestCIavgspikerateclick,...
            MMCgroupstats.ttestStatsavgspikerateclick] = ttest2(MMCsortedavg.RCavgspikerateclick,...
            MMCsortedavg.RefRCavgspikerateclick);
        
        
        [MMCgroupstats.ttestHavgspikeratepostclick,...
            MMCgroupstats.ttestpanovaavgspikeratepostclick,...
            MMCgroupstats.ttestCIavgspikeratepostclick,...
            MMCgroupstats.ttestStatsavgspikeratepostclick] = ttest2(MMCsortedavg.RCavgspikeratepostclick,...
            MMCsortedavg.RefRCavgspikeratepostclick);
        
        [MMCgroupstats.ttestHavgspikerateclickrep1,...
            MMCgroupstats.ttestpanovaavgspikerateclickrep1,...
            MMCgroupstats.ttestCIavgspikerateclickrep1,...
            MMCgroupstats.ttestStatsavgspikerateclickrep1] = ttest2(MMCsortedavg.RCavgspikerateclickrep1,...
            MMCsortedavg.RefRCavgspikerateclickrep1);
        
        [MMCgroupstats.ttestHavgspikerateclickrep2,...
            MMCgroupstats.ttestpanovaavgspikerateclickrep2,...
            MMCgroupstats.ttestCIavgspikerateclickrep2,...
            MMCgroupstats.ttestStatsavgspikerateclickrep2] = ttest2(MMCsortedavg.RCavgspikerateclickrep2,...
            MMCsortedavg.RefRCavgspikerateclickrep2);
        
        [MMCgroupstats.ttestHavgspikerateclickrep3,...
            MMCgroupstats.ttestpanovaavgspikerateclickrep3,...
            MMCgroupstats.ttestCIavgspikerateclickrep3,...
            MMCgroupstats.ttestStatsavgspikerateclickrep3] = ttest2(MMCsortedavg.RCavgspikerateclickrep3,...
            MMCsortedavg.RefRCavgspikerateclickrep3);
        
        [MMCgroupstats.ttestHclickavgspikerate,...
            MMCgroupstats.ttestpanovaclickavgspikerate,...
            MMCgroupstats.ttestCIclickavgspikerate,...
            MMCgroupstats.ttestStatsclickavgspikerate] = ttest2(MMCsortedavg.RCclickavgspikerate,...
            MMCsortedavg.RefRCclickavgspikerate);
        
        [MMCgroupstats.ttestHclickavgspikerateperclickfirst5,...
            MMCgroupstats.ttestpanovaclickavgspikerateperclickfirst5,...
            MMCgroupstats.ttestCIclickavgspikerateperclickfirst5,...
            MMCgroupstats.ttestStatsclickavgspikerateperclickfirst5] = ttest2(MMCsortedavg.RCclickavgspikerateperclickfirst5,...
            MMCsortedavg.RefRCclickavgspikerateperclickfirst5);
        
        [MMCgroupstats.ttestHclickavgspikerateperclickfirst,...
            MMCgroupstats.ttestpanovaclickavgspikerateperclickfirst,...
            MMCgroupstats.ttestCIclickavgspikerateperclickfirst,...
            MMCgroupstats.ttestStatsclickavgspikerateperclickfirst] = ttest2(MMCsortedavg.RCclickavgspikerateperclickfirst,...
            MMCsortedavg.RefRCclickavgspikerateperclickfirst);
        
        
    end
    
    
    %% save group data
    
    if(savegroupdataflag)
        save([MMCsorted.savepath 'Results/GroupResults' filesep 'GroupResults.mat'],'MMCsortedall','MMCsortedavg','MMCgroupstats');
    end
    
end
