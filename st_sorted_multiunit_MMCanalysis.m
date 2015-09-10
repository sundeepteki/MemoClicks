function [MMCsorted,y] = st_sorted_multiunit_MMCanalysis(varargin)
% clear all; close all; clc;

%{
function to analyze MMC sorted multi unit data

ST. v1:            07.07.15
Last modified:     07.07.15
Backup:            st_analyze_multiunit_MMCsorted_copy
Last backup saved: 06.07.15

%}

%% Set paths

MMCsorted.datapath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#multiunit/';
MMCsorted.savepath = [MMCsorted.datapath '#analyzed/'];
addpath(MMCsorted.datapath);
addpath(MMCsorted.savepath);

%% analysis plots

flag = 0;
loaddata                    = ~flag;
loadMMCunsorteddata         = ~flag;
getclicktimes               = flag;
getspiketimes               = flag;
plotrasterstim              = flag;
plotspikeratestim           = flag;
plotspikeraterepeat         = flag;
plotpsth                    = flag;
plotclickpsth               = flag;
plotclickspikeratepertrial  = flag;
plotclickspikerateperclick  = flag;
plotclickintervalpsth       = flag;
savedataflag                = flag;


%% load sorted data

tmpdir = dir(MMCsorted.datapath);

for ij = 1:(length(tmpdir)-4)
    
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
        y.RefRC{ij} = MMCsorted.RefRCavgspikerateclick;
        y.RC{ij} = MMCsorted.RCavgspikerateclick;
        y.C{ij} = MMCsorted.Cavgspikerateclick;
        
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
                saveas(fig1,[MMCsorted.savepath filesep 'Raster/Sorted_1_Raster_' MMCsorted.filename(10:end) '.tiff']);
        
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
                saveas(fig2,[MMCsorted.savepath filesep 'PSTH/Sorted_2_PSTH_' MMCsorted.filename(10:end) '.tiff']);
        
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
        MMCsorted.Cclickavgspikerateperclicksecond = [];
        MMCsorted.Cclickavgspikerateperclicksecondsem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            tmpCavg = [];
            tmpCavgsem = [];
            
            % all clicks
            for i=1:length(MMCsorted.Ctrials)
                tmpCavg    = [tmpCavg MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclick    = [MMCsorted.Cclickavgspikerateperclick; mean(tmpCavg)];
            MMCsorted.Cclickavgspikerateperclicksem = [MMCsorted.Cclickavgspikerateperclicksem; std(tmpCavg)/sqrt(length(tmpCavg))];
            
            tmpCavgfirst5    = [];
            tmpCavgsemfirst5 = [];
            
            % first 5 clicks
            for i=1:5
                tmpCavgfirst5    = [tmpCavgfirst5 MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclickfirst5    = [MMCsorted.Cclickavgspikerateperclickfirst5; mean(tmpCavgfirst5)];
            MMCsorted.Cclickavgspikerateperclickfirst5sem = [MMCsorted.Cclickavgspikerateperclickfirst5sem; std(tmpCavgfirst5)/sqrt(length(tmpCavgfirst5))];
            
            tmpCavgsecond    = [];
            tmpCavgsemsecond = [];
            
            % second click
            for i=2
                tmpCavgsecond    = [tmpCavgsecond MMCsorted.Cclickspikerate{i}{j}];
            end
            MMCsorted.Cclickavgspikerateperclicksecond    = [MMCsorted.Cclickavgspikerateperclicksecond; mean(tmpCavgsecond)];
            MMCsorted.Cclickavgspikerateperclicksecondsem = [MMCsorted.Cclickavgspikerateperclicksecondsem; std(tmpCavgsecond)/sqrt(length(tmpCavgsecond))];
        end
        
        
        % RC
        MMCsorted.RCclickavgspikerateperclick = [];
        MMCsorted.RCclickavgspikerateperclicksem = [];
        MMCsorted.RCclickavgspikerateperclickfirst5 = [];
        MMCsorted.RCclickavgspikerateperclickfirst5sem = [];
        MMCsorted.RCclickavgspikerateperclicksecond = [];
        MMCsorted.RCclickavgspikerateperclicksecondsem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            tmpRCavg1 = [];
            
            for i=1:length(MMCsorted.RCtrials)
                tmpRCavg1 = [tmpRCavg1 MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclick    = [MMCsorted.RCclickavgspikerateperclick; mean(tmpRCavg1)];
            MMCsorted.RCclickavgspikerateperclicksem = [MMCsorted.RCclickavgspikerateperclicksem; std(tmpRCavg1)/sqrt(length(tmpRCavg1))];
            
            tmpRCavgfirst5    = [];
            tmpRCavgsemfirst5 = [];
            
            for i=1:5
                tmpRCavgfirst5    = [tmpRCavgfirst5 MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclickfirst5    = [MMCsorted.RCclickavgspikerateperclickfirst5; mean(tmpRCavgfirst5)];
            MMCsorted.RCclickavgspikerateperclickfirst5sem = [MMCsorted.RCclickavgspikerateperclickfirst5sem; std(tmpRCavgfirst5)/sqrt(length(tmpRCavgfirst5))];
            
            tmpRCavgsecond    = [];
            tmpRCavgsemsecond = [];
            
            for i=2
                tmpRCavgsecond    = [tmpRCavgsecond MMCsorted.RCclickspikerate{i}{j}];
            end
            MMCsorted.RCclickavgspikerateperclicksecond    = [MMCsorted.RCclickavgspikerateperclicksecond; mean(tmpRCavgsecond)];
            MMCsorted.RCclickavgspikerateperclicksecondsem = [MMCsorted.RCclickavgspikerateperclicksecondsem; std(tmpRCavgsecond)/sqrt(length(tmpRCavgsecond))];
        end
        
        
        % RefRC
        
        MMCsorted.RefRCclickavgspikerateperclick = [];
        MMCsorted.RefRCclickavgspikerateperclicksem = [];
        MMCsorted.RefRCclickavgspikerateperclickfirst5 = [];
        MMCsorted.RefRCclickavgspikerateperclickfirst5sem = [];
        MMCsorted.RefRCclickavgspikerateperclicksecond = [];
        MMCsorted.RefRCclickavgspikerateperclicksecondsem = [];
        
        for j = 1:MMCsorted.num_common_clicks
            tmpRefRCavg = [];
            
            for i=1:length(MMCsorted.RefRCtrials)
                tmpRefRCavg = [tmpRefRCavg MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclick    = [MMCsorted.RefRCclickavgspikerateperclick; mean(tmpRefRCavg)];
            MMCsorted.RefRCclickavgspikerateperclicksem = [MMCsorted.RefRCclickavgspikerateperclicksem; std(tmpRefRCavg)/sqrt(length(tmpRefRCavg))];
            
            tmpRefRCavgfirst5    = [];
            tmpRefRCavgsemfirst5 = [];
            
            for i=1:5
                tmpRefRCavgfirst5    = [tmpRefRCavgfirst5 MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclickfirst5    = [MMCsorted.RefRCclickavgspikerateperclickfirst5; mean(tmpRefRCavgfirst5)];
            MMCsorted.RefRCclickavgspikerateperclickfirst5sem = [MMCsorted.RefRCclickavgspikerateperclickfirst5sem; std(tmpRefRCavgfirst5)/sqrt(length(tmpRefRCavgfirst5))];
            
            tmpRefRCavgsecond    = [];
            tmpRefRCavgsemsecond = [];
            
            for i=2
                tmpRefRCavgsecond    = [tmpRefRCavgsecond MMCsorted.RefRCclickspikerate{i}{j}];
            end
            MMCsorted.RefRCclickavgspikerateperclicksecond    = [MMCsorted.RefRCclickavgspikerateperclicksecond; mean(tmpRefRCavgsecond)];
            MMCsorted.RefRCclickavgspikerateperclicksecondsem = [MMCsorted.RefRCclickavgspikerateperclicksecondsem; std(tmpRefRCavgsecond)/sqrt(length(tmpRefRCavgsecond))];
        end
        
        % plot
        fig7 = figure;
        ax1 = subplot(313); hold on; % all clicks
        errorbar(MMCsorted.Cclickavgspikerateperclick,MMCsorted.Cclickavgspikerateperclicksem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclick,MMCsorted.RCclickavgspikerateperclicksem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclick,MMCsorted.RefRCclickavgspikerateperclicksem,'ro-','LineWidth',2, 'MarkerSize',8);
        xlim([0 MMCsorted.num_common_clicks+1]);
        xlabel('Click number','FontWeight','bold');
        
        ax2 = subplot(312); hold on; % first 5 clicks
        errorbar(MMCsorted.Cclickavgspikerateperclickfirst5,MMCsorted.Cclickavgspikerateperclickfirst5sem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclickfirst5,MMCsorted.RCclickavgspikerateperclickfirst5sem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclickfirst5,MMCsorted.RefRCclickavgspikerateperclickfirst5sem,'ro-','LineWidth',2, 'MarkerSize',8);
        ylabel('Click-evoked mean spike rate per click (Hz)','FontSize',12,'FontWeight','bold');
        xlim([0 MMCsorted.num_common_clicks+1]);
        
        ax3 = subplot(311); hold on; % second click
        errorbar(MMCsorted.Cclickavgspikerateperclicksecond,MMCsorted.Cclickavgspikerateperclicksecondsem,'ko-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RCclickavgspikerateperclicksecond,MMCsorted.RCclickavgspikerateperclicksecondsem,'bo-','LineWidth',2, 'MarkerSize',8);
        errorbar(MMCsorted.RefRCclickavgspikerateperclicksecond,MMCsorted.RefRCclickavgspikerateperclicksecondsem,'ro-','LineWidth',2, 'MarkerSize',8);
        xlim([0 MMCsorted.num_common_clicks+1]);
        legend('C','RC','RefRC','Location','Best');
        
        tmpval = [MMCsorted.Cclickavgspikerateperclick MMCsorted.Cclickavgspikerateperclickfirst5 MMCsorted.Cclickavgspikerateperclicksecond,...
            MMCsorted.RCclickavgspikerateperclick MMCsorted.RCclickavgspikerateperclickfirst5 MMCsorted.RCclickavgspikerateperclicksecond,...
            MMCsorted.RefRCclickavgspikerateperclick MMCsorted.RefRCclickavgspikerateperclickfirst5 MMCsorted.RefRCclickavgspikerateperclicksecond];
        linkaxes([ax1,ax2,ax3],'y');
        ax3.YLim = ([0 max(max(tmpval))]);
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
        
        C
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
    
    
    %% Save data
    
    if(savedataflag)
        save([MMCsorted.savepath 'Results' filesep MMCsorted.filename(10:end) '_' num2str(MMCsorted.maxgap*1000,'%0.2dmsMaxGap.mat')],'MMCsorted','bin');
    end
    
    
end
