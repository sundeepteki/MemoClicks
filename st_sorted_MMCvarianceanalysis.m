function [resultC,resultRC,resultRefRC,stats] = st_sorted_MMCvarianceanalysis(varargin)
% function [result] = st_sorted_MMCvarianceanalysis(varargin)

clear all; close all; clc;

%{
function to analyze fano factor of MMC sorted data
usage: [x,y,z] = st_sorted_MMCvarianceanalysis();
(cf. Churchland et al 2010, Nat Neurosci)

ST. v1:            02.07.15
Last modified:     07.07.15
Backup:            st_sorted_MMCvarianceanalysis_copy
Last backup saved: 07.07.15

% Data format
See VarVsMean.m in Variance Toolbox
alignTime = 0; % align to stimulus onset

1st 11 datasets: click onset at 0.4s; rest at 0.6s

% to do
- for 1st 5 vs. last 5 trials

%}

%% Set paths

MMCsingleunit.path = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/SpikeTimes/';
VarVsMeanpath      = '/Users/sundeepteki/Documents/ENS/Resources/Electrophysiology/Cunningham&Yu/Variance_toolbox/';
addpath(MMCsingleunit.path);
addpath(VarVsMeanpath);

%% Load data

tmpdir      = dir(MMCsingleunit.path);
x.onsets    = zeros(1,length(tmpdir)-2);
field       = 'spikes';
values      = {}; valuesC     = {}; valuesRC    = {}; valuesRefRC = {};
% trialsR     = [21:25]; % for last 5 trial analysis
% trialsC     = [46:50];
% tmp = [12 37 55 63 70 72:79 85 91:93 107 119 124 125 127];

for i1 = 1:length(tmpdir)-2
    
    xdata{i1} =  [tmpdir(i1+2).name];
    load(xdata{i1});
    x.onsets(i1) = mean([clickonset.start])*1000;
    disp(['Loading file - ' tmpdir(i1+2).name(1:end-4)]);
    
%             values{i1}      = spikes(:,:);
    valuesC{i1}     = Cspikes(:,:);
    valuesRC{i1}    = RCspikes(:,:);
    valuesRefRC{i1} = RefRCspikes(:,:);
    
%             mmcdata         = struct(field,values);
    Cdata           = struct(field,valuesC);
    RCdata          = struct(field,valuesRC);
    RefRCdata       = struct(field,valuesRefRC);
    
    % for finding tuned cells
    %     Cdata           = struct(field,valuesC{i1});
    %     RCdata          = struct(field,valuesRC{i1});
    %     RefRCdata       = struct(field,valuesRefRC{i1});
    
    %     Alldata(i1).cond(1).spikes = Cdata.spikes;
    %     Alldata(i1).cond(2).spikes = RCdata.spikes;
    %     Alldata(i1).cond(3).spikes = RefRCdata.spikes;
    
end

%% Find tuned cells

findtunedcells = 0;
if(findtunedcells)
    baselineInt = [100 600];
    stimulusInt = [600 1100]; % [600 2100];
    minresponseindex = 1;
    mintuningrange = 0.1;
    
    [goodData, goodCells, cellSummary] = TunedCells(Alldata, baselineInt, stimulusInt,minresponseindex,mintuningrange);
end

%% Perform VarVsMean analysis

computeVar = 1;

if (computeVar)
    timebin = 10;
    timeon  = 100;
    timeoff = 2300;
    times = timeon:timebin:timeoff; % save 100:25/50:2300
    % drop in FF during second repeat only for RefRC not for RC or C
    
%     result = VarVsMean(mmcdata, times);
    resultC     = VarVsMean(Cdata, times);
    resultRC    = VarVsMean(RCdata, times);
    resultRefRC = VarVsMean(RefRCdata, times);
end

%% Plot FF

plotFF = 1;
if(plotFF)
    
%         plotFano(result);
    plotFano(resultC);
    plotFano(resultRC);
    plotFano(resultRefRC);
    
end

%% Plot scatter

plotscatter = false;
if(plotscatter)
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultC,-100,scatterParams);
    text(2.5,7,'100ms before target - C','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultC,100,scatterParams);
    text(2.5,7,'100ms after target - C','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultC,500,scatterParams);
    text(2.5,7,'500ms after target - C','hori','center');
    
    %% 2RC
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRC,-100,scatterParams);
    text(2.5,7,'100ms before target - RC','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRC,100,scatterParams);
    text(2.5,7,'100ms after target - RC','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRC,500,scatterParams);
    text(2.5,7,'500ms after target - RC','hori','center');
    
    %% 2RefRC
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRefRC,-100,scatterParams);
    text(2.5,7,'100ms before target - RefRC','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRefRC,100,scatterParams);
    text(2.5,7,'100ms after target - RefRC','hori','center');
    
    scatterParams.axLim1 = 'auto';
    scatterParams.axLen  = 5;
    plotScatter(resultRefRC,500,scatterParams);
    text(2.5,7,'500ms after target - RefRC','hori','center');
    
end

%% 3 Movie

plotmovie = 0;
if(plotmovie)
    
    ScatterMovie(resultC);     pause(5);
    ScatterMovie(resultRC);    pause(5);
    ScatterMovie(resultRefRC); pause(5);
    
end


%% Check that FF is not an artifact of non-Poisson spiking

checkFF=0;
if(checkFF)
    
    resultC=VarVsMean(Fakerize(Cdata,'gamma'),times);
    resultRC=VarVsMean(Fakerize(RCdata,'gamma'),times);
    resultRefRC=VarVsMean(Fakerize(RefRCdata,'gamma'),times);
    plotFanoParams.plotRawF=1;
    
    plotFano(resultC,plotFanoParams);
    plotFano(resultRC,plotFanoParams);
    plotFano(resultRefRC,plotFanoParams);
    
end

%%

computestats = 1;
% compare FFs in 500ms windows - pre-stim and repeats 1,2,3 (timebin: 10ms)

if(computestats)
    
    %% compare FF in 500ms window during pre-stim vs. onset
    
    z=length([timeon:timebin:timeon+500]);
    
    % RefRC - significant: p = 0.0025; also significant for pre-stim vs. complete stimulus period of 1.5s
    t1RefRC = resultRefRC.FanoFactor(1:z);
    t2RefRC = resultRefRC.FanoFactor(z+1:4*z);
    [stats.hRefRConset,stats.pRefRConset,stats.CIRefRConset,stats.statsRefRConset] = ttest2(t1RefRC,t2RefRC);
    
    % RC - significant: p = 9.67e-14; also significant for pre-stim vs. complete stimulus period of 1.5s
    t1RC = resultRC.FanoFactor(1:z);
    t2RC = resultRC.FanoFactor(z+1:4*z);
    [stats.hRConset,stats.pRConset,stats.CIRConset,stats.statsRConset] = ttest2(t1RC,t2RC);
    
    % C - - significant: p = 3.05e-13; also significant for pre-stim vs. complete stimulus period of 1.5s
    t1C = resultC.FanoFactor(1:z);
    t2C = resultC.FanoFactor(z+1:4*z);
    [stats.hConset,stats.pConset,stats.CIConset,stats.statsConset] = ttest2(t1C,t2C);
    
    
    %% compare FF in 500ms window during repeat1 vs. 2
    
    z=length([timeon:timebin:timeon+500]);
    
    % RefRC - n.s. at p = 0.91; t = 0.1157; df = 100
    t1RefRC = resultRefRC.FanoFactor(z+1:2*z);
    t2RefRC = resultRefRC.FanoFactor(2*z+1:3*z);
    [stats.hRefRCrep1vs2,stats.pRefRCrep1vs2,stats.CIRefRCrep1vs2,stats.statsRefRCrep1vs2] = ttest2(t1RefRC,t2RefRC);
    
    % RC - significant at p = 0.0057; t = 2.83; df = 100
    t1RC = resultRC.FanoFactor(z+1:2*z);
    t2RC = resultRC.FanoFactor(2*z+1:3*z);
    [stats.hRCrep1vs2,stats.pRCrep1vs2,stats.CIRCrep1vs2,stats.statsRCrep1vs2] = ttest2(t1RC,t2RC);
    
    % C - significant at p = 0.0015; t = 3.26, df = 100
    t1C = resultC.FanoFactor(z+1:2*z);
    t2C = resultC.FanoFactor(2*z+1:3*z);
    [stats.hCrep1vs2,stats.pCrep1vs2,stats.CICrep1vs2,stats.statsCrep1vs2] = ttest2(t1C,t2C);
    
     %% compare FF in 500ms window during repeat2 vs. 3
    
    z=length([timeon:timebin:timeon+500]);
    
    % RefRC - n.s. at p = 0.12, t = -1.55, df = 100
    t1RefRC = resultRefRC.FanoFactor(2*z+1:3*z);
    t2RefRC = resultRefRC.FanoFactor(3*z+1:4*z);
    [stats.hRefRCrep2vs3,stats.pRefRCrep2vs3,stats.CIRefRCrep2vs3,stats.statsRefRCrep2vs3] = ttest2(t1RefRC,t2RefRC);
    
    % RC - significant at p = 6.999e-4, t = -3.5, df = 100
    t1RC = resultRC.FanoFactor(2*z+1:3*z);
    t2RC = resultRC.FanoFactor(3*z+1:4*z);
    [stats.hRCrep2vs3,stats.pRCrep2vs3,stats.CIRCrep2vs3,stats.statsRCrep2vs3] = ttest2(t1RC,t2RC);
    
    % C - significant at p = 0.019, t = -2.38, df = 100
    t1C = resultC.FanoFactor(2*z+1:3*z);
    t2C = resultC.FanoFactor(3*z+1:4*z);
    [stats.hCrep2vs3,stats.pCrep2vs3,stats.CICrep2vs3,stats.statsCrep2vs3] = ttest2(t1C,t2C);
    
    %% compare FF in 500ms window during stimulus RefRC vs. RC
    
    z=length([timeon:timebin:timeon+500]);
    
    % RefRC vs. RC during entire stimulus - significant at p = 0.0153, t = 2.44, df = 304
    t1RefRC = resultRefRC.FanoFactor(3*z+1:4*z);
    t1RC = resultRC.FanoFactor(3*z+1:4*z);
    [stats.hRefRCvsRC,stats.pRefRCvsRC,stats.CIRefRCvsRC,stats.statsRefRCvsRC] = ttest2(t1RefRC,t1RC);
    
    % RefRC vs. C during entire stimulus - significant at p = 0.011, t = -2.56, df = 304
    t1RefRC = resultRefRC.FanoFactor(3*z+1:4*z);
    t1C = resultC.FanoFactor(3*z+1:4*z);
    [stats.hRefRCvsC,stats.pRefRCvsC,stats.CIRefRCvsC,stats.statsRefRCvsC] = ttest2(t1RefRC,t1C);
    
    % RC vs. C during entire stimulus - significant at p = 1.5e-8, t = -5.82, df = 304
    t1RC = resultRC.FanoFactor(3*z+1:4*z);
    t1C = resultC.FanoFactor(3*z+1:4*z);
    [stats.hRCvsC,stats.pRCvsC,stats.CIRCvsC,stats.statsRCvsC] = ttest2(t1RC,t1C);
        
    % test: during repeat 1
    % RefRC vs. RC: n.s. at p = 0.51, t = 0.66, df = 100
    % RefRC vs. C: sign. at p = 0.0035, t = -2.99, df = 100
    % RC vs. C: sign. at p = 3.12e-5, t = -4.36, df = 100
    
    % test: during repeat 2
    % RefRC vs. RC: sign. at p = 0.012, t = 2.55, df = 100
    % RefRC vs. C: n.s. at p = 0.44, t = -0.78, df = 100
    % RC vs. C: sign. at p = 9.38e-7, t = -5.22, df = 100
    
    % test: during repeat 3
    % RefRC vs. RC: n.s. at p = 0.31, t = 1.02, df = 100
    % RefRC vs. C: n.s. at p = 0.46, t = -0.75, df = 100
    % RC vs. C: n.s. at p = 0.0798, t = -1.77, df = 100
    
    
end
