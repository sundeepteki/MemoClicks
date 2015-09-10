% function [] = st_make_figure(varargin)

%{
   ST. 11.07.15
   Make group analysis plots for MMC data
   Created:  11.07.2015
   Backup:   11.07.2015
   
%}

plotMMCgroup = 0;
plotMMCcorr  = 0;
plotMMCshuffledcorr = 0;
plotMMCcorrbtw = 1;
plotMMCshuffledcorrbtw = 1;

%% 1

if(plotMMCgroup)
    
    close all;
    fig = figure; hold on;
    savepath = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/Results/GroupResults/';
    xticks = [0.5 1 1.5];
    
    errorbar(xticks,[nanmean(z.Cclickavgspikerateperclickfirst5) nanmean(z.RCclickavgspikerateperclickfirst5) nanmean(z.RefRCclickavgspikerateperclickfirst5)],...
        [nanstd(z.Cclickavgspikerateperclickfirst5)/sqrt(length(z.Cclickavgspikerateperclickfirst5)-1),...
        nanstd(z.RCclickavgspikerateperclickfirst5)/sqrt(length(z.RCclickavgspikerateperclickfirst5)-1),...
        nanstd(z.RefRCclickavgspikerateperclickfirst5)/sqrt(length(z.RefRCclickavgspikerateperclickfirst5)-1)],...
        'ko','LineWidth',2,'MarkerSize',8);
    
    set(gca,'XTick', xticks);
    set(gca,'XTickLabel',[{'C'} {'RC'} {'RefRC'}],'FontSize',12,'FontWeight','Bold');
    ylim([5.5 8.5]);
    ylabel('Average Spike rate (Hz)','FontSize',12,'FontWeight','bold');
    title(['Average click-evoked spike rate during first 5 trials'],'FontSize',12,'FontWeight','Bold');
    saveas(fig,[savepath 'ClickSpikerateFirst5.tiff']);
    
end

%% 2a

if (plotMMCcorr) % within-pair
    
    fig = figure; hold on;
    savepathcorr = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/';
    xticks = [0.5 1 1.5];
    
    errorbar(xticks,[nanmean(ravg.C) nanmean(ravg.RC) nanmean(ravg.RefRC)],...
        [nanstd(ravg.C)/sqrt(length(ravg.C)-1),...
        nanstd(ravg.RC)/sqrt(length(ravg.RC)-1),...
        nanstd(ravg.RefRC)/sqrt(length(ravg.RefRC)-1)],...
        'ko','LineWidth',2,'MarkerSize',8);
    
    set(gca,'XTick', xticks);
    set(gca,'XTickLabel',[{'C'} {'RC'} {'RefRC'}],'FontSize',12,'FontWeight','Bold');
    ylabel('Average Pearson correlation coefficient','FontSize',12,'FontWeight','bold');
    ylim([-5 12]*0.0001);
    title(['Average within-pair noise correlations'],'FontSize',12,'FontWeight','Bold');
    saveas(fig,[savepathcorr 'AvgWithinPairCorr.tiff']);
    
end


%% 2b

if (plotMMCshuffledcorr) % within-pair
    
    fig = figure; hold on;
    savepathcorr = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/';
    xticks = [0.5 1 1.5];
    
    errorbar(xticks,[nanmean(ravg.shC) nanmean(ravg.shRC) nanmean(ravg.shRefRC)],...
        [nanstd(ravg.shC)/sqrt(length(ravg.shC)-1),...
        nanstd(ravg.shRC)/sqrt(length(ravg.shRC)-1),...
        nanstd(ravg.shRefRC)/sqrt(length(ravg.shRefRC)-1)],...
        'ko','LineWidth',2,'MarkerSize',8);
    
    set(gca,'XTick', xticks);
    set(gca,'XTickLabel',[{'C'} {'RC'} {'RefRC'}],'FontSize',12,'FontWeight','Bold');
    ylabel('Average Pearson correlation coefficient','FontSize',12,'FontWeight','bold');
    ylim([-5 12]*0.0001);
    title(['Average within-pair shuffled noise correlations'],'FontSize',12,'FontWeight','Bold');
    saveas(fig,[savepathcorr 'AvgWithinPairShuffledCorr.tiff']);
    
end

%% 3a

if (plotMMCcorrbtw) % btw pop
    
    fig = figure; hold on;
    savepathcorr = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/';
    xticks = [0.5 1 1.5];
    
    errorbar(xticks,[nanmean(ravgbtw.C) nanmean(ravgbtw.RC) nanmean(ravgbtw.RefRC)],...
        [nanstd(ravgbtw.C)/sqrt(length(ravgbtw.C)-1),...
        nanstd(ravgbtw.RC)/sqrt(length(ravgbtw.RC)-1),...
        nanstd(ravgbtw.RefRC)/sqrt(length(ravgbtw.RefRC)-1)],...
        'ko','LineWidth',2,'MarkerSize',8);
    
    set(gca,'XTick', xticks);
    set(gca,'XTickLabel',[{'C'} {'RC'} {'RefRC'}],'FontSize',12,'FontWeight','Bold');
    ylabel('Average Pearson correlation coefficient','FontSize',12,'FontWeight','bold');
    ylim([-5 12]*0.0001);
    title(['Average between-population noise correlations'],'FontSize',12,'FontWeight','Bold');
    saveas(fig,[savepathcorr 'AvgBtwPairCorr.tiff']);
    
end


%% 3b

if (plotMMCshuffledcorrbtw) % btw pop
    
    fig1 = figure; hold on;
    savepathcorr = '/Users/sundeepteki/Dropbox (Equipe Audition)/Work/Shihab/#DATA/Maroille/#SortedData/#MMC/#analyzed/#singleunit/#analyzed/NoiseCorrelation/Results/';
    xticks = [0.5 1 1.5];
    
    errorbar(xticks,[nanmean(ravgbtw.shC) nanmean(ravgbtw.shRC) nanmean(ravgbtw.shRefRC)],...
        [nanstd(ravgbtw.shC)/sqrt(length(ravgbtw.shC)-1),...
        nanstd(ravgbtw.shRC)/sqrt(length(ravgbtw.shRC)-1),...
        nanstd(ravgbtw.shRefRC)/sqrt(length(ravgbtw.shRefRC)-1)],...
        'ko','LineWidth',2,'MarkerSize',8);
    
    set(gca,'XTick', xticks);
    set(gca,'XTickLabel',[{'C'} {'RC'} {'RefRC'}],'FontSize',12,'FontWeight','Bold');
    ylabel('Average Pearson correlation coefficient','FontSize',12,'FontWeight','bold');
    ylim([-5 12]*0.0001);
    title(['Average between-population shuffled noise correlations'],'FontSize',12,'FontWeight','Bold');
    saveas(fig1,[savepathcorr 'AvgBtwPairShuffledCorr.tiff']);
    
end
