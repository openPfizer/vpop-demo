% Main script for response to Duffull and Gulati.
% This script does the following:
% 1. Generates simulated data for matching by the model by adding data to an "acceptable" profile.
% 2. Launches a Metropolis-Hastings search to find plausible patients that match the distribution from (1).
% 3. Simulate all the retained plausible patients and plot the final figure.

% Basic Matlab setup:
clear;clc;close all;
rng('shuffle');

% Set the number of plausible patients desired in the final cohort:
npp = 2500;

%% Generate simulated clinical data:
[t,bg_data,bi_data] = sim_clin_data(100,0.1,0.5,1); 
close all; % Close any unneeded intermediate figures.

%% Create a criteria for selection:
t = t(:);
mudata = mean([bg_data' bi_data']);
covdata = cov([bg_data' bi_data']);

%% Generate PPs, using M-H:
% Parameter bounds, this is PROBLEM specific, here we choose a wide search range
% to demonstrate a high degreee of prior knowledge is not necessary (although helpful):
p_bnds = [0.1 10;
    0.1 10;
    0.4 40;
    0.01 1;
    0.1 10];

[pps,pp_yield] = mh_generate_pps(npp,mudata,covdata,0.1,t,p_bnds); % Main plausible patient search

%% Commentary Figure 1:

% Fetch a certain percentile of the simulated clinical data:
bg_pctile = prctile(bg_data',[0 100]);
bi_pctile = prctile(bi_data',[0 100]);

f1 = figure('Name','Summary of Simulations','Units','Inches','Position',[1 1 10.64 5.36]);

subplot(1,2,1);

for ii = 1:size(pps,2)
    [s(ii), t_sim, bg_sim(:,ii), bi_sim(:,ii)] = sim_clin_data_for_mh(pps(:,ii),t,mudata,covdata);
end

eb_bg = prctile(bg_sim',[0 100])';

patch([t;flipud(t)],[eb_bg(:,1);flipud(eb_bg(:,2))],[0.6 0.6 0.6],'facealpha',0.3,'edgecolor','none');
hold on;

f11(1) = plot(t,mean(bg_data'),'-o','LineWidth',3,'Color',[0.85 0.325 0.098],'MarkerSize',12,'MarkerFaceColor',[0.85 0.325 0.098]);
for ii = 1:numel(t)
    plot([t(ii) t(ii)],bg_pctile(:,ii),'-','LineWidth',3,'Color',[0.85 0.325 0.098]);
    plot([t(ii)-0.3333 t(ii)+0.3333],[bg_pctile(1,ii) bg_pctile(1,ii)],'-','LineWidth',3,'Color',[0.85 0.325 0.098]);
    plot([t(ii)-0.3333 t(ii)+0.3333],[bg_pctile(2,ii) bg_pctile(2,ii)],'-','LineWidth',3,'Color',[0.85 0.325 0.098]);
end

set(gca,'LineWidth',2,'FontSize',18,'box','on');
xlabel('Time');
ylabel('A_1');
ylim([0 12]);
xlim([0 12]);
%title('A_1 - Simulated Data Summary and Generated PPs');

subplot(1,2,2);

eb_bi = prctile(bi_sim',[0 100])';
patch([t;flipud(t)],[eb_bi(:,1);flipud(eb_bi(:,2))],[0.6 0.6 0.6],'facealpha',0.3,'edgecolor','none');
hold on;

f12(1) = plot(t,mean(bi_data'),'-s','LineWidth',3,'Color',[0 0.4470 0.741],'MarkerSize',12,'MarkerFaceColor',[0 0.4470 0.741]);
for ii = 1:numel(t)
    plot([t(ii) t(ii)],bi_pctile(:,ii),'-','LineWidth',3,'Color',[0 0.4470 0.741]);
    plot([t(ii)-0.3333 t(ii)+0.3333],[bi_pctile(1,ii) bi_pctile(1,ii)],'-','LineWidth',3,'Color',[0 0.4470 0.741]);
    plot([t(ii)-0.3333 t(ii)+0.3333],[bi_pctile(2,ii) bi_pctile(2,ii)],'-','LineWidth',3,'Color',[0 0.4470 0.741]);
end

set(gca,'LineWidth',2,'FontSize',18,'box','on');
xlabel('Time');
ylabel('B_1');
xlim([0 12]);
ylim([0 3]);
%title('B_1 - Simulated Data Summary and Generated PPs');

saveas(f1,'Figure1.png');
