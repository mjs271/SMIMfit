%===============================================================================
%===============================================================================
%
% This script runs the conservative and reactive fit for the data,
% 'ExampleExperimentalBTC.mat', found in the directory 'ExampleData'.
%
% NOTE: the addpath() and rmpath() commands adds the entire SMIMfit directory to
% the path and removes is after running this script. If the script does not
% finish, these directories may remain on your path.
%
%===============================================================================
%===============================================================================

clear variables

%% turn the profiler on

% profile on

%% add SMIM code to path

addpath(genpath('./'), '-end');

%% load data

load('./ExampleData/ExampleExperimentalBTC.mat');

% tc_exp - vector of measurement times for the conservative tracer
% cc_exp - vector of conservative concentrations (mg/L)
% Equivalent tr_exp and cr_exp for reactive tracer (mg/L)
% xBTC - sensor location downstream of injection point (48.5 m)
% mc - mass of injected conservative tracer
% mr - mass of injected reactive tracer

%% conservative fit

% normalize BTC
[ccNorm, ~, Qdg] = cNorm(tc_exp, cc_exp, 'c', 'cMass', mc);

Mcons = SMIMfit(tc_exp, ccNorm, 'c', 'L', xBTC)

%% reactive fit

[crNorm, fmr] = cNorm(tr_exp, cr_exp, 'r', 'rMass', mr, 'Q', Qdg);

Mreact = SMIMfit(tr_exp, crNorm, 'r', 'L', xBTC, 'M', Mcons)

%% remove SMIM code from path

rmpath(genpath('./'));

%% plots

colors = lines(7);

fig = 1;
figure(fig)
clf
% conservative fit
subplot(1, 2, 1)
plot(tc_exp, ccNorm, 'linewidth', 5, 'color', colors(2, :))
hold on
plot(Mcons.tcfit, Mcons.ccfit, 'o', 'markersize', 10, 'linewidth', 3,...
     'color', colors(1, :))
ax = gca;
ax.FontSize = 16; 
legend('\textbf{Experimental}', '\textbf{Fitted}', 'Interpreter', 'latex',...
       'FontSize', 32, 'Location', 'northeast')
title('\textbf{Conservative}', 'Interpreter', 'latex', 'FontSize', 32)
xlabel('\textbf{Time}', 'Interpreter', 'latex', 'FontSize', 32)
ylabel('\textbf{Normalized Concentration [mol L{\boldmath$^{-1}$}]}',...
       'Interpreter', 'latex', 'FontSize', 32)

axis tight
subplot(1, 2, 2)
% reactive fit
plot(tr_exp, crNorm, 'linewidth', 5, 'color', colors(2, :))
hold on
plot(Mreact.kfits.kboth.trfit, Mreact.kfits.kboth.crfit, 'o',...
     'markersize', 10, 'linewidth', 3, 'color', colors(1, :))
ax = gca;
ax.FontSize = 16; 
legend('\textbf{Experimental}', '\textbf{Fitted}', 'Interpreter', 'latex',...
       'FontSize', 32, 'Location', 'northeast')
title('\textbf{Reactive}', 'Interpreter', 'latex', 'FontSize', 32)
xlabel('\textbf{Time}', 'Interpreter', 'latex', 'FontSize', 32)
ylabel('\textbf{Normalized Concentration [mol L{\boldmath$^{-1}$}]}',...
       'Interpreter', 'latex', 'FontSize', 32)
axis tight
figure(fig)
set(gcf, 'Position', [0, 100, 1800, 900])

%% save/load the profiler data and view in the explorer

% % save
% p = profile('info');
% save myprofiledata p

% % load
% load myprofiledata
% profview(0, p)
