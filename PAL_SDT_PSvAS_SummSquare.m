function [all_results] = PAL_SDT_PSvAS_SummSquare(StimLevels,PC,CMps,CMas,NFPps,NFPas,Qps,Qas,filename,RF_1,RF_2)
%%% PAL_SDT_PSvAS_SummSquare
%%% Frederick Kingdom, McGill University, Montreal, Canada, 2017
% Demonstrates use of Palamedes functions to determine whether the 
% data from a 5-PF (psychometric function) summation square 
% experiment, for the detection of two stimuli in the target interval, 
% accords more with probability summation (PS) or additive summation (AS) 
% under the assumptions of signal-detection-theory (SDT) and assuming
% that the observer is monitoring both channels sensitive to the two
% stimuli
%
% Note that the Bootstrap and Goodness-of-fit routines in this script will
% likely take some time to execute depending on the number of specified 
% simulations and the speed of the computer.  An indication of the time 
% to execute on a MacBook Pro running under OSX 10.9.5 is given when 
% the script is executed
%
%Demonstrates usage of Palamedes SDT summation multiple psychometric 
% function (PF) fitting routines:
%-PAL_SDT_Summ_MultiplePFML_Fit
%-PAL_SDT_Summ_MultiplePFML_BootstrapParametric
%-PAL_SDT_Summ_MultiplePFML_GoodnessOfFit
%
%Demonstrates usage of Palamedes SDT PS (probability) and AS (additive)
% summation routines:
%-PAL_SDT_PS_uneqSLtoPC
%-PAL_SDT_AS_uneqSLtoPC
%
%Demonstrates usage of Palamedes PF fitting routines:
%-PAL_PFML_Fit;
%
%Demonstrates usage of PF routines, e.g.
%-PAL_Logistic
%
%More information on any of these functions may be found by typing
% help followed by the name of the function, 
% e.g., help PAL_SDT_AS_uneqSLtoPC
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

% clear all;

current_Dir = pwd;
figure_folder = [current_Dir '/Figures'];




% message = sprintf('Number of simulations to determine standard errors: ');
% Bse = input(message);
Bse = 200;
% message = sprintf('Number of simulations to determine Goodness-of-Fit p-values: ');
% Bmc = input(message);
Bmc= 200;
% MacProTimePerSim=0.5;
% timeToExecute=round((Bse+Bmc)*MacProTimePerSim);
% message=sprintf('\rApprox. time to execute simulations on a MacBook Pro OSX 10.9.5 is %4d mins',timeToExecute); 
% disp(message);

% EXAMPLE
% Input stimulus amplitudes and PC data
% StimLevels(1,1,:)=[0 0 0 0 0 0]; % stim A 
% StimLevels(1,2,:)=[0.002 0.0032 0.005 0.0079 0.0126 0.0199]; % stim B 
%  
% StimLevels(2,1,:)=[0.0006 0.00095 0.0015 0.002375 0.003775 0.005975]; % stim A 
% StimLevels(2,2,:)=[0.0018 0.00285 0.0045 0.007125 0.011325 0.017925]; % stim B 
%  
% StimLevels(3,1,:)=[0.0014 0.0022 0.0035 0.00555 0.0088 0.01395]; % stim A 
% StimLevels(3,2,:)=[0.0014 0.0022 0.0035 0.00555 0.0088 0.01395]; % stim B
%  
% StimLevels(4,1,:)=[0.0024 0.00375 0.006 0.009525 0.015075 0.02385]; % stim A 
% StimLevels(4,2,:)=[0.0008 0.00125 0.002 0.003175 0.005025 0.00795]; % stim B
%  
% StimLevels(5,1,:)=[0.0048 0.0076 0.012 0.019 0.0301 0.0478]; % stim A 
% StimLevels(5,2,:)=[0 0 0 0 0 0]; % stim B 
%  
% PC(1,:)=[0.57 0.67 0.83 0.83 0.97 0.97];
% PC(2,:)=[0.67 0.67 0.9 0.9 1 0.97];
% PC(3,:)=[0.6 0.7 0.67 0.87 0.97 1];
% PC(4,:)=[0.67 0.67 0.73 0.8 0.87 1];
% PC(5,:)=[0.57 0.6 0.73 0.77 1 1];


% Reorganise the component condition stim leveles and PC data
tempPC=PC;
PC(1,:)=tempPC(5,:);
PC(5,:)=tempPC(1,:);

tempStim=StimLevels;
StimLevels(1,:,:)=tempStim(5,:,:);
StimLevels(5,:,:)=tempStim(1,:,:);

OutOfNum=ones(5,6).*30; % Note 30 trials per data point - Gunnar is this correct??
NumPos=PC.*OutOfNum;


M=2; % Number of alternatives in forced-choice task

% -------------------------- do some data plotting -----------------
% The following creates
% in order to estimate thresholds that can be put onto the graph.  This is 
% for illustrative purposes only - the fitted parameters here are NOT used 
% in the modeling of probability and additiv summation

searchGrid.alpha = 0:0.001:1; % range of possible threshold values
searchGrid.beta = logspace(.01,0.1,100); % range of possible slope values
searchGrid.gamma = .5;  %guessing rate
searchGrid.lambda = 0.0;  %lapse rate
paramsFree=[1 1 0 0];

%Fit Logistic functions
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant

[nrows, ncols, nnums]=size(StimLevels);

% Fit each data PF with a standard model in order to obtain thresholds
% to show on the graph
threshX=zeros(1,nrows);
threshY=zeros(1,nrows);

for i=1:nrows
    if (mean(StimLevels(i,1,:))==0)
        threshX(i)=0;
    else
    Stim(1,:)=StimLevels(i,1,:);
    [params, trash, trash, trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshX(i)=params(1);
    end
    
    if (mean(StimLevels(i,2,:))==0)
        threshY(i)=0;
    else
    Stim(1,:)=StimLevels(i,2,:);
    [params, trash, trash, trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshY(i)=params(1);
    end
end

% Plot data thresholds
fig=figure('name','Individual Psychometric Functions','units','pixels',...
    'position',[50 50 400 400]); 

a=plot(threshX,threshY,'ko','markersize',8,'markerfacecolor','k');
hold on
axis([0 .015 0 .015]);
set(gca,'YTick',0:0.005:0.015);
set(gca,'XTick',0:0.005:0.015);
axis square
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 20)
set(gca, 'LineWidth', 2)
xlabel(['Amplitude RF',num2str(RF_1)],'FontSize', 24);
ylabel(['Amplitude RF',num2str(RF_2)],'FontSize', 24);
title(filename, 'FontSize',16)
grid on





%--------- Fit PFs with probability summation (PS) model ----------

% Fit PS model to stim A, B and the combinations of A+B in the 5 PFs, 
% simultaneously, in order to obtain parameters gA, gB, pA and pB, 
% where g is stimulus gain and p the transducer exponent.  The fitting 
% procedure also provides the negative log likeleihood of the fit which
% is usd to compare the PS with the AS model

% Initial guesses for gA, gB, pA and pB
PSgParams=[200 200];
PSpParams=[0.8 0.8];


SummFunc=@PAL_SDT_PS_uneqSLtoPC;


% Fit all five data PFs simultaneously with PS model
[PSgParams, PSpParams, PSnegLL, trash, trash] = Test_PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,'PS',M,Qps,CMps);


%----------- Fit PFs with additive summation (AS) model ------------------


% Initial guesses for gA, gB, pA and pB
ASgParams=[200 200];
ASpParams=[1.0 1.0];

% Fit all PF data simultaneously with AS summation model
[ASgParams, ASpParams, ASnegLL, exitflag, output] = Test_PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,'AS',M,Qas,CMas);



% ------ Generate new PFs from fitted summation model and fit these with the same model as for the original data
% This is for illustrative purposes in the graph and is not part of the modelling itself -----

% ---- First for PS model---
StimLevelsBaseline=StimLevels(1,1,:);
meanOutNum=mean(mean(OutOfNum));
newOutNum=ones(1,nnums).*meanOutNum;
newPC=zeros(1,nnums);

numModelPoints=100;
minLogR=-2.5;
maxLogR=2.5;
incLogR=(maxLogR-minLogR)/numModelPoints;

PSthresh=zeros(2,numModelPoints);

for k=1:numModelPoints
    r=minLogR+k*incLogR;
    r=10.^r;
    x1=(r.*StimLevelsBaseline)./(1+r);
    x2=StimLevelsBaseline-x1;
   
    for i=1:nnums
        newPC(i)=PAL_SDT_PS_uneqSLtoPC([x1(i) x2(i)],PSgParams,PSpParams,M,Qps(2));
    end
    
    Pos=newPC.*meanOutNum;
    
    [params, trash, trash, trash] = PAL_PFML_Fit(x1,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    PSthresh(1,k)=params(1);
    [params, trash, trash, trash] = PAL_PFML_Fit(x2,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    PSthresh(2,k)=params(1);
    
end


% -----Second for the AS model -----

ASthresh=zeros(2,numModelPoints);

for k=1:numModelPoints
    r=minLogR+k*incLogR;
    r=10.^r;
    x1=(r.*StimLevelsBaseline)./(1+r);
    x2=StimLevelsBaseline-x1;
    
    for i=1:nnums
        newPC(i)=PAL_SDT_AS_uneqSLtoPC([x1(i) x2(i)],ASgParams,ASpParams,M,Qas(2));
    end
    
    Pos=newPC.*meanOutNum;
    
    [params, trash, trash, trash] = PAL_PFML_Fit(x1,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(1,k)=params(1);
    [params, trash, trash, trash] = PAL_PFML_Fit(x2,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(2,k)=params(1);
   
end


% Plot model predictions on graph
b=plot(PSthresh(1,:),PSthresh(2,:),'-','linewidth',2,'color',[0 .9 0]);
hold on;
c=plot(ASthresh(1,:),ASthresh(2,:),'-','linewidth',2,'color',[.9 0 0]);

% Add in legend
h = legend([a b c],'Data','PS Model','AS model');
set(h,'Interpreter','none','fontsize',16,'Location','NorthEast');
drawnow




%---------Determine bootstrap errors on fitted g and p parameters--------

fprintf('\rConducting bootstrap simulations for Prob. Summ. (PS) model:   ')
[PSgSE, PSpSE] = Test_PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,PSgParams,PSpParams,OutOfNum,'PS',M,Qps,Bse,CMps);

fprintf('\nConducting bootstrap simulations for Add. Summ. (AS) model:   ')
[ASgSE, ASpSE] = Test_PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,ASgParams,ASpParams,OutOfNum,'AS',M,Qas,Bse,CMas);


% ------------ Type out parameter estimates for both models-----------

message=sprintf('\rEstimated PS model params for stimulus A: gA: %6.4f se: %6.4f,  pA: %6.4f se: %6.4f',PSgParams(1),PSgSE(1),PSpParams(1),PSpSE(1)); 
disp(message);
message=sprintf('Estimated PS model params for stimulus B: gB: %6.4f se: %6.4f,  pB: %6.4f se: %6.4f',PSgParams(2),PSgSE(2),PSpParams(2),PSpSE(2)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus A: gA: %6.4f sd: %6.4f,  pA: %6.4f se: %6.4f',ASgParams(1),ASgSE(1),ASpParams(1),ASpSE(1)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus B: gB: %6.4f sd: %6.4f,  pB: %6.4f se: %6.4f',ASgParams(2),ASgSE(2),ASpParams(2),ASpSE(2)); 
disp(message);



%---------------------- Model comparisons ------------------------------

% AS-PS difference between Akaike's AIC.  Note that a positive value 
% means that PS model is better, a negative value that AS model is better
PSaic = -2*(-PSnegLL)+2*NFPps; % Note NFP eual number of free parameters
ASaic = -2*(-ASnegLL)+2*NFPas;
deltaAIC = ASaic - PSaic;
message = sprintf('\rAkaike Information Criterion (AIC) is: PS=%8.4f, AS=%8.4f, AS-PS AIC difference=%8.4f',PSaic,ASaic,deltaAIC);
disp(message);
message = sprintf('Note that positive value = PS better, negative = AS better');
disp(message);


%-----------------------Goodness of fits-------------------------------

fprintf('Conducting Goodness-of-fit (GOF) simulations for PS model:   ')
[trash, pDev, trash, trash] = Test_PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,'PS',M,Qps,Bmc,CMps);
PSpDev=pDev;

fprintf('\nConducting Goodness-of-fit (GOF) simulations for AS model:   ')
[trash, pDev, DevSim, converged] = Test_PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,'AS',M,Qas,Bmc,CMas);
ASpDev=pDev;

message=sprintf('\rGOF p-values (e.g. p<0.05 to reject model): PS p=%6.4f;  AS p=%6.4f',PSpDev,ASpDev); 
disp(message);

all_results = [PSgParams(1) PSgParams(2) PSpParams(1) PSpParams(2) ASgParams(1) ASgParams(2) ASpParams(1) ASpParams(2) PSaic ASaic deltaAIC];


cd(figure_folder)
savefig(fig,[filename,'.fig'])
cd(current_Dir)

end