%
% Gunnar_SummSquareSingleChannel
% single channel AS model
% 

clear all;



RF=input('Enter the two RFs as a vector e.g. [3 5]: ');
phase=input('Enter their phase relationship 1 = in phase, 2 = out of phase, 3 = intermediate phase: ');
detVdisc=input('Detection or discrimination data 1 or 0: ');
if detVdisc==0
    detVdiscFlag=input('Treat discrimination data as discrimination "y" or detection "n" ? ','s');
end
Observer=input('Enter observer ','s');

Bse=input('Number of simulations to determine standard errors: ');
Bmc=input('Number of simulations to determine Goodness-of-Fit p-values: ');

% get some data
[StimLevels, PC] = RF_compound_data(detVdisc,Observer,RF(1),RF(2),phase);

% Code phase in terms of radians
P(1)=0.0;
if (phase==1) 
    P(2)=0.0;
end
if (phase==2) 
    P(2)=pi;
end
if (phase==3) 
    P(2)=pi/2;
end

% Set pedestal amplitude
pedAmp=0.0; %default
if ((detVdisc==0) && (detVdiscFlag=='y'))
    pedAmp=0.05;
end


%Number of trials for A, B and A+B. May contain zeros. 
OutOfNum(1,:) = [30 30 30 30 30 30];
OutOfNum(2,:) = [30 30 30 30 30 30];
OutOfNum(3,:) = [30 30 30 30 30 30];
OutOfNum(4,:) = [30 30 30 30 30 30];
OutOfNum(5,:) = [30 30 30 30 30 30];


%Plot raw data, proportion correct PC against stimulus level
NumPos=PC.*OutOfNum; 


%------- Fit each PF with a standard model, such as a Logistic or Wibull,
% in order to estimate thresholds that can be put onto the graph.  This is 
% for illustrative purposes only - the fitted parameters here are NOT used 
% in the modeling of probability and additiv summation

[nrows, ncols, nnums]=size(StimLevels);

searchGrid.alpha = 0:.001:1; % range of possible threshold values
searchGrid.beta = logspace(0.1,0.1,100); % range of possible slope values
searchGrid.gamma = .5;  %guessing rate
searchGrid.lambda = 0.0;  %lapse rate

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
 
%Fit Logistic functions
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant


% Fit each data PF with a standard model in order to obtain thresholds
% to show on the graph
for i=1:nrows
    if (mean(StimLevels(i,1,:))==0)
        threshX(i)=0;
    else
    Stim(1,:)=StimLevels(i,1,:);
    [params trash trash trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshX(i)=params(1);
    end
    
    if (mean(StimLevels(i,2,:))==0)
        threshY(i)=0;
    else
    Stim(1,:)=StimLevels(i,2,:);
    [params trash trash trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshY(i)=params(1);
    end
end

% Plot data thresholds
figure('name','Individual Psychometric Functions','units','pixels',...
    'position',[50 50 400 400]); 

a=plot(threshX,threshY,'ko','markersize',8,'markerfacecolor','k');
hold on;
set(gca, 'fontsize',18);
set(gca, 'Xtick',0:0.005:0.015);
set(gca, 'Ytick',0:0.005:0.015);
axis([0 0.015 0 0.015]); axis square;
xlabel('Stimulus A');
ylabel('Stimulus B');
title('Summation square');
drawnow



Q=1; % Number of monitored channels
M=2; % Number of alternatives in forced-choice task



%----------- Fit PFs with additive summation (AS) model ------------------

% Same as above for addiitve summation model

SummFunc=@PAL_SDT_AS_uneqSLtoPC;

% Initial guesses for gA, gB, pA and pB
ASgParams=[120.0];
ASpParams=[1.0];

% Fit all PF data simultaneously with AS summation model
[ASgParams, ASpParams, ASnegLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q);


%---------Generate individual PFs from AS model but fitted with the standard PF fitting model as for the original data
% This is for illustrative purposes in the graph and is not part of the modelling itself -----

% first generate model PCs
ASmodelPC=zeros(nrows,nnums);
for i=1:nrows
    for k=1:nnums
        ASmodelPC(i,k)=PAL_SDT_AS_uneqSLtoPC([StimLevels(i,1,k) StimLevels(i,2,k)],ASgParams,ASpParams,M,Q);
    end
end

% next fit with standard PF model
NumPosASmodel=ASmodelPC.*OutOfNum;


% Fit each model PF with a standard PF model to obtain thresholds
for i=1:nrows
    if (mean(StimLevels(i,1,:))==0)
        threshX(i)=0;
    else
    Stim(1,:)=StimLevels(i,1,:);
    [paramsPFout trash trash trash] = PAL_PFML_Fit(Stim(1,:),NumPosASmodel(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshX(i)=paramsPFout(1);
    end
    
    if (mean(StimLevels(i,2,:))==0)
        threshY(i)=0;
    else
    Stim(1,:)=StimLevels(i,2,:);
    [paramsPFout trash trash trash] = PAL_PFML_Fit(Stim(1,:),NumPosASmodel(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshY(i)=paramsPFout(1);
    end
end

% Plot model thresholds  
b=plot(threshX,threshY,'ko','markersize',8,'markerfacecolor','r');
hold on;
set(gca, 'fontsize',18);
set(gca, 'Xtick',0:0.005:0.015);
set(gca, 'Ytick',0:0.005:0.015);
axis([0 0.015 0 0.015]); axis square;
xlabel('Stimulus A');
ylabel('Stimulus B');
title('Summation square');
drawnow



% ------ Generate continuous line plots of thresholds from fitted summation model and fit these with the same model as for the original data
% This is for illustrative purposes in the graph and is not part of the modelling itself -----

StimLevelsBaseline=StimLevels(5,1,:);

meanOutNum=mean(mean(OutOfNum));
newOutNum=ones(1,nnums).*meanOutNum;
numModelPoints=100;
minLogR=-2.5;
maxLogR=2.5;
incLogR=(maxLogR-minLogR)/numModelPoints;


ASthresh=zeros(2,numModelPoints);

for k=1:numModelPoints
    r=minLogR+k*incLogR;
    r=10.^r;
    x1=(r.*StimLevelsBaseline)./(1+r);
    x2=StimLevelsBaseline-x1;
    
    for i=1:nnums
        PCnew(i)=PAL_SDT_AS_uneqSLtoPC([x1(i) x2(i)],ASgParams,ASpParams,M,Q);
    end
    
    Pos=PCnew.*meanOutNum;
    
    [params LL exitflag output] = PAL_PFML_Fit(x1,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(1,k)=params(1);
    [params LL exitflag output] = PAL_PFML_Fit(x2,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(2,k)=params(1);
   
end


% Plot model predictions on graph
c=plot(ASthresh(1,:),ASthresh(2,:),'r-','linewidth',2);

% Add in legend
h = legend([a c],'Data','AS model');
set(h,'Interpreter','none','fontsize',16,'Location','NorthEast');
drawnow




%---------Determine bootstrap errors on fitted g and p parameters--------

message = sprintf('\rPerforming bootstrap simulations for standard errors...');
disp(message);

SummFunc=@PAL_SDT_AS_uneqSLtoPC;
fprintf('\nCounting simulations for Add. Summ. model:   ')
[ASgSE, ASpSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,ASgParams,ASpParams,OutOfNum,SummFunc,M,Q,Bse);


% ------------ Type out parameter estimates for both models-----------

message=sprintf('\rEstimated Add. Summ. model params: gA: %6.4f sd: %6.4f,  pA: %6.4f se: %6.4f',ASgParams(1),ASgSE(1),ASpParams(1),ASpSE(1)); 
disp(message);



%---------------------- Model comparisons ------------------------------

% AS-PS difference between Akaike's AIC.  Note that a positive value 
% means that PS model is better, a negative value that AS model is better
ASaic = -2*(-ASnegLL)+2*2; % Note 2 free parameters

message=sprintf('AIC AS: %8.4f',ASaic); 
disp(message);

%-----------------------Goodness of fits-------------------------------
message = sprintf('\rPerforming Goodness of Fits (GOF) of each model...');
disp(message);

SummFunc=@PAL_SDT_AS_uneqSLtoPC;
fprintf('\nCounting GOF simulations for Add. Summ. model:   ')
[trash, pDev, DevSim, converged] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q,Bmc);
ASpDev=pDev;

message=sprintf('\rp-value for GOF (e.g. p<0.05 to reject): AS p=%6.4f',ASpDev); 
disp(message);

