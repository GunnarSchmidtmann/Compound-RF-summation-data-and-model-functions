
%%% Analyse_Data.m
%%% Gunnar Schmidtmann, University of Plymouth, UK - January 2018
%%% Frederick A.A. Kingdom, McGill Vision Research, McGill University,Canada,January - 2018


%%% Requirements: Palamedes Toolbox 1.8.2 (Prins,N. & Kingdom,F.A.A., (2009), Palamedes:Matlab routines for analyzing psychophysical data.  http://www.palamedestoolbox.org)

%%% The folders "Figures" and "Data" in the current directory are essential %%%

%%% IMPORTANT
%%% phase_arrangement 1,2,3 only for RF_1 = 3 & RF_2 = 5 / for
%%% compounds RF_1 = 3, RF_2 = 8 & RF_1 = 4 & RF_2 = 7 only
%%% phase_arrangment 1 & 3....there is no out of phase for those
%%% observer = 'GS';  % 'GS', 'GK' , 'RL'
%%% phase_arrangment = 1;   % 1 = in phase, 2 = out of phase, 3 = intermediate phase

%%% RF_1 = 3;
%%% RF_2 = 5;
%%% RF_1 = 3;
%%% RF_2 = 8;
%%% RF_1 = 4;
%%% RF_2 = 7;
%%% detection = 1; % enter 1 for detection data / 0 for discrimination data

clear all
close all
clc
commandwindow


current_Dir = pwd;
data_folder = [current_Dir '/Data'];
figure_folder = [current_Dir '/Figures'];

filename = 'conditions.xlsx';

[num,txt,raw] = xlsread(filename);


sequence = raw;
all_data=zeros(length(sequence),11);


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
[StimLevels, PC] = RF_compound_data(detVdisc,Observer,RF_1,RF_2,phase);

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





for i=1:length(sequence)

observer = char(sequence(i,1));
RF_1 = cell2mat(sequence(i,2));
RF_2 = cell2mat(sequence(i,3));
detection = cell2mat(sequence(i,4));
phase_arrangment = cell2mat(sequence(i,5));



if detection == 1
    condi_det = 'A=0';
else
    condi_det = 'A=0.05';
end


if phase_arrangment == 1
    condi_phase = 'in phase';
elseif phase_arrangment == 2
    condi_phase = 'out of phase';
elseif phase_arrangment == 3
    condi_phase = 'intermediate phase';
end



file_name = [observer,' ','RF1=',num2str(RF_1),' ','RF2=',num2str(RF_2),' ',condi_phase,' ',condi_det];

[ StimLevels, PC ] = RF_compound_data( detection, observer, RF_1, RF_2, phase_arrangment );




[all_results] = PAL_SDT_PSvAS_SummSquare(StimLevels,PC,CMps,CMas, NFPps,NFPas,Qps,Qas,file_name,RF_1,RF_2);

all_data(i,:)=all_results;

end

cd(data_folder)

%%% save the file in a folder called Data

%xlswrite('All_T_T_MAW_1.xlsx',all_data);

cd(current_Dir);