%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 6/6/2016
%
%   Description: 3) Measure background noise level in the tank.
%
%   Inputs:     n/a
%
%   Outputs:	n/a
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

figConst;
duration = 70; %recording duration in seconds
fS = 192000; %hz

%% Stop all currently-running data acquisition objects
if(~isempty(daqfind))
    stop(daqfind)
end

%% Collect info on Data Translation i/o board
daqHardwareInfo=daqhwinfo;
daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;

daqdtolInfo = daqhwinfo('dtol');

%% Create analog input and output objects and add channels

ai = analoginput('dtol',0);
addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
ChIn0 = ai.Channel;
ActualRateIn = setverify(ai,'SampleRate',floor(fS));
ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');

if(ActualRateIn~=fS)
    error('Input or Output Sample Rate Not Set Properly')
end

%% Define Input Sample Period and Wait Period
sampleRate = ai.SampleRate;
requiredSamples = floor(sampleRate*duration);
ai.SamplesPerTrigger = requiredSamples;
% 
waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input


%% Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
ai.TriggerType = 'Manual';

start(ai)
trigger(ai)

[data, time] = getdata(ai);

yLimits = [min([data]) max([data])];

figure;

plot(time,data);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
title(['Input Noise Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
xlim([0 duration]);
ylim(yLimits);
grid on;

temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
save(['fishDaq3_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);