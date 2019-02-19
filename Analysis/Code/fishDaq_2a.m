%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 6/6/2016
%
%   Description: 2a) Calibrate calculated tag source waveform so that
%   received and isolated sources have approximately the same peak power.
%
%   Inputs:     n/a
%
%   Outputs:	n/a
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc
figConst

load('fishDaq2_','sEst','tsEst','data','time','fS'); %data = recorded tag waveform

% load('fishDaq2a_','sEst','tsEst','data','time','fS'); %data = recorded tag waveform, modified 
duration = 7; %recording duration in seconds
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

ao = analogoutput('dtol',0);
addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
ChOut0 = ao.Channel;
ActualRateOut = setverify(ao,'SampleRate',floor(fS));
ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');

if((ActualRateIn~=fS)|(ActualRateOut~=fS))
    error('Input or Output Sample Rate Not Set Properly')
end

%% Define Input Sample Period and Wait Period
sampleRate = ai.SampleRate;
requiredSamples = floor(sampleRate*duration);
ai.SamplesPerTrigger = requiredSamples;
% 
waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input


%% Prepare and Queue Output Data
simCase = 'Tag Estimate: Raw';
putdata(ao,sEst);

    
%% Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
ai.TriggerType = 'Manual';
ao.TriggerType = 'Manual';

start([ai ao])
trigger(ai)
trigger(ao);

[dataSest, timeSest] = getdata(ai);

yLimits = [min([data;dataSest]) max([data;dataSest])];

figure;

subplot(2,1,1);
plot(timeSest,dataSest);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
title('Observed Data, Estimated Tag Source Waveform','FontSize',TITLE_FONTSIZE);
xlim([0 duration]);
ylim(yLimits);
grid on;

subplot(2,1,2);
plot(time,data);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
title(['Observed Data, Measured Tag Source Waveform for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
xlim([0 duration]);
ylim(yLimits);
grid on;

[sEstMax,isEstMax] = max(dataSest);
[dataMax,idataMax] = max(data);

sEst = sEst * (dataMax/sEstMax);

temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
save(['fishDaq2a_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);