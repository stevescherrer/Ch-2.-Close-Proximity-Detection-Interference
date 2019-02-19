%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 6/6/2016
%
%   Description: 2) Record tag waveforms, remove tank IR to estimate 
%   tag source waveforms.
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
duration = 45; %recording duration in seconds
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
addchannel(ai,2); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
ChIn0 = ai.Channel;
ActualRateIn = setverify(ai,'SampleRate',floor(fS));
ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');

if(ActualRateIn~=fS)
    error('Output Sample Rate Not Set Properly')
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

recordStart = datetime;
[data, time] = getdata(ai);

figure;

plot(time,data);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
xlim([0 duration]);
% ylim(yLimits);
grid on;

%% Remove Tank IR from Received Waveform to Estimate Tag Source Waveform
%%% Trim leading zeroes in data and y

temp = find(data~=0);
iDataStart = min(temp);
iDataEnd = max(temp);
if(iDataEnd<length(y))
    data((iDataEnd+1):end) = [];
    time((iDataEnd+1):end) = [];
end
data(1:(iDataStart-1)) = [];
time(1:(iDataStart-1)) = [];

%%% Input tank impulse response
load('fishDaq1_','hEst');

sEst = fdeconv(data,hEst);
tsEst = (1:length(sEst))/fS;
 
yLimits = [min([data;sEst]) max([data;sEst])];

figure;
subplot(2,1,1);
plot(time,data);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
title(['Input Tag Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
xlim([0 duration]);
ylim(yLimits);
grid on;

subplot(2,1,2);
plot(tsEst,sEst);
xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
title('Estimated Tag Source Waveform','FontSize',TITLE_FONTSIZE);
xlim([0 duration]);
ylim(yLimits);
grid on;


%% Manually isolate start and end indices for tag waveform in recording
iStart = [];
iEnd = [];
tagSource = sEst(iStart:iEnd);

temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
save(['TankTest_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);