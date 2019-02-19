%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 1/6/2016
%
%   Description: 1) Produce IR-characterization chirp, record response, 
%       calculate IR. Analog input/output using a Data Translation 
%       Simultaneous A/D, USB Series, DT9832-4-2-BNC-IO.
%
%   Inputs:     n/a
%
%   Outputs:	n/a
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

plotOn=true;
figConst;
duration = 4.5; %recording duration in seconds
fS = 192000; %hz

SourceChirp = false;
SourceTag = true;
SourceTone = false;
SourceBlank = false;

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
if(SourceChirp)
    %%% Wide Bandwidth Chirp
    f0 = 5000;
    f1 = 8000;
    tMaxData = 6;
    t = (0:1/fS:tMaxData);
    y_orig = 2*chirp(t,f0,t(end),f1,'linear');
    h = hann(length(y_orig))';
    y = (y_orig.*h)';
    y = y/5;
%     z = zeros(1,length(y));
    simCase = 'Chirp';
    putdata(ao,[y]);
elseif(SourceTag)
    load('sEst, Trial','sEstTrial'); %sync tag source waveform
    y = 6*sEstTrial;
    t = (0:(length(y)-1))/fS;
   
    simCase = 'Tank Tag Recording';
    putdata(ao,y);
    
elseif(SourceTone)
    f = 69000;
    tMaxData = 1;
    t_in = (0:1/fS:tMaxData);
    y = sin(2*pi*f*t_in)';
    y = [y;zeros(192000,1);y];
    t = (1:length(y))/fS;
    simCase = 'Tone';
    putdata(ao,y);
elseif(SourceBlank)
    putdata(ao,[1]);
end
    
%% Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
ai.TriggerType = 'Manual';
ao.TriggerType = 'Manual';

start([ai ao])
trigger(ai)
trigger(ao);

[data, time] = getdata(ai);

yLimits = [min([data;y]) max([data;y])];

if(plotOn)
    figure;

    subplot(2,1,1);
    plot(t,y);
    xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
    ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
    title('Output Data','FontSize',TITLE_FONTSIZE);
    xlim([0 duration]);
    ylim(yLimits);
    grid on;

    subplot(2,1,2);
    plot(time,data);
    xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
    ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
    title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
    xlim([0 duration]);
    ylim([-3 3]);
    grid on;
end

%% Deconvolve Chirp Source from Received Waveform to Estimate Tank IR
%%% Trim leading zeroes in data and y
temp = find(y~=0);
iYStart = min(temp);
iYEnd = max(temp);
if(iYEnd<length(y))
    y((iYEnd+1):end) = [];
    t((iYEnd+1):end) = [];
end
y(1:(iYStart-1)) = [];
t(1:(iYStart-1)) = [];

temp = find(data~=0);
iDataStart = min(temp);
iDataEnd = max(temp);
if(iDataEnd<length(y))
    data((iDataEnd+1):end) = [];
    time((iDataEnd+1):end) = [];
end
data(1:(iDataStart-1)) = [];
time(1:(iDataStart-1)) = [];


hEst = fdeconv(data,y); %% deconvolve y out of data using FFT division

temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
save(['fishDaq1_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);