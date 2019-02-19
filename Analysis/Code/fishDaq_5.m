%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 6/6/2016
%
%   Description: 5) Convolve the estimated source waveform with each of 
%   the calculated IR vectors.  Project this waveform into the tank and 
%   record the result. 
%
%   Inputs:   n/a  
%
%   Outputs:	n/a
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

%% Test and Case Switch
tDM = true;
tDO = false;

c1A = false;
c1B = false;
c2A = false;
c2B = false;
c2C = false;
c3A = false;
c3B = false;
c4 = false;

cSG1 = false;
cSG2 = false;
cSG3 = false;

cMG1 = false;
cMG2 = false;
cMG3 = false;
cMG4 = false;
cMG5 = false;
cMG6 = false;

cDG1 = false;
cDG2 = true;
cDG3 = false;
cDG4 = false;
cDG5 = false;
cDG6 = false;

%% Constants and Data Structures

figConst;
duration = 5; %recording duration in seconds for each trial
nRepeat = 3; %number of times to run each simulation condition
pauseDuration = 5;
fS = 192000; %hz
plotOn = false;

load('sEst, Sync','sEstSync'); %sync tag source waveform 
nSync = length(sEstSync); 
syncPause = zeros(duration*fS - nSync,1);
syncTrain = 6*[sEstSync; syncPause];

load('sEst, Trial','sEstTrial'); %sync tag source waveform
% sEstTrial = sEstTrial * 0.6 / 1.4;
sEstTrial = sEstTrial*6;

% Sync storage data structures
% sync5 = cell(5*2,2);        iSync5 = 1; %before and after tests
sync3 = cell(3*2,1);        iSync3 = 1; %before and after cases
% sync1 = cell(1*12*nRepeat,1); iSync1 = 1; %before and after trials 

% %% 5-Fold Sync Pulse Train
% for i=1:5
%     [sync5{iSync5,1}, sync5{iSync5,2}] = sync(fS, duration, syncTrain);
%     iSync5 = iSync5 + 1;
% end

%% Test:Case:Trial Execution
if(tDM)
    if(c1A)
        load('fishDaq4_1A_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM1A'
    elseif(c1B)
        load('fishDaq4_1B_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM1B'
    elseif(c2A)
        load('fishDaq4_2A_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM2A'
    elseif(c2B)
        load('fishDaq4_2B_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM2B'
    elseif(c2C)
        load('fishDaq4_2C_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM2C'
    elseif(c3A)
        load('fishDaq4_3A_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM3A'
    elseif(c3B)
        load('fishDaq4_3B_2016_7_11__16_43_40','IRstore','nSim') 
        caseTest = 'Test DM3B'
    elseif(c4)
        load('fishDaq4_4_2016_7_11__16_43_40','IRstore','nSim')
        caseTest = 'Test DM4'
    elseif(cSG1)
        load('fishDaq4_Shallow_rS1TO101_2016_7_12__22_32_49');
        caseTest = 'Test SG1';
    elseif(cSG2)
        load('fishDaq4_Shallow_rS121TO201_2016_7_12__22_33_22');
        caseTest = 'Test SG2';
    elseif(cSG3)
        load('fishDaq4_Shallow_rS221TO301_2016_7_12__22_33_38');
        caseTest = 'Test SG3';
    elseif(cMG1)
        load('fishDaq4_Mid_rS1TO81_2016_7_12__22_36_22');
        caseTest = 'Test MG1';
    elseif(cMG2)
        load('fishDaq4_Mid_rS101TO181_2016_7_12__22_36_40');
        caseTest = 'Test MG2';
    elseif(cMG3)
        load('fishDaq4_Mid_rS201TO281_2016_7_12__22_36_55');
        caseTest = 'Test MG3';
    elseif(cMG4)
        load('fishDaq4_Mid_rS301TO381_2016_7_12__22_37_18');
        caseTest = 'Test MG4';
    elseif(cMG5)
        load('fishDaq4_Mid_rS401TO481_2016_7_12__22_37_25');
        caseTest = 'Test MG5';
    elseif(cMG6)
        load('fishDaq4_Mid_rS501TO521_2016_7_12__22_37_48');
        caseTest = 'Test MG6';
    elseif(cDG1)
        load('fishDaq4_Deep_rS1TO81_2016_7_13__1_46_28');
        caseTest = 'Test DG1';
    elseif(cDG2)
        load('fishDaq4_Deep_rS101TO181_2016_7_13__1_46_39');
        caseTest = 'Test DG2';
    elseif(cDG3)
        load('fishDaq4_Deep_rS201TO281_2016_7_13__1_46_52');
        caseTest = 'Test DG3';
    elseif(cDG4)
        load('fishDaq4_Deep_rS301TO381_2016_7_13__1_47_5');
        caseTest = 'Test DG4';
    elseif(cDG5)
        load('fishDaq4_Deep_rS401TO481_2016_7_13__1_47_30');
        caseTest = 'Test DG5';
    elseif(cDG6)
        load('fishDaq4_Deep_rS501TO581_2016_7_13__1_48_4');
        caseTest = 'Test DG6';
    else
        error('Choose Case Number');
    end
elseif(tDO)
    if(c1A)
        load('fishDaq4_1A_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO1A'
    elseif(c1B)
        load('fishDaq4_1B_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO1B'
    elseif(c2A)
        load('fishDaq4_2A_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO2A'
    elseif(c2B)
        load('fishDaq4_2B_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO2B'
    elseif(c2C)
        load('fishDaq4_2C_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO2C'
    elseif(c3A)
        load('fishDaq4_3A_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO3A'
    elseif(c3B)
        load('fishDaq4_3B_2016_7_11__16_49_19','IRstore','nSim') 
        caseTest = 'Test DO3B'
    elseif(c4)
        load('fishDaq4_4_2016_7_11__16_49_19','IRstore','nSim')
        caseTest = 'Test DO4'
    else
        error('Choose Case Number');
    end
else
    error('Choose Test Number');
end

%% Test 1A (June 2014, Deep Receiver)
% load('fishDaq4_1A_2016_6_23__16_42_55','IRstore','nSim') %shallow case impulse responses
% a = 'Test 1A'
% dataStore = cell(nSim,nRepeat);
% trialTime = cell(nSim*nRepeat);
% 
% for i=1:3
%     [sync3{iSync3,1}, sync3{iSync3,2}] = sync(fS, duration, syncTrain);
%     iSync3 = iSync3 + 1;
% end
% 
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         [sync1{iSync1,1}, sync1{iSync1,2}] = sync(fS, duration, syncTrain);
%         iSync1 = iSync1 + 1;
%                 
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
%         
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
%         
%         trialTime{(j-1)*nSim+i} = datetime;
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, ~] = getdata(ai);
%         
%         dataStore{i,j} = data;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% [sync1{iSync1,1}, sync1{iSync1,2}] = sync(fS, duration, syncTrain);
% iSync1 = iSync1 + 1;
% 
% for i=1:3
%     [sync3{iSync3,1}, sync3{iSync3,2}] = sync(fS, duration, syncTrain);
%     iSync3 = iSync3 + 1;
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save([fishDaq5_1A,'_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore','trialTime','sync3','sync1');
% clear IRstore nSim dataStore trialTime sync3 sync1

%% Test 1B (June 2014, Shallow Receiver)
% load('fishDaq4_1B_2016_6_23__16_42_55','IRstore','nSim') %shallow case impulse responses
% a = 'Test 1B'
% dataStore = cell(nSim,nRepeat);
% trialTime = cell(nSim*nRepeat);
% 
% for i=1:3
%     [sync3{iSync3,1}, sync3{iSync3,2}] = sync(fS, duration, syncTrain);
%     iSync3 = iSync3 + 1;
% end
% 
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         [sync1{iSync1,1}, sync1{iSync1,2}] = sync(fS, duration, syncTrain);
%         iSync1 = iSync1 + 1;
%         
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
%         
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         trialTime{(j-1)*nSim+i} = datetime;
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, ~] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% [sync1{iSync1,1}, sync1{iSync1,2}] = sync(fS, duration, syncTrain);
% iSync1 = iSync1 + 1;
% 
% for i=1:3
%     [sync3{iSync3,1}, sync3{iSync3,2}] = sync(fS, duration, syncTrain);
%     iSync3 = iSync3 + 1;
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_1B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 2A (November 2014, Deep Receiver)
% load('fishDaq4_2A_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 2A'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_2A_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 2B (November 2014, Mid-Depth Receiver)
% load('fishDaq4_2B_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 2B'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_2B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 2C (November 2014, Shallow)
% load('fishDaq4_2C_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 2C'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_2C_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 3A (March 2015, Shallow Water Depth)
% load('fishDaq4_3A_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 3A'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_3A_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 3B (March 2015, Deep Water Depth)
% load('fishDaq4_3B_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 3B'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_3B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Test 4 (June 2015)
% load('fishDaq4_4_2016_6_23__16_42_56','IRstore','nSim') %shallow case impulse responses
% a = 'Test 4A'
% dataStore = cell(nSim,nRepeat);
% % timeStore = dataStore;
% for j=1:nRepeat
%     j
%     for i=1:nSim
%         h=IRstore{i};
% 
%         % Stop all currently-running data acquisition objects
%         if(~isempty(daqfind))
%             stop(daqfind)
%         end
% 
%         % Collect info on Data Translation i/o board
%         daqHardwareInfo=daqhwinfo;
%         daqInstalledAdaptors = daqHardwareInfo.InstalledAdaptors;
% 
%         daqdtolInfo = daqhwinfo('dtol');
% 
%         % Create analog input and output objects and add channels
% 
%         ai = analoginput('dtol',0);
%         addchannel(ai,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChIn0 = ai.Channel;
%         ActualRateIn = setverify(ai,'SampleRate',floor(fS));
%         ActualBufferingModeIn = setverify(ai,'BufferingMode','Auto');
% 
%         ao = analogoutput('dtol',0);
%         addchannel(ao,0); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
%         ChOut0 = ao.Channel;
%         ActualRateOut = setverify(ao,'SampleRate',floor(fS));
%         ActualBufferingModeOut = setverify(ao,'BufferingMode','Auto');
% 
%         if((ActualRateIn~=fS)|(ActualRateOut~=fS))
%             error('Input or Output Sample Rate Not Set Properly')
%         end
% 
%         % Define Input Sample Period and Wait Period
%         sampleRate = ai.SampleRate;
%         requiredSamples = floor(sampleRate*duration);
%         ai.SamplesPerTrigger = requiredSamples;
% 
%         waitTime = duration*1.1+0.5; %ensure that the Matlab script doesn't overrun the i/o device during input
% 
%         % Prepare and Queue Output Data
%         y = conv(h,sEst);
%         putdata(ao,y);
% 
%         % Configure data production and recovery triggers, Collect Data, Bring Into Matlab Workspace, and Plot
%         ai.TriggerType = 'Manual';
%         ao.TriggerType = 'Manual';
% 
%         start([ai ao])
%         trigger(ai)
%         trigger(ao);
% 
%         [data, time] = getdata(ai);
%         
%         dataStore{i,j} = data;
% %         timeStore{i,j} = time;
% 
%         if(plotOn)
%             yLimits = [min([data;y]) max([data;y])];
% 
%             figure;
% 
%             subplot(2,1,1);
%             plot(t,y);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE);
%             title('Output Data','FontSize',TITLE_FONTSIZE);
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
% 
%             subplot(2,1,2);
%             plot(time,data);
%             xlabel('Time [s]','FontSize',LABEL_FONTSIZE);  % Setting up the xlabel
%             ylabel('Signal [V]','FontSize',LABEL_FONTSIZE); % Setting up the ylabel
%             title(['Input Data Acquired using Data Translation DT-9832-4-2-BNC for ',num2str(duration), ' seconds'],'FontSize',TITLE_FONTSIZE); % Setting up the title
%             xlim([0 duration]);
%             ylim(yLimits);
%             grid on;
%         end
%     end
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_4_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)],'dataStore','ActualRateOut','IRstore');
% clear IRstore nSim dataStore

%% Grid - Shallow Case

% dW = 50;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 48];               %Test 2 Receiver, Closest to Sea Bed   
% 
% rSres = 20;
% zSres = 7.5;
% rSstore = 241:rSres:301; nrS = length(rSstore);
% zSstore = 1:zSres:49.5; nzS = length(zSstore);
% nSim = nrS*nzS;
% 
% Sstore = NaN(nSim,2);
% sampArrStore = NaN(nSim,sum(arrivePath));
% tArrStore = sampArrStore;
% RLstore = sampArrStore;
% 
% for i=1:nrS
%     rS = rSstore(i);
%     for j=1:nzS
%         offset = (i-1)*nzS+j;
%         zS = zSstore(j);
%         Sstore(offset,:) = [rS, zS];
%         
%         m = [Sstore(offset,1),0,Sstore(offset,2),tS,tSync,R(1),R(2),R(3),dW,cW,cW];
%         nArr = sum(arrivePath);
%         [tArr,lArr] = MOI_XYZ(m,arrivePath); tArrStore(offset,:) = tArr'; lArr = lArr';
%         sampArr = round(tArr*fS) + 1; %time=0 corresponds to the first sample
%         sampArrStore(offset,:) = sampArr - min(sampArr) + 1; %each channel starts at first sample
%         RL = 10.^(-1*log10(lArr)); %assumes omnidirectional source level of 0dB, no reflection losses/absorption/scattering
%         RL = RL_PhaseInvert(arrivePath,RL);
%         RLstore(offset,:) = RL;
%     end
% end
% 
% IRstore = cell(nSim,1);
% 
% for i=1:nSim
%     hTemp = zeros(1,max(sampArrStore(i,:)));
%     hTemp(sampArrStore(i,:)) = RLstore(i,:);
%     IRstore{i} = hTemp;
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_Shallow_rS241TO301_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);

%% Grid - Mid-Depth Case

% dW = 212;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 210];               %Test 2 Receiver, Closest to Sea Bed   
% 
% rSres = 20;
% zSres = 40;
% rSstore = 421:rSres:501; nrS = length(rSstore);
% zSstore = 1:zSres:210; nzS = length(zSstore);
% nSim = nrS*nzS;
% 
% Sstore = NaN(nSim,2);
% sampArrStore = NaN(nSim,sum(arrivePath));
% tArrStore = sampArrStore;
% RLstore = sampArrStore;
% 
% for i=1:nrS
%     rS = rSstore(i);
%     for j=1:nzS
%         offset = (i-1)*nzS+j;
%         zS = zSstore(j);
%         Sstore(offset,:) = [rS, zS];
%         
%         m = [Sstore(offset,1),0,Sstore(offset,2),tS,tSync,R(1),R(2),R(3),dW,cW,cW];
%         nArr = sum(arrivePath);
%         [tArr,lArr] = MOI_XYZ(m,arrivePath); tArrStore(offset,:) = tArr'; lArr = lArr';
%         sampArr = round(tArr*fS) + 1; %time=0 corresponds to the first sample
%         sampArrStore(offset,:) = sampArr - min(sampArr) + 1; %each channel starts at first sample
%         RL = 10.^(-1*log10(lArr)); %assumes omnidirectional source level of 0dB, no reflection losses/absorption/scattering
%         RL = RL_PhaseInvert(arrivePath,RL);
%         RLstore(offset,:) = RL;
%     end
% end
% 
% IRstore = cell(nSim,1);
% 
% for i=1:nSim
%     hTemp = zeros(1,max(sampArrStore(i,:)));
%     hTemp(sampArrStore(i,:)) = RLstore(i,:);
%     IRstore{i} = hTemp;
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_Mid_rS421TO501_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);

%% Grid - Deep Case

% dW = 500;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 498];               %Test 2 Receiver, Closest to Sea Bed   
% 
% rSres = 20;
% zSres = 90;
% rSstore = 321:rSres:401; nrS = length(rSstore);
% zSstore = 1:zSres:498; nzS = length(zSstore);
% nSim = nrS*nzS;
% 
% Sstore = NaN(nSim,2);
% sampArrStore = NaN(nSim,sum(arrivePath));
% tArrStore = sampArrStore;
% RLstore = sampArrStore;
% 
% for i=1:nrS
%     rS = rSstore(i);
%     for j=1:nzS
%         offset = (i-1)*nzS+j;
%         zS = zSstore(j);
%         Sstore(offset,:) = [rS, zS];
%         
%         m = [Sstore(offset,1),0,Sstore(offset,2),tS,tSync,R(1),R(2),R(3),dW,cW,cW];
%         nArr = sum(arrivePath);
%         [tArr,lArr] = MOI_XYZ(m,arrivePath); tArrStore(offset,:) = tArr'; lArr = lArr';
%         sampArr = round(tArr*fS) + 1; %time=0 corresponds to the first sample
%         sampArrStore(offset,:) = sampArr - min(sampArr) + 1; %each channel starts at first sample
%         RL = 10.^(-1*log10(lArr)); %assumes omnidirectional source level of 0dB, no reflection losses/absorption/scattering
%         RL = RL_PhaseInvert(arrivePath,RL);
%         RLstore(offset,:) = RL;
%     end
% end
% 
% IRstore = cell(nSim,1);
% 
% for i=1:nSim
%     hTemp = zeros(1,max(sampArrStore(i,:)));
%     hTemp(sampArrStore(i,:)) = RLstore(i,:);
%     IRstore{i} = hTemp;
% end
% 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq5_Deep_rS321TO401_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);