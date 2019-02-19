%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (21396083)
%
%   Course: Scherrer - Tank Experiment
%
%   Date: 11/7/2016
%
%   Description: Project a sync pulse into the tank.
%
%   Inputs:     
%
%   Outputs:	
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [syncTime] = sync(fS, duration, y)

% Stop all currently-running data acquisition objects
if(~isempty(daqfind))
    stop(daqfind)
end

ai = analoginput('dtol',0);
addchannel(ai,1); %channel number is the same as the number on the analog channel numbers on the DAQ board, starting with channel 0
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

% Define Input Sample Period and Wait Period
sampleRate = ai.SampleRate;
requiredSamples = floor(sampleRate*duration);
ai.SamplesPerTrigger = requiredSamples;

% Queue output data
putdata(ao,y);

ai.TriggerType = 'Manual';
ao.TriggerType = 'Manual';

syncTime = datetime;
start([ai ao])
trigger(ai)
trigger(ao);

[~, ~] = getdata(ai);
q=1;