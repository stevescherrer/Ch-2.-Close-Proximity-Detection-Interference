%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI Modeling
%
%   Date: 25/06/2017
%
%   Description: Calculate checksum invalidation success/failure at each
%   combination of rS:zS:zP.  Assume dW, cW, tBlanking, max detection
%   range, target fish depth interval given.
%
%   Inputs:     n/a
%
%   Outputs:	n/a
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
clc

LABEL_FONTSIZE = 20;
TITLE_FONTSIZE = 26;
AXIS_FONTSIZE = 22;
MARKER_FONTSIZE = 15;
LINE_WIDTH = 1.5;

%% Simulation Conditions
%%% FIXED
dW = 300; %m, water depth
cW = 1530; %m/s, uniform sound speed
tBlank = 0.26; %s, blanking period length
maxDetectRange = 847; %maximum detectable acoustic path length
dTargFish = [1 dW-1]; %m, min and max depth fish of interest found at
rP = 0; %m, phone horizontal range from the coordinate systme origin
tS = 0; %time signal is created, leave at zero
tSync = 0; %constant required by the propagation model
arrivePath =  [1 1 1 1  1 1 1 1  0 0 0 0  0 0 0 0  0 0 0 0]; %path lengths to simulate from the first 20, refer to MOI_XYZ for more details

%%% SEARCH RANGES
rSstore = 1:1:900; nrSstore = length(rSstore); %set of candidate fish horizontal ranges to test
zSstore = dTargFish(1):1:dTargFish(2); nzSstore = length(zSstore); %set of candidate fish depths (measured positive downward from the surface) to test
zPstore = [dW-2]; nzPstore = length(zPstore); %receiver depth
if(max(zSstore)>dW)
    error('Candidate source depth range exceeds water depth');
end

%% Checksum Status Calculation
checkSumPass = NaN(length(zSstore),length(rSstore),length(zPstore));
maxTarr = checkSumPass;
lArrStore = checkSumPass;
m = [NaN,0,NaN,tS,tSync,rP,0,NaN,dW,cW,cW];
for i=1:nrSstore
    i
    m(1) = rSstore(i);
    for j=1:nzSstore
        m(3) = zSstore(j);
        for k=1:nzPstore
            m(8) = zPstore(k);
            [tArr,lArr] = MOI_XYZ(m,arrivePath);
            tArr = tArr - tArr(1);
            tempIndL = find(lArr<maxDetectRange); %find all paths whose path length is less than the maximum path length detectable by the Vemco receiver
            [tempVal, tempInd] = max(tArr(tempIndL)); %find the longest path smaller than the maximum detectable path length
            if(isempty(tempVal)) %if no path lengths 
                maxTarr(j,i,k) = NaN;
                lArrStore(j,i,k) = NaN;
            else            
                maxTarr(j,i,k) = tempVal;
                lArrStore(j,i,k) = lArr(tempIndL(tempInd));
            end
            if(tempVal>tBlank) %if the longest path smaller than the maximum detectable path length arrives after the end of the blanking period
                checkSumPass(j,i,k) = true;
            else
                checkSumPass(j,i,k) = false;
            end
        end
    end
end

%% Create figure for checkSumPass matrix
figure(1)
pcolor(rSstore, zSstore,checkSumPass(:,:,n));
shading flat
ylabel('Source Depth [m]','FontSize',LABEL_FONTSIZE);
xlabel('Source Range [m]','FontSize',LABEL_FONTSIZE);
title(['Checksum Success (0) or Failure (1), Receiver Depth = ', num2str(zPstore(n))],'FontSize',TITLE_FONTSIZE);
set(gca,'FontSize',MARKER_FONTSIZE);
set(gca,'YDir','reverse');
colorbar
caxis([0 1]);