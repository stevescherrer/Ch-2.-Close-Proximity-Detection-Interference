%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI - Tank Experiment
%
%   Date: 6/6/2016
%
%   Description: 4) Calculate IR for each candidate source and receiver
%   combination.  Assume no interface or volumetric scattering, and first
%   20 multipaths.
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

fS = 192000;
arrivePath = [1 1 1 1   1 1 1 1   1 1 1 1   1 1 1 1   1 1 1 1];

%% Test 1A (June 2014, Deep Receiver)

% dW = 300;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 298];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1 199 399 578 766 989];
% zSstore = [270 297];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_1A_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% 
% %% Test 1B (June 2014, Shallow Receiver)
% 
% dW = 300;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 270];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1 199 399 578 766 989];
% zSstore = [270 297];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_1B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 2A (Nov 2014, Deep Receiver)
% 
% dW = 25;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 24];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1 75 150 300 600 1200];
% zSstore = [17.5];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_2A_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 2B (Nov 2014, Mid-Depth Receiver)
% 
% dW = 25;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 17.5];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1 75 150 300 600 1200];
% zSstore = [17.5];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_2B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 2C (Nov 2014, Shallow Receiver)
% 
% dW = 25;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 10];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1 75 150 300 600 1200];
% zSstore = [17.5];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_2C_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 3A (March 2015, Shallow Water Depth)
% 
% dW = 50;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 48];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1];
% zSstore = [48];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_3A_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 3B (March 2015, Deep Water Depth)
% 
% dW = 212;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 210];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [1];
% zSstore = [210];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_3B_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);
% 
% %% Test 4 (June 2015)
% 
% dW = 300;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 298];               %Test 2 Receiver, Closest to Sea Bed   
% 
% % rSstore = 1:rSres:300; 
% % zSstore = 1:zSres:24.5; 
% rSstore = [60 508];
% zSstore = [298];
% nrS = length(rSstore);
% nzS = length(zSstore);
% 
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
% % 
% temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
% save(['fishDaq4_4_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);

%% Shallow Case

% dW = 50;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 48];               %Test 2 Receiver, Closest to Sea Bed   
% 
% rSres = 20;
% zSres = 10;
% rSmin = 221; rSmax = 301;
% rSstore = rSmin:rSres:rSmax; nrS = length(rSstore);
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
% save(['fishDaq4_Shallow_rS',num2str(rSmin),'TO',num2str(rSmax),'_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);

%% Mid-Depth Case

% dW = 212;
% cW = 1530;
% tS = 0;
% tSync = 0;
% 
% R = [0, 0, 210];               %Test 2 Receiver, Closest to Sea Bed   
% 
% rSres = 20;
% zSres = 40;
% rSmin = 501; rSmax = 521;
% rSstore = rSmin:rSres:rSmax; nrS = length(rSstore);
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
% save(['fishDaq4_Mid_rS',num2str(rSmin),'TO',num2str(rSmax),'_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);

%% Deep Case

dW = 500;
cW = 1530;
tS = 0;
tSync = 0;

R = [0, 0, 498];               %Test 2 Receiver, Closest to Sea Bed   

rSres = 20;
zSres = 100;
rSmin = 501; rSmax = 581;
rSstore = rSmin:rSres:rSmax; nrS = length(rSstore);
zSstore = 1:zSres:498; nzS = length(zSstore);
nSim = nrS*nzS;

Sstore = NaN(nSim,2);
sampArrStore = NaN(nSim,sum(arrivePath));
tArrStore = sampArrStore;
RLstore = sampArrStore;

for i=1:nrS
    rS = rSstore(i);
    for j=1:nzS
        offset = (i-1)*nzS+j;
        zS = zSstore(j);
        Sstore(offset,:) = [rS, zS];
        
        m = [Sstore(offset,1),0,Sstore(offset,2),tS,tSync,R(1),R(2),R(3),dW,cW,cW];
        nArr = sum(arrivePath);
        [tArr,lArr] = MOI_XYZ(m,arrivePath); tArrStore(offset,:) = tArr'; lArr = lArr';
        sampArr = round(tArr*fS) + 1; %time=0 corresponds to the first sample
        sampArrStore(offset,:) = sampArr - min(sampArr) + 1; %each channel starts at first sample
        RL = 10.^(-1*log10(lArr)); %assumes omnidirectional source level of 0dB, no reflection losses/absorption/scattering
        RL = RL_PhaseInvert(arrivePath,RL);
        RLstore(offset,:) = RL;
    end
end

IRstore = cell(nSim,1);

for i=1:nSim
    hTemp = zeros(1,max(sampArrStore(i,:)));
    hTemp(sampArrStore(i,:)) = RLstore(i,:);
    IRstore{i} = hTemp;
end

temp = clock; y=temp(1); mo=temp(2); d=temp(3); h=temp(4); mi=temp(5); sec=floor(temp(6));
save(['fishDaq4_Deep_rS',num2str(rSmin),'TO',num2str(rSmax),'_',num2str(y),'_',num2str(mo),'_',num2str(d),'__',num2str(h),'_',num2str(mi),'_',num2str(sec)]);