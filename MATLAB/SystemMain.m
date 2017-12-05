function SystemMain()
close all;
%SystemMain Run the expierment.
%To Do:
% *Learning on other spurs
% *Learning on main carrier.
% *Case for if CCs are different BWs
% *Add suppport for TX/RX on different freqs

%% Simulations 

% Setup the data source
myLTE = LTE(5,'QPSK','uplink',200,-8e6);            % 5 MHz QPSK SCFDMA 200 symbols at -6 MHz in BB
myLTE = newcomponentcarrier(myLTE,2,5,'QPSK',8e6);  % 5 MHz QPSK SCFDMA 200 symbols at 6 MHz in BB
myLTE.sampleArray = myLTE.normalizeSignal(myLTE.sampleArray,1); %Normalize to be within [-0.7.0.7]
myLTE.signalWithDPD = myLTE.sampleArray;              % Initialize the with DPD signal to be the original signal.

% Setup the PA
myPA = PA(9);     

% Broadcast double pre signal
out = broadcast(myPA,myLTE.sampleArray);

% Plot the frequency domain of the signal
LTE.plot_freqdomain(out,myLTE.CCs.CC1.systemFs,'','No DPD',100);

% Double Precision. Setup DPD and perform learning
myDPD = SubBandDPD(myPA,myLTE,'IM3+',5,[0.5;0.5],1); % Use the myPA model on myLTE signal to perform 3rd order learning on IM3+ spur with \mu = 1

%Apply learned DPD to signal
myLTE.signalWithDPD = applyDPDtoSignal(myDPD,myLTE);

%Broadcast through PA
myLTE.signalWithDPD = broadcast(myPA,myLTE.signalWithDPD);

% Plot the frequency domain of the signal
LTE.plot_freqdomain(myLTE.signalWithDPD ,myLTE.CCs.CC1.systemFs,'','Double',100);

% % 8 bits
% myDAC = DAC(8,7);
% myFrontend = Frontend(myPA,myDAC);
% myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[0.25;2],1);
% DPDout1 = applyDPDtoSignal(myDPD,myLTE);
% signalWithDPD = broadcast(myPA,DPDout1);
% LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','8 bits',100);
%
 % 6 Bits
 %myDAC = DAC(6,5);
 %myFrontend = Frontend(myPA,myDAC);
 %myDPD = SubBandDPD(myPA,myLTE,'IM3+',5,[0.5;5],1); 
 %DPDout1 = applyDPDtoSignal(myDPD,myLTE);
 %signalWithDPD = broadcast(myPA,DPDout1);
 %LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','6 bits',100);

%  %% 4 bits
% myDAC = DAC(4,3);
% myFrontend = Frontend(myPA,myDAC);
% myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[0.25/2;1],1);
% DPDout1 = applyDPDtoSignal(myDPD,myLTE);
% signalWithDPD = broadcast(myPA,DPDout1);
% LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','4 bits',100);
%
%  %% 2 bits
%  myDAC = DAC(2,1);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','2 bits',100);

%
myDAC = DAC(1,0);
myFrontend = Frontend(myPA,myDAC);
myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[0.5;0.5],1);
DPDout1 = applyDPDtoSignal(myDPD,myLTE);
signalWithDPD = broadcast(myPA,DPDout1);
LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','1 bits',100);
legend('show')


%% WARP BOARD TEST
myLTE = LTE(1.4,'QPSK','uplink',200,-3e6);            % 5 MHz QPSK SCFDMA 200 symbols at -6 MHz in BB
myLTE = newcomponentcarrier(myLTE,2,1.4,'QPSK',3e6);  % 5 MHz QPSK SCFDMA 200 symbols at 6 MHz in BB
myLTE.sampleArray = myLTE.normalizeSignal(myLTE.sampleArray,0.6); %Normalize to be within [-0.7.0.7]
myLTE.signalWithDPD = myLTE.sampleArray;              % Initialize the with DPD signal to be the original signal.

% Setup the PA
myPA = WARP(1); %Set up WARP board

% Broadcast double pre signal
out = broadcast(myPA,myLTE.sampleArray);

% Plot the frequency domain of the signal
LTE.plot_freqdomain(out,myLTE.CCs.CC1.systemFs,'','No DPD',200);

% Double Precision. Setup DPD and perform learning
myDPD = SubBandDPD(myPA,myLTE,'IM3+',5,[1;100],0.5403 - 0.8415*i); % Use the myPA model on myLTE signal to perform 3rd order learning on IM3+ spur with \mu = 1

%Apply learned DPD to signal
myLTE.signalWithDPD = applyDPDtoSignal(myDPD,myLTE);

%Broadcast through PA
myLTE.signalWithDPD = broadcast(myPA,myLTE.signalWithDPD);

% Plot the frequency domain of the signal
LTE.plot_freqdomain(myLTE.signalWithDPD ,myLTE.CCs.CC1.systemFs,'','Double',200);

% 8 bits
%myDAC = DAC(8,7);
%myFrontend = Frontend(myPA,myDAC);
%myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[0.25;2],0.5403 - 0.8415*i);
%DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%signalWithDPD = broadcast(myPA,DPDout1);
%LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','8 bits',200);
%
%   %% 6 Bits
%   myDAC = DAC(6,5);
%   myFrontend = Frontend(myPA,myDAC);
%   myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[1;100],0.5403 - 0.8415*i);
%   DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%   signalWithDPD = broadcast(myPA,DPDout1);
%   LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','6 bits',200);

%  %% 4 bits
%myDAC = DAC(4,3);
%myFrontend = Frontend(myPA,myDAC);
%myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[0.25/2;1],0.5403 - 0.8415*i);
%DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%signalWithDPD = broadcast(myPA,DPDout1);
%LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','4 bits',200);
%
%  %% 2 bits
%  myDAC = DAC(2,1);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','2 bits',200);

%
myDAC = DAC(1,0);
myFrontend = Frontend(myPA,myDAC);
myDPD = SubBandDPD(myFrontend,myLTE,'IM3+',5,[2;400],0.5403 - 0.8415*i);
DPDout1 = applyDPDtoSignal(myDPD,myLTE);
signalWithDPD = broadcast(myPA,DPDout1);
LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','1 bits',200);
legend('show')

end

