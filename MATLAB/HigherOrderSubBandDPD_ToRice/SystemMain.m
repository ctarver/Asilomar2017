clear;
clc;
close all;

rng(0);
set(0,'defaultlinelinewidth',1.5)
Sparsity_Indx = 10;

%% Simulation Paramters %%
TxScenario                          =   3;     % 1- Contigious Single CC
                                               % 2- Multi-Cluster Single CC
                                               % 3- Intra-Band CA Two CC 
N_symbols                           =   160;    % No of OFDM Symbols for Simulation
TxSignalType                        =   2;     % 1- OFDM for DL Transmission
                                               % 2- SCFDMA for UL Transmission
ModulationType                      =   1;     % 1- QPSK
                                               % 2- 16-QAM  
                                               % 3- 64-QAM
POWER_PLOT_1MHZ                     =   1;     % Use 1 MHz power plot

% RF and PA paramters
Pout_dBm                            =   23;    % Desired Tx Output Power after the PA in dBm    

MemoryLessPA                        =   1;     % 1- MemoryLess PA model
                                               % 0- PH PA model with memory
MemoryLessDPD                       =   1;     % 1- MemoryLess DPD
                                               % 0- MemoryDPD

%% Baseband Tx signal paramters                                               
LTE_Bandwidth = 1.4; % Carrier BW of the LTE Tx Signal
NRB1 = 6;            % Number of RB allocated to first CC
NRB2 = 6;            % Number of RB allocated to second CC  
CarrierSpacing = 12; % Spacing in MHz between the 2 CC
IM3_Freq = 3*(CarrierSpacing/2);
IM5_Freq = 5*(CarrierSpacing/2);
IM7_Freq = 7*(CarrierSpacing/2);
IM9_Freq = 9*(CarrierSpacing/2);
Signal_Bandwidth = NRB1*0.18;

%% Power Amplifier Model
% Estimated PA parameters at 26 dBm output power, 1.88 GHz, using a 10 MHz signal with 120 MHz Fs
PA_Power_Measured = 26; 
MemorylessPA_Paramters = [0.9512 - 0.0946i;
                          4*0.0239 + 0.1632i;
                          4*0.0082 - 0.0727i;
                         -0.0016 + 0.0147i;
                         -0.0001 - 0.0011i];
% Memory depth = 4 per non-lineaity order
Memory_depth = 4;
MemoryPA_Paramters =  [0.5823 - 0.0608i;
                       1.1417 - 0.1198i;
                      -1.1964 + 0.1273i;
                       0.4264 - 0.0407i;
                       0.0074 + 0.1609i;
                       0.0292 + 0.0037i;
                      -0.0185 + 0.0002i;
                       0.0032 - 0.0016i;
                       0.0096 - 0.0727i;
                       0.0034 + 0.0010i;
                      -0.0060 - 0.0026i;
                       0.0028 + 0.0012i;
                      -0.0017 + 0.0149i;
                      -0.0011 - 0.0008i;
                       0.0014 + 0.0012i;
                      -0.0007 - 0.0005i;
                      -0.0001 - 0.0012i;
                       0.0001 + 0.0001i;
                      -0.0001 - 0.0001i;
                       0.0001 + 0.0001i];  


%% Baseband Equivilent LTE Signal Transmitter
[LTE_Signal, CC1, CC2, SystemFs, UpsamplingFactor] = LTE_Transmitter(LTE_Bandwidth,CarrierSpacing,NRB1,NRB2,...
                                                  N_symbols,TxSignalType,ModulationType,TxScenario);
% Scale the Baseband generated signal to have a unit RMS power
TX_PowerScale = sqrt(10^((PA_Power_Measured-Pout_dBm)/10));
RMS_PAin = sqrt(mean(abs(LTE_Signal).^2));
ScalingForPA = 1/(RMS_PAin*TX_PowerScale);                                              
PA_InputSignal = LTE_Signal*ScalingForPA; 
CC1 = CC1*ScalingForPA; 
CC2 = CC2*ScalingForPA;

%% Power Amplifier Model 
if MemoryLessPA
    PA_OutputSignal = MemoryLess_PA(PA_InputSignal,MemorylessPA_Paramters);
else
    PA_OutputSignal = MemoryPH_PA(PA_InputSignal,MemoryPA_Paramters,Memory_depth);
end

%% Decorrelating DPD Basis Functions Generation
% Generate the positive IM3 third and higher order basis functions
IM3_BasisThirdOrder   = conj(CC2).*(CC1.^2); 
IM3_BasisFifthOrder   = IM3_BasisThirdOrder.*(2*(abs(CC1).^2) + 3*(abs(CC2).^2));
IM3_BasisSeventhOrder = IM3_BasisThirdOrder.*(3*(abs(CC1).^4) + 6*(abs(CC2).^4) + ...
                                             12*(abs(CC1).^2).*(abs(CC2).^2));
IM3_BasisNinthOrder   = IM3_BasisThirdOrder.*(4*(abs(CC1).^6) + 10*(abs(CC2).^6) + ...
                                             30*(abs(CC1).^4).*(abs(CC2).^2) + ...
                                             40*(abs(CC1).^2).*(abs(CC2).^4)); 
% Orthogonalize the IM3 basis functions for better convergence in MBF Decorr DPD
IM3_Basis_5th = [IM3_BasisThirdOrder IM3_BasisFifthOrder];
IM3_Basis_7th = [IM3_BasisThirdOrder IM3_BasisFifthOrder IM3_BasisSeventhOrder];
IM3_Basis_9th = [IM3_BasisThirdOrder IM3_BasisFifthOrder IM3_BasisSeventhOrder IM3_BasisNinthOrder];
[Q,~] = qr(IM3_Basis_5th,0);
IM3_Basis_Orth_5th = Q.';
[Q,~] = qr(IM3_Basis_7th,0);
IM3_Basis_Orth_7th = Q.';
[Q,~] = qr(IM3_Basis_9th,0);
IM3_Basis_Orth_9th = Q.';

%% Negative IM3 basis functions 
IM3_Basis_3rdOrder  = conj(CC1).*(CC2.^2);
IM3_Basis_5thOrder  = IM3_Basis_3rdOrder.*(2*(abs(CC2).^2) + 3*(abs(CC1).^2));
IM3_Basis_7thOrder  = IM3_Basis_3rdOrder.*(3*(abs(CC2).^4) + 6*(abs(CC1).^4) + ...
                                             12*(abs(CC2).^2).*(abs(CC1).^2));
IM3_Basis_9thOrder  = IM3_Basis_3rdOrder.*(4*(abs(CC2).^6) + 10*(abs(CC1).^6) + ...
                                             30*(abs(CC2).^4).*(abs(CC1).^2) + ...
                                             40*(abs(CC2).^2).*(abs(CC1).^4)); 
IM3_Basis_9th_Neg = [IM3_Basis_3rdOrder IM3_Basis_5thOrder IM3_Basis_7thOrder IM3_Basis_9thOrder];
[Q,~] = qr(IM3_Basis_9th_Neg,0);
IM3_Basis_9th_Neg = Q.';

%% Negative IM5 basis functions 
IM5_Basis_5thOrder  = (conj(CC1).^2).*(CC2.^3); 
IM5_Basis_7thOrder  = IM5_Basis_5thOrder.*(4*(abs(CC1).^2) + 3*(abs(CC2).^2));
IM5_Basis_9thOrder  = IM5_Basis_5thOrder.*(10*(abs(CC1).^4) + 6*(abs(CC2).^4) + ...
                                          20*(abs(CC1).^2).*(abs(CC2).^2));
IM5_Basis_9th = [IM5_Basis_5thOrder IM5_Basis_7thOrder IM5_Basis_9thOrder];
[Q,~] = qr(IM5_Basis_9th,0);
IM5_Basis = Q.';

%% Negative IM7 and IM9 basis functions
IM7_Basis_7thOrder = (conj(CC1).^3).*(CC2.^4);
IM7_Basis_9thOrder  = IM7_Basis_7thOrder.*(5*(abs(CC1).^2) + 4*(abs(CC2).^2));

IM7_Basis_9th = [IM7_Basis_7thOrder IM7_Basis_9thOrder];
[Q,~] = qr(IM7_Basis_9th,0);
IM7_Basis = Q.';

IM9_Basis = (conj(CC1).^4).*(CC2.^5);

%% DPD Loop Delay Estimation
[LoopDelay, FeedBackFilter] = DPD_LoopDelayEst(PA_InputSignal, IM3_BasisThirdOrder, MemoryLessPA, ...
    SystemFs,Signal_Bandwidth,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,0);
                                           
%% IM3+ Block Decorrelating DPD Estimation
% Third Order IM3+ DPD learning
if MemoryLessDPD
    DPD_LearningRate = 0.05;           % Decorrelating DPD Learning Rate
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
else 
    DPD_LearningRate = 0.006;         % Decorrelating DPD Learning Rate
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
end
IM3_ThirdOrder_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM3_BasisThirdOrder.',MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay,....
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
disp('DPD Filtering BlockSize in msec = ');
disp(1e3*DPD_FilteringBlockSize/SystemFs); 

% Fifth Order IM3+ DPD learning
if MemoryLessDPD
    DPD_LearningRate = 40;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 800000;      % Number of samples used for learning
else 
    DPD_LearningRate = 80;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 800000;      % Number of samples used for learning
end

IM3_FifthOrder_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM3_Basis_Orth_5th,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM3_FifthOrder_Coeffs = IM3_FifthOrder_Coeffs.';

% Seventh Order IM3+ DPD learning
if MemoryLessDPD
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
else 
    DPD_LearningRate = 80;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 800000;      % Number of samples used for learning
end
IM3_SeventhOrder_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM3_Basis_Orth_7th,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM3_SeventhOrder_Coeffs = IM3_SeventhOrder_Coeffs.';

% Ninth Order IM3+ DPD learning
if MemoryLessDPD
    DPD_LearningRate = 10;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
else 
    DPD_LearningRate = 40;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
end
IM3_NinthOrder_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM3_Basis_Orth_9th,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM3_NinthOrder_Coeffs = IM3_NinthOrder_Coeffs.';

%% Applying IM3+ Decorrelating DPD
if  MemoryLessDPD
    % Applying Third Order IM3+ Decorr DPD
    AdaptiveFilterDelay  = length(IM3_ThirdOrder_Coeffs) - 1;
    IM3_DPD_Signal = conv(IM3_BasisThirdOrder,IM3_ThirdOrder_Coeffs);
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_ThirdOrderDPD   = PA_InputSignal + IM3_DPD_Signal;

    % Applying Fifth Order IM3+ Decorr DPD 
    AdaptiveFilterDelay = length(IM3_FifthOrder_Coeffs(1:2:2)) - 1;
    H53 = IM3_FifthOrder_Coeffs(1:2:end);
    IM3_DPD_Signal = conv(IM3_Basis_Orth_5th(1,:),IM3_FifthOrder_Coeffs(1:2:2)) + ...
                     conv(IM3_Basis_Orth_5th(2,:),IM3_FifthOrder_Coeffs(2:2:2));
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_FifthOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));

    % Applying Seventh Order IM3+ Decorr DPD 
    AdaptiveFilterDelay  = length(IM3_SeventhOrder_Coeffs(1:3:3)) - 1;
    IM3_DPD_Signal = conv(IM3_Basis_Orth_7th(1,:),IM3_SeventhOrder_Coeffs(1:3:3)) + ...
                     conv(IM3_Basis_Orth_7th(2,:),IM3_SeventhOrder_Coeffs(2:3:3)) + ...
                     conv(IM3_Basis_Orth_7th(3,:),IM3_SeventhOrder_Coeffs(3:3:3));
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_SeventhOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));

    % Applying Ninth Order IM3+ Decorr DPD 
    AdaptiveFilterDelay  = length(IM3_NinthOrder_Coeffs(1:4:4)) - 1;
    IM3_DPD_Signal = conv(IM3_Basis_Orth_9th(1,:),IM3_NinthOrder_Coeffs(1:4:4)) + ...
                     conv(IM3_Basis_Orth_9th(2,:),IM3_NinthOrder_Coeffs(2:4:4)) + ...
                     conv(IM3_Basis_Orth_9th(3,:),IM3_NinthOrder_Coeffs(3:4:4)) + ...
                     conv(IM3_Basis_Orth_9th(4,:),IM3_NinthOrder_Coeffs(4:4:4));
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_NinthOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));
else
    % Applying Third Order IM3+ Decorr DPD
    H33 = [IM3_ThirdOrder_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM3_ThirdOrder_Coeffs(2); zeros(Sparsity_Indx-1,1)];
    AdaptiveFilterDelay  = length(H33) - 1;
    IM3_DPD_Signal = conv(IM3_BasisThirdOrder,H33);
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_ThirdOrderDPD   = PA_InputSignal + IM3_DPD_Signal;

    % Applying Fifth Order IM3+ Decorr DPD 
    H53 = [IM3_FifthOrder_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM3_FifthOrder_Coeffs(3); zeros(Sparsity_Indx-1,1)];
    H55 = [IM3_FifthOrder_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM3_FifthOrder_Coeffs(4); zeros(Sparsity_Indx-1,1)];
    AdaptiveFilterDelay = length(H55) - 1;
    IM3_DPD_Signal = conv(IM3_Basis_Orth_5th(1,:),H53) + ...
                     conv(IM3_Basis_Orth_5th(2,:),H55);
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_FifthOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));

    % Applying Seventh Order IM3+ Decorr DPD 
    H73 = [IM3_SeventhOrder_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM3_SeventhOrder_Coeffs(4); zeros(Sparsity_Indx-1,1)];
    H75 = [IM3_SeventhOrder_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM3_SeventhOrder_Coeffs(5); zeros(Sparsity_Indx-1,1)];
    H77 = [IM3_SeventhOrder_Coeffs(3); zeros(Sparsity_Indx-1,1); ...
           IM3_SeventhOrder_Coeffs(6); zeros(Sparsity_Indx-1,1)];
    AdaptiveFilterDelay = length(H77) - 1;
    IM3_DPD_Signal = conv(IM3_Basis_Orth_7th(1,:),H73) + ...
                     conv(IM3_Basis_Orth_7th(2,:),H75) + ...
                     conv(IM3_Basis_Orth_7th(3,:),H77);
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_SeventhOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));

    % Applying Ninth Order IM3+ Decorr DPD 
    H93 = [IM3_NinthOrder_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM3_NinthOrder_Coeffs(5); zeros(Sparsity_Indx-1,1)];
    H95 = [IM3_NinthOrder_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM3_NinthOrder_Coeffs(6); zeros(Sparsity_Indx-1,1)];
    H97 = [IM3_NinthOrder_Coeffs(3); zeros(Sparsity_Indx-1,1); ...
           IM3_NinthOrder_Coeffs(7); zeros(Sparsity_Indx-1,1)];
    H99 = [IM3_NinthOrder_Coeffs(4); zeros(Sparsity_Indx-1,1); ...
           IM3_NinthOrder_Coeffs(8); zeros(Sparsity_Indx-1,1)];   
    AdaptiveFilterDelay  = length(H99) - 1;
    IM3_DPD_Signal = conv(IM3_Basis_Orth_9th(1,:),H93) + ...
                     conv(IM3_Basis_Orth_9th(2,:),H95) + ...
                     conv(IM3_Basis_Orth_9th(3,:),H97) + ...
                     conv(IM3_Basis_Orth_9th(4,:),H99);
    IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_IM3_NinthOrderDPD   = PA_InputSignal + IM3_DPD_Signal(1:length(PA_InputSignal));
end
    
%% Negative IM3 Block Decorrelating DPD Estimation and Filtering
% Ninth Order IM3- DPD learning
if MemoryLessDPD
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 400000;      % Number of samples used for learning
else 
    DPD_LearningRate = 40;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
end
IM3_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM3_Basis_9th_Neg,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,-IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM3_Coeffs = IM3_Coeffs.';

% Applying Ninth Order IM3- Decorr DPD 
if  MemoryLessDPD
    AdaptiveFilterDelay  = length(IM3_NinthOrder_Coeffs(1:4:end)) - 1;
    IM3_DPD_Signal_Neg = conv(IM3_Basis_9th_Neg(1,:),IM3_Coeffs(1:4:end)) + ...
                     conv(IM3_Basis_9th_Neg(2,:),IM3_Coeffs(2:4:end)) + ...
                     conv(IM3_Basis_9th_Neg(3,:),IM3_Coeffs(3:4:end)) + ...
                     conv(IM3_Basis_9th_Neg(4,:),IM3_Coeffs(4:4:end));
    IM3_DPD_Signal_Neg = IM3_DPD_Signal_Neg(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal_Neg = IM3_DPD_Signal_Neg.*exp(-2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_NegIM3_DPD   = PA_InputSignal + IM3_DPD_Signal_Neg(1:length(PA_InputSignal));
else
    H93 = [IM3_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM3_Coeffs(5); zeros(Sparsity_Indx-1,1)];
    H95 = [IM3_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM3_Coeffs(6); zeros(Sparsity_Indx-1,1)];
    H97 = [IM3_Coeffs(3); zeros(Sparsity_Indx-1,1); ...
           IM3_Coeffs(7); zeros(Sparsity_Indx-1,1)];
    H99 = [IM3_Coeffs(4); zeros(Sparsity_Indx-1,1); ...
           IM3_Coeffs(8); zeros(Sparsity_Indx-1,1)];   
    AdaptiveFilterDelay  = length(H99) - 1;
    IM3_DPD_Signal_Neg = conv(IM3_Basis_9th_Neg(1,:),H93) + ...
                     conv(IM3_Basis_9th_Neg(2,:),H95) + ...
                     conv(IM3_Basis_9th_Neg(3,:),H97) + ...
                     conv(IM3_Basis_9th_Neg(4,:),H99);
    IM3_DPD_Signal_Neg = IM3_DPD_Signal_Neg(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM3_DPD_Signal)).';
    IM3_DPD_Signal_Neg = IM3_DPD_Signal_Neg.*exp(-2*pi*1i*n*IM3_Freq*1e6/SystemFs);
    PAin_NegIM3_DPD   = PA_InputSignal + IM3_DPD_Signal_Neg(1:length(PA_InputSignal));             
end

if 1
%% Negative IM5 Block Decorrelating DPD Estimation and Filtering
% Ninth Order IM5- DPD learning
if MemoryLessDPD
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 400000;      % Number of samples used for learning
else 
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
end
IM5_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM5_Basis,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,-IM5_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM5_Coeffs = IM5_Coeffs.';

% Applying Ninth Order IM5- Decorr DPD 
if  MemoryLessDPD
    AdaptiveFilterDelay  = length(IM5_Coeffs(1:3:end)) - 1;
    IM5_DPD_Signal = conv(IM5_Basis(1,:),IM5_Coeffs(1:3:end)) + ...
                     conv(IM5_Basis(2,:),IM5_Coeffs(2:3:end)) + ...
                     conv(IM5_Basis(3,:),IM5_Coeffs(3:3:end));
    IM5_DPD_Signal = IM5_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM5_DPD_Signal)).';
    IM5_DPD_Signal = IM5_DPD_Signal.*exp(-2*pi*1i*n*IM5_Freq*1e6/SystemFs);
    PAin_IM5_DPD   = PA_InputSignal + IM5_DPD_Signal(1:length(PA_InputSignal));
else
    H95 = [IM5_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM5_Coeffs(4); zeros(Sparsity_Indx-1,1)];
    H97 = [IM5_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM5_Coeffs(5); zeros(Sparsity_Indx-1,1)];
    H99 = [IM5_Coeffs(3); zeros(Sparsity_Indx-1,1); ...
           IM5_Coeffs(6); zeros(Sparsity_Indx-1,1)]; 
    AdaptiveFilterDelay  = length(H99) - 1;
    IM5_DPD_Signal = conv(IM5_Basis(1,:),H95) + ...
                     conv(IM5_Basis(2,:),H97) + ...
                     conv(IM5_Basis(3,:),H99);
    IM5_DPD_Signal = IM5_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM5_DPD_Signal)).';
    IM5_DPD_Signal = IM5_DPD_Signal.*exp(-2*pi*1i*n*IM5_Freq*1e6/SystemFs);
    PAin_IM5_DPD   = PA_InputSignal + IM5_DPD_Signal(1:length(PA_InputSignal));
end
%% Negative IM7 Block Decorrelating DPD Estimation and Filtering
% Ninth Order IM7- DPD learning
if MemoryLessDPD
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 400000;      % Number of samples used for learning
else 
    DPD_LearningRate = 20;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 1200000;      % Number of samples used for learning
end
IM7_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM7_Basis,MemoryLessPA,MemoryLessDPD, ...
      SystemFs,-IM7_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);
IM7_Coeffs = IM7_Coeffs.';

% Applying Ninth Order IM7- Decorr DPD 
if  MemoryLessDPD
    AdaptiveFilterDelay  = length(IM7_Coeffs(1:2:end)) - 1;
    IM7_DPD_Signal = conv(IM7_Basis(1,:),IM7_Coeffs(1:2:end)) + ...
                     conv(IM7_Basis(2,:),IM7_Coeffs(2:2:end));
    IM7_DPD_Signal = IM7_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM7_DPD_Signal)).';
    IM7_DPD_Signal = IM7_DPD_Signal.*exp(-2*pi*1i*n*IM7_Freq*1e6/SystemFs);
    PAin_IM7_DPD   = PA_InputSignal + IM7_DPD_Signal(1:length(PA_InputSignal));
else
    H97 = [IM7_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM7_Coeffs(3); zeros(Sparsity_Indx-1,1)];
    H99 = [IM7_Coeffs(2); zeros(Sparsity_Indx-1,1); ...
           IM7_Coeffs(4); zeros(Sparsity_Indx-1,1)];  
    AdaptiveFilterDelay  = length(H99) - 1;
    IM7_DPD_Signal = conv(IM7_Basis(1,:),H97) + ...
                     conv(IM7_Basis(2,:),H99);
    IM7_DPD_Signal = IM7_DPD_Signal(1:end-AdaptiveFilterDelay).';
    n = (1:length(IM7_DPD_Signal)).';
    IM7_DPD_Signal = IM7_DPD_Signal.*exp(-2*pi*1i*n*IM7_Freq*1e6/SystemFs);
    PAin_IM7_DPD   = PA_InputSignal + IM7_DPD_Signal(1:length(PA_InputSignal));
end
if 0
%% Negative IM9 Block Decorrelating DPD Estimation and Filtering
% Ninth Order IM9- DPD learning
if MemoryLessDPD
    DPD_LearningRate = 0.1;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 400000;      % Number of samples used for learning
else 
    DPD_LearningRate = 0.001;         % Decorrelating DPD Learning Rate 
    DPD_LearningBlockSize  = 1000;    % Decorrelating DPD Learning Block size
    DPD_FilteringBlockSize = 1000;    % Decorrelating DPD Filtering Block size
    NumLearningSamples = 800000;      % Number of samples used for learning
end
IM9_Coeffs = BlockDecorrDPD_MEM(PA_InputSignal,IM9_Basis.',MemoryLessPA,MemoryLessDPD, ...
      SystemFs,-IM9_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay, ...
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,0,Sparsity_Indx);

% Applying Ninth Order IM9- Decorr DPD
if  MemoryLessDPD
    AdaptiveFilterDelay  = length(IM9_Coeffs(1:2:end)) - 1;
    IM9_DPD_Signal = conv(IM9_Basis,IM9_Coeffs);
    n = (1:length(IM9_DPD_Signal)).';
    IM9_DPD_Signal = IM9_DPD_Signal.*exp(-2*pi*1i*n*IM9_Freq*1e6/SystemFs);
    PAin_IM9_DPD   = PA_InputSignal + IM9_DPD_Signal(1:length(PA_InputSignal));
else
    H99 = [IM9_Coeffs(1); zeros(Sparsity_Indx-1,1); ...
           IM7_Coeffs(3); zeros(Sparsity_Indx-1,1)];  
    AdaptiveFilterDelay  = length(H99) - 1;
    IM9_DPD_Signal = conv(IM9_Basis,H99);
    n = (1:length(IM9_DPD_Signal)).';
    IM9_DPD_Signal = IM9_DPD_Signal.*exp(-2*pi*1i*n*IM9_Freq*1e6/SystemFs);
    PAin_IM9_DPD   = PA_InputSignal + IM9_DPD_Signal(1:length(PA_InputSignal));             
end
end
end

%% Power Amplifier Model  
if MemoryLessPA
    PA_Output_IM3_ThirdOrderDPD   = MemoryLess_PA(PAin_IM3_ThirdOrderDPD,MemorylessPA_Paramters);
    PA_Output_IM3_FifthOrderDPD   = MemoryLess_PA(PAin_IM3_FifthOrderDPD,MemorylessPA_Paramters);
    PA_Output_IM3_SeventhOrderDPD = MemoryLess_PA(PAin_IM3_SeventhOrderDPD,MemorylessPA_Paramters);
    PA_Output_IM3_NinthOrderDPD   = MemoryLess_PA(PAin_IM3_NinthOrderDPD,MemorylessPA_Paramters);
    
    PA_Output_IM3_DPD   = MemoryLess_PA(PAin_NegIM3_DPD,MemorylessPA_Paramters);
    PA_Output_IM5_DPD   = MemoryLess_PA(PAin_IM5_DPD,MemorylessPA_Paramters);
    PA_Output_IM7_DPD   = MemoryLess_PA(PAin_IM7_DPD,MemorylessPA_Paramters);
    % PA_Output_IM9_DPD   = MemoryLess_PA(PAin_IM9_DPD,MemorylessPA_Paramters);

else
    PA_Output_IM3_ThirdOrderDPD   = MemoryPH_PA(PAin_IM3_ThirdOrderDPD,MemoryPA_Paramters,Memory_depth);
    PA_Output_IM3_FifthOrderDPD   = MemoryPH_PA(PAin_IM3_FifthOrderDPD,MemoryPA_Paramters,Memory_depth);
    PA_Output_IM3_SeventhOrderDPD = MemoryPH_PA(PAin_IM3_SeventhOrderDPD,MemoryPA_Paramters,Memory_depth);
    PA_Output_IM3_NinthOrderDPD   = MemoryPH_PA(PAin_IM3_NinthOrderDPD,MemoryPA_Paramters,Memory_depth);
    
    PA_Output_IM3_DPD   = MemoryPH_PA(PAin_NegIM3_DPD,MemoryPA_Paramters,Memory_depth);
    PA_Output_IM5_DPD   = MemoryPH_PA(PAin_IM5_DPD,MemoryPA_Paramters,Memory_depth);
    PA_Output_IM7_DPD   = MemoryPH_PA(PAin_IM7_DPD,MemoryPA_Paramters,Memory_depth);
    % PA_Output_IM9_DPD   = MemoryPH_PA(PAin_IM9_DPD,MemoryPA_Paramters,Memory_depth);
end

%% Spectral plots
if 0
figure(100);
plot_freqdomain(PA_OutputSignal,SystemFs,'','r',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM3_ThirdOrderDPD,SystemFs,'','b',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM3_FifthOrderDPD,SystemFs,'','g',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM3_SeventhOrderDPD,SystemFs,'','m',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM3_NinthOrderDPD,SystemFs,'','k',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
grid on;
L2= legend('No DPD','Third order $IM3+$ DPD','Fifth order $IM3+$ DPD',...
                    'Seventh order $IM3+$ DPD','Ninth order $IM3+$ DPD');
set(L2,'Interpreter','Latex');
axis([-60 60 -Inf 25])
setStyle(gcf, 'ieee_column', false)

figure(200);
plot_freqdomain(PA_OutputSignal,SystemFs,'','r',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM3_DPD,SystemFs,'','b',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM5_DPD,SystemFs,'','g',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
hold on;
plot_freqdomain(PA_Output_IM7_DPD,SystemFs,'','k',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
% hold on;
% plot_freqdomain(PA_Output_IM9_DPD,SystemFs,'','k',UpsamplingFactor,POWER_PLOT_1MHZ,PA_Power_Measured);
grid on;
L2= legend('No DPD','Ninth order $IM3-$ DPD','Ninth order $IM5-$ DPD',...
                    'Ninth order $IM7-$ DPD');
set(L2,'Interpreter','Latex');
axis([-60 60 -Inf 25])
setStyle(gcf, 'ieee_column', false)
end

%% Plot with markers
figure(106);
linewidth = 1.5;
MarkerSize = 8;
Y_Shift = 100;    % Shift the Y-axis to plot the figure with marker for legend
Step = 200e4;     % Step between markers
Step_1 = 1.36e6;  % Initial step

L = length(PA_OutputSignal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of Rx_Signal
SubcarrierSpacing = SystemFs/NFFT;
Bins_1MHz = round(1e6/SubcarrierSpacing);
f = (SystemFs/2)*linspace(-1,1,NFFT);
Power_Rx_Signal = 10*log10(mean(abs(PA_OutputSignal).^2));
PwrScaleFactor = 10^((Pout_dBm - Power_Rx_Signal)*0.1);

SignalToPlot = PA_OutputSignal;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_1 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM3_ThirdOrderDPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_2 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM3_FifthOrderDPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_3 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM3_SeventhOrderDPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_4 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM3_NinthOrderDPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_5 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));


plot(f(1:Step:end)/1e6,PSD_1MHz_1(1:Step:end)-Y_Shift,'-bs','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_2(1:Step:end)-Y_Shift,'-ro','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_3(1:Step:end)-Y_Shift,'-mx','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_4(1:Step:end)-Y_Shift,'-gv','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_5(1:Step:end)-Y_Shift,'-ko','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
L2 = legend('No DPD','$3^{rd}$ order IM3+ DPD','$5^{th}$ order IM3+ DPD','$7^{th}$ order IM3+ DPD','$9^{th}$ order IM3+ DPD');
set(L2,'Interpreter','Latex');
xlabel('Frequency (MHz)')
ylabel('Power (dBm/MHz)')
plot(f/1e6,PSD_1MHz_1,'-b','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_2,'-r','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_3,'-m','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_4,'-g','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_5,'-k','LineWidth',linewidth);hold on;
plot(f(Step_1:Step:end)/1e6,PSD_1MHz_1(Step_1:Step:end),'bs','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1:Step:end)/1e6,PSD_1MHz_2(Step_1:Step:end),'ro','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1:Step:end)/1e6,PSD_1MHz_3(Step_1:Step:end),'mx','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1:Step:end)/1e6,PSD_1MHz_4(Step_1:Step:end),'gv','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1:Step:end)/1e6,PSD_1MHz_5(Step_1:Step:end),'ko','MarkerSize',MarkerSize,'LineWidth',linewidth);hold off;
setStyle(gcf, 'ieee_column', false)
grid on;
axis([-Inf Inf -80 25])

figure(206);
Step = 21.2e4;     % Step between markers
Step_1 = 10.3e4;  % Initial step

SignalToPlot = PA_Output_IM3_DPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_2 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM5_DPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_3 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

SignalToPlot = PA_Output_IM7_DPD;
SignalToPlot = SignalToPlot*sqrt(PwrScaleFactor);
Rx_Signal_FD = fftshift(fft(SignalToPlot,NFFT)/L);
Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
PSD_1MHz_4 = 10*log10(Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving'));

plot(f(1:Step:end)/1e6,PSD_1MHz_1(1:Step:end)-Y_Shift,'-bs','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_2(1:Step:end)-Y_Shift,'-ro','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_3(1:Step:end)-Y_Shift,'-mx','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(1:Step:end)/1e6,PSD_1MHz_4(1:Step:end)-Y_Shift,'-ko','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
L2 = legend('No DPD','$9^{th}$ order IM3- DPD','$9^{th}$ order IM5- DPD','$9^{th}$ order IM7- DPD');
set(L2,'Interpreter','Latex');
xlabel('Frequency (MHz)')
ylabel('Power (dBm/MHz)')
plot(f/1e6,PSD_1MHz_1,'-b','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_2,'-r','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_3,'-m','LineWidth',linewidth);hold on;
plot(f/1e6,PSD_1MHz_4,'-k','LineWidth',linewidth);hold on;
plot(f(Step_1:Step:Step_1+3*Step)/1e6,PSD_1MHz_1(Step_1:Step:Step_1+3*Step),'bs','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1+3*Step)/1e6,PSD_1MHz_2(Step_1+3*Step),'ro','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1+2*Step)/1e6,PSD_1MHz_3(Step_1+2*Step),'mx','MarkerSize',MarkerSize,'LineWidth',linewidth);hold on;
plot(f(Step_1+Step)/1e6,PSD_1MHz_4(Step_1+Step),'ko','MarkerSize',MarkerSize,'LineWidth',linewidth);hold off;
setStyle(gcf, 'ieee_column', false)
grid on;
axis([-Inf Inf -80 25])


