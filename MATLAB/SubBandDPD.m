classdef SubBandDPD
   %SubBandDPD Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      learningBlockLength
      filteringBlockLength
      learningRateMu
      nLearningSamples
      spur
      dpdFreq
      feedbackFilter
      lmsSignal
      alpha
      loopDelay
      normalizedFreqShift
      DPD_Coeff
      correlation
      IM3powers
      phaseshiftDPD
   end
   
   methods
      function obj = SubBandDPD(pa,signal,spur)
         % Copy spur to object
         obj.spur = spur;
         
         % Use signal and spur to calculate the freq of the spur
         obj = getDPDFreq(obj,signal);
         
         % We'll set it up to loop from spur to maximum order. Change the
         % basis function depending on what has been done.
         
         obj.learningBlockLength  = 1024*4;
         obj.filteringBlockLength = 1024*4;
         obj.learningRateMu       = 0.25;
         obj.nLearningSamples     = 180000;
         obj.alpha                = 0;
         obj.phaseshiftDPD        = 1;%(0.5403 - 0.8415*i);  %Need for WARP;
         
         % Set up filter
         obj = makeFeedbackFilter(obj,signal);
         
         %Calculate the LMS basis function
         obj = getLMSBasis(obj,signal);
         
         % Go perform learning.
         obj = performDpdLearning(obj,pa,signal);
         %obj.alpha = 0;
         %obj = performDpdLearning3(obj,pa,signal);
         %Go perform fixed point learning.
         %obj = performFixedPointDpdLearning(obj,pa,signal,myDAC,12,11);
         
      end
      function obj = getDPDFreq(obj,signal)
         %Calculate Carrier Spacing
         carrierSpacing = signal.CCs.CC2.centerFreq - ...
            signal.CCs.CC1.centerFreq;
         switch obj.spur
            case 'IM3+'
               obj.dpdFreq  =  (3/2) * carrierSpacing;
            case 'IM3-'
               obj.dpdFreq  = -(3/2) * carrierSpacing;
            case 'IM5+'
               obj.dpdFreq  =  (5/2) * carrierSpacing;
            case 'IM5-'
               obj.dpdFreq  = -(5/2) * carrierSpacing;
         end
      end
      function obj = getLMSBasis(obj,signal)
         CC1 = signal.CCs.CC1.signalArray;
         CC2 = signal.CCs.CC2.signalArray;
         
         switch obj.spur
            case 'IM3+'
               obj.lmsSignal =  conj(CC1).*(CC2.^2);
            case 'IM3-'
               obj.lmsSignal  =  conj(CC1).*(CC2.^2);
            case 'IM5+'
               %Need to fill this in
            case 'IM5-'
               %Need to fill this in
         end
      end
      
      function obj = makeFeedbackFilter(obj,signal)
         % Third order IMD extraction filter, all frequency values are in Hz.
         signalBW = signal.CCs.CC1.signalBandwidth;
         Fs    = signal.CCs.CC1.systemFs;  % Sampling Frequency
         N     = 200;  % Order
         Fpass = (3*signalBW)/2; % Passband Frequency
         Fstop = (5*signalBW)/2; % Stopband Frequency
         Wpass = 1;    % Passband Weight
         Wstop = 5;    % Stopband Weight
         dens  = 20;   % Density Factor
         % Calculate the coefficients using the FIRPM function.
         obj.feedbackFilter = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
            {dens});
         
         %Estimate the loopdelay of the filter
         obj.loopDelay = mean(grpdelay(obj.feedbackFilter));
      end
      function obj = performDpdLearning(obj,pa,signal)
         %Set some things up
         dpdBlockIndx = 0;
         obj.DPD_Coeff(1) = obj.alpha;
         obj.normalizedFreqShift = (obj.dpdFreq)/signal.CCs.CC1.systemFs;
         
         if obj.nLearningSamples > length(obj.lmsSignal)
            obj.nLearningSamples = length(obj.lmsSignal)-obj.filteringBlockLength;
         end
         
         for sample = 1:obj.filteringBlockLength:obj.nLearningSamples
           dpdBlockIndx    = dpdBlockIndx + 1;
           obj = processBlock(obj,signal,pa,sample,dpdBlockIndx);
         end
%          figure();
%          plot(0:length(obj.DPD_Coeff)-1,real(obj.DPD_Coeff));
%          hold on
%          plot(0:length(obj.DPD_Coeff)-1,imag(obj.DPD_Coeff));
%          
%          figure()
%          plot(real(obj.correlation))
%          hold on
%          plot(imag(obj.correlation))         
   
      end
      function obj = performNewtonDpdLearning(obj,pa,signal)
         %Set some things up
         dpdBlockIndx = 0;
         obj.DPD_Coeff(1) = obj.alpha;
         obj.normalizedFreqShift = (obj.dpdFreq)/signal.CCs.CC1.systemFs;
         
         if obj.nLearningSamples > length(obj.lmsSignal)
            obj.nLearningSamples = length(obj.lmsSignal)-obj.filteringBlockLength;
         end
         
         for sample = 1:obj.filteringBlockLength:obj.nLearningSamples
           dpdBlockIndx    = dpdBlockIndx + 1;
           obj = processBlock(obj,signal,pa,sample,dpdBlockIndx);
         end
         figure();
         plot(0:length(obj.DPD_Coeff)-1,real(obj.DPD_Coeff));
         hold on
         plot(0:length(obj.DPD_Coeff)-1,imag(obj.DPD_Coeff));
   
      end      
      function [obj,dpdBlockIndx]  = processBlock(obj,signal,pa,sample,dpdBlockIndx)
            IM3Filter    = obj.feedbackFilter; %filter function prefers this
            PA_In_StartIndx = sample;
            PA_In_EndIndx   = sample + ...
               obj.filteringBlockLength + obj.loopDelay - 1;
            
            % Apply LMS block filtering
            LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
            dpdSignal = (obj.alpha*LMS_FilteringBlock);
            
            % Add the LMS filter output to the PA input
            signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
            PA_InBlock = signalInput ...
               + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
            
            %Brodcast block
            PA_OutBlock = broadcast(pa,PA_InBlock);
            obj.IM3powers(dpdBlockIndx,:) = IM3Power(obj,PA_OutBlock,signal.CCs.CC1.systemFs);
            
            % Shift the PA output such that the IM3 frequency is at baseband
            PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
            
            % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
            IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
            
            % Extract the IM3 block from the PA output in the feedback receiver
            ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
            
            %% Perform Correlation Calculation
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = sample;
            LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
            LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
            MeanCorrelation = MeanCorrelation * obj.phaseshiftDPD;  %Constant shift for this freq on WARP
            
            obj.correlation(dpdBlockIndx,:) = MeanCorrelation; 
            
            % Block LMS filter update
            W = obj.alpha' - (obj.learningRateMu).*MeanCorrelation;
            
            % Store DPD filter coeff.
            obj.DPD_Coeff(dpdBlockIndx+1,:) = W';  %+1 because we started on the 0th block with alpha = 0;
            obj.alpha = W';
      end
      
      function obj = performDpdLearning2(obj,pa,signal)
         %Set some things up
         dpdBlockIndx = 0;
         IM3Filter    = obj.feedbackFilter; %filter function prefers this
         W = 0;
         
         obj.normalizedFreqShift = obj.dpdFreq/signal.CCs.CC1.systemFs;
         
         if obj.nLearningSamples > length(obj.lmsSignal)
            obj.nLearningSamples = length(obj.lmsSignal)-obj.filteringBlockLength;
         end
         
         %% Part 1 of Newton
         sample = 1;
         PA_In_StartIndx = sample;
         PA_In_EndIndx   = sample + ...
            obj.filteringBlockLength + obj.loopDelay - 1;
         dpdBlockIndx    = dpdBlockIndx + 1;
         
         % Apply LMS block filtering
         LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
         dpdSignal = (obj.alpha*LMS_FilteringBlock);
         
         % Add the LMS filter output to the PA input
         signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
         PA_InBlock = signalInput ...
            + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
         
         %Brodcast block
         PA_OutBlock = broadcast(pa,PA_InBlock);
         
         % Shift the PA output such that the IM3 frequency is at baseband
         PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
         
         % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
         IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
         
         % Extract the IM3 block from the PA output in the feedback receiver
         ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
         
         % Perform Correlation Calculation
         % Prepare block of data to be used for learning
         LMS_In_StartIndx = sample;
         LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
         LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
         
         % Correlation between filter input and error block
         MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
         correlation(dpdBlockIndx,:) = MeanCorrelation;
         
         % Block LMS filter update
         W = W - (obj.learningRateMu).*MeanCorrelation;
         change1 = - (obj.learningRateMu).*MeanCorrelation';
         
         % Store DPD filter coeff.
         DPD_Coeff(dpdBlockIndx,:) = W';
         obj.alpha = double(W');
         sample = sample + obj.filteringBlockLength;
         
         %% Part 2 of Newton
         PA_In_StartIndx = sample;
         PA_In_EndIndx   = sample + ...
            obj.filteringBlockLength + obj.loopDelay - 1;
         dpdBlockIndx    = dpdBlockIndx + 1;
         
         % Apply LMS block filtering
         LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
         dpdSignal = (obj.alpha*LMS_FilteringBlock);
         
         % Add the LMS filter output to the PA input
         signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
         PA_InBlock = signalInput ...
            + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
         
         %Brodcast block
         PA_OutBlock = broadcast(pa,PA_InBlock);
         
         % Shift the PA output such that the IM3 frequency is at baseband
         PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
         
         % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
         IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
         
         % Extract the IM3 block from the PA output in the feedback receiver
         ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
         
         % Perform Correlation Calculation
         % Prepare block of data to be used for learning
         LMS_In_StartIndx = sample;
         LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
         LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
         
         % Correlation between filter input and error block
         MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
         correlation(dpdBlockIndx,:) = MeanCorrelation;
         % Block LMS filter update
         
         change2 = - (obj.learningRateMu).*MeanCorrelation';
         
         %Newtons Method
         m = change2-change1;
         real_x_intercept = -real(change1)/real(m);
         integral_real = 0.5*real_x_intercept*real(change1);
         
         imag_x_intercept = -imag(change1)/imag(m);
         integral_imag = 0.5*imag_x_intercept*imag(change1);
         
         W = change1' + change2' + (integral_real - integral_imag*i);
         %W = W + (integral_real - integral_imag*i);
         DPD_Coeff(dpdBlockIndx,:) = W';
         obj.alpha = W';
         sample = sample + obj.filteringBlockLength;
         
         for sample = sample:obj.filteringBlockLength:obj.nLearningSamples
            PA_In_StartIndx = sample;
            PA_In_EndIndx   = sample + ...
               obj.filteringBlockLength + obj.loopDelay - 1;
            dpdBlockIndx    = dpdBlockIndx + 1;
            
            % Apply LMS block filtering
            LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
            %            LMS_FilteringBlock = fi(LMS_FilteringBlock,1,12,11);
            dpdSignal = (obj.alpha*LMS_FilteringBlock);
            
            % Add the LMS filter output to the PA input
            signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
            %             signalInput = fi(signalInput,1,12,11);
            PA_InBlock = signalInput ...
               + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
            
            %Brodcast block
            PA_OutBlock = broadcast(pa,PA_InBlock);
            
            % Shift the PA output such that the IM3 frequency is at baseband
            PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
            
            % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
            IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
            %             IM3FilteredBlock = fi(IM3FilteredBlock,1,12,11);
            
            % Extract the IM3 block from the PA output in the feedback receiver
            ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
            
            %% Perform Correlation Calculation
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = sample;
            LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
            LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
            correlation(dpdBlockIndx,:) = MeanCorrelation;
            
            % Block LMS filter update
            W = W - (obj.learningRateMu).*MeanCorrelation;
            
            % Store DPD filter coeff.
            DPD_Coeff(dpdBlockIndx,:) = W';
            obj.alpha = double(W');
         end
         
         
         figure();
         plot(real(DPD_Coeff));
         hold on;
         plot(imag(DPD_Coeff));
         figure();
         plot(real(correlation));
         hold on
         plot(imag(correlation));
         
         
      end
      function obj = performDpdLearning3(obj,pa,signal)
         %Set some things up
         dpdBlockIndx = 0;
         IM3Filter    = obj.feedbackFilter; %filter function prefers this
         W = 0;
         
         obj.normalizedFreqShift = obj.dpdFreq/signal.CCs.CC1.systemFs;
         
         if obj.nLearningSamples > length(obj.lmsSignal)
            obj.nLearningSamples = length(obj.lmsSignal)-obj.filteringBlockLength;
         end
         
         for sample = 1:2*obj.filteringBlockLength:obj.nLearningSamples
            %% Part 1 of Newton
            PA_In_StartIndx = sample;
            PA_In_EndIndx   = sample + ...
               obj.filteringBlockLength + obj.loopDelay - 1;
            dpdBlockIndx    = dpdBlockIndx + 1;
            
            % Apply LMS block filtering
            LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
            %            LMS_FilteringBlock = fi(LMS_FilteringBlock,1,12,11);
            dpdSignal = (obj.alpha*LMS_FilteringBlock);
            
            % Add the LMS filter output to the PA input
            signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
            %             signalInput = fi(signalInput,1,12,11);
            PA_InBlock = signalInput ...
               + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
            
            %Brodcast block
            PA_OutBlock = broadcast(pa,PA_InBlock);
            power_IM3(dpdBlockIndx,:) = IM3Power(obj,PA_OutBlock);
            
            % Shift the PA output such that the IM3 frequency is at baseband
            PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
            
            % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
            IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
            %             IM3FilteredBlock = fi(IM3FilteredBlock,1,12,11);
            
            % Extract the IM3 block from the PA output in the feedback receiver
            ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
            
            % Perform Correlation Calculation
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = sample;
            LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
            LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
            correlation(dpdBlockIndx,:) = MeanCorrelation;
            
            % Block LMS filter update
            W = W - (obj.learningRateMu).*MeanCorrelation;
            change1 = - (obj.learningRateMu).*MeanCorrelation';
            
            % Store DPD filter coeff.
            DPD_Coeff(dpdBlockIndx,:) = W';
            obj.alpha = double(W');
            sample = sample + obj.filteringBlockLength;
            
            %% Part 2 of Newton
            PA_In_StartIndx = sample;
            PA_In_EndIndx   = sample + ...
               obj.filteringBlockLength + obj.loopDelay - 1;
            dpdBlockIndx    = dpdBlockIndx + 1;
            
            % Apply LMS block filtering
            LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
            %            LMS_FilteringBlock = fi(LMS_FilteringBlock,1,12,11);
            dpdSignal = (obj.alpha*LMS_FilteringBlock);
            
            % Add the LMS filter output to the PA input
            signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
            %             signalInput = fi(signalInput,1,12,11);
            PA_InBlock = signalInput ...
               + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
            
            %Brodcast block
            PA_OutBlock = broadcast(pa,PA_InBlock);
            power_IM3(dpdBlockIndx,:) = IM3Power(obj,PA_OutBlock);
            
            % Shift the PA output such that the IM3 frequency is at baseband
            PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
            
            % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
            IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
            %             IM3FilteredBlock = fi(IM3FilteredBlock,1,12,11);
            
            % Extract the IM3 block from the PA output in the feedback receiver
            ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
            
            % Perform Correlation Calculation
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = sample;
            LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
            LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
            %             MeanCorrelation= fi(MeanCorrelation,1,12,11);
            correlation(dpdBlockIndx,:) = MeanCorrelation;
            % Block LMS filter update
            
            change2 = - (obj.learningRateMu).*MeanCorrelation';
            
            %Newtons Method
            m = change2-change1;
            real_x_intercept = -real(change1)/real(m);
            integral_real = 0.5*real_x_intercept*real(change1);
            
            imag_x_intercept = -imag(change1)/imag(m);
            integral_imag = 0.5*imag_x_intercept*imag(change1);
            
            W = W + (integral_real - integral_imag*i);
            DPD_Coeff(dpdBlockIndx,:) = W';
            obj.alpha = W';
         end
         
         figure();
         yyaxis left
         plot(real(DPD_Coeff));
         hold on
         plot(imag(DPD_Coeff));
         
         %figure();
         %plot(real(correlation));
         %hold on
         %plot(imag(correlation));
         
         %figure();
         yyaxis right
         plot(power_IM3);
         
         
      end
      function obj = performFixedPointDpdLearning(obj,pa,signal,myDAC,mathWordLength,mathFracLength)
         %Set some things up
         dpdBlockIndx = 0;
         IM3Filter    = obj.feedbackFilter; %filter function prefers this
         W = 0;
         obj.normalizedFreqShift = obj.dpdFreq/signal.CCs.CC1.systemFs;
         
         if obj.nLearningSamples > length(obj.lmsSignal)
            obj.nLearningSamples = length(obj.lmsSignal)-obj.filteringBlockLength;
         end
         
         for sample = 1:obj.filteringBlockLength:obj.nLearningSamples
            PA_In_StartIndx = sample;
            PA_In_EndIndx   = sample + ...
               obj.filteringBlockLength + obj.loopDelay - 1;
            dpdBlockIndx    = dpdBlockIndx + 1;
            
            % Apply LMS block filtering
            LMS_FilteringBlock = obj.lmsSignal(PA_In_StartIndx:PA_In_EndIndx);
            %            LMS_FilteringBlock = fi(LMS_FilteringBlock,1,12,11);
            dpdSignal = (obj.alpha*LMS_FilteringBlock);
            
            % Add the LMS filter output to the PA input
            signalInput = signal.sampleArray(PA_In_StartIndx:PA_In_EndIndx);
            %             signalInput = fi(signalInput,1,12,11);
            PA_InBlock = signalInput ...
               + dpdSignal.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'* obj.normalizedFreqShift);
            
            %Brodcast block
            PA_OutBlock = broadcast(pa,PA_InBlock);
            PA_OutBlock = quantize(myDAC,PA_OutBlock);
            
            % Shift the PA output such that the IM3 frequency is at baseband
            PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*obj.normalizedFreqShift);
            
            % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
            IM3FilteredBlock = filter(IM3Filter,1,double(PA_OutBlockShifted));
            %             IM3FilteredBlock = fi(IM3FilteredBlock,1,12,11);
            
            % Extract the IM3 block from the PA output in the feedback receiver
            ErrorBlock = IM3FilteredBlock(obj.loopDelay+1:obj.loopDelay+obj.learningBlockLength);
            
            %% Perform Correlation Calculation
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = sample;
            LMS_In_EndIndx = LMS_In_StartIndx + obj.learningBlockLength - 1;
            LMS_InputBlock = obj.lmsSignal(LMS_In_StartIndx:LMS_In_EndIndx).';
            %             LMS_InputBlock = fi(LMS_InputBlock,1,12,11);
            %             LMS_InputBlock.fimath.ProductMode = 'SpecifyPrecision';
            %             LMS_InputBlock.fimath.ProductWordLength     = 12;
            %             LMS_InputBlock.fimath.ProductFractionLength = 11;
            %             LMS_InputBlock.fimath.SumMode = 'SpecifyPrecision';
            %             LMS_InputBlock.fimath.SumWordLength     = 12;
            %             LMS_InputBlock.fimath.SumFractionLength = 11;
            %
            %             ErrorBlock.fimath.ProductMode = 'SpecifyPrecision';
            %             ErrorBlock.fimath.ProductWordLength     = 12;
            %             ErrorBlock.fimath.ProductFractionLength = 11;
            %             ErrorBlock.fimath.SumMode = 'SpecifyPrecision';
            %             ErrorBlock.fimath.SumWordLength     = 12;
            %             ErrorBlock.fimath.SumFractionLength = 11;
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);
            %             MeanCorrelation= fi(MeanCorrelation,1,12,11);
            
            
            
            % Block LMS filter update
            W = W - (obj.learningRateMu).*MeanCorrelation;
            
            % Store DPD filter coeff.
            %DPD_Coeff(dpdBlockIndx,:) = W';
            obj.alpha = double(W');
         end
      end
      function out = applyDPDtoSignal(obj,in)
         dpdSignal = conv(obj.alpha,obj.lmsSignal); %Baseband centered DPD
         dpdSignal = dpdSignal .*exp(2*pi*1i*(1:length(dpdSignal)).'* obj.normalizedFreqShift);
         out = in.sampleArray + dpdSignal;
      end
      function out = IM3Power(obj,Rx_Signal,fs)
         PA_Power_Measured = 23;
         Power_Rx_Signal = 10*log10(mean(abs(Rx_Signal).^2));
         PwrScaleFactor = 10^((PA_Power_Measured - Power_Rx_Signal)*0.1);
         Rx_Signal = Rx_Signal*sqrt(PwrScaleFactor);
         [pxx,f] = pwelch(Rx_Signal,500,300,500,fs,'centered','power');
         log_version = 10*log10(pxx);
         f_9MHz = 359;
         out = log_version(f_9MHz);
      end
      
      
   end
end
