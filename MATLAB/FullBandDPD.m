classdef FullBandDPD
   %FullBandDPD Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      poly_type
      NumIterations % Number of ILA iterations
      P_M_1         % Model order for the PH main branch
      P_C_1         % Model order for the PH conjugate branch
      BlockSize     % Block size used for DPD paramter estimation
      PH_FilterLenghts_M_1  % Model filter lengths for the PH main branch
      PH_FilterLenghts_C_1  % Model filter lengths for the PH conjugate branch
      dpd_est1      % DPD coeffs
      
   end
   
   methods
      function obj = FullBandDPD(myPA, mySignal)
         % ILA DPD Paramters
         obj.poly_type = 'conv';
         obj.NumIterations = 50;
         obj.P_M_1 = 5;
         obj.P_C_1 = 3;
         obj.BlockSize = 30970;
         
         MemoryOrder_1 = 5; % Memory depth per branch
         obj.PH_FilterLenghts_M_1 = [MemoryOrder_1 MemoryOrder_1 MemoryOrder_1];
         obj.PH_FilterLenghts_C_1 = [MemoryOrder_1 MemoryOrder_1];
         
         obj = dpdLearning(obj,mySignal,myPA);     
      end
      
      function obj = dpdLearning(obj,mySignal,myPA)         
         DPDout1 = mySignal.sampleArray;
         
         %Broadcast through PA
         DPD_Estimator_In1 = broadcast(myPA,DPDout1);
         
         for Iteration = 1:obj.NumIterations
            % DPD estimation:
            obj.dpd_est1 = DPD_APH_Estimation(obj,DPD_Estimator_In1,DPDout1);
            
            % Predistort a new block of data:  FOR NOW,JUST DO WHOLE THING
            DPDin = mySignal.sampleArray;
            DPDout1 = DPD_APH(obj,DPDin);
            
            %Broadcast through PA
            PA_OutputSignal_DPD = broadcast(myPA,DPDout1);
            
            % Dividing by the PA gain and the additional attenuation factor
            DPD_Estimator_In1 = PA_OutputSignal_DPD;
         end
      end
      
      function [h_est] = DPD_APH_Estimation(obj,IN,OUT)
         % Usage: [h_est]=estimateAPH(IN,OUT,P,Mph,N)
         % Estimate parameters of an Augmented Parallel Hammerstein (APH) nonlinear
         % system with polynomial nonlinearities and FIR filters.
         % INPUTS:
         %    IN   -  input signal to the APH nonlinearity
         %    OUT  -  output signal of the APH nonlinearity ("reference signal")
         %    P    -  polynomial order of the non-conjugate branch (it is assumed
         %            that only odd orders are used)
         %    Q    -  polynomial order of the conjugate branch
         %    Mp   -  filter lengths of the non-conjugate branch FIR filters (vector with (P+1)/2
         %             elements)
         %    Mq   -  filter lengths of the conjugate branch FIR filters (vector with (Q+1)/2
         %             elements)
         %    lo_est - flag which indicates whether LO leakage compensator coefficient
         %             is estimated (1 == estimated; 0 == not estimated)
         %    N    -  estimation block length (default value=min(length(IN),length(REF)))
         %    poly_type - polynomial type: either 'conv' for the conventional monomial
         %                form, or 'orth', which uses a class of orthogonal polynomials
         %                designed for uniformly distributed data
         % OUTPUTS:
         %    h_est -  estimated APH parameters
         %
         % Lauri Anttila / TUT, April 2013
         
         IN  = IN(:);
         OUT = OUT(:);
         Mp = obj.PH_FilterLenghts_M_1(:);
         Mq = obj.PH_FilterLenghts_C_1(:);
         Mmax = max([Mp;Mq]);  % length of the longest PD filter
         P = obj.P_M_1;
         Q = obj.P_C_1;
         N = obj.BlockSize;
         lo_est = 0;
         
         %% SYSTEM IDENTIFICATION
         PHI = BasisMatrix(obj,IN(1:N),P,Mp,obj.poly_type); % build polynomial basis matrices
         % PHI = BasisMatrixWeighted(IN(1:N),P,Mp,poly_type,b); % build weighted polynomial basis matrix
         
         R_win = PHI(Mmax:N,:);  % Covariance windowing method
         if sum(Mq)
            CPHI = BasisMatrix(obj,conj(IN(1:N)),Q,Mq,obj.poly_type);
            % CPHI = BasisMatrixWeighted(conj(IN(1:N)),Q,Mq,poly_type,b); % build polynomial basis matrix
            CPHI_win = CPHI(Mmax:N,:);  % Covariance windowing method
            R_win = [R_win CPHI_win];
         end
         R_win = [R_win ones(length(R_win),lo_est)]; % append with a column of all ones
         % to account for LO leakage
         z = OUT(Mmax:N);  % Covariance windowing method
         h_est = R_win\z;  % Least Squares estimation
         
      end
      
      function out = DPD_APH(obj,IN)
         
         % Usage: [OUT]=predistorterAPH(IN,param,P,Q,Mp,Mq,poly_type)
         % Parallel Hammerstein (PH) predistorter with polynomial nonlinearities
         % and FIR filters.
         % INPUTS:
         %    IN   -  input signal to the PH predistorter
         %   param -  parameters of the PH predistorter in vector form (the first Mph(1)
         %             elements contain the linear filter parameters, the next Mph(2)
         %             elements contain the 3rd order polynomial filter parameters, etc.)
         %    P    -  polynomial order for the non-conjugate branch (it is assumed that
         %            only odd orders are used)
         %    Q    -  polynomial order for the conjugate branch
         %    Mp   -  filter lengths of the branch FIR filters (vector with (P+1)/2
         %             elements)
         %    Mq   -  filter lengths of the branch FIR filters (vector with (Q+1)/2
         %             elements)
         %    lo_comp - flag which indicates whether LO leakage is compensated
         %              (1 == compensated; 0 == not compensated)
         %    poly_type - polynomial type: either 'conv' for the conventional monomial
         %                form, or 'orth', which uses a class of orthogonal polynomials
         %                designed for uniformly distributed data
         %
         % OUTPUTS:
         %    OUT  -  predistorter output signal
         %
         % Lauri Anttila / TUT, April 2013
         IN  = IN(:);
         Mp = obj.PH_FilterLenghts_M_1(:);
         Mq = obj.PH_FilterLenghts_C_1(:);
         Mpmax = max(Mp); Mqmax = max(Mq);
         Mmax = max([Mpmax Mqmax]);  % length of the longest PD filter
         
         P = obj.P_M_1;
         Q = obj.P_C_1;
         lo_comp = 0;
         param = obj.dpd_est1;
         
         %% PREDISTORTION:
         PHI = BasisMatrix(obj,IN,P,Mp,obj.poly_type); % build the polynomial basis matrix
         PHI = [PHI;zeros(Mmax-Mpmax,sum(Mp))];
         CPHI = BasisMatrix(obj,conj(IN),Q,Mq,obj.poly_type); % build the polynomial basis matrix
         CPHI = [CPHI;zeros(Mmax-Mqmax,sum(Mq))];
         LO = ones(length(PHI),lo_comp);
         out = [PHI CPHI LO]*param;
      end
      
      function PHI = BasisMatrix(obj,IN,P,Mph,poly_type)
         
         % Usage: PHI = basisMatrix(IN,P,Mph,poly_type)
         % Generates the polynomial basis matrix to be used in estimation of the
         % parameters of a Parallel Hammerstein/Memory Polynomial, or a memoryless
         % complex-valued polynomial (when Mph = all ones), nonlinearity.
         %
         % INPUTS:
         %        IN  -  input signal
         %         P  -  Polynomial order
         %       Mph  -  vector of filter lengths
         % poly_type  -  'conv' for conventional (default), 'orth' for orthogonal
         %
         % Lauri Anttila / TUT, April 2013
         
         Mmax = max(Mph);  % length of the longest filter
         param_tot = sum(Mph); % total number of parameters
         N = length(IN);
         
         % Create the "elementary" basis functions:
         X = zeros(N,(P+1)/2);
         
         switch poly_type
            case 'orth'
               % [X,~,~] = OrthoPolyUniform(IN,P,'odd'); % Orthogonal polynomials for
               % uniformly distributed signals with abs(IN)<=1
               % X = X(1:N,:);
            otherwise
               % Conventional polynomials (monomials):
               for k=1:2:P,
                  X(:,(k+1)/2) = IN.*abs(IN).^(k-1);
               end
         end
         
         % Build the polynomial basis matrix PHI from X:
         PHI = zeros(N+Mmax-1,param_tot);
         for Ip = 1:length(Mph)
            for ll = 1:Mph(Ip),
               PHI(ll:N+ll-1,sum(Mph(1:Ip-1))+ll) = X(:,Ip);
            end
         end
      end
   end
end

