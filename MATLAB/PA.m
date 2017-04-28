classdef PA
   %PA Class that implements a power amplifier
   %   Detailed explanation goes here
   
   properties
      order
      paParameters
      memory
      memorydepth
   end
   
   methods
      function obj = PA(order)
         obj.memory = 0;
         obj.order = order;
         if(obj.memory)
            obj.paParameters = [0.5823 - 0.0608i;
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
            obj.memorydepth = 4;
         else
            obj.paParameters = [0.9512 - 0.0946i;
               2*(0.239 + 0.1632i);
               2*(0.082 - 0.0727i);
               -0.016 + 0.0147i;
               -0.0001 - 0.0011i];
         end
         obj.paParameters = obj.paParameters*2;
      end
      function out = broadcast(obj,in)
         if(obj.memory)
            PH_f1  = obj.paParameters(1:obj.memorydepth);
            PH_f3 = obj.paParameters(obj.memorydepth + 1:2*obj.memorydepth);
            PH_f5 = obj.paParameters(2*obj.memorydepth + 1:3*obj.memorydepth);
            PH_f7 = obj.paParameters(3*obj.memorydepth + 1:4*obj.memorydepth);
            PH_f9 = obj.paParameters(4*obj.memorydepth + 1:end);
            
            Epsi_1 = in;
            PH_Branch_1 = filter(PH_f1,1,Epsi_1);
            Epsi_3 = in.*abs(in).^2;
            PH_Branch_3 = filter(PH_f3,1,Epsi_3);
            Epsi_5 = in.*abs(in).^4;
            PH_Branch_5 = filter(PH_f5,1,Epsi_5);
            Epsi_7 = in.*abs(in).^6;
            PH_Branch_7 = filter(PH_f7,1,Epsi_7);
            Epsi_9 = in.*abs(in).^8;
            PH_Branch_9 = filter(PH_f9,1,Epsi_9);
            
            out = PH_Branch_1 + PH_Branch_3 + PH_Branch_5 + PH_Branch_7 + PH_Branch_9;
         else
            out = zeros(length(in),1);
            j = 1;
            for i = 1:2:obj.order
               out = out + obj.paParameters(j)*PA.epsilonBranch(in,i);
               j = j + 1;
            end
         end
      end
      
      
   end
   methods(Static)
      function epsilon = epsilonBranch(in, order)
         epsilon = in.*abs(in).^(order-1);
      end
   end
   
end

