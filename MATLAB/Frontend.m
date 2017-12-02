classdef Frontend
   %UNTITLED2 Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      PA
      DAC
   end
   
   methods
      function obj = Frontend(PA,DAC)
         obj.PA = PA;
         obj.DAC =DAC;
      end
      function out = broadcast(obj,in)
         out = zeros(length(in),1);
         out = broadcast(obj.PA, in);
         %maxreal = max(real(out));
         %maximag = max(imag(out));
         %MAX     = max(maxreal,maximag);
         %out = out / MAX; %Normalize to best use resolution
         
         %out = quantize(obj.DAC,out);
      end
      function out = quantize(obj,in)
         maxreal = max(abs(real(in)));
         maximag = max(abs(imag(in)));
         MAX     = max(maxreal,maximag);
         out = in / MAX; %Normalize to best use resolution          
         out = quantize(obj.DAC,out);
      end
   end
   
end

