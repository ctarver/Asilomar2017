classdef DAC
   %DAC Implements finite precision sampling for DAC and ADC
   %   Direction I don't think matters so much.
   % bits: integer for total number of bits
   % fractional: integer for number of fractional bits. should be less than
   % the total number of bits
   
   properties
      bits
      fractional
   end
   
   methods
      function obj = DAC(bits,fractional)
         obj.bits = bits;
         obj.fractional = fractional;
      end
      function out = quantize(obj,signal)
         if obj.bits == 1
           real_signal = real(signal);
           imag_signal = imag(signal);
           %real
            real_signal(real_signal>=0) = 0.5;
            real_signal(real_signal<0)  = -0.5;
            %imag
            imag_signal(imag_signal>=0) = 0.5;
            imag_signal(imag_signal<0)  = -0.5;
            x = real_signal+imag_signal*1i;
         else
            x = fi(signal,1,obj.bits,obj.fractional);
         end
         out = double(x);
      end
   end
end

