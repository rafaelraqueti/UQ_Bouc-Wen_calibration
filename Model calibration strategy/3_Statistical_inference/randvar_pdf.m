function [bins,freq,area] = randvar_pdf(data,numbins)

   Ns = length(data);    
        
   data_max = max(data);
   data_min = min(data);
   binwidth = (data_max-data_min)/(numbins-1);
   
   bins     = (data_min:binwidth:data_max);
   freq     = histc(data,bins);
   freq     = freq/(Ns*binwidth);
   area     = binwidth*sum(freq);

end