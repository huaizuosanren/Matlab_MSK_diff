function out = delay(data,n,sample_number) 
 %data:延迟的数据 
%n:延迟码元个数 
%sample_number:码元采样个数 
out = zeros(1,length(data)); 
 out(n*sample_number+1:length(data)) = data(1:length(data)-n*sample_number);