function out = delay(data,n,sample_number) 
 %data:�ӳٵ����� 
%n:�ӳ���Ԫ���� 
%sample_number:��Ԫ�������� 
out = zeros(1,length(data)); 
 out(n*sample_number+1:length(data)) = data(1:length(data)-n*sample_number);