 function [signal_out,I_out,Q_out] = mod_msk(data,data_len,sample_number,Rb) 
 %MSK�������� 
%************************************************************************** 
 % data              �����ź� 
% data_len          ��Ԫ���� 
% sample_number     ÿ����Ԫ�������� 
% Rb                ��Ԫ���� 
% signal_out        ����������� 
% I_out             I·��� 
% Q_out             Q·��� 
%************************************************************************** 
  
 % data_len = 10;                %��Ԫ���� 
% sample_number = 8;            %�������� 
% Rb = 16000;                   %��Ԫ���� 
% data1 = randint(1,data_len); 
 % data = 2*data1-1;             %��������� 
 
Tb = 1/Rb;                      %��Ԫʱ�� 
fs = Rb*sample_number;          %�������� 
 
 %-------------------------------------------------------------------------- 
%��ֱ��� 
[data_diff] = difference(data); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %����ת������ʱ 
I(1) = 1;             %fai0 = 0,cos(fai0) = 1 
 for i = 1:2:data_len 
     Q(i) = data_diff(i); 
     Q(i+1) = data_diff(i); 
 end 
 for i = 2:2:data_len 
     I(i+1) = data_diff(i); 
     I(i) = data_diff(i); 
 end 
  
 for i = 1:sample_number 
     I1(i:sample_number:data_len*sample_number) = I(1:data_len); 
     Q1(i:sample_number:data_len*sample_number) = Q(1:data_len); 
 end 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %�˼�Ȩ���� 
t=1/fs:1/fs:data_len*Tb; 
 I_out = I1 .* cos(pi*t/2/Tb); 
 Q_out = Q1 .* sin(pi*t/2/Tb); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %�����źŲ��� 
signal_out = I_out + j*Q_out; 
 %************************************************************************** 
  
 % %-------------------------------------------------------------------------- 
 % %��ͼ 
% subplot(221) 
 % plot(data,'.-');title('MSK���������');xlabel('ʱ��');ylabel('����') 
 % subplot(222) 
 % plot(data_diff,'.-');title('��ֺ������');xlabel('ʱ��');ylabel('����') 
 % subplot(223) 
 % plot(I1,'.-');title('��ȨǰI·');xlabel('ʱ��');ylabel('����'); 
 % subplot(224) 
 % plot(Q1,'.-');title('��ȨǰQ·');xlabel('ʱ��');ylabel('����'); 
 %  
 % figure(2) 
 % subplot(221) 
 % plot(cos(pi*t/2/Tb),'.-');title('��Ȩ����cos(��t/(2Tb))');xlabel('ʱ��');ylabel('����') 
 % subplot(222) 
 % plot(sin(pi*t/2/Tb),'.-');title('��Ȩ����sin(��t/(2Tb))');xlabel('ʱ��');ylabel('����') 
 % subplot(223) 
 % plot(I_out,'.-');title('��Ȩ��I·');xlabel('ʱ��');ylabel('����'); 
 % subplot(224) 
 % plot(Q_out,'.-');title('��Ȩ��Q·');xlabel('ʱ��');ylabel('����'); 
 % %************************************************************************** 
