 function [signal_out,I_out,Q_out] = mod_msk(data,data_len,sample_number,Rb) 
 %MSK基带调制 
%************************************************************************** 
 % data              调制信号 
% data_len          码元个数 
% sample_number     每个码元采样点数 
% Rb                码元速率 
% signal_out        基带调制输出 
% I_out             I路输出 
% Q_out             Q路输出 
%************************************************************************** 
  
 % data_len = 10;                %码元个数 
% sample_number = 8;            %采样点数 
% Rb = 16000;                   %码元速率 
% data1 = randint(1,data_len); 
 % data = 2*data1-1;             %传输的序列 
 
Tb = 1/Rb;                      %码元时间 
fs = Rb*sample_number;          %采样速率 
 
 %-------------------------------------------------------------------------- 
%差分编码 
[data_diff] = difference(data); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %并串转换，延时 
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
 %乘加权函数 
t=1/fs:1/fs:data_len*Tb; 
 I_out = I1 .* cos(pi*t/2/Tb); 
 Q_out = Q1 .* sin(pi*t/2/Tb); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %调制信号产生 
signal_out = I_out + j*Q_out; 
 %************************************************************************** 
  
 % %-------------------------------------------------------------------------- 
 % %画图 
% subplot(221) 
 % plot(data,'.-');title('MSK传输的数据');xlabel('时间');ylabel('幅度') 
 % subplot(222) 
 % plot(data_diff,'.-');title('差分后的数据');xlabel('时间');ylabel('幅度') 
 % subplot(223) 
 % plot(I1,'.-');title('加权前I路');xlabel('时间');ylabel('幅度'); 
 % subplot(224) 
 % plot(Q1,'.-');title('加权前Q路');xlabel('时间');ylabel('幅度'); 
 %  
 % figure(2) 
 % subplot(221) 
 % plot(cos(pi*t/2/Tb),'.-');title('加权函数cos(πt/(2Tb))');xlabel('时间');ylabel('幅度') 
 % subplot(222) 
 % plot(sin(pi*t/2/Tb),'.-');title('加权函数sin(πt/(2Tb))');xlabel('时间');ylabel('幅度') 
 % subplot(223) 
 % plot(I_out,'.-');title('加权后I路');xlabel('时间');ylabel('幅度'); 
 % subplot(224) 
 % plot(Q_out,'.-');title('加权后Q路');xlabel('时间');ylabel('幅度'); 
 % %************************************************************************** 
