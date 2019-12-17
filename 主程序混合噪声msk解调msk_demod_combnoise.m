%function [ber]=mskdemod_NBPF(SNR,alpha)
%MSK调制，差分解调方法一  
clear all  
close all 
clc
%--------------------------------------------------------------------------  %参数设置  
%GSNR= [-16 -14 -12 -10 -8 -6 -4 -2 0]    %[-8 -7 -6 -5 -4 -3 -2 -1 0]    %[0 1 2 3 4 5 6 7 8]
%EbN0_db = GSNR./2
EbN0_db =[-8 -7 -6 -5 -4 -3 -2 -1 0]   %  [-12 -11 -10 -9]
GSNR = EbN0_db+3
EbN0=10.^(EbN0_db/10);
alpha=1.5  %1.2
[m n]=size(GSNR);
%Hd=tianxian;
data_len = 1000;           %码元个数  
sample_number = 8;          %采样个数  
Rb = 50;                 %码元速率 
 fc = 20e3; %96000;                 %载波频率 
 yita = 0.2  % the ratio of impulsive noise
%**************************************************************************    
%--------------------------------------------------------------------------  %随机产生传输信号  
data=rand_binary(data_len);  
%**************************************************************************    
%--------------------------------------------------------------------------  %MSK基带调制  
[signal_out,I_out,Q_out] = mod_msk2(data,data_len,sample_number,Rb);  
%**************************************************************************  
%--------------------------------------------------------------------------  %中频搬移  
multi = fc/Rb;  
I_temp=interp(I_out,multi); 
 Q_temp=interp(Q_out,multi); 
Fs=fc*sample_number; 
 t=1/Fs:1/Fs:length(I_temp)*1/Fs;  
signal_i=I_temp.*cos(2*pi*fc*t);  
signal_q=Q_temp.*sin(2*pi*fc*t);  
signal_mod=I_temp.*cos(2*pi*fc*t)-Q_temp.*sin(2*pi*fc*t);  %**************************************************************************  
p_signal=0.5

 %*******************Pass Through The tianxian***************************
%  signal_mod=filter(Hd,signal_mod);
BER = zeros(1,n)
%--------------------------------------------------------------------------  %加噪声  
for SNRi = 1:n  
    num = 0;  %initial the number of error
    num1 = 0;
    SNRpbit=10.^(GSNR(SNRi)/10);
      p_noise = p_signal/SNRpbit;        
    a_noise = sqrt(p_noise);
%    alpha=1.5;
 %   alpha=1.2;
%SNRpbit=10.^(dB/10);                % Eb/No conversion from dB to decimal
    R=1;% R=k/n1;                     % code rate 编码率
    No=1./SNRpbit;
    %由信噪比计算得到噪声产生的参数值
        S0=1/(sqrt(7.12*SNRpbit*R));        %3.56
        gamma=((1.78*S0)^alpha)/1.78;
        scale=gamma^(1/alpha);
        %Nsymbols = time/20;
        symbol_sample_stream_num= (400*Fs/fc)*data_len;
        pulsenoise = yita*stblrnd(alpha,0,scale,0,[1,symbol_sample_stream_num]);
 %-------calculate the comb noise*----------------

    Ntotal=0.5/EbN0(SNRi);%总噪声功率
    N1=Ntotal*yita;%脉冲
    N2=Ntotal-N1;%梳状
    delta=sqrt(N1)/sqrt(2);
    noise2=sqrt(N2)*randn(1,symbol_sample_stream_num).*cos(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs)-...
       sqrt(N2)*randn(1,symbol_sample_stream_num).*sin(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs);
%---------------add noise----------------------
        signal_mod1 = signal_mod +pulsenoise+noise2 ;   %awgn(signal_mod,SNR);      %  signal_mod +pulsenoise 
      %******calcaulate the gaosi noise***********
       % noise=a_noise*randn(1,symbol_sample_stream_num);
       % signal_mod1 = signal_mod +noise;
       
    
       
           %****************xianfulvbo****************************
       for mm=1:length (signal_mod1)
           if signal_mod1(mm)>=1.5
               signal_mod1(mm)=1;
           elseif signal_mod1(mm)<=-1.5
                 signal_mod1(mm)=-1;
           end
       end
 
%for ii=1:length(signal_mod1)
%   if signal_mod1(ii)>1.5+5
%    signal_mod1(ii)= 1+5;
%   elseif signal_mod1(ii)<-1.5+5
%    signal_mod1(ii)=-1+5;
%   end
%end
maxerronum = 50;
maxloopnum = 100;
loopnum = 0;
    while((num1<maxerronum)&&(loopnum<maxloopnum))
%--------------------------------------------------------------------------      %去载波  
    N=300;                                              % 滤波器的阶数为(N+1)        
    F=[0,fc-100,fc+100,Fs/2]*2/Fs;    
    A=[1,1,0,0];      
    lpf=firls(N,F,A);    
    [amp_lpf,w]=freqz(lpf);        
    I_dem=signal_mod1.*cos(2*pi*fc*t)*2; 
    I_dem=conv(I_dem,lpf);  
    I_dem=I_dem(N/2+1:N/2+length(I_temp));      
    Q_dem=signal_mod1.*sin(2*pi*fc*t)*2;      
    Q_dem=conv(Q_dem,lpf);  
    Q_dem=-Q_dem(N/2+1:N/2+length(I_temp));        
    I_dem_out=zeros(1,length(I_dem)/multi);         % 抽取      
    Q_dem_out=zeros(1,length(Q_dem)/multi);      
    for i=1:length(I_dem_out)  
       I_dem_out(i)=I_dem(multi*(i-1)+1);         
       Q_dem_out(i)=Q_dem(multi*(i-1)+1);      
    end;  
    %************************************************************************** 
%--------------------------------------------------------------------------      %差分解调  
    demod_data = zeros(1,data_len);  
    demod_data(1) = Q_dem_out(sample_number);      
    for i = 2:data_len          
    demod_data(i) = Q_dem_out(i*sample_number)*I_dem_out((i-1)*sample_number) - I_dem_out(i*sample_number)*Q_dem_out((i-1)*sample_number);      
    end  
    %************************************************************************** 
%--------------------------------------------------------------------------      %判决  
    demod_data = demod_data>0;      
    demod_data = 2*demod_data-1;  
    %**************************************************************************    
    %--------------------------------------------------------------------------      %计算误码率  
    [num,ber]=symerr(demod_data,data);  
    num1 = num1+num;
    loopnum = loopnum + 1;
    %**************************************************************************  
    end  
    BER(SNRi) = num1/(data_len*loopnum)
end

%************************************************************************** 
%误码率曲线  
figure (101)
semilogy(EbN0_db,BER,'r*-');  
%**************************************************************************    
%--------------------------------------------------------------------------  %误码率理论值  
snr = 0:0.1:8;  
for i = 1:n  
    snr1(1,i) = 10^(EbN0_db(1,i)/10);     
 ps(1,i) = 1/2 * erfc(sqrt(snr1(1,i)));     
 pe(1,i) = 2 * ps(1,i); 
 end 
 hold on  
semilogy(EbN0_db,pe);
legend('fact','calculate');
xlabel('Eb/N0(dB)');
ylabel('BER');
title ('MSK');
hold on;
save BER_02_12_9.txt -ascii BER 

%% ****************************change the yita=0.5********
 yita = 0.5  % the ratio of impulsive noise
%**************************************************************************    
%--------------------------------------------------------------------------  %随机产生传输信号  
data=rand_binary(data_len);  
%**************************************************************************    
%--------------------------------------------------------------------------  %MSK基带调制  
[signal_out,I_out,Q_out] = mod_msk2(data,data_len,sample_number,Rb);  
%**************************************************************************  
%--------------------------------------------------------------------------  %中频搬移  
multi = fc/Rb;  
I_temp=interp(I_out,multi); 
 Q_temp=interp(Q_out,multi); 
Fs=fc*sample_number; 
 t=1/Fs:1/Fs:length(I_temp)*1/Fs;  
signal_i=I_temp.*cos(2*pi*fc*t);  
signal_q=Q_temp.*sin(2*pi*fc*t);  
signal_mod=I_temp.*cos(2*pi*fc*t)-Q_temp.*sin(2*pi*fc*t);  %**************************************************************************  
p_signal=0.5

 %*******************Pass Through The tianxian***************************
%  signal_mod=filter(Hd,signal_mod);
BER = zeros(1,n)
%--------------------------------------------------------------------------  %加噪声  
for SNRi = 1:n  
    num = 0;  %initial the number of error
    num1 = 0;
    SNRpbit=10.^(GSNR(SNRi)/10);
      p_noise = p_signal/SNRpbit;        
    a_noise = sqrt(p_noise);
%    alpha=1.5;
 %   alpha=1.2;
%SNRpbit=10.^(dB/10);                % Eb/No conversion from dB to decimal
    R=1;% R=k/n1;                     % code rate 编码率
    No=1./SNRpbit;
    %由信噪比计算得到噪声产生的参数值
        S0=1/(sqrt(7.12*SNRpbit*R));        %3.56
        gamma=((1.78*S0)^alpha)/1.78;
        scale=gamma^(1/alpha);
        %Nsymbols = time/20;
        symbol_sample_stream_num= (400*Fs/fc)*data_len;
        pulsenoise = yita*stblrnd(alpha,0,scale,0,[1,symbol_sample_stream_num]);
 %-------calculate the comb noise*----------------

    Ntotal=0.5/EbN0(SNRi);%总噪声功率
    N1=Ntotal*yita;%脉冲
    N2=Ntotal-N1;%梳状
    delta=sqrt(N1)/sqrt(2);
    noise2=sqrt(N2)*randn(1,symbol_sample_stream_num).*cos(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs)-...
       sqrt(N2)*randn(1,symbol_sample_stream_num).*sin(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs);
%---------------add noise----------------------
        signal_mod1 = signal_mod +pulsenoise+noise2 ;   %awgn(signal_mod,SNR);      %  signal_mod +pulsenoise 
      %******calcaulate the gaosi noise***********
       % noise=a_noise*randn(1,symbol_sample_stream_num);
       % signal_mod1 = signal_mod +noise;
       
    
       
           %****************xianfulvbo****************************
       for mm=1:length (signal_mod1)
           if signal_mod1(mm)>=1.5
               signal_mod1(mm)=1;
           elseif signal_mod1(mm)<=-1.5
                 signal_mod1(mm)=-1;
           end
       end
 
%for ii=1:length(signal_mod1)
%   if signal_mod1(ii)>1.5+5
%    signal_mod1(ii)= 1+5;
%   elseif signal_mod1(ii)<-1.5+5
%    signal_mod1(ii)=-1+5;
%   end
%end
maxerronum = 50;
maxloopnum = 100;
loopnum = 0;
    while((num1<maxerronum)&&(loopnum<maxloopnum))
%--------------------------------------------------------------------------      %去载波  
    N=300;                                              % 滤波器的阶数为(N+1)        
    F=[0,fc-100,fc+100,Fs/2]*2/Fs;    
    A=[1,1,0,0];      
    lpf=firls(N,F,A);    
    [amp_lpf,w]=freqz(lpf);        
    I_dem=signal_mod1.*cos(2*pi*fc*t)*2; 
    I_dem=conv(I_dem,lpf);  
    I_dem=I_dem(N/2+1:N/2+length(I_temp));      
    Q_dem=signal_mod1.*sin(2*pi*fc*t)*2;      
    Q_dem=conv(Q_dem,lpf);  
    Q_dem=-Q_dem(N/2+1:N/2+length(I_temp));        
    I_dem_out=zeros(1,length(I_dem)/multi);         % 抽取      
    Q_dem_out=zeros(1,length(Q_dem)/multi);      
    for i=1:length(I_dem_out)  
       I_dem_out(i)=I_dem(multi*(i-1)+1);         
       Q_dem_out(i)=Q_dem(multi*(i-1)+1);      
    end;  
    %************************************************************************** 
%--------------------------------------------------------------------------      %差分解调  
    demod_data = zeros(1,data_len);  
    demod_data(1) = Q_dem_out(sample_number);      
    for i = 2:data_len          
    demod_data(i) = Q_dem_out(i*sample_number)*I_dem_out((i-1)*sample_number) - I_dem_out(i*sample_number)*Q_dem_out((i-1)*sample_number);      
    end  
    %************************************************************************** 
%--------------------------------------------------------------------------      %判决  
    demod_data = demod_data>0;      
    demod_data = 2*demod_data-1;  
    %**************************************************************************    
    %--------------------------------------------------------------------------      %计算误码率  
    [num,ber]=symerr(demod_data,data);  
    num1 = num1+num;
    loopnum = loopnum + 1;
    %**************************************************************************  
    end  
    BER(SNRi) = num1/(data_len*loopnum)
end
semilogy(EbN0_db,BER,'r^-');  

save BER_05_12_9.txt -ascii BER


%% ****************************change the yita=0.8********
 yita = 0.8  % the ratio of impulsive noise
%**************************************************************************    
%--------------------------------------------------------------------------  %随机产生传输信号  
data=rand_binary(data_len);  
%**************************************************************************    
%--------------------------------------------------------------------------  %MSK基带调制  
[signal_out,I_out,Q_out] = mod_msk2(data,data_len,sample_number,Rb);  
%**************************************************************************  
%--------------------------------------------------------------------------  %中频搬移  
multi = fc/Rb;  
I_temp=interp(I_out,multi); 
 Q_temp=interp(Q_out,multi); 
Fs=fc*sample_number; 
 t=1/Fs:1/Fs:length(I_temp)*1/Fs;  
signal_i=I_temp.*cos(2*pi*fc*t);  
signal_q=Q_temp.*sin(2*pi*fc*t);  
signal_mod=I_temp.*cos(2*pi*fc*t)-Q_temp.*sin(2*pi*fc*t);  %**************************************************************************  
p_signal=0.5

 %*******************Pass Through The tianxian***************************
%  signal_mod=filter(Hd,signal_mod);
BER = zeros(1,n)
%--------------------------------------------------------------------------  %加噪声  
for SNRi = 1:n  
    num = 0;  %initial the number of error
    num1 = 0;
    SNRpbit=10.^(GSNR(SNRi)/10);
      p_noise = p_signal/SNRpbit;        
    a_noise = sqrt(p_noise);
%    alpha=1.5;
 %   alpha=1.2;
%SNRpbit=10.^(dB/10);                % Eb/No conversion from dB to decimal
    R=1;% R=k/n1;                     % code rate 编码率
    No=1./SNRpbit;
    %由信噪比计算得到噪声产生的参数值
        S0=1/(sqrt(7.12*SNRpbit*R));        %3.56
        gamma=((1.78*S0)^alpha)/1.78;
        scale=gamma^(1/alpha);
        %Nsymbols = time/20;
        symbol_sample_stream_num= (400*Fs/fc)*data_len;
        pulsenoise = yita*stblrnd(alpha,0,scale,0,[1,symbol_sample_stream_num]);
 %-------calculate the comb noise*----------------

    Ntotal=0.5/EbN0(SNRi);%总噪声功率
    N1=Ntotal*yita;%脉冲
    N2=Ntotal-N1;%梳状
    delta=sqrt(N1)/sqrt(2);
    noise2=sqrt(N2)*randn(1,symbol_sample_stream_num).*cos(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs)-...
       sqrt(N2)*randn(1,symbol_sample_stream_num).*sin(2*pi*(fc)*(1:symbol_sample_stream_num)/Fs);
%---------------add noise----------------------
        signal_mod1 = signal_mod +pulsenoise+noise2 ;   %awgn(signal_mod,SNR);      %  signal_mod +pulsenoise 
      %******calcaulate the gaosi noise***********
       % noise=a_noise*randn(1,symbol_sample_stream_num);
       % signal_mod1 = signal_mod +noise;
       
    
       
           %****************xianfulvbo****************************
       for mm=1:length (signal_mod1)
           if signal_mod1(mm)>=1.5
               signal_mod1(mm)=1;
           elseif signal_mod1(mm)<=-1.5
                 signal_mod1(mm)=-1;
           end
       end
 
%for ii=1:length(signal_mod1)
%   if signal_mod1(ii)>1.5+5
%    signal_mod1(ii)= 1+5;
%   elseif signal_mod1(ii)<-1.5+5
%    signal_mod1(ii)=-1+5;
%   end
%end
maxerronum = 50;
maxloopnum = 100;
loopnum = 0;
    while((num1<maxerronum)&&(loopnum<maxloopnum))
%--------------------------------------------------------------------------      %去载波  
    N=300;                                              % 滤波器的阶数为(N+1)        
    F=[0,fc-100,fc+100,Fs/2]*2/Fs;    
    A=[1,1,0,0];      
    lpf=firls(N,F,A);    
    [amp_lpf,w]=freqz(lpf);        
    I_dem=signal_mod1.*cos(2*pi*fc*t)*2; 
    I_dem=conv(I_dem,lpf);  
    I_dem=I_dem(N/2+1:N/2+length(I_temp));      
    Q_dem=signal_mod1.*sin(2*pi*fc*t)*2;      
    Q_dem=conv(Q_dem,lpf);  
    Q_dem=-Q_dem(N/2+1:N/2+length(I_temp));        
    I_dem_out=zeros(1,length(I_dem)/multi);         % 抽取      
    Q_dem_out=zeros(1,length(Q_dem)/multi);      
    for i=1:length(I_dem_out)  
       I_dem_out(i)=I_dem(multi*(i-1)+1);         
       Q_dem_out(i)=Q_dem(multi*(i-1)+1);      
    end;  
    %************************************************************************** 
%--------------------------------------------------------------------------      %差分解调  
    demod_data = zeros(1,data_len);  
    demod_data(1) = Q_dem_out(sample_number);      
    for i = 2:data_len          
    demod_data(i) = Q_dem_out(i*sample_number)*I_dem_out((i-1)*sample_number) - I_dem_out(i*sample_number)*Q_dem_out((i-1)*sample_number);      
    end  
    %************************************************************************** 
%--------------------------------------------------------------------------      %判决  
    demod_data = demod_data>0;      
    demod_data = 2*demod_data-1;  
    %**************************************************************************    
    %--------------------------------------------------------------------------      %计算误码率  
    [num,ber]=symerr(demod_data,data);  
    num1 = num1+num;
    loopnum = loopnum + 1;
    %**************************************************************************  
    end  
    BER(SNRi) = num1/(data_len*loopnum)
end
semilogy(EbN0_db,BER,'rs-');  

save BER_08_12_9.txt -ascii BER
