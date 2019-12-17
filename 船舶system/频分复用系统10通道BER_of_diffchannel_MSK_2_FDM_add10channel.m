 %MSK调制，差分解调方法二 
clear all 
 close all 
  
 %-------------------------------------------------------------------------- 
 %参数设置 
data_len = 1000;           %码元个数 
sample_number = 8;          %采样个数 
Rb = 200;                 %码元速率 
%fc = 20e3-200;                 %载波频率 
%fc1 = 20e3+200;
fc = 20e3;  %96000;                 %载波频率 
df=400;
channel_num=10;
f=zeros(1,channel_num);
f(1)=fc-(df/2);
%f(2)=f(1)-df;
for i= 2:(channel_num/2)
    f(i)=f(i-1)-df;
end
f((channel_num/2)+1)=fc+(df/2);
for i= ((channel_num/2)+2):(channel_num)
    f(i)=f(i-1)+df;
end

%f(3)=fc+(df/2);
%f(4)=f(3)+df;
f
%f1=f(1);
%************************************************************************** 
  
 %-------------------------------------------------------------------------- 
%随机产生传输信号 
data=rand_binary(data_len); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %MSK基带调制 
[signal_out,I_out,Q_out] = mod_msk(data,data_len,sample_number,Rb); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %中频搬移 
multi = fc/Rb; 
 I_temp=interp(I_out,multi); 
 Q_temp=interp(Q_out,multi); 
  
 Fs=fc*sample_number; 
 t=1/Fs:1/Fs:length(I_temp)*1/Fs; 
% signal_i=I_temp.*cos(2*pi*fc*t); 
% signal_q=Q_temp.*sin(2*pi*fc*t); 
signal_mod=zeros(channel_num,length(I_temp));
 for i =1 :channel_num
     signal_mod(i,:)=I_temp.*cos(2*pi*f(i)*t)-Q_temp.*sin(2*pi*f(i)*t);
 end


 %************************************************************************** 
%--------------------add all the channels----------------
s=zeros(1,length(I_temp));
s=signal_mod(1,:);
for i= 2:channel_num
    s=s+signal_mod(i,:);
end
%s=signal_mod(1,:)+signal_mod(2,:)+signal_mod(3,:)+signal_mod(4,:);
%s=signal_mod(1,:)+signal_mod(2,:);
%---------------设计Butterworth低通滤波器-----------------
%---------------parameters of butter filter------------
fs=Fs;
n=4;
Wn=[(f(1)-200)/(fs/2) (f(1)+200)/(fs/2)]
[a,b]=butter(n,Wn);
a
b
[h,f_filt]=freqz(a,b,'whole',fs);        %求数字低通滤波器的频率响应
f_filt=(0:length(f_filt)-1*fs/length(f_filt));     %进行对应的频率转换
figure;
plot(f_filt(1:length(f_filt)/2),abs(h(1:length(f_filt)/2)));       %绘制幅频响应图
title('巴特沃斯带通滤波器');xlabel('频率/Hz');ylabel('幅度');
grid;
sF=filter(a,b,s);                   %叠加函数s经过低通滤波器以后的新函数
figure;
subplot(121);
plot(t,sF);                         %绘制叠加函数s经过低通后时域图形
title('输出信号');xlabel('t/s');ylabel('幅度');
SF=fft(sF);
subplot(122);
plot((1:length(SF)/2)*fs/length(SF),2*abs(SF(1:length(SF)/2))/length(SF));
title('带通滤波后频谱');xlabel('频率/Hz');ylabel('幅度');
%-------------------------------------------------------------------------- 
ber=zeros(channel_num,9)
for channel= 1:channel_num
    %---------------设计Butterworth低通滤波器-----------------
%---------------parameters of butter filter------------
fs=Fs;
n=4;
Wn=[(f(channel)-200)/(fs/2) (f(channel)+200)/(fs/2)]
[a,b]=butter(n,Wn);
[h,f_filt]=freqz(a,b,'whole',fs);        %求数字低通滤波器的频率响应
f_filt=(0:length(f_filt)-1*fs/length(f_filt));     %进行对应的频率转换
%figure;
%plot(f_filt(1:length(f_filt)/2),abs(h(1:length(f_filt)/2)));       %绘制幅频响应图
%title('巴特沃斯带通滤波器');xlabel('频率/Hz');ylabel('幅度');
%grid;
sF=filter(a,b,s);                   %叠加函数s经过低通滤波器以后的新函数
%figure;
%subplot(121);
%plot(t,sF);                         %绘制叠加函数s经过低通后时域图形
%title('输出信号');xlabel('t/s');ylabel('幅度');
%SF=fft(sF);
%subplot(122);
%plot((1:length(SF)/2)*fs/length(SF),2*abs(SF(1:length(SF)/2))/length(SF));
%title('带通滤波后频谱');xlabel('频率/Hz');ylabel('幅度');
 %加噪声 
for SNR = 0:8 
 signal_mod1 = awgn(sF,SNR); 
 cycle_num = 1000;
 num=zeros(1,cycle_num);
 for cycle= 1:cycle_num
     %-------------------------------------------------------------------------- 
     %去载波 
    N=300;                                              % 滤波器的阶数为(N+1)   
     F=[0,f(channel)-400,f(channel)+400,Fs/2]*2/Fs;  
     A=[1,1,0,0]; 
     lpf=firls(N,F,A); 
     [amp_lpf,w]=freqz(lpf); 
      
     I_dem=signal_mod1.*cos(2*pi*f(channel)*t)*2; 
     I_dem=conv(I_dem,lpf); 
     I_dem=I_dem(N/2+1:N/2+length(I_temp)); 
     Q_dem=signal_mod1.*sin(2*pi*f(channel)*t)*2; 
     Q_dem=conv(Q_dem,lpf); 
     Q_dem=-Q_dem(N/2+1:N/2+length(I_temp)); 
      
     I_dem_out=zeros(1,length(I_dem)/multi);         % 抽取 
    Q_dem_out=zeros(1,length(Q_dem)/multi); 
     for i=1:length(I_dem_out) 
        I_dem_out(i)=I_dem(multi*(i-1)+1); 
        Q_dem_out(i)=Q_dem(multi*(i-1)+1); 
     end; 
     %************************only plot when breakpoints on************************************************** 
   %  figure;
   %   plot(I_dem_out);
   %  plot(I_dem_out);
     %-------------------------------------------------------------------------- 
     %差分解调 
    demod_data = zeros(1,data_len); 
     demod_data(1) = Q_dem_out(sample_number); 
     for i = 2:2:data_len 
         demod_data(i) = -I_dem_out(i*sample_number)*Q_dem_out((i-1)*sample_number); 
     end 
     for i = 3:2:data_len 
         demod_data(i) = Q_dem_out(i*sample_number)*I_dem_out((i-1)*sample_number); 
     end 
     %************************************************************************** 
      
     %-------------------------------------------------------------------------- 
     %判决 
    demod_data = demod_data>0; 
     demod_data = 2*demod_data-1; 
     %************************************************************************** 
     [num(cycle),ber(channel,SNR+1)]=symerr(demod_data,data); 
 end %end for cycle
     %-------------------------------------------------------------------------- 
     %计算误码率 
   % [num,ber(SNR+1)]=symerr(demod_data,data); 
     %************************************************************************** 
    errnum=sum(num);
    ber(channel,SNR+1)=errnum/(data_len*cycle_num)
    disp(int2str(SNR));
    disp('dB have done');
 end 
 %************************************************************************** 
disp(int2str(channel));
disp('channel have done');
end %end for channel
 %-------------------------------------------------------------------------- 
 %误码率曲线 
 figure;
semilogy([0:8],ber(1,:),'r*'); 
hold on;
semilogy([0:8],ber(2,:),'r-'); 
hold on;
semilogy([0:8],ber(3,:),'r-.'); 
hold on;
semilogy([0:8],ber(4,:),'r:'); 
hold on;
semilogy([0:8],ber(5,:),'rs'); 
hold on;
semilogy([0:8],ber(6,:),'rv'); 
 %************************************************************************** 
  
 %-------------------------------------------------------------------------- 
 %误码率理论值 
snr = 0:0.1:8; 
 for i = 1:length(snr) 
     snr1(1,i) = 10^(snr(1,i)/10); 
     ps(1,i) = 1/2 * erfc(sqrt(snr1(1,i))); 
     pe(1,i) = 2 * ps(1,i); 
 end 
 hold on 
 semilogy([0:.1:8],pe); 
 %**************************************************************************
 %save 'ber_for_10channel.mat' ber;