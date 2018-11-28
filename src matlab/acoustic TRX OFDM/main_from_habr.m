%DVB-T 2K Transmission
%��������� ������ 8 MHz
%2K ��� ��������� ��������
clear all;
close all;
%DVB-T ���������
Tu=224e-6; %�������� ������ OFDM �������
T=Tu/2048; %������������ ������
G=1/4; %���������� 1/4, 1/8, 1/16, � 1/32
delta=G*Tu; %�������� ��������
Ts=delta+Tu; %������ ������ OFDM �������
Kmax=1705; % ������������ ���������� ����������
Kmin=0;
FS=4096; %IFFT/FFT �����
q=10; %������ ����������
fc=q*1/T; %������� �������
Rs=4*fc; %������ ���������
t=0:1/Rs:Tu;
 
%��������� ������
M=Kmax+1;
rand('state',0);
a=-1+2*round(rand(M,1)).'+i*(-1+2*round(rand(M,1))).';
A=length(a);
info=zeros(FS,1);
plot(info);
info(1:(A/2)) = [ a(1:(A/2)).'];
info((FS-((A/2)-1)):FS) = [ a(((A/2)+1):A).'];
 
%��������� ����������
carriers=FS.*ifft(info,FS);
tt=0:T/2:Tu;
figure(1);
subplot(211);
stem(tt(1:20),real(carriers(1:20)));%�������� ����� ��������� �������������� �����
subplot(212);
stem(tt(1:20),imag(carriers(1:20)));%������ ����� ��������� �������������� �����
figure(2);
f=(2/T)*(1:(FS))/(FS);
subplot(211);
plot(f,abs(fft(carriers,FS))/FS);
subplot(212);
pwelch(carriers,[],[],[],2/T);
 
% ��������� ���
L = length(carriers);
chips = [ carriers.';zeros((2*q)-1,L)]; %����
p=1/Rs:1/Rs:T/2;
g=ones(length(p),1);
figure(3);
stem(p,g); 
 
dummy=conv(g,chips(:)); %������
u=[dummy(1:length(t))]; % 
figure(4);
subplot(211);
plot(t(1:400),real(u(1:400)));
subplot(212);
plot(t(1:400),imag(u(1:400)));
 
figure(5);
ff=(Rs)*(1:(q*FS))/(q*FS);
subplot(211);
plot(ff,abs(fft(u,q*FS))/FS);
subplot(212);
pwelch(u,[],[],[],Rs);
 
[b,a] = butter(13,1/20); %������ ������
[H,F] = FREQZ(b,a,FS,Rs);
figure(6);
plot(F,20*log10(abs(H)));
uoft = filter(b,a,u); %��������� ������ 
 
figure(7);
subplot(211);
plot(t(80:480),real(uoft(80:480)));
subplot(212);
plot(t(80:480),imag(uoft(80:480)));
 
figure(8);
subplot(211);
plot(ff,abs(fft(uoft,q*FS))/FS);
subplot(212);
pwelch(uoft,[],[],[],Rs);
 
%Upconverter
s_tilde=(uoft.').*exp(1i*2*pi*fc*t);
s=real(s_tilde);
figure(9);
plot(t(80:480),s(80:480));
figure(10);
subplot(211);
%plot(ff,abs(fft(((real(uoft).').*cos(2*pi*fc*t)),q*FS))/FS);
%plot(ff,abs(fft(((imag(uoft).').*sin(2*pi*fc*t)),q*FS))/FS);
plot(ff,abs(fft(s,q*FS))/FS);
subplot(212);
%pwelch(((real(uoft).').*cos(2*pi*fc*t)),[],[],[],Rs);
%pwelch(((imag(uoft).').*sin(2*pi*fc*t)),[],[],[],Rs);
pwelch(s,[],[],[],Rs);