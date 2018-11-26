%2016-01-12
%Sound transfer
close all,clc,clear all;
%******************************
%*******Transmitter************
%******************************
%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)  
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)  
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)
%2016-01-12
%Sound transfer
close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)  
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)  
%2016-01-12
%Sound transfer

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%carrier signal forming(start)
%Fs = 8192;
Fs = 44100;     %sample rate
F = 44100/100;  %frequency of signal, 200<F<Fs/2, [Hz]
t = 1;          %duration of signal, [seconds]
Td = 2*pi/Fs;   %sampling interval
Ft = F*t;       %Number of bits/symbols/waves

x=0:F*Td:Ft*2*pi;
yC = sin(x);     %carrier signal
disp(['Freq of signal = ',num2str(F),' Hz']);
disp(['Number of transmitted bits = ',num2str(Ft),' [bits]']);
disp(['Troughput = ',num2str(F),' [bits/s]']);
%figure,plot(y(1:ceil(Fs/F)+1));
%carrier signal forming(stop)

%modulation(start)
signal = randi([0,1],Ft,1); %information signal
Samples = ceil(Fs/F);       %!!!! number of samples per one wave/symbol
yM = yC;
for n = 1:Ft
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if signal(n) == 0
        yM(iA:iB) = - yM(iA:iB);
    end
end
signal(1:5)
%modulation(stop)