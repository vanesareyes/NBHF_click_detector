clear all;
directory=pwd;
files=dir('*.wav');

%Define parameters
fs=500e3;
% Bandpass filter detector
lf=15e3; %lower freq of the filter
hf=50e3; %higher freq of the filter
[b,a]=fir1(10,[lf/(fs/2),hf/(fs/2)]);
% Bandpass filter espectros
[x,y] = butter(10,[3000 240000]/(fs/2),'bandpass');
nfft=512;
funtran=xlsread('I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\función de transferencia','Hoja2','e2:e257');

%Preallocation
clicks=[];
datos=[];
espectros=[];
espectrosdB=[];
parametros=[];

for i=1:size(files,1) 
    disp(['File no: ' num2str(i)])
    file=files(i).name(1:8);
    siz=wavread([files(i).name],'size');
    w=1;
    while w<=siz(1)-100
        sonido=wavread([files(i).name],[w w+99]);% 100 samples corresponden a 200 microseg
        sonidof=filter(b,a,sonido);
% Defino señal y ruido
        signal=sonidof(51:100);
        noise=sonidof(1:50);
% Calculo RMSsignal y SNR
        SNR=rms(signal)/rms(noise);
% Condición click
        if SNR>=4  % SNR>=12dB                                                                                                         
           c=w-50;
           d=w+249;
           fragmento=wavread([files(i).name],[c d]);
           env=abs(hilbert(fragmento));
           env_prom=env;
%            env_prom=smooth(env,2);
           thr=0.1;
           if max(env_prom)>thr
              int=find(env_prom>rms(fragmento));
              mi=min(int);
            if mi-10 >= 1 && mi+139 <= size(fragmento,1)
               inicio=mi-10;
               fin=mi+139;  
               click=fragmento(inicio:fin); %800 micros en total
                if max(click)<1 & min(click)>-1
                   clicks(:,size(clicks,2)+1)=click; % matriz con click por columna
                   si=c+inicio; %muestra inicial
                   sf=c+fin; %muestra final
                   ti_click=si/fs;
                   tf_click=sf/fs;
                   datos(end+1,:)=[w SNR si sf ti_click tf_click];
                   
                   % Calcula espectro del click
                   dc=mean(click);
                   click=click-dc; %correct for DC offset
                   click=click.*(2^15);%wav en counts2
                   filtsound=filter(x,y,click);
                   fsigrec=filtsound.*rectwin(length(click));
                   fsigrec=fft(fsigrec(:),nfft);
                   FSrec=abs(fsigrec(1:length(fsigrec)/2)); %me quedo con la mitad que es relevante
                   FSreccal=20*log10(FSrec)+funtran;
                   espectros(:,size(espectros,2)+1)=FSrec;
                   espectrosdB(:,size(espectrosdB,2)+1)=FSreccal;
                   
                   %Parámetros acústicos
                   [fpeak, fcenter , dur10dB, bw3dB, bw10dB] = acoustic_params(click, fs, nfft, FSrec, FSreccal);
                   parametros(size(parametros,1)+1,:)=[fpeak fcenter dur10dB bw3dB bw10dB];
                end
            end
           end
        end
    w=w+50;
    clear click si sf SNR ti_click tf_click
    end
    
end
m=[]              
save 'I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\output\clicks.mat' 
save 'I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\output\datos.mat'
save 'I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\output\espectros.mat'
save 'I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\output\espectrosdB.mat'
save 'I:\Las Grutas 2017\Análisis Acústica tursiops 2014\probando\output\parametros.mat'

% Vanesa Reyes 01/2019

