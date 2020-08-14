function [fpeak, fcenter , dur10dB, bw3dB, bw10dB] = acoustic_params(click, fs, nfft, FSrec, FSreccal)

% click: Time domain click signal  
% fs: Sample rate
% NFFT: number of DFT points
% magnitude: click spectrum
% magnitudedB: click spectrum in dB
%
% Calculates   pf: peak frequency 
%              cf: centroid frequency
%              dur: 10dB duration  
%              3dbbw: 3dB bandwidth  
%              10dbbw: 10dB bandwidth

persistent freq_kHz 
freq_kHz = fs/2*linspace(0,1,nfft/2)*1e-3;
%
% Peak Frequency
[y_max index] = max(FSreccal(2:end)); % 2:end Exclude DC
pos_max = index + 1;
fpeak = freq_kHz(pos_max);
% Centroid Frequency
espcalli=10.^(FSreccal/20);%paso a escala lineal
fcenter=sum((espcalli(:).^2).*(freq_kHz'))/sum(espcalli(:).^2);
% 3dB Bandwidth
magnitudedBflip = flipud(FSreccal); 
bw3dB_first = index - (find(magnitudedBflip(end-index+1:end)<=(y_max-3), 1, 'first')) + 1;
bw3dB_last = index + (find(FSreccal(index+1:end)<=(y_max-3), 1, 'first'));
bw3dB = freq_kHz(bw3dB_last)-freq_kHz(bw3dB_first);

%-10dB bandwidth
 low=y_max-10;
    %walk along spectrogram until low is reached on either side
    slopeup=fliplr(FSreccal(1:pos_max));
    slopedown=FSreccal(pos_max:end);

    for e10dB=1:length(slopeup)
       if slopeup(e10dB)<low %stop at value < -3dB: point of lowest frequency
           break
       end
    end

    for o10dB=1:length(slopedown)
       if slopedown(o10dB)<low %stop at value < -3dB: point of highest frequency
           break
       end
    end
    high10db=freq_kHz(pos_max+o10dB); %-10dB highest frequency in kHz
    low10db=freq_kHz(pos_max-e10dB); %-10dB lowest frequency in kHz
    bw10db=high10db-low10db;
    
%-3dB bandwidth   
%calculation of -3dB bandwidth - amplitude associated with the halfpower points of a pressure pulse (see Au 1993, p.118); 
    low=y_max-3; 
    %walk along spectrogram until low is reached on either side
    slopeup=fliplr(FSreccal(1:pos_max));
    slopedown=FSreccal(pos_max:saz(1));

    for e3dB=1:length(slopeup)
       if slopeup(e3dB)<low %stop at value < -3dB: point of lowest frequency
           break
       end
    end

    for o3dB=1:length(slopedown)
       if slopedown(o3dB)<low %stop at value < -3dB: point of highest frequency
           break
       end
    end
   %calculation from spectrogram in 256 steps (FFT=512)
    high3db=fax(pos_max+o3dB); %-3dB highest frequency in kHz
    low3db=fax(pos_max-e3dB); %-3dB lowest frequency in kHz
    bw3db=high3db-low3db;
% 10dB Duration
% Find the envelope of the signal
    env = abs(hilbert(click));
    [e_max pos_t]=max(env);
    % Obtain the -10dB time interval relative to the envelope of the signal
    thr10db = max(env)/(10^(10/20)); % 20*log (max(env)/limite10db)=10dB
%     int10db=find(env>thr10db);
%     dur10db=length(int10db)/fs*1e6;
%     
% clickenv = abs(hilbert(click));   % Envolvente de señal
% [max_t,pos_t] = max(clickenv);
clickenvflip = flipud(env);
dur10dB_first = pos_t - (find(clickenvflip(end-pos_t+1:end)<=(thr10db), 1, 'first')) + 1; 
dur10dB_last = pos_t + (find(clickenv(pos_t+1:end)<=(thr10db), 1, 'first'));
if isempty(dur10dB_last) == 1; 
    dur10dB_last = length(click);
end
if isempty(dur10dB_first) == 1; 
    dur10dB_first = 1;
end
dur10dB = (dur10dB_last - dur10dB_first)/fs*1e6; 

end

%Vanesa Reyes 01/2019