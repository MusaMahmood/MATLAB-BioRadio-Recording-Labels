%% Analyse Data
clear;clc;close all;
load('output_data\SSVEP_DATA_2017_11_28__12.36.22.mat');
% Plot raw:
Fs = 250; %% implicit
[~,b,a] = customFilt(zeros(32,1), Fs, [1 40], 3);
wlen = 1024; h=64; nfft = 4096; K = sum(hamming(wlen, 'periodic'))/wlen; winLim = [1 40];
for i = 1:8
    figure(1);subplot(2,4,i); plot(relevant_data(:,i)); title(['Ch' num2str(i)])
    [PSD, f] = welch_psd(relevant_data(:,i), Fs, hann(1024));
    figure(2);subplot(4,2,i); plot(f,PSD); title(['Ch' num2str(i)]), xlim([5 40])
    [S1,f1,t1]=stft(filtfilt(b,a,relevant_data(:,i)), wlen, h, nfft, Fs); S2 = 20*log10(abs(S1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6);
    figure(3); subplot(2,4,i); imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),S2),xlim([min(t1) max(t1)]),ylim(winLim); set(gca,'YDir','normal');xlabel('Time, s');ylabel('Frequency, Hz');colormap(jet)
    cb = colorbar;ylabel(cb, 'Power (db)'); title(['Ch' num2str(i)]);
end
