%% Analyse Data
%{
% 15s cleanup
clear;clc;close all
nmc = 10;
classes = 2;
dirname = ['output_data\r5_S002\' num2str(nmc) 's\'];
files = dir([dirname '*.mat']);
[~,b,a] = customFilt(zeros(32,1), 250.0, [5, 100], 3);
for f = 1:length(files)
    load([dirname files(f).name]);
    relevant_data = (relevant_data(1:250*classes*nmc, :));
    for ch = 1:8
        relevant_data(:,ch) = filtfilt(b,a,relevant_data(:,ch));
    end
    dn = [dirname(1:end-4), 'filtered_2class_' +num2str(nmc) 's\'];
    mkdir(dn);
    save([dn files(f).name], 'relevant_data');
end
%}
% %{
clear;clc;close all;
% dirr = dir('output_data\r5_S001\10sfiltered_5class');
load('output_data\r5_S002\SSVEP_DATA_2017_12_05__17.17.32.mat');
% load('output_data\SSVEP_DATA_2017_12_05__17.15.57.mat');
% Plot raw:
NUMBER_CHANNELS = 2;
Fs = 250; %% implicit
[~,b,a] = customFilt(zeros(32,1), Fs, [5 40], 3);
wlen = 1024; h=64; nfft = 1024; K = sum(hamming(wlen, 'periodic'))/wlen; winLim = [0 125];
for i = 1:NUMBER_CHANNELS
%     figure(1);subplot(2,4,i); plot(relevant_data(:,i)); title(['Ch' num2str(i)])
    [PSD, f] = welch_psd(relevant_data(:,i), Fs, hann(1024));
    figure(2);subplot(4,2,i); plot(f,PSD); title(['Ch' num2str(i)]), xlim([5 40])
    [S1,f1,t1]=stft(filtfilt(b,a,relevant_data(:,i)), wlen, h, nfft, Fs); S2 = 20*log10(abs(S1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6);
    figure(3); subplot(2,4,i); imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),S2),xlim([min(t1) max(t1)]),ylim(winLim); set(gca,'YDir','normal');xlabel('Time, s');ylabel('Frequency, Hz');colormap(jet)
    cb = colorbar;ylabel(cb, 'Power (db)'); title(['Ch' num2str(i)]);
end
%}