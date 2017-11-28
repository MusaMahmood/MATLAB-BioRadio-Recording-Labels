function varargout = ssvep_gui_orig(varargin)
% SSVEPGUI MATLAB code for SSVEPGUI.fig
%      SSVEPGUI, by itself, creates a new SSVEPGUI or raises the existing
%      singleton*.
%
%      H = SSVEPGUI returns the handle to a new SSVEPGUI or the handle to
%      the existing singleton*.
%
%      SSVEPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSVEPGUI.M with the given input arguments.
%
%      SSVEPGUI('Property','Value',...) creates a new SSVEPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the SSVEPGUI before SSVEPGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSVEPGUI_OpeningFcn via varargin.
%
%      *See SSVEPGUI Options on GUIDE's Tools menu.  Choose "SSVEPGUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help SSVEPGUI
% Last Modified by GUIDE v2.5 27-Feb-2017 13:16:51
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSVEPGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SSVEPGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
% --- Executes just before SSVEPGUI is made visible.
function SSVEPGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSVEPGUI (see VARARGIN)
 
% global countFar countMiddle countClose
% Choose default command line output for SSVEPGUI
handles.output = hObject;
% global trainingData
% trainingData = cell(1,1);
% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes SSVEPGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
 
% --- Outputs from this function are returned to the command line.
function varargout = SSVEPGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
 
 
% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myDevice;
global deviceName;
count = 0;
 
if get(hObject,'Value') && count == 0
% load the BioRadio API using a MATLAB's .NET interface
if isequal(getenv('computername'),'EGR-YEO-12')
    [ deviceManager , flag ] = load_API('C:\Users\mahmoodms\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\BioRadioSDK.dll');
elseif isequal(getenv('computername'),'MusaMahmood-PC')
    [ deviceManager , flag ] = load_API('C:\Users\Musa Mahmood\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\BioRadioSDK.dll');
end
% input = full path to api dll file
% outputs = deviceManager object, success flag
 
if ~flag % if API not successfully loaded, do not continue
    return
end
 
% search for available sensors and select one
%
[ deviceName , macID , ok ] = BioRadio_Find( deviceManager );
% input = deviceManager object
% outputs = device name, macid, and flag if selection was canceled out
%
 
if ~ok %if no sensors selected, do not continue
    errordlg('Please select a BioRadio.')
    return   
% initialize BioRadio object
end
 
 
[ myDevice, flag ] = BioRadio_Connect ( deviceManager , macID , deviceName );
% input = deviceManager object, 64-bit mac address of BioRadio, and name of
% BioRadio
% outputs = BioRadio object, success flag for connection
%global myDevice;
if ~flag %if connection failed, do not continue
    return
end
count = 1;
 
else
    BioRadio_Disconnect( myDevice )
    count = 0;
     
end
% Hint: get(hObject,'Value') returns toggle state of togglebutton1
 
 
% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myDevice Idx 
clear trainingData
global totalCount which_pc OUTPUT
 
totalCount = cell(2);
totalCount{1} = 0;
totalCount{2} = 0;
BioRadio_Name = 'EEG-SSVEP';
numEnabledBPChannels = double(myDevice.BioPotentialSignals.Count);
 
if numEnabledBPChannels == 0
    myDevice.Disconnect;
    RawBioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end
 
Fs = double(myDevice.BioPotentialSignals.SamplesPerSecond);
%Preallocating and setting up which area in the SSVEPGUI the plot will go into
    % First two are raw data
    % Next two are data analysis features. 
numAxes = 14; 
axis_handles = zeros(1,numAxes);
for ch = 1:numAxes
    axis_handles(ch) = handles.(['axes',num2str(ch)]);
end
 
%Preallocating BPSignals
BioPotentialSignals = cell(1,numEnabledBPChannels);
filteredEEG{ch} = cell(1,numEnabledBPChannels);
Idx = cell(1,numEnabledBPChannels);
if get(hObject,'Value') == 1
myDevice.StartAcquisition;
end
plotWindow = 5;
plotGain_BP = 1;
fft_len = plotWindow*Fs;
% fft_len = 2^(nextpow2(plotWindow*Fs)); 
% USE WITH fft(X,fft_len_pow2)
f0 = [8.75 18.5];
xl = [5 20];
N = 5;
spect_1 = handles.axes5;
spect_2 = handles.axes6;
t = cell(numEnabledBPChannels, 1);
fp1_data_unfilt = zeros(plotWindow*Fs,1);
fp2_data_unfilt = zeros(plotWindow*Fs,1);
eog3_data_unfilt = zeros(plotWindow*Fs,1); 
eog4_data_unfilt = zeros(plotWindow*Fs,1);
%LOAD tXtY:: (NOW EMBEDDED IN CODE AS STATIC VAR)
% load('allEOGtD.mat');
% SSVEP Training Data:
% load('tXtY_SSVEP.mat');
 
startLen = 499;
ln = startLen;
cIdx = cell(numEnabledBPChannels, 1);
cIdx1 = cell(numEnabledBPChannels, 1);
cIdx2 = cell(numEnabledBPChannels, 1);
for i = 1:numEnabledBPChannels
    cIdx{i} = 1;
    cIdx1{i} = 1;
    cIdx2{i} = 1;
end
 
YEOG = cell(1,1);
% OUTPUT = cell(1,1);
W = cell(4,1); % # channels = # windows
op=1;
waitPlot = 250;
waitEOG = 60;
WAITDEFAULT = 60;
wait = WAITDEFAULT;
cnt = 1;
ocnt = 1;
maxTimeout = 1750;
r = startLen:WAITDEFAULT:maxTimeout;
% HISTORY = zeros(1,10);
plotData = true;
OUTPUT = zeros(1,50);
OUT = zeros(1,size(r,2));
while get(hObject,'Value') == 1
    pause(0.05)
    for ch = 1:numEnabledBPChannels
        BioPotentialSignals{ch} = [BioPotentialSignals{ch}; myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
        Idx{ch} = 1:length(BioPotentialSignals{ch});            
        %Plot the Axes in the SSVEPGUI
        if length(BioPotentialSignals{ch}) <= plotWindow*Fs
            t{ch} = (0:(length(BioPotentialSignals{ch})-1))*(1/Fs);
            if ch==1 || ch==2
                if length(BioPotentialSignals{ch}) >= 31
                    filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}, Fs, f0, N);
                    plot(axis_handles(ch),t{ch},filteredEEG{ch});
                end
                set(handles.(['axes',num2str(ch)]),'XLim',[0 plotWindow]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
            end
            if ch==1
                set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 1 Data')
                if length(BioPotentialSignals{ch})>Fs
                    fp1_data_filtered = eeg_h_custom(BioPotentialSignals{ch}, Fs, f0, N);
                    [f, P1] = get_fft_data(fp1_data_filtered, Fs);
                    plot(axis_handles(3),f,P1);
                    set(handles.(['axes',num2str(3)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                end
            elseif ch==2
                set(get(handles.(['axes',num2str(ch)]),'Title'),'String','Channel 2 Data')
                if length(BioPotentialSignals{ch})>Fs
                    fp2_data_filtered = eeg_h_custom(BioPotentialSignals{ch}, Fs, f0, N);
                    [f, P1] = get_fft_data(fp2_data_filtered, Fs);
                    plot(axis_handles(4),f,P1);
                    set(handles.(['axes',num2str(4)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                end
            end
        else %once plot window is exceeded:
            if ch==1
                 t{1} = ((length(BioPotentialSignals{ch})-(plotWindow*Fs-1)):length(BioPotentialSignals{ch}))*(1/Fs);
                 filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*Fs+1:end), Fs, f0, N);
                 plot(axis_handles(ch),t{1},filteredEEG{ch})
                 set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 1 Filtered')
            elseif ch==2
                 t{2} = ((length(BioPotentialSignals{ch})-(plotWindow*Fs-1)):length(BioPotentialSignals{ch}))*(1/Fs);
                 filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*Fs+1:end), Fs, f0, N);
                 plot(axis_handles(ch),t{2},filteredEEG{ch})
                 set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 2 Filtered')
            elseif ch==3
                 t{3} = ((length(BioPotentialSignals{ch})-(plotWindow*Fs-1)):length(BioPotentialSignals{ch}))*(1/Fs);
            elseif ch==4
                 t{4} = ((length(BioPotentialSignals{ch})-(plotWindow*Fs-1)):length(BioPotentialSignals{ch}))*(1/Fs);
            end
            if ch == 1 || ch == 2
                set(handles.(['axes',num2str(ch)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                set(handles.(['axes',num2str(ch)]),'YLim',[-2E-4 2E-4]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
            end
            if ch==1
                % Window and Filter
                fp1_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*Fs+1:end);
                fp1_data_filtered = eeg_h_custom(fp1_data_unfilt, Fs, f0, N);
                [f, P1] = get_fft_data(fp1_data_filtered, Fs);
                plot(axis_handles(3),f,P1);
                set(handles.(['axes',num2str(3)]),'XLim',xl);
                set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                % Spectrogram (STFT):
                if which_pc == 0
%                         [S, Fspect, T, P] = spectrogram(fp1_data_filtered, 5*Fs,4*Fs,10*Fs,Fs);
%                         imagesc(spect_1, T, Fspect(Fspect<30 & Fspect>1), 10*log10(P(Fspect<30 & Fspect>1,:)));
%                         set(spect_1,'YDir','normal')
%                         cb = colorbar(spect_1);
%                         ylabel(cb, 'Power (db)')
%                         colormap(spect_1,jet)
%                         set(get(handles.(['axes',num2str(5)]), 'XLabel'), 'String', 'Time (s)')
%                         set(get(handles.(['axes',num2str(5)]), 'YLabel'), 'String', 'Frequency (Hz)')
%                         set(get(handles.(['axes',num2str(5)]), 'Title'), 'String', 'Spectrogram (Fp1)')
                end
                % Pwelch
%                 if mod(length(BioPotentialSignals{ch}),2) == 0
%                     hW = hannWin(length(BioPotentialSignals{ch}));
%                     temp25926 = eegcfilt(fp1_data_unfilt);
%                     [PSD, fPSD] = welch_psd(temp25926,Fs,hW);
%                     plot(axis_handles(5), fPSD, PSD); 
%                     set(handles.(['axes',num2str(5)]),'XLim',xl);
%                     set(get(handles.(['axes',num2str(5)]), 'XLabel'), 'String', 'Frequency (Hz)')
%                     set(get(handles.(['axes',num2str(5)]), 'YLabel'), 'String', 'Power (dB)')
%                     set(get(handles.(['axes',num2str(5)]), 'Title'), 'String', 'Pwelch (Fp1)')
%                 end
%                 [Pxx, F] = pwelch(fp1_data_filtered,[],[],250);
                %EOG filt:
                eog_data_1 = eog_h_fcn(fp1_data_unfilt, Fs);
                plot(axis_handles(9),t{ch},eog_data_1);
                set(handles.(['axes',num2str(9)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                set(handles.(['axes',num2str(9)]),'YLim',[-2.5E-4 2.5E-4]);
                set(get(handles.(['axes',num2str(9)]), 'XLabel'), 'String', 'Time (s)')
                set(get(handles.(['axes',num2str(9)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(9)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
            elseif ch==2
                fp2_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*Fs+1:end);
                fp2_data_filtered = eeg_h_custom(fp2_data_unfilt, Fs, f0, N);
                [f, P1] = get_fft_data(fp2_data_filtered, Fs);
                plot(axis_handles(4),f,P1);
                set(handles.(['axes',num2str(4)]),'XLim',xl);
                set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P2(f)|')
                set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                % Spect:
                if which_pc == 0
%                         [S, Fspect, T, P] = spectrogram(fp2_data_filtered, 5*Fs,4*Fs,10*Fs,Fs);
%                         imagesc(spect_2, T, Fspect(Fspect<30 & Fspect>1), 10*log10(P(Fspect<30 & Fspect>1,:)));
%                         set(spect_2,'YDir','normal') %gca
%                         cb2 = colorbar(spect_2);
%                         ylabel(cb2, 'Power (db)')
%                         colormap(spect_2,jet)
%                         set(get(handles.(['axes',num2str(6)]), 'XLabel'), 'String', 'Time (s)')
%                         set(get(handles.(['axes',num2str(6)]), 'YLabel'), 'String', 'Frequency (Hz)')
%                         set(get(handles.(['axes',num2str(6)]), 'Title'), 'String', 'Spectrogram (Fp2)')
                end
                % P welch
%                 if mod(length(BioPotentialSignals{ch}),2) == 0
%                     hW = hannWin(length(BioPotentialSignals{ch}));
%                     [PSD, fPSD] = welch_psd(eegcfilt(fp1_data_unfilt),Fs,hW);
%                     plot(axis_handles(6), fPSD, PSD); 
%                     set(handles.(['axes',num2str(6)]),'XLim',xl);
%                     set(get(handles.(['axes',num2str(6)]), 'XLabel'), 'String', 'Frequency (Hz)')
%                     set(get(handles.(['axes',num2str(6)]), 'YLabel'), 'String', 'Power (dB)')
%                     set(get(handles.(['axes',num2str(6)]), 'Title'), 'String', 'Pwelch (Fp2)')
%                 end
%                 [Pxx, F] = pwelch(fp2_data_filtered,[],[],250);
                %EOG Filt
                eog_data_1 = eog_h_fcn(fp2_data_unfilt, Fs);
                plot(axis_handles(10),t{ch},eog_data_1);
                set(handles.(['axes',num2str(10)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                set(handles.(['axes',num2str(10)]),'YLim',[-2.5E-4 2.5E-4]);
                set(get(handles.(['axes',num2str(10)]), 'XLabel'), 'String', 'Time (s)')
                set(get(handles.(['axes',num2str(10)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(10)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
            end
            if ch==3
                eog3_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*Fs+1:end);
                eog_data_1 = eog_h_fcn(eog3_data_unfilt, Fs);
                plot(axis_handles(7),t{ch},eog_data_1);
                set(handles.(['axes',num2str(7)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                set(handles.(['axes',num2str(7)]),'YLim',[-4E-4 4E-4]);
                set(get(handles.(['axes',num2str(7)]), 'XLabel'), 'String', 'Time (s)')
                set(get(handles.(['axes',num2str(7)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(7)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                %%FFT CH3:
                CH3_data_filtered = eeg_h_custom(eog3_data_unfilt, Fs, f0, N);
                [f, P1] = get_fft_data(CH3_data_filtered, Fs);
                plot(axis_handles(11),f,P1);
                set(handles.(['axes',num2str(11)]),'XLim',xl);
                set(get(handles.(['axes',num2str(11)]), 'XLabel'), 'String', 'f (Hz)')
                set(get(handles.(['axes',num2str(11)]), 'YLabel'), 'String', '|P1(f)|')
                set(get(handles.(['axes',num2str(11)]), 'Title'), 'String', 'FFT(Ch3, Fpz)')
                %PSD CH3:
%                 [Pxx, F] = pwelch(CH3_data_filtered,[],[],250);
%                 plot(axis_handles(12), Pxx); 
%                 set(handles.(['axes',num2str(12)]),'XLim',xl);
%                 set(get(handles.(['axes',num2str(12)]), 'XLabel'), 'String', 'Frequency (Hz)')
%                 set(get(handles.(['axes',num2str(12)]), 'YLabel'), 'String', 'Power (dB)')
%                 set(get(handles.(['axes',num2str(12)]), 'Title'), 'String', 'Pwelch (Ch3)')
            elseif ch==4
                eog4_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*Fs+1:end);
                eog_data_2 = eog_h_fcn(eog4_data_unfilt, Fs);
                plot(axis_handles(8),t{ch},eog_data_2);
                set(handles.(['axes',num2str(8)]),'XLim',[t{2}(end)-plotWindow t{2}(end)]);
                set(handles.(['axes',num2str(8)]),'YLim',[-4E-4 4E-4]);
                set(get(handles.(['axes',num2str(8)]), 'XLabel'), 'String', 'Time (s)')
                set(get(handles.(['axes',num2str(8)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(8)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                %ffT
                CH4_data_filtered = eeg_h_custom(eog4_data_unfilt, Fs, f0, N);
                [f, P1] = get_fft_data(CH4_data_filtered, Fs);
                plot(axis_handles(13),f,P1);
                set(handles.(['axes',num2str(13)]),'XLim',xl);
                set(get(handles.(['axes',num2str(13)]), 'XLabel'), 'String', 'f (Hz)')
                set(get(handles.(['axes',num2str(13)]), 'YLabel'), 'String', '|P1(f)|')
                set(get(handles.(['axes',num2str(13)]), 'Title'), 'String', 'FFT(Ch4, Fpz)')
                %PSD CH4:
%                 [Pxx, F] = pwelch(CH4_data_filtered,[],[],250); 
%                 plot(axis_handles(14), Pxx); 
%                 set(handles.(['axes',num2str(14)]),'XLim',xl);
%                 set(get(handles.(['axes',num2str(14)]), 'XLabel'), 'String', 'Frequency (Hz)')
%                 set(get(handles.(['axes',num2str(14)]), 'YLabel'), 'String', 'Power (dB)')
%                 set(get(handles.(['axes',num2str(14)]), 'Title'), 'String', 'Pwelch (Ch4)')
            end
            %% %%%%%% FULL HYBRID CLASSIFIER FUNCTION CALLS: %%%%%% %%
            %CLASSIFY EOG DATA:
               %take each ch unfilt (5s) & cut into 5 pieces, & pass
               %thru fullHybridClassifier: Put in columns:
               %{
            for i = 1:numEnabledBPChannels
                W_EOG(i,:) = BioPotentialSignals{i}(end-249:end); %last second of data; always
                b1(i) = length(BioPotentialSignals{i}) > cIdx{i}+wait; % has 1 second passed?
                b2(i) = length(BioPotentialSignals{i}) > cIdx1{i}+waitEOG;
                b3(i) = length(BioPotentialSignals{i}) > cIdx2{i}+waitPlot;
            end
            if sum(b2) == numEnabledBPChannels
                [YEOG{1}] = fHC(W_EOG(1,:), W_EOG(2,:), W_EOG(3,:), W_EOG(4,:), Fs, true);
                if(YEOG{1}(1)~=0)
                    RESULT_EOG = YEOG{1}(1)
                end
                for i = 1:numEnabledBPChannels
                    cIdx1{i} = length(BioPotentialSignals{i});       %Assign new current indx
                end
            end
            if (YEOG{1}(1) == 2) || (YEOG{1}(1) == 1) %Double Blink Detected: Reset Command (Stop and reset):
                %Reset cell (ON ANDROID, RESET ARRAYS TO ZERO).
                for i = 1:numEnabledBPChannels
                    cIdx{i} = length(BioPotentialSignals{i});       %Assign new current indx, forced to wait 1s
                end
                OUT = zeros(1,size(r,2));
                ln = startLen;
                cnt = 1;
                wait = 750; % changed from 250
            end
            if sum(b1) == numEnabledBPChannels
                fprintf('len = %d\r\n',ln);
                if ln<=maxTimeout %12s
                    for i = 1:numEnabledBPChannels
                        cIdx{i} = length(BioPotentialSignals{i});       %Assign new current indx
                        W{i} = BioPotentialSignals{i}(end-ln:end);
                    end
                    if sum(b3)==numEnabledBPChannels
                        plotData = true;
                        for i = 1:numEnabledBPChannels
                            cIdx2{i} = length(BioPotentialSignals{i});       %Assign new current indx
                        end
                    else
                        plotData = false;
                    end
                    [~, FEA0] = fHC(W{1}, W{2}, W{3}, W{4}, Fs, false, plotData);
                    [History(cnt,:), OUT(cnt)] = featureAnalysis(FEA0,ln);
                    meanH = mean(History(1:cnt,:))
                    if OUT(cnt)~=0
                        countH(cnt) = countOccurrences(OUT(:,1:cnt), OUT(cnt));
                    else
                        countH(cnt) = 0;
                    end
                    if countH(cnt)>3 %(max(meanH)>7) &&
                        OUTPROPER(cnt) = OUT(cnt); %This is the command.;
                        OUTPUT(ocnt) = OUTPROPER(cnt);
                        ocnt = ocnt+1;
                    else
                        OUTPROPER(cnt) = 0;
                    end
                    OUT
                    OUTI = OUT(cnt)
                     
                    % Mean of 4 ch:
                    ln = ln+WAITDEFAULT;
                    cnt = cnt+1;
                    wait = WAITDEFAULT;
                else
                    OUT = zeros(1,size(OUT,2));
                    ln = startLen;
                    cnt = 1;
                    wait = WAITDEFAULT;
                end
            end
               %}
            %{
            YEOG{1} = fullHybridClassifier(W_EOG(1,:), W_EOG(2,:), W_EOG(3,:), W_EOG(4,:), Fs, true);
            if (YEOG{1}(1) == 1)%Double Blink Detected: Reset Command (Stop and reset):
                Y = cell(7,1); %Reset cell (ON ANDROID, RESET ARRAYS TO ZERO).
                for i = 1:numEnabledBPChannels
                    cIdx{i} = length(BioPotentialSignals{i});       %Assign new current indx, forced to wait 1s
                end
                ln = 249;
                wait = 500;
                fprintf('\n >>>>OUTPUT = DOUBLE BLINK :: RESET COMMAND \n \n');
            end
            if sum(b1) == numEnabledBPChannels
                for i = 1:numEnabledBPChannels
                    cIdx{i} = length(BioPotentialSignals{i});       %Assign new current indx
                    W{i} = BioPotentialSignals{i}(end-ln:end);
                end
                if ln<500%ln==249
                    Y{1} = fullHybridClassifier(W{1}, W{2}, W{3}, W{4}, Fs, YEOG{1}(1)==1)';
                    fprintf('Y: \n'); disp(Y{1})
                    ln = ln+250;
                    wait = 250;
                elseif ln >= 500 && ln < 2000
                    Y{2} = fullHybridClassifier(W{1}, W{2}, W{3}, W{4}, Fs, YEOG{1}(1)==1)';
                    fprintf('Y: \n'); disp(Y{2})
                    B1 = Y{2}(1) == 1; %6==7
                    B2 = (Y{2}(4)) && (Y{2}(5) == 1);
                    if B1 && B2
                        if Y{2}(6) == Y{2}(2)
                            %Almost certainly correct: (for any time window
                            %> 500samples (2s)
                            OUTPUT{1}(op) = Y{2}(6);        %%% ACCEPT %%% 
                            ln = 249;
                        else
                            OUTPUT{1}(op) = 0;              %%% REJECT %%%
                            ln = ln+250;
                        end
                    else
                        if ln > 1000 %TODO: NOT SURE ABOUT THIS?!
                            if B1
                                OUTPUT{1}(op) = Y{2}(6);    %%% ACCEPT %%%
                                ln=249;
                            elseif B2
                                OUTPUT{1}(op) = Y{2}(2);    %%% ACCEPT %%%
                                ln=249;
                            else
                                OUTPUT{1}(op) = 0;          %%% REJECT %%%
                                ln=ln+250;
                            end
                        else
                            OUTPUT{1}(op) = 0;              %%% REJECT %%%
                            ln=ln+250;
                        end
                    end
%                     if ~B1 && ~B2
%                         OUTPUT{1}(op) = 0;
%                         ln=ln+250;
%                     end
                    if OUTPUT{1}(op) ~= 0
                        fprintf('\n >>>>OUTPUT = %d\n',OUTPUT{1}(op));
                    end
                    op=op+1;
                    wait = 250;
                elseif ln >= 2000 % RESET:
                    Y = cell(2,1);
                    ln = 249;
                    fprintf('\n >>>>Timer exceeded, RESETTING\n');
                end     %/if ln<499
            end     %/if sum(b1) == numEnabledBPChannels
            %}
            %% %%%%%% NEW HYBRID CLASSIFIER FUNCTION CALLS: %%%%%% %%
        end     %/if length(BioPotentialSignals{ch}) <= plotWindow*Fs
    end     %/for ch = 1:numEnabledBPChannels
end     %/while connected==1
 
if get(hObject,'Value') == 0
    myDevice.StopAcquisition;  
    Trial = cell(1,numEnabledBPChannels);
    for i=1:numEnabledBPChannels
        Trial{1,i} = BioPotentialSignals{i};
    end
    assignin('base','NumberOfChannels',numEnabledBPChannels);
    assignin('base','Trial',Trial)
    SamplingRate = Fs;
    assignin('base','SamplingRate',SamplingRate);
    %% Change into button function.
    OP_1{1} = OUTPUT;
    H_Notes = [handles.edit2 handles.edit_ChannelLocations, handles.editImpedanceValues, ...
        handles.editElectrodeType, handles.editStimulusSource, handles.editMiscInfo];
    RecordingNotes = cell(length(H_Notes)+1, 2);
    RecordingNotes{1,1} = 'Filename (.mat)'; RecordingNotes{2,1} = 'Channel Locations';
    RecordingNotes{3,1} = 'Electrode Impedance'; RecordingNotes{4,1} = 'Electrode Type';
    RecordingNotes{5,1} = 'Source of Stimulus'; RecordingNotes{6,1} = 'Misc Notes';
    for i=1:length(H_Notes)
        RecordingNotes{i,2} = get(H_Notes(i),'String');
        if isempty(RecordingNotes{i,2})
            RecordingNotes{i,2} = 'No Notes Recorded for This Session'; 
        end
    end
    RecordingNotes{length(H_Notes)+1,1} = 'Sampling Rate';
    RecordingNotes{length(H_Notes)+1,2} = num2str(SamplingRate);
    assignin('base','RecordingNotes',RecordingNotes);
    assignin('base','OUTPUT',OP_1);
    tD = trainingData;
    assignin('base','tD',tD);
    %% Auto-save variables
    filename = get(handles.edit2,'String');
    if isempty(filename)
        c=clock;
        filename = ['RawBioRadioData_',num2str(c(2)),'-',num2str(c(3)),'-',num2str(c(1)),'_',num2str(c(4)),'.',num2str(c(5)),',',num2str(c(6)),'.mat'];
    end
    save(filename,'Trial','SamplingRate','RecordingNotes','trainingData','OUTPUT');
end
 
% --- Executes on button press in pushbuttonDoubleBlink.
function pushbuttonDoubleBlink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDoubleBlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx doubleBlinkIdx doubleBlinkCount trainingData totalCount
if doubleBlinkCount==1
    doubleBlinkIdx = zeros(1,1);
     
    doubleBlinkIdx(1,1) = size(Idx{1,1},2)
     
    doubleBlinkCount = doubleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = doubleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 2; %class
     
    assignin('base','tD',trainingData);
else
    doubleBlinkIdx(1,1) = size(Idx{1,1},2)
     
    doubleBlinkCount = doubleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = doubleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 2; %class
    assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbuttonSingleBlink.
function pushbuttonSingleBlink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSingleBlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx singleBlinkIdx singleBlinkCount trainingData totalCount
if singleBlinkCount==1
    singleBlinkIdx = zeros(1,1);
     
    singleBlinkIdx(1,1) = size(Idx{1,1},2)
     
    singleBlinkCount = singleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = singleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 1; %class
     
%     assignin('base','tD',trainingData);
else
    singleBlinkIdx(1,1) = size(Idx{1,1},2)
     
    singleBlinkCount = singleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = singleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 1; %class
    
%     assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbuttonUp.
function pushbuttonUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveUpIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveUpIdx = zeros(1,1);
     
    eyeMoveUpIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveUpIdx; %index
    trainingData{1}(totalCount{1},2) = 3; %class
%     assignin('base','tD',trainingData);
else
    eyeMoveUpIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveUpIdx; %index
    trainingData{1}(totalCount{1},2) = 3; %class
     
%     assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbuttonDown.
function pushbuttonDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveDownIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveDownIdx = zeros(1,1);
     
    eyeMoveDownIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveDownIdx; %index
    trainingData{1}(totalCount{1},2) = 4; %class
%     assignin('base','tD',trainingData);
else
    eyeMoveDownIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveDownIdx; %index
    trainingData{1}(totalCount{1},2) = 4; %class
     
%     assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbuttonLeft.
function pushbuttonLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveLeftIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveLeftIdx = zeros(1,1);
     
    eyeMoveLeftIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveLeftIdx; %index
    trainingData{1}(totalCount{1},2) = 5; %class
%     assignin('base','tD',trainingData);
else
    eyeMoveLeftIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveLeftIdx; %index
    trainingData{1}(totalCount{1},2) = 5; %class
     
%     assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbuttonRight.
function pushbuttonRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveRightIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveRightIdx = zeros(1,1);
     
    eyeMoveRightIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveRightIdx; %index
    trainingData{1}(totalCount{1},2) = 6; %class
%     assignin('base','tD',trainingData);
else
    eyeMoveRightIdx(1,1) = size(Idx{1,1},2)
     
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
     
    trainingData{1}(totalCount{1},1) = eyeMoveRightIdx; %index
    trainingData{1}(totalCount{1},2) = 6; %class
     
%     assignin('base','tD',trainingData);
end
 
 
% --- Executes on button press in pushbutton11. (10Hz)
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq10Idx freqCount trainingData totalCount
if freqCount==1
    freq10Idx = zeros(1,1);
    freq10Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq10Idx;
    trainingData{1}(totalCount{2},2) = 7; %class
%     assignin('base','tD',trainingData);
else
    freq10Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq10Idx;
    trainingData{1}(totalCount{2},2) = 7; %class
%     assignin('base','tD',trainingData);
end
 
 
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq12Idx freqCount trainingData totalCount
if freqCount==1
    freq12Idx = zeros(1,1);
    freq12Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq12Idx;
    trainingData{1}(totalCount{2},2) = 8; %class
%     assignin('base','tD',trainingData);
else
    freq12Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq12Idx;
    trainingData{1}(totalCount{2},2) = 8; %class
%     assignin('base','tD',trainingData);
end
 
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq15Idx freqCount trainingData totalCount
if freqCount==1
    freq15Idx = zeros(1,1);
    freq15Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq15Idx;
    trainingData{1}(totalCount{2},2) = 9; %class
%     assignin('base','tD',trainingData);
else
    freq15Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq15Idx;
    trainingData{1}(totalCount{2},2) = 9; %class
%     assignin('base','tD',trainingData);
end
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq16Idx freqCount trainingData totalCount
 
if freqCount==1
    freq16Idx = zeros(1,1);
    freq16Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq16Idx;
    trainingData{1}(totalCount{2},2) = 10; %class
%     assignin('base','tD',trainingData);
else
    freq16Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
     
    trainingData{1}(totalCount{2},1) = freq16Idx;
    trainingData{1}(totalCount{2},2) = 10; %class
%     assignin('base','tD',trainingData);
end
 
 
%------------------------ %%%%%%%%%%%% ?IGNORE? %%%%%%%%%%% ---------------------------%
 
function edit_ChannelLocations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelLocations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit_ChannelLocations as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelLocations as a double
 
 
% --- Executes during object creation, after setting all properties.
function edit_ChannelLocations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelLocations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
 
 
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
 
 
function editImpedanceValues_Callback(hObject, eventdata, handles)
% hObject    handle to editImpedanceValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of editImpedanceValues as text
%        str2double(get(hObject,'String')) returns contents of editImpedanceValues as a double
 
 
% --- Executes during object creation, after setting all properties.
function editImpedanceValues_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImpedanceValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function editElectrodeType_Callback(hObject, eventdata, handles)
% hObject    handle to editElectrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of editElectrodeType as text
%        str2double(get(hObject,'String')) returns contents of editElectrodeType as a double
 
 
% --- Executes during object creation, after setting all properties.
function editElectrodeType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editElectrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function editStimulusSource_Callback(hObject, eventdata, handles)
% hObject    handle to editStimulusSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of editStimulusSource as text
%        str2double(get(hObject,'String')) returns contents of editStimulusSource as a double
 
 
% --- Executes during object creation, after setting all properties.
function editStimulusSource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStimulusSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function editMiscInfo_Callback(hObject, eventdata, handles)
% hObject    handle to editMiscInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of editMiscInfo as text
%        str2double(get(hObject,'String')) returns contents of editMiscInfo as a double
 
 
% --- Executes during object creation, after setting all properties.
function editMiscInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMiscInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end