function varargout = SSVEPGUI(varargin)
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
% Last Modified by GUIDE v2.5 28-Nov-2017 12:25:35
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
if isequal(getenv('computername'),'ME11W0207GRD01')
    [ deviceManager , flag ] = load_API('C:\Users\mmahmood31\Dropbox (GaTech)\YeoLab\_SSVEP\MATLAB-BioRadio-Recording-Labels\BioRadioSDK.dll');
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
global myDevice 
clear trainingData
global OUTPUT

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
numAxes = 9; 
axis_handles = zeros(1,numAxes);
for ch = 1:numAxes
    axis_handles(ch) = handles.(['axes',num2str(ch)]);
end

%Preallocating BPSignals
BioPotentialSignals = cell(1,numEnabledBPChannels+1);
if get(hObject,'Value') == 1
myDevice.StartAcquisition;
end
plotWindow = 5;
t = cell(numEnabledBPChannels, 1);
OUTPUT = zeros(1,50);
% Filter coeficients:
[Y,Sound_Fs] = audioread('beep_01a.mp3');
beep = audioplayer(Y,Sound_Fs);
% play(beep)
[~,b,a] = customFilt(zeros(1000,1), Fs, [1.0 40.0], 3);
INTERVAL = 10; % SECONDS
CURRENT_CLASS = 0;
CLASS = cell(1,1);
while get(hObject,'Value') == 1
    pause(0.05)
    second = length(BioPotentialSignals{1})/Fs
    if (second < INTERVAL)
        CURRENT_CLASS = 0;
    elseif second >= INTERVAL*1 && second < INTERVAL*2
        if CURRENT_CLASS ~= 1
            CURRENT_CLASS = 1;
            play(beep);
        end
    elseif second >= INTERVAL*2 && second < INTERVAL*3
        if CURRENT_CLASS ~= 2
            CURRENT_CLASS = 2;
            play(beep);
        end
    elseif second >= INTERVAL*3 && second < INTERVAL*4
        if CURRENT_CLASS ~= 3
            CURRENT_CLASS = 3;
            play(beep);
        end
    elseif second >= INTERVAL*4 && second < INTERVAL*5
        if CURRENT_CLASS ~= 4
            CURRENT_CLASS = 4;
            play(beep);
        end
    elseif second >= INTERVAL*5 && second < INTERVAL*6
        if CURRENT_CLASS ~= 5
            CURRENT_CLASS = 5;
            play(beep);
        end
    elseif second >= INTERVAL*6 && second < INTERVAL*7
        if CURRENT_CLASS ~= 0
            CURRENT_CLASS = 0;
            play(beep);
        end
    end
    CURRENT_CLASS
    cla(axis_handles(9), 'reset');
    hold( axis_handles(9), 'on' )
    for ch = 1:numEnabledBPChannels
        NewData = myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double';
        BioPotentialSignals{ch} = [BioPotentialSignals{ch}; NewData];
        if (ch == 1)
            BioPotentialSignals{1+numEnabledBPChannels} = [BioPotentialSignals{1+numEnabledBPChannels}; assignall(length(NewData), CURRENT_CLASS)];
        end
        if length(BioPotentialSignals{ch}) <= plotWindow*Fs
            t{ch} = (0:(length(BioPotentialSignals{ch})-1))*(1/Fs);
            if (ch<=8)
                if length(BioPotentialSignals{ch}) >= 31
                    filtered_data = filtfilt(b, a, BioPotentialSignals{ch} );
                    plot(axis_handles(ch), t{ch},filtered_data );
                end
                set(handles.(['axes',num2str(ch)]),'XLim',[0 plotWindow]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
            end
        else %once plot window is exceeded:
            if (ch<=8)
                t{ch} = ((length(BioPotentialSignals{ch})-(plotWindow*Fs-1)):length(BioPotentialSignals{ch}))*(1/Fs);
                filtered_data = filtfilt(b, a, BioPotentialSignals{ch}(end-plotWindow*Fs+1:end) );
                plot(axis_handles(ch),t{ch},filtered_data)
                set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', ['Channel ' num2str(ch)])
                set(handles.(['axes',num2str(ch)]),'XLim',[t{ch}(end)-plotWindow t{ch}(end)]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                % Plot PSD:
                [PSD(:, ch), f] = welch_psd(filtered_data, Fs, hann(1024));
                plot(axis_handles(9), f, PSD(:,ch));
                set(handles.(['axes',num2str(9)]),'XLim',[f(1) f(165)]);
            end
        end     %/if length(BioPotentialSignals{ch}) <= plotWindow*Fs
    end     %/for ch = 1:numEnabledBPChannels
    % PLOT PSD DATA:
end     %/while connected==1

if get(hObject,'Value') == 0
    mkdir('output_data')
    filename = ['output_data/SSVEP_DATA_' datestr(now,'yyyy_mm_dd__HH.MM.SS') '.mat']
    myDevice.StopAcquisition;  
    Trial = cell(1,numEnabledBPChannels);
    for i=1:numEnabledBPChannels + 1
        Trial{1,i} = BioPotentialSignals{i};
    end
    relevant_data = cell2mat(Trial);
    SamplingRate = Fs;
    
    assignin('base','NumberOfChannels',numEnabledBPChannels);
    assignin('base','TrialData',Trial)
    assignin('base','SamplingRate',SamplingRate);
    assignin('base','relevant_data',relevant_data);
    if(~isempty(relevant_data))
        save(filename, 'relevant_data');
    end
end
