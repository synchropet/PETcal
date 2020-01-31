function varargout = DacStepper(varargin)
% DACSTEPPER MATLAB code for DacStepper.fig
%      DACSTEPPER, by itself, creates a new DACSTEPPER or raises the existing
%      singleton*.
%
%      H = DACSTEPPER returns the handle to a new DACSTEPPER or the handle to
%      the existing singleton*.
%
%      DACSTEPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DACSTEPPER.M with the given input arguments.
%
%      DACSTEPPER('Property','Value',...) creates a new DACSTEPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DacStepper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DacStepper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DacStepper

% Last Modified by GUIDE v2.5 12-May-2019 16:32:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DacStepper_OpeningFcn, ...
    'gui_OutputFcn',  @DacStepper_OutputFcn, ...
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


% --- Executes just before DacStepper is made visible.
function DacStepper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DacStepper (see VARARGIN)

% Choose default command line output for DacStepper
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DacStepper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DacStepper_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_start_scan.
function pushbutton1_start_scan_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_start_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton1_start_scan
get(handles.pushbutton1_start_scan,'value')
if get(handles.pushbutton1_start_scan,'value')
    set(handles.pushbutton1_start_scan,'BackgroundColor',[0 1 0],'string','Working')
    disp('*')
    if get(handles.pushbutton1_start_scan,'value')
        set(handles.pushbutton1_start_scan,'BackgroundColor',[0.9400 0.9400 0.9400],'string','Start','value',0)
    end
else
    disp('val 0')
    set(handles.pushbutton1_start_scan,'BackgroundColor',[1 0 0],'string','Aborting')
    pause(1)
    set(handles.pushbutton1_start_scan,'BackgroundColor',[0.9400 0.9400 0.9400],'string','Start','value',0)
end


function data_folder_Callback(hObject, eventdata, handles)
% hObject    handle to data_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_folder as text
%        str2double(get(hObject,'String')) returns contents of data_folder as a double


% --- Executes during object creation, after setting all properties.
function data_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
