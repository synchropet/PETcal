function varargout = PET_spi_editor(varargin)
% PET_SPI_EDITOR M-file for PET_spi_editor.fig
%      PET_SPI_EDITOR, by itself, creates a new PET_SPI_EDITOR or raises the existing
%      singleton*.
%
%      H = PET_SPI_EDITOR returns the handle to a new PET_SPI_EDITOR or the handle to
%      the existing singleton*.
%
%      PET_SPI_EDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PET_SPI_EDITOR.M with the given input arguments.
%
%      PET_SPI_EDITOR('Property','Value',...) creates a new PET_SPI_EDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PET_spi_editor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PET_spi_editor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PET_spi_editor

% Last Modified by GUIDE v2.5 20-Aug-2015 10:47:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PET_spi_editor_OpeningFcn, ...
                   'gui_OutputFcn',  @PET_spi_editor_OutputFcn, ...
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


% --- Executes just before PET_spi_editor is made visible.
function PET_spi_editor_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PET_spi_editor (see VARARGIN)

% Choose default command line output for PET_spi_editor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PET_spi_editor wait for user response (see UIRESUME)
% uiwait(handles.PET_spi_editor);
CurFig = hObject; pause(0.1);
setappdata(0,'hPET_spi_editor',CurFig);
Cfg = disp_init_info();
Cfg.NC = 32; % number of channels
Cfg.NA = 12; % number of asics
Cfg.data = zeros(Cfg.NC*Cfg.NA,4);
setappdata(CurFig,'Cfg',Cfg);
setappdata(CurFig,'fun_PET_spi_editor',  @fun_PET_spi_editor);

% --- Outputs from this function are returned to the command line.
function varargout = PET_spi_editor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Cfg = disp_init_info()
try
    Cfg.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    Cfg.hostname = ''; 
end
Cfg.folder=pwd;
[~, scriptname] =fileparts(mfilename('fullpath'));
Cfg.modulename = scriptname;
disp(sprintf('%s started %s on %s runing %s',scriptname,datestr(now),Cfg.hostname,computer()))

% --- Executes on button press in pushbutton1_load.
function pushbutton1_load_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton1_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

[File2read, Path2read,tmp] = uigetfile([Cfg.folder,'/*.mat'],'Select gain file to read');
if tmp==0, % cancel
    return
else
    Cfg.folder = Path2read;
    Cfg.file = File2read;
    try
        q=load([Path2read,File2read]);
        Cfg.data = q.data;
        set(handles.filename,'string',[Cfg.folder,Cfg.file])
    catch
        errordlg(sprintf('Cannot read file %s',File2read))
    end
    setappdata(CurFig,'Cfg',Cfg)
end

% --- Executes on selection change in asic_number.
function asic_number_Callback(hObject, eventdata, handles)
% hObject    handle to asic_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns asic_number contents as cell array
%        contents{get(hObject,'Value')} returns selected item from asic_number
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

anum = get(handles.asic_number,'Value');

curtbl = Cfg.data( (1+(anum-1)*Cfg.NC):(anum*Cfg.NC),: );
set(handles.gain_table,'Data',curtbl);



% --- Executes during object creation, after setting all properties.
function asic_number_CreateFcn(hObject, ~, ~)
% hObject    handle to asic_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filename_Callback(~, ~, ~)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, ~, ~)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2_save.
function pushbutton2_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

[File2read, Path2read,tmp] = uiputfile([Cfg.folder,'/*.mat'],'Select gain file to write');
if tmp==0, % cancel
    return
else
    Cfg.folder = Path2read;
    Cfg.file = File2read;
    try
        data = Cfg.data; 
        save([Path2read,File2read], 'data');
        set(handles.filename,'string',[Cfg.folder,Cfg.file])
    catch
        errordlg(sprintf('Cannot write file %s',File2read))
    end
    setappdata(CurFig,'Cfg',Cfg)
end


% --- Executes when entered data in editable cell(s) in gain_table.
function gain_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to gain_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and set_selection indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

val = eventdata.NewData;
n = eventdata.Indices(1);
m = eventdata.Indices(2);

anum = get(handles.asic_number,'Value');

Cfg.data((anum-1)*Cfg.NC+n,m) = val;
setappdata(CurFig,'Cfg',Cfg);


% --- Executes on button press in pushbutton3_set.
function pushbutton3_set_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

curtbl = get(handles.gain_table,'Data');
value = str2double(get(handles.value_set,'String'));
anum = get(handles.asic_number,'Value');

column = get(handles.set_selection,'value');

if column <4,
    value = (value>0);
    set(handles.value_set,'string',sprintf('%d',value));
else
    value = mod(value,32);
    set(handles.value_set,'string',sprintf('%d',value));
end

curtbl(:,column) = value;
set(handles.gain_table,'Data',curtbl);

Cfg.data( (1+(anum-1)*Cfg.NC):(anum*Cfg.NC),: ) = curtbl;
setappdata(CurFig,'Cfg',Cfg);


function value_set_Callback(~, ~, ~)
% hObject    handle to value_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value_set as text
%        str2double(get(hObject,'String')) returns contents of value_set as a double


% --- Executes during object creation, after setting all properties.
function value_set_CreateFcn(hObject, ~, ~)
% hObject    handle to value_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4_set_all.
function pushbutton4_set_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_set_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

curtbl = get(handles.gain_table,'Data');
value = str2double(get(handles.value_set,'String'));

column = get(handles.set_selection,'value');

if column <4,
    value = (value>0);
    set(handles.value_set,'string',sprintf('%d',value));
else
    value = mod(value,32);
    set(handles.value_set,'string',sprintf('%d',value));
end

curtbl(:,column) = value;
set(handles.gain_table,'Data',curtbl);

Cfg.data(:,column) = value;
setappdata(CurFig,'Cfg',Cfg);


% --- Executes on button press in pushbutton5_prepare_408.
function w = pushbutton5_prepare_408_Callback(~, ~, ~)
% hObject    handle to pushbutton5_prepare_408 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

d = Cfg.data;

DEBUG_PRINT = 0;


for a=1:Cfg.NA
    
    offset = (a-1)*(Cfg.NC+2);
    
    ir = (1+(a-1)*Cfg.NC):(a*Cfg.NC);
    
    amon = (d(ir,1)>0); % only zero or one allowed
    vth = (d(ir,2)>0); % only zero or one allowed
    tp = (d(ir,3)>0); %#ok<*NASGU> % only zero or one allowed
    gains = d(ir,4);
    for k=1:32,
        if k<=16,
            w1(2*k-1) = vth(k); %#ok<*AGROW>
            w1(2*k) = amon(k);
        else
            w2(2*k-1-32) = vth(k);
            w2(2*k-32) = amon(k);
        end
    end
    if DEBUG_PRINT, fprintf('### asic %d\n',a), end %#ok<*UNRCH>
        
    sw1 = (sprintf('%1d',fliplr(w1)));
    curw1 = uint32(bin2dec(sw1));
    if DEBUG_PRINT, fprintf('%s %08x\n',sw1, curw1), end
    w(offset+1) = curw1;
    
    sw2 = (sprintf('%1d',fliplr(w2)));
    curw2 = uint32(bin2dec(sw2));
    if DEBUG_PRINT, disp(sprintf('%s %08x',sw2, curw2)), end
    w(offset+2) = curw2;

    for c=1:Cfg.NC
        bits = zeros(1,32);
        curgain = gains(c);
        if curgain>0 && curgain<32,
            bits(32-(1:gains(c)))=1;
        end
        if tp(c), bits(32)=1; end
        sw = (sprintf('%1d',bits));
        curw = uint32(bin2dec(sw));
        if DEBUG_PRINT, disp(sprintf('%s %08x %d',sw, curw, curw)), end
        w(offset+2+Cfg.NC-c+1) = curw; % fix 8/25/15 for proper channel order
    end
end
tmp = [w(205:408),w(1:204)]; w=tmp; % fix 8/19/15 for proper asic order
% tmp = w(1:34); w = [tmp, w]; % test one extra 8/20/15
% Cfg.data
% reshape(w,[],12)
setappdata(CurFig,'spi408words',w);

%%
% interface function
function w = fun_PET_spi_editor()
% disp('This is PET gain control function')
CurFig = getappdata(0,'hPET_spi_editor');
w=getappdata(CurFig,'spi408words');


% --- Executes on selection change in set_selection.
function set_selection_Callback(~, ~, ~)
% hObject    handle to set_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns set_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from set_selection


% --- Executes during object creation, after setting all properties.
function set_selection_CreateFcn(hObject, ~, ~)
% hObject    handle to set_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6_add.
function pushbutton6_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

curtbl = get(handles.gain_table,'Data');
value = str2double(get(handles.value_set,'String'));
anum = get(handles.asic_number,'Value');

column = get(handles.set_selection,'value');
value = curtbl(:,column) + value;

if column <4,
    value = (value>0);
else
    ix = find(value<0);
    if ~isempty(ix), value(ix) = 0; end
    ix = find(value>31);
    if ~isempty(ix), value(ix) = 31; end
end

curtbl(:,column) = value;
set(handles.gain_table,'Data',curtbl);

Cfg.data( (1+(anum-1)*Cfg.NC):(anum*Cfg.NC),: ) = curtbl;
setappdata(CurFig,'Cfg',Cfg);

% --- Executes on button press in pushbutton7_add_all.
function pushbutton7_add_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPET_spi_editor');
Cfg = getappdata(CurFig,'Cfg');

curtbl = get(handles.gain_table,'Data');
anum = get(handles.asic_number,'Value');
value = str2double(get(handles.value_set,'String'));
column = get(handles.set_selection,'value');

value = Cfg.data(:,column) + value;

if column <4,
    value = (value>0);
else
    ix = find(value<0);
    if ~isempty(ix), value(ix) = 0; end
    ix = find(value>31);
    if ~isempty(ix), value(ix) = 31; end
end

Cfg.data(:,column) = value;

curtbl(:,column) = Cfg.data( (1+(anum-1)*Cfg.NC):(anum*Cfg.NC), column );

set(handles.gain_table,'Data',curtbl);
setappdata(CurFig,'Cfg',Cfg);
