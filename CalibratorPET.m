function varargout = CalibratorPET(varargin)
% CALIBRATORPET M-file for CalibratorPET.fig
%      CALIBRATORPET, by itself, creates a new CALIBRATORPET or raises the existing
%      singleton
%
%      H = CALIBRATORPET returns the handle to a new CALIBRATORPET or the handle to
%      the existing singleton*.
%
%      CALIBRATORPET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATORPET.M with the given input arguments.
%
%      CALIBRATORPET('Property','Value',...) creates a new CALIBRATORPET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CalibratorPET_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CalibratorPET_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CalibratorPET

% Last Modified by GUIDE v2.5 20-Jan-2017 15:53:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CalibratorPET_OpeningFcn, ...
    'gui_OutputFcn',  @CalibratorPET_OutputFcn, ...
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
set_default_figure_properties()

% --- Executes just before CalibratorPET is made visible.
function CalibratorPET_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CalibratorPET (see VARARGIN)

% Choose default command line output for CalibratorPET
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CalibratorPET wait for user response (see UIRESUME)
% uiwait(handles.CalibratorPET);
CurFig = hObject; pause(0.1);
setappdata(0,'hPETmonitor_new',CurFig);

initalize_device_names_selection(CurFig,handles)

Cfg = initialize_configuration(CurFig,handles,[]);
setappdata(CurFig,'Cfg',Cfg);

%Create tab group
tmp=version; matver = str2num(tmp(1));
if matver<8,
    handles.tgroup = uitabgroup('v0','Parent', hObject,'TabLocation', 'top');
else
    handles.tgroup = uitabgroup('Parent', hObject,'TabLocation', 'top');
end
pos = get(handles.P1,'position');
uni = get(handles.P1,'units');
set(handles.tgroup,'units',uni)
set(handles.tgroup,'position',pos);
% set(handles.tgroup,'BackgroundColor',[0 0 0]);
tgrpos = get(handles.tgroup,'position');
if matver<8,
    handles.tab1 = uitab('v0','Parent', handles.tgroup, 'Title', 'Display');
    handles.tab2 = uitab('v0','Parent', handles.tgroup, 'Title', 'Image');
    handles.tab3 = uitab('v0','Parent', handles.tgroup, 'Title', 'Control');
    handles.tab4 = uitab('v0','Parent', handles.tgroup, 'Title', 'Batch');
    handles.tab5 = uitab('v0','Parent', handles.tgroup, 'Title', 'Calibr');
    handles.tab6 = uitab('v0','Parent', handles.tgroup, 'Title', 'IO');
else
    handles.tab1 = uitab('Parent', handles.tgroup, 'Title', 'Display');
    handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Image');
    handles.tab3 = uitab('Parent', handles.tgroup, 'Title', 'Control');
    handles.tab4 = uitab('Parent', handles.tgroup, 'Title', 'Batch');
    handles.tab5 = uitab('Parent', handles.tgroup, 'Title', 'Calibr');
    handles.tab6 = uitab('Parent', handles.tgroup, 'Title', 'IO');
end
%Place panels into each tab
set(handles.P1,'Parent',handles.tab1)
set(handles.P2,'Parent',handles.tab2)
set(handles.P3,'Parent',handles.tab3)
set(handles.P4,'Parent',handles.tab4)
set(handles.P5,'Parent',handles.tab5)
set(handles.P6,'Parent',handles.tab6)
%Reposition each panel to same location as panel 1
corpos = [tgrpos(1) tgrpos(2) 0 0];
set(handles.P1,'position',get(handles.P1,'position')-corpos);
set(handles.P2,'position',get(handles.P1,'position'));
set(handles.P3,'position',get(handles.P1,'position'));
set(handles.P4,'position',get(handles.P1,'position'));
set(handles.P5,'position',get(handles.P1,'position'));
set(handles.P6,'position',get(handles.P1,'position'));

set(handles.calibr_root,'string',pwd)
set(handles.batch_root,'string',pwd)

% set( findall(CurFig, '-property', 'Units' ), 'Units', 'Normalized' )
set( CurFig, 'Units', 'Normalized' )
set( CurFig, 'Resize', 'on' )
h = findobj( CurFig, '-property', 'Units' );
set( h, 'Units', 'Normalized' )

%%
% Initialize device types according to description in PETdevice.cfg
function initalize_device_names_selection(CurFig,handles)

for k=1:10,
    [device_name,nring,nblk,nasic,bmap,cmap] = read_device_type_address_map(k);
    if ~isempty(device_name),
        dnms{k} = device_name;
    end
end
set(handles.device_num,'string',dnms)
drawnow

%%
% initialize configuration
function Cfg = initialize_configuration(CurFig, handles, device_num)
if ~isempty(device_num)
    disp(sprintf('Device %d is selected', device_num))
end

Cfg = disp_init_info();
Cfg.defaults = read_petmonitor_defaults(device_num);
Cfg.coinc_times = Cfg.defaults.Coincidence_time_window(1):1:Cfg.defaults.Coincidence_time_window(2);
Cfg.start = 0;
Cfg.readmode = 1;
Cfg.eventfile_r = 'petmon.evt';
Cfg.coincfile_r = 'petmon.coinc';
Cfg.writemode = 1;
Cfg.eventfile_w = 'petmon.evt';
Cfg.coincfile_w = 'petmon.coinc';
Cfg.gaindata=[];
Cfg.tspm.delay = 0.2;
% Cfg.tspm.delay = determine_tspm_delay(0.2);
[nax, nco] = initalize_image_views(Cfg, handles,1);
Cfg.image.NumberSlices = [nax, nco, nco];
Cfg.sindex = CalculateSinogramIndex(Cfg);
Cfg.effcorr = read_efficiency_corrections(Cfg,'');

setappdata(CurFig,'Cfg',Cfg);
setappdata(CurFig,'Calibr',[]);

%%
% initialize image view panel
function [Nax, Nco] = initalize_image_views(Cfg, handles, nentry)

Nax = 2*Cfg.defaults.Device.nring*Cfg.defaults.Device.bmap(1)-1; % Number of axial slices
N = Cfg.defaults.Device.nblk*Cfg.defaults.Device.nasic*Cfg.defaults.Device.bmap(2);

switch nentry
    case 1
        Nco = 2*floor(N/(2*sqrt(2))); % based on iradon description
    case 2
        Nco = round(str2double(get(handles.image_zoom,'string')));
    otherwise
        disp('Error: unrecognized choise for image initialization'), return;
end
Nsa = Nco;

first_time = ~iscell(get(handles.image_axial_slice,'string'));

if first_time,
    cur_a_s = ceil(Nax/2);
    cur_c_s = ceil(Nco/2);
    cur_s_s = ceil(Nsa/2);
    cur_a_onoff = 1;
else
    cur_a_s = get(handles.image_axial_slice, 'value');
    cur_c_s = get(handles.image_coronal_slice, 'value');
    cur_s_s = get(handles.image_sagital_slice, 'value');
    cur_a_onoff = get(handles.control_asic_onoff, 'value');
end

tmp=[]; for k=1:Nax, tmp{k}=sprintf('%2d',k); end
set(handles.image_axial_slice,'string',tmp);
tmp=[]; for k=1:Nco, tmp{k}=sprintf('%2d',k); end
set(handles.image_coronal_slice,'string',tmp);
tmp=[]; for k=1:Nsa, tmp{k}=sprintf('%2d',k); end
set(handles.image_sagital_slice,'string',tmp);

set(handles.image_view,'value',1);
set(handles.image_axial_slice,'enable','on','value', cur_a_s);
set(handles.image_coronal_slice,'enable','off','value', cur_c_s);
set(handles.image_sagital_slice,'enable','off','value', cur_s_s);

NA = Cfg.defaults.Device.nring*Cfg.defaults.Device.nblk*Cfg.defaults.Device.nasic;
tmp=[]; for k=1:NA, tmp{k}=sprintf('%2d',k); end
set(handles.control_asic_onoff,'string',tmp,'value',cur_a_onoff);

set(handles.image_zoom,'string',sprintf('%2d',Nco));

%%
% Detrmine TSPM connectivity delay
function tdly = determine_tspm_delay(recommended_delay)
fprintf('TSPM command delay is set to ');
tdly=0; rtn=-1; cnt=0;
while rtn~=256 & cnt < 10
    a = tic;
    TSPMcontrol('write',0,256);
    rtn = TSPMcontrol('read',0);
    if isempty(rtn), rtn=-1; end
    tdly = toc(a);
    pause(0.1);
    cnt=cnt+1;
end
if tdly<recommended_delay, tdly=recommended_delay; end
fprintf('%5.3f seconds\n',tdly);

%%
% Read efficiency corrections
function effcorr = read_efficiency_corrections(Cfg,file2load)

XX = Cfg.defaults.Device.bmap;
DZ = XX(1); % detectors along Z axis in each asic
DT = XX(2); % detectors along circumference in each asic

NC = prod(Cfg.defaults.Device.bmap);
NA = Cfg.defaults.Device.nring*Cfg.defaults.Device.nblk*Cfg.defaults.Device.nasic;
NTALL = DT*NA; % total number of detectors in single z virtual ring
NZALL = DZ*Cfg.defaults.Device.nring; % total number of detectors in z directions

effcorr = double(ones(NTALL,NTALL,2*NZALL-1));
try
    load(file2load);
    for nfi1=0:NTALL-1;
        for nz1=0:NZALL-1;
            for nfi2=0:NTALL-1;
                for nz2=0:NZALL-1;
                    iringdiff=nz2-nz1;
                    if abs(iringdiff)<=1
                        n1=1+nfi1;
                        n2=1+mod(nfi1-nfi2,NTALL);
                        n3=2*(nz1+1)+iringdiff-1;
                        effco_sgm(n1,n2,n3)=eHM(nz1+1,nfi1+1)*eHM(nz1+1,nfi2+1);
                    end
                end
            end
        end
    end
    effcorr = effco_sgm;
catch
    disp('Efficiency corrections were not loaded')
end


% --- Outputs from this function are returned to the command line.
function varargout = CalibratorPET_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function rtn = clean_up_network_conn()
rtn=0;
tmp = instrfind('type','udp');
% tmp = instrfind('type','udp','localport',32003);
if ~isempty(tmp),
    try fclose(tmp); catch, disp('Error closing udp'), end
    try delete(tmp); catch, disp('Error deleting udp'), end
end

% --- Executes on button press in pushbutton1_start.
function pushbutton1_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

% disp(sprintf('>>> start entry with %s',Cfg.eventfile_w))
Cfg.start = ~Cfg.start;
if Cfg.start, set(hObject,'String','Stop')
else, set(hObject,'String','Start'), end

[nax, nco] = initalize_image_views(Cfg, handles,2);
Cfg.image.NumberSlices = [nax, nco, nco];
setappdata(CurFig,'Cfg',Cfg);
drawnow;

MACHINE_FORMAT = 'ieee-be';
DEFAULT_READ_PACKET_SIZE = 520; % in 2 bytes words; UDP will read 520
DEFAULT_READ_BUFFER_SIZE = Cfg.defaults.UDP_buffer_size;
DEFAULT_WRITE_BUFFER_SIZE = DEFAULT_READ_PACKET_SIZE*100;
NumAngularProjections = Cfg.defaults.Device.nblk*Cfg.defaults.Device.nasic*Cfg.defaults.Device.bmap(2);

% disp(sprintf('### start %d ',Cfg.start))
if ~Cfg.start, return; end

NumAverages = str2double(get(handles.averages_indicator,'String'));
readmode = get(handles.readmode,'Value');
monfile_r = get(handles.monfile_r,'string');
writemode = get(handles.writemode,'Value');
monfile_w = get(handles.monfile_w,'string');
ShowMode = get(handles.showplot,'Value');
MonitorAxes = handles.monitoraxes;
NetworkBufferAxes = handles.axesbuffer;

if strcmp(monfile_r,monfile_w) & length(monfile_r)>1 & length(monfile_w)>1
    errordlg('Input and output files must not be the same')
    return
end

WorkFlow.PacketHeaders.Calc = 1;
WorkFlow.PacketHeaders.Show = 0;
WorkFlow.EventMatrix.Calc = 1;
WorkFlow.Display.Show = ShowMode;
WorkFlow.PacketIntegrity.Calc = 1;
WorkFlow.PacketsNumber2Read = str2double(get(handles.packets_indicator,'String'));
WorkFlow.ProcessPacketHistory = 1;
WorkFlow.ProcessEventHistory = 1;
WorkFlow.ProcessCoinTimeHistogram = 1;
WorkFlow.UpdateIndicators = 1;
WorkFlow.DisplayMode = 1;

coinctimevector = Cfg.defaults.Coincidence_time_window(1):2:Cfg.defaults.Coincidence_time_window(2); % in ns
coinc_in_index = find(coinctimevector>Cfg.defaults.In_time_window(1) & coinctimevector<Cfg.defaults.In_time_window(2));
coinc_out_index = find(coinctimevector>Cfg.defaults.Out_time_window(1) & coinctimevector<Cfg.defaults.Out_time_window(2));

eventtimevector = 0:(Cfg.defaults.Event_show_time/1000):Cfg.defaults.Event_show_time;
historytimemax = Cfg.defaults.History_show_time*1e-9; % in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = instrfind('type','udp','localport',32003);
if ~isempty(tmp),
    try fclose(tmp); catch, disp('Error closing high speed port'), end
    try delete(tmp); catch, disp('Error deleting high speed port'), end
end

% try
switch readmode
    case 1 % Monitor only
        WorkFlow.ReadNetwork = 1;
        WorkFlow.ReadEventFile = 0;
        WorkFlow.ReadCoincFile = 0;
        WorkFlow.NetworkBuffer.Show = 1;
    case 2 % Read event file
        WorkFlow.ReadNetwork = 0;
        WorkFlow.ReadEventFile = 1;
        WorkFlow.ReadCoincFile = 0;
        WorkFlow.NetworkBuffer.Show = 0;
    case 3 % Read coinc file
        WorkFlow.ReadNetwork = 0;
        WorkFlow.ReadEventFile = 0;
        WorkFlow.ReadCoincFile = 1;
        WorkFlow.NetworkBuffer.Show = 0;
    otherwise
        errordlg('Error: unrecognized acquizition mode')
        return
end

switch writemode
    case 1 % Nothing to write
        WorkFlow.WriteEventFile = 0;
        WorkFlow.WriteCoincFile = 0;
    case 2 % Write event file
        WorkFlow.WriteEventFile = 1;
        WorkFlow.WriteCoincFile = 0;
    case 3 % Write coinc file
        WorkFlow.WriteEventFile = 0;
        WorkFlow.WriteCoincFile = 1;
    otherwise
        errordlg('Error: unrecognized output mode')
        return
end

if WorkFlow.WriteEventFile
    if isempty(Cfg.eventfile_w),
        uiwait(warndlg({'Event file is not valid','Disabling output'})),
        WorkFlow.WriteEventFile = 0;
        set(handles.writemode,'Value',1);
        set(handles.monfile_w,'string','')
    else
        monfile_w = get(handles.monfile_w,'string');
        [a b c] = fileparts(Cfg.eventfile_w);
        int_monfile_w = [b,c]; %clear a b c
        if ~strcmp(int_monfile_w, monfile_w)
%             uiwait(warndlg({'Output file naming discrepancy','Disabling output'}))
            WorkFlow.WriteEventFile = 0;
            set(handles.writemode,'Value',1);
            set(handles.monfile_w,'string','')
        else
            [fidw, msg] = fopen(Cfg.eventfile_w,'w',MACHINE_FORMAT);
            if ~isempty(msg), disp(sprintf('Error: %s',msg)), return; end
            DataStream=uint16([]);
%             disp(['Opened ',Cfg.eventfile_w,' to write data']);
        end
    end
end

if WorkFlow.WriteCoincFile
    uiwait(errordlg('Not developed function'))
    return;
%     if isempty(Cfg.coincfile_w),
%         uiwait(errordlg({'Coinc file is not valid','Changing to no output'})),
%         WorkFlow.WriteCoincFile = 0;
%         set(handles.writemode,'Value',1);
%     else
%         [fidw, msg] = fopen(Cfg.coincfile_w,'w',MACHINE_FORMAT);
%         if ~isempty(msg), disp(sprintf('Error: %s',msg)), return; end
%         CoincStream=uint16([]);
%         disp(['Opened ',Cfg.coincfile_w,' to write data']);
%     end
end

if WorkFlow.ReadEventFile
    curfilename = fullfile(Cfg.eventfile_r);
    [fidr, msg] = fopen(curfilename,'r',MACHINE_FORMAT);
    if ~isempty(msg), disp(sprintf('Error: %s',msg)), return; end
    curfilesize = subsref(dir(curfilename), substruct('.','bytes'));
end

if WorkFlow.ReadCoincFile
    curfilename = fullfile(Cfg.coincfile_r);
    [fidr, msg] = fopen(curfilename,'r',MACHINE_FORMAT);
    if ~isempty(msg), disp(sprintf('Error: %s',msg)), return; end
    curfilesize = subsref(dir(curfilename), substruct('.','bytes'));
end

if WorkFlow.ReadNetwork
    clean_up_network_conn();
    u = udp('192.168.120.1',...
        'localport', Cfg.defaults.UDP_port,...
        'ByteOrder', 'bigEndian',...
        'Timeout',0,...
        'Terminator','',...
        'InputBufferSize', DEFAULT_READ_BUFFER_SIZE);
    fopen(u);
    set(u,'DatagramTerminateMode','on')
    flushinput(u);
    flushoutput(u);
end

% WorkFlow.ProcessPacketHistory
PacketHistory = [];
% WorkFlow.ProcessEventHistory
EventHistory = [];
% WorkFlow.ProcessCoinTimeHistogram
dt=NaN;
detime=[];
detime_hist=[];
ldt=0;
% WorkFlow.EventMatrix.Calc
evr = InitializeEventRate(Cfg); evr_hist= evr;
sgr=0; zgr=0; t0corr=0;

startTime = tic;
proc_rate = [];
singleReadTimer = tic; % just once here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while Cfg.start,
    
    %%
    % communicate with gui
    ShowMode = get(handles.showplot,'Value');
    if WorkFlow.Display.Show ~= ShowMode,
        history_e = [];
        history_p = [];
        WorkFlow.Display.Show = ShowMode;
    end
    NumAverages = str2double(get(handles.averages_indicator,'String'));
    WorkFlow.PacketsNumber2Read = str2double(get(handles.packets_indicator,'String'));
    DisplayMode = get(handles.display_mode,'Value');
    Cfg = getappdata(CurFig,'Cfg');
    
    if WorkFlow.DisplayMode~=DisplayMode
        WorkFlow.DisplayMode=DisplayMode;
        switch DisplayMode
            case 1
                WorkFlow.NetworkBuffer.Show = 1;
                WorkFlow.PacketIntegrity.Calc = 1;
                WorkFlow.EventMatrix.Calc = 1;
                WorkFlow.ProcessPacketHistory = 1;
                WorkFlow.ProcessEventHistory = 1;
                WorkFlow.ProcessCoinTimeHistogram = 1;
                WorkFlow.UpdateIndicators = 1;
            case 2
                WorkFlow.NetworkBuffer.Show = 1;
                WorkFlow.PacketIntegrity.Calc = 1;
                WorkFlow.EventMatrix.Calc = 1;
                WorkFlow.ProcessPacketHistory = 0;
                WorkFlow.ProcessEventHistory = 0;
                WorkFlow.ProcessCoinTimeHistogram = 1;
                WorkFlow.UpdateIndicators = 1;
            case 3
                WorkFlow.NetworkBuffer.Show = 0;
                WorkFlow.PacketIntegrity.Calc = 0;
                WorkFlow.EventMatrix.Calc = 0;
                WorkFlow.ProcessPacketHistory = 0;
                WorkFlow.ProcessEventHistory = 0;
                WorkFlow.ProcessCoinTimeHistogram = 0;
                WorkFlow.UpdateIndicators = 1;
            case 4
                WorkFlow.NetworkBuffer.Show = 0;
                WorkFlow.PacketIntegrity.Calc = 0;
                WorkFlow.EventMatrix.Calc = 0;
                WorkFlow.ProcessPacketHistory = 0;
                WorkFlow.ProcessEventHistory = 0;
                WorkFlow.ProcessCoinTimeHistogram = 0;
                WorkFlow.UpdateIndicators = 0;
        end
    end
    
    Cfg = getappdata(CurFig,'Cfg');
    sfigure(CurFig);
    drawnow;
    
    %%
    % Acquire data
    x=uint16([]); y=uint16([]); DataLength=0;
    for pck_cnt=1:WorkFlow.PacketsNumber2Read
        if WorkFlow.ReadNetwork
            if u.BytesAvailable>DEFAULT_READ_PACKET_SIZE
                [cx, DataLength, ~, ~, ~] = fread(u,DEFAULT_READ_PACKET_SIZE,'uint16');
            else
                break
            end
        end
        if WorkFlow.WriteEventFile & DataLength>0
            DataStream = [DataStream;cx(1:end-8)];
            % minus 8 bytes to comply with BNL coincidence processor
            % expectations, also assuming full packet length
        end
        if WorkFlow.ReadEventFile
            [cx, DataLength] = fread(fidr, DEFAULT_READ_PACKET_SIZE,'uint16',MACHINE_FORMAT);
            curfilepos = ftell(fidr);
            tmp =sprintf('%s: %5.2f%%',datestr(now,13),100*curfilepos/curfilesize);
            set(handles.io_file_status,'string',tmp);
            if feof(fidr), Cfg.start=0; setappdata(CurFig,'Cfg',Cfg); disp('EOF'); 
%                 save('PETmonitor_dump.mat'), 
                break, end
        end
        
        if DataLength,
            [pix, qix] = find_packets(cx); % get packet starts
            npix = length(pix);
            if npix % at least one packet header found
                if pix(end) > DataLength-7,
                    disp('Partial header detected: skipping')
                    %                     pix
                    break
                end
                for k=1:npix
                    if mod(pix(k)-1,4)
                        disp('Byte loss detected: skipping')
                        %                         pix, DataLength, y, cx(1:16)
                        break
                    end
                    if k==1 & pix(1)>1 % see if there are preceeding events from previous packet
                        x=[x;cx(1:pix(1)-1)];
                    end
                    %                     if k<npix  % events in packet
                    %                         x=[x;cx(pix(k)+8:pix(k+1)-1)];
                    %                     elseif k==npix
                    %                         x=[x;cx(pix(k)+8:end)];
                    %                     end
                    y = [y; cx(pix(k):pix(k)+7)]; % get first 8 words for packet info
                    x=[x; cx(pix(k)+8:pix(k)+qix(k)-1)];
                    
                end
                %             print_packet_content(cx);
            else
                x=[x;cx];
            end
        else
            break; % stop acquisition and go to display if there is nothing in the stream
        end
    end
    
    if WorkFlow.WriteEventFile
        lz = length(DataStream);
        %         disp(sprintf('Write buffered %d words out of %d', lz,DEFAULT_WRITE_BUFFER_SIZE))
        if Cfg.start==0 | lz>=DEFAULT_WRITE_BUFFER_SIZE
            cnt = fwrite(fidw,DataStream,'uint16',0,MACHINE_FORMAT);
            fprintf('Wrote to file %d out of %d words\n',cnt, lz); pause(0.1)
            DataStream=uint16([]);
        end
    end
    
    
    DataLength=length(x);
    
    if DataLength==0
        %         disp(sprintf('%s: no data',datestr(now,31)))
    elseif mod(DataLength,4)
        disp('### Sqewed event data. Something wrong ###')
        DataLength
        return
    else
        if WorkFlow.NetworkBuffer.Show
            tmp = u.BytesAvailable/DEFAULT_READ_BUFFER_SIZE;
            if tmp>1, tmp=1; end
            bar(NetworkBufferAxes, tmp), axis tight
            set(NetworkBufferAxes,'xlim',[0.5,1.5],'ylim',[0 1], 'xtick',[],'ytick',[])
            %             set(NetworkBufferAxes,'ylim',[0 1], 'xtick',[], 'ytick',[])
        end
        
        [psn, era] = extract_packet_headers(y);
        [etime, easic, echan] = extract_events(x);
        
        NumberEvents = length(etime);
        NumberPackets = length(psn);
        
        if WorkFlow.PacketIntegrity.Calc & NumberPackets
            pchk = (diff(psn)~=1);
            %             fprintf('### %s %05.1f events in %3d paskets with %3d dropped\n',...
            %                 datestr(now,13),NumberEvents/NumberPackets,NumberPackets,psn(end)-psn(1)+1-NumberPackets),
        end
        
        if WorkFlow.EventMatrix.Calc
            [NC, NA] = size(evr);
            evr = zeros([NC, NA]);
            for k=1:NumberEvents
                na = 1+easic(k);
                nc = 1+echan(k);
                if na <= NA & nc <= NC,
                    evr(nc,na) = evr(nc,na)+1;
                else
                    disp(sprintf('Error calc with asic %d and channel %d length of data %d',na, nc, DataLength))
                    %                     return
                end
            end
            evr_hist = (evr_hist*(NumAverages-1) + evr)/NumAverages;
            
%             nevr = numel(evr(:));
%             me = mean(evr(:)); ms = std(evr(:));
%             meh = mean(evr_hist(:)); msh = std(evr_hist(:));
%             
%             dead_elem_hist = numel(find(evr_hist(:)<(meh-msh)))/nevr*100;
%             dead_elem = numel(find(evr(:)<1))/nevr*100;
%             disp(sprintf('Dead elements %2.0f %2.0f',dead_elem,dead_elem_hist))
            
        end
        
        if WorkFlow.ProcessPacketHistory & NumberPackets
            curtime = toc(startTime); lpsn = double(psn(end)-psn(1)+1);
            PacketHistory = [PacketHistory; [curtime, [100*(lpsn-NumberPackets), NumberEvents]/lpsn]];
            elim = find(curtime-PacketHistory(:,1)>historytimemax);
            if ~isempty(elim), PacketHistory(elim,:)=[]; end
        end
        
        if WorkFlow.ProcessEventHistory & NumberEvents
            EventHistory = double([EventHistory; etime]);
            mih = min(EventHistory);
            mah = max(EventHistory);
            if WorkFlow.Display.Show == 5
                est = max(eventtimevector);
            elseif WorkFlow.Display.Show == 8
                est = historytimemax*1e9;
            else
                est=0;
            end
            if est
                if (mah-mih)> est
                    elim = find(EventHistory<(mah-est));
                    EventHistory(elim)=[];
                end
            end
        end
        
        if WorkFlow.ProcessCoinTimeHistogram
            if ~isempty(etime)
                detime=diff(double(etime));
                ldt = length(detime);
                if NumAverages>1 & ~isempty(detime_hist),
                    if ldt>1
                        try
                            detime_hist = (detime_hist*(NumAverages-1) + histc(double(detime),coinctimevector))/NumAverages;
                            ldt_hist = (ldt_hist*(NumAverages-1) + ldt)/NumAverages;
                        catch
                            disp('Coinc historgram callc problem')
                        end
                    end
                else
                    detime_hist = histc(double(detime),coinctimevector);
                    ldt_hist = ldt;
                end
            end
            
        end
        
        if WorkFlow.WriteCoincFile
            %             if DataLength,
            %                 CoincStream = [CoincStream; x];
            %             end
            if ~WorkFlow.ProcessCoinTimeHistogram
                detime=diff(double(etime));
            end
            if ~isempty(detime)
                tmp = find(detime>Cfg.defaults.Coincidence_time_window(1) & detime<Cfg.defaults.Coincidence_time_window(2));
                tmp = [tmp,tmp+1]';
                icoinc = tmp(:)';
                if ~isempty(icoinc)
                    tmp = [x(icoinc*4-3),x(icoinc*4-2),x(icoinc*4-1),x(icoinc*4-0)]';
                    xcoinc = tmp(:);
                    CoincStream=[CoincStream;xcoinc];
                    lx = length(CoincStream);
                    %         disp(sprintf('Write buffered %d words out of %d', lz,DEFAULT_WRITE_BUFFER_SIZE))
                    if Cfg.start==0 | lx>=DEFAULT_WRITE_BUFFER_SIZE
                        cnt = fwrite(fidw,CoincStream,'uint16',0,MACHINE_FORMAT);
                        fprintf('Wrote to file %d out of %d coinc words\n',cnt, lx);
                        %                         reshape(CoincStream,4,[])'
                        CoincStream=uint16([]);
                    end
                end
            end
        end
        
        
        %%
        % Show graphics
        if WorkFlow.Display.Show
            if ishandle(CurFig), set(0, 'CurrentFigure', CurFig);
            else CurFig = figure(CurFig); end
            set(CurFig,'CurrentAxes',MonitorAxes),
            switch WorkFlow.Display.Show
                case 1 % Events 2D map
                    if WorkFlow.EventMatrix.Calc
                        imagesc(evr_hist)
                        set(MonitorAxes, 'ydir', 'normal')
                        axis tight
                        xlabel('Asic'), ylabel('Channel');
                    end
                case 2 % Events Asic bar
                    if WorkFlow.EventMatrix.Calc
                        bar(sum(evr_hist))
                        axis tight
                        xlabel('Asic'), ylabel('Counts');
                    end
                case 3 % Events Channel bar
                    if WorkFlow.EventMatrix.Calc
                        bar(sum(evr_hist'))
                        axis tight
                        xlabel('Channel'), ylabel('Counts');
                    end
                case 4 % Coincidence time
                    if WorkFlow.ProcessCoinTimeHistogram
                        if ~isempty(detime_hist)
                            plot(coinctimevector,detime_hist);
                            axis tight
                            xlabel('Time (ns)')
                        end
                    end
                case 5 % Events density
                    if WorkFlow.ProcessEventHistory
                        if ~isempty(EventHistory)
                            ets = max(eventtimevector)*1e-9;
                            evdens = histc(EventHistory-EventHistory(1),eventtimevector);
                            edstat(1) = length(find(evdens==1))/ets;
                            edstat(2) = length(find(evdens==2))/ets;
                            edstat(3) = length(find(evdens>2))/ets;
                            disp(sprintf('%6.2f  \t%6.2f  \t%6.2f',...
                                edstat(1), edstat(2), edstat(3)));
                            plot(eventtimevector/1e6, evdens),
                            set(MonitorAxes,'ylim',[0 4]), grid
                            xlabel('Time (ms)');
                        end
                    end
                case 6 % Packets drop out
                    if WorkFlow.ProcessPacketHistory
                        if ~isempty(PacketHistory)
                            plot(PacketHistory(:,1)-PacketHistory(1,1),PacketHistory(:,2))
                            axis tight
                            xlabel('Time (s)');
                        end
                    end
                case 7 % Events per packet
                    if WorkFlow.ProcessPacketHistory
                        if ~isempty(PacketHistory)
                            plot(PacketHistory(:,1)-PacketHistory(1,1),PacketHistory(:,3))
                            axis tight
                            xlabel('Time (s)');
                        end
                    end
                case 8 % Events times
                    if WorkFlow.ProcessEventHistory
                        if ~isempty(EventHistory)
                            plot((EventHistory-EventHistory(1))/1e9), view(-90, 90),
                            set(MonitorAxes, 'ydir', 'reverse'), axis tight, grid
                            ylabel('Time (s)');
                        end
                    end
                case 9 % Sinogram
                    slice = get(handles.image_axial_slice,'value');
                    [sgr, zgr] = CalculateSinogram(Cfg, etime, easic, echan, sgr, zgr, NumAverages);
                    if numel(sgr(:))>1, imagesc(sgr(:,:,slice)); end
                case 10 % Zinogram
                    [sgr, zgr] = CalculateSinogram(Cfg, etime, easic, echan, sgr, zgr, NumAverages);
                    if numel(zgr(:))>1, imagesc(zgr); end
                case 11 % Image
                    img3d=[];
                    image_view = get(handles.image_view,'value');
                    image_reconstruction = get(handles.image_reconstruction,'value');
                    %                     sum(Cfg.effcorr(:))
                    [sgr, zgr] = CalculateSinogram(Cfg, etime, easic, echan, sgr, zgr, NumAverages);
                    dtheta = 360/NumAngularProjections;
                    if numel(sgr(:))>1,
                        try
                            for k=1:Cfg.image.NumberSlices(1),
                                switch image_reconstruction
                                    case 1, % inverse radon
                                        img3d(:,:,k) = iradon(sgr(:,:,k)',dtheta,'linear','Cosine',Cfg.image.NumberSlices(2));
                                    case 2, % simple back projection
                                        img3d(:,:,k) = bpg_unfiltered(sgr(:,:,k)',dtheta);
                                    case 3, % filtered back projection in spatial domain
                                        img3d(:,:,k) = bpg_filtered_spatial(sgr(:,:,k)',dtheta);
                                end
                            end
                        catch
                            disp('Error: backward transformation is not ready')
                        end
                    end
                    if ~isempty(img3d)
                        
                        [Nc,Ns,Na]=size(img3d);
                        switch image_view
                            case 1 % axial view
                                slice = get(handles.image_axial_slice,'value');
                                if slice>Na, img = []; disp('Error: axial slice is out of range')
                                else, img = img3d(:,:,slice); end
                            case 2 % coronal view
                                slice = get(handles.image_coronal_slice,'value');
                                if slice>Nc, img = []; disp('Error: coronal slice is out of range')
                                else, img = shiftdim(img3d(slice,:,:),1); end
                            case 3 % sagital view
                                slice = get(handles.image_sagital_slice,'value');
                                if slice>Nc, img = []; disp('Error: sagital slice is out of range')
                                else, img = shiftdim(img3d(:,slice,:),2); end
                            otherwise
                                img=[];
                        end
                        
                        % 'linear','Ram-Lak');
                        % 'linear','Shepp-Logan');
                        % 'linear','Cosine');
                        % 'linear','Hamming');
                        % 'linear','Hann');
                        % 'linear','None');
                        
                        image_scaling = get(handles.image_scaling,'value');
                        switch image_scaling
                            case 1 % global
                                img_max = max(img3d(:));
                            case 2 % local
                                img_max = max(img(:));
                            case 3 % fixed
                                img_max = str2double(get(handles.image_max,'string'));
                        end
                        
                        if ~isempty(img) & img_max,
                            imagesc(img,[0 img_max]);
                            set(handles.image_status,'string',sprintf('Mean %6.4f (%6.4f) Max %6.4f',mean(img(:)),std(img(:)),max(img(:))));
                        end
                    end
                case 12 % t0 shift
                    [t0corr, tshift] = Calculate_t0_correction(Cfg, etime, easic, echan, t0corr, NumAverages, coinctimevector);
                    %                     mt0=reshape(shiftdim(mean(t0corr,[],1),1),12*32,12*32);
                    %                     imagesc(mt0), colorbar;
                    %                     mt0=reshape(t0corr,length(coinctimevector),12*32*12*32);
                    %                     imagesc(mt0), colorbar;
                    imagesc(tshift), colorbar;
                case 13 % t0 signal
                case 14 % t0 noise
                otherwise
                    ;
            end
            %%
            % update information boxes
            if WorkFlow.UpdateIndicators
                if ~isempty(era)
                    set(handles.event_rate,'String',sprintf('%d',era(end)));
                    
                end
                sret = toc(singleReadTimer);
                if ~isempty(proc_rate)
                    proc_rate = round((proc_rate*(NumAverages-1) + round(NumberEvents/sret))/NumAverages);
                else
                    proc_rate = round(NumberEvents/sret);
                end
                set(handles.proc_rate,'String',sprintf('%d',proc_rate));
                singleReadTimer = tic;
                if ~isempty(PacketHistory)
                    dra = mean(PacketHistory(:,2));
                    set(handles.drop_rate,'String',sprintf('%4.1f',dra));
                end
                if ~isempty(detime_hist)
                    incoinc = sum(detime_hist(coinc_in_index));
                    outcoinc = sum(detime_hist(coinc_out_index));
                    alcoinc = ldt_hist;
                    if alcoinc>0, cra = 100*(incoinc-outcoinc)/alcoinc;
                    else, cra=0; end
                    set(handles.coinc_rate,'String',sprintf('%4.2f',cra));
                end
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if WorkFlow.WriteEventFile,
    fclose(fidw);
end
if WorkFlow.WriteCoincFile,
    fclose(fidw);
end
if WorkFlow.ReadEventFile,
    fclose(fidr);
end
if WorkFlow.ReadCoincFile,
    fclose(fidr);
end
if WorkFlow.ReadNetwork,
    fclose(u);
end
set(hObject,'String','Start')

% catch
%     disp('Error reading from high speed port')
%     Cfg.start = 0;
%     setappdata(CurFig,'Cfg',Cfg);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save('PETmonitor_new_dump.mat')


% --- Executes on selection change in readmode.
function readmode_Callback(hObject, eventdata, handles)
% hObject    handle to readmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns readmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from readmode
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

if Cfg.start,
    tmp = questdlg({'Data acquisition is in progress.','Do you want to stop?'},'Warning');
    if strcmpi(tmp,'yes')
        Cfg.start = 0;
    else
        set(hObject,'Value',Cfg.readmode); drawnow;
        return
    end
end

tmp = get(hObject,'String');
if isempty(tmp), return, end
contents = cellstr(tmp);
selection = get(hObject,'Value');
str = lower(contents{selection});
if isempty(Cfg.folder), Cfg.folder = pwd; end
switch selection
    case 1 % Read from network
        Cfg.readmode = 1;
        set(handles.monfile_r,'string','')
    case 2 % Read event file
        [fname,paths] = uigetfile(fullfile(Cfg.folder,'*.evt'),'Specify input event file');
        Cfg.eventfile_r = fullfile(paths,fname);
        Cfg.readmode = 2;
        set(handles.monfile_r,'string',fname)
        Cfg.folder = paths;
    case 3 % Read coinc file
        [fname,paths] = uigetfile(fullfile(Cfg.folder,'*.coinc'),'Specify input coinc file');
        Cfg.coincfile_r = fullfile(paths,fname);
        Cfg.readmode = 3;
        set(handles.monfile_r,'string',fname)
        Cfg.folder = paths;
end
setappdata(CurFig,'Cfg',Cfg);

% --- Executes during object creation, after setting all properties.
function readmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to readmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in showplot.
function showplot_Callback(hObject, eventdata, handles)
% hObject    handle to showplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns showplot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from showplot


% --- Executes during object creation, after setting all properties.
function showplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function averages_indicator_Callback(hObject, eventdata, handles)
% hObject    handle to averages_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of averages_indicator as text
%        str2double(get(hObject,'String')) returns contents of averages_indicator as a double
val = str2double(get(hObject,'String'));
gas = log2(val)/10;
set(handles.averages_control,'Value',gas)
drawnow

% --- Executes during object creation, after setting all properties.
function averages_indicator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to averages_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% Read monitor configuration file
function ParamList = read_petmonitor_defaults(requested_device_num)

try
    [localpath, ~] =fileparts(mfilename('fullpath'));
    fprintf('Loading monitor configuration form %s\n',fullfile(localpath,'PETmonitor.cfg'))
    fid = fopen(fullfile(localpath,'PETmonitor.cfg'));
    if fid<0
        fprintf('Error: cannot open monitor configuration %s\n',fullfile(localpath,'PETmonitor.cfg'))
    end
    tline = fgets(fid);
    while ischar(tline)
        if ~strcmp(tline(1),'%')
            s = strrep(tline,sprintf('\n'),'');
            s(min(strfind(s,'%')):end)=[]; % remove comments
            p = char(regexp(s,'^\w*','match'));
            q = regexp(s,'-?\d+\.?\d*e?\+?-?\d*','match');
            for n=1:length(q),
                eval(['ParamList.',p,'(n) = ',char(q{n}),';']);
            end
            if strcmp(p, 'Device_type')
                devtype = str2double(char(q{1}));
                if ~isempty(requested_device_num)
                    devtype = requested_device_num;
                    ParamList.Device_type = requested_device_num;
                end
                [device_name,nring,nblk,nasic,bmap,cmap] = read_device_type_address_map(devtype);
                if isempty(cmap)
                    disp('Warning: Device type configuration file is missing');
                    cmap = [28 23 19 15 8 4 31 24 20 16 12 7 32 27 11 3 5 10 14 18 25 29 2 9 13 17 21 26 1 6 22 30];
                    devtype=1;nring=1;nblk=12;nasic=1;bmap=[8 4];
                    fprintf('Assuming device type %d\n',devtype)
                else
                    fprintf('Device type %d information loaded\n',devtype)
                end
                eval(['ParamList.Device.name  = ',sprintf('''%s''',device_name),';']);
                eval(['ParamList.Device.type  = ',sprintf('%d ',devtype),';']);
                eval(['ParamList.Device.nring = ',sprintf('%d ',nring),';']);
                eval(['ParamList.Device.nblk  = ',sprintf('%d ',nblk),';']);
                eval(['ParamList.Device.nasic = ',sprintf('%d ',nasic),';']);
                eval(['ParamList.Device.bmap = [',sprintf('%d ',bmap),'];']);
                eval(['ParamList.Device.cmap = [',sprintf('%d ',cmap),'];']);
            end
        end
        tline = fgets(fid);
    end
    fclose(fid);
    
catch
    ParamList=[];
    uiwait(errordlg('Error reading default parameters'));
end

%%
% Read device configuration file
function [device_name,nring,nblk,nasic,bmap,cmap] = read_device_type_address_map(dt)
device_name=[];nring=[];nblk=[];nasic=[];bmap=[];cmap=[];
try
    [localpath , ~] =fileparts(mfilename('fullpath'));
    fid = fopen(fullfile(localpath,'PETdevice.cfg'));
    tline = fgets(fid);
    while ischar(tline)
        if ~strcmp(tline(1),'%')
            s = strrep(tline,sprintf('\n'),'');
            s(min(strfind(s,'%')):end)=[]; % remove comments
            ix = strfind(s,':');
            cdt = sscanf(s(ix+1:end),'%d')';
            if cdt(1)==dt,
                device_name = s(1:ix-1);
                nring = int16(cdt(2));
                nblk = int16(cdt(3));
                nasic = int16(cdt(4));
                bmap = int16(cdt(5:6));
                cmap = int16(cdt(7:end));
            end
        end
        tline = fgets(fid);
    end
    fclose(fid);
catch
    
end

%%
% Display initial information
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

%%
% Extract packet headers
function [sn, rate] = extract_packet_headers(y)
sn = uint32(y(2:8:end-6)) + bitshift( uint32(y(1:8:end-7)),16 );
rate = uint32(y(8:8:end)) + bitshift( uint32(mod(y(7:8:end-1),2^12)),16 );

%%
% Extarct events
function [etime, easic, echan] = extract_events(p)
Tf = uint32(p(1:4:end-3))+bitshift(uint32(p(2:4:end-2)),16);
Tc = uint32(mod(p(3:4:end-1),2^14));
% etime=int64(bitor(uint64(Tf),bitshift(uint64(Tc),32))); % in clock units
etime=int64(bitshift(bitor(uint64(Tf),bitshift(uint64(Tc),32)),1)); % in ns units

HighBits = bitshift(uint32(p(3:4:end-1)),-14)+bitshift(uint32(p(4:4:end)),2);
% ev.counter=uint32(mod(HighBits,2^6));
% ev.gate=uint32(mod(bitshift(HighBits,-16),2^2));
TenBits=uint32(mod(bitshift(HighBits,-6),2^10));
easic=int8(bitshift(TenBits,-5));
echan=int8(mod(TenBits,2^5));

%%
% Intialize event matrix to the device configuration
function evr = InitializeEventRate(Cfg)
Na = Cfg.defaults.Device.nring * Cfg.defaults.Device.nblk * Cfg.defaults.Device.nasic;
Nc = prod(Cfg.defaults.Device.bmap);
evr = zeros(Nc,Na);

%%
% prepare figure to be updated in the background
function sfigure(h)
if nargin>=1
    if ishandle(h), set(0, 'CurrentFigure', h);
    else h = figure(h); end %#ok<*NASGU>
else h = figure;
end

function set_default_figure_properties()
% set(0, 'DefaultFigureColor', 'White',...
%     'DefaultFigurePaperType', 'a4letter',...
%     'DefaultAxesColor', 'white',...
%     'DefaultAxesFontUnits', 'points',...
%     'DefaultAxesFontSize', 12,...
%     'DefaultAxesFontAngle', 'normal',...
%     'DefaultAxesFontName', 'Times',...
%     'DefaultAxesGridLineStyle', ':',...
%     'DefaultAxesInterruptible', 'on',...
%     'DefaultAxesLayer', 'Bottom',...
%     'DefaultAxesNextPlot', 'replace',...
%     'DefaultAxesUnits', 'normalized',...
%     'DefaultAxesXcolor', [0, 0, 0],...
%     'DefaultAxesYcolor', [0, 0, 0],...
%     'DefaultAxesZcolor', [0, 0, 0],...
%     'DefaultAxesVisible', 'on',...
%     'DefaultAxesDrawmode', 'fast',...
%     'DefaultLineColor', 'Red',...
%     'DefaultLineLineStyle', '-',...
%     'DefaultLineLineWidth', 1,...
%     'DefaultLineMarker', 'none',...
%     'DefaultLineMarkerSize', 8,...
%     'DefaultTextColor', [0, 0, 0],...
%     'DefaultTextFontUnits', 'Points',...
%     'DefaultTextFontSize', 12,...
%     'DefaultTextFontName', 'Times',...
%     'DefaultTextVerticalAlignment', 'middle',...
%     'DefaultTextHorizontalAlignment', 'left');
% w = warning('off','all');
% rmpath('folderthatisnotinpath')
% warning(w)


% --- Executes on slider movement.
function packets_control_Callback(hObject, eventdata, handles)
% hObject    handle to packets_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gas = get(hObject,'Value');
val = round(2^(10*gas));
set(handles.packets_indicator,'String',sprintf('%d',val))


% --- Executes during object creation, after setting all properties.
function packets_control_CreateFcn(hObject, eventdata, handles)
% hObject    handle to packets_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function packets_indicator_Callback(hObject, eventdata, handles)
% hObject    handle to packets_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of packets_indicator as text
%        str2double(get(hObject,'String')) returns contents of packets_indicator as a double
val = str2double(get(hObject,'String'));
gas = log2(val)/10;
set(handles.packets_control,'Value',gas)
drawnow

% --- Executes during object creation, after setting all properties.
function packets_indicator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to packets_indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% find indexes of packets with events
function [pix, qix] = find_packets(x)
%
% pix - index of packet header start
% qix - packet length
%
pix=[];
qix=[];
if isempty(x), return, end
xb = [1; (x>0)];
cxb = cumsum(xb);
ix4 = find( cxb(5:end) == cxb(1:end-4) ); % find indixes to 4 zeros
pix = ix4-2;
if isempty(pix), return, end
elim=[];
lx = length(x);
lp = length(pix);
for k=2:lp
    % handle the case of the lower (second) byte of the packet serial is
    % zero, which creates 5 consequetive zeros, and is being interpreted
    % and two pix indexes separated by one
    if pix(k-1)==(pix(k)-1) & x(pix(k)+1)==0,
        elim = [elim, k-1];
    end
    % handle the case of empty packets starting from 1 to lp-1
    if pix(k-1)+8<lx
        try
            if x(pix(k-1))==x(pix(k-1)+8)
                elim = [elim, k-1];
            end
        catch
            %             reshape(x,4,[])', k, elim, pix
            disp('Warning: empty packets detected')
        end
    end
    
end
% what if the last one is empty
if pix(lp)+8<lx
    if x(pix(lp))==x(pix(lp)+8)
        elim = [elim, lp];
        empty_packets_flag = 1;
    end
end
elim=unique(elim);

if ~isempty(elim)
    pix(elim)=[];
end

if isempty(pix), return, end
%%
% check for real packet length
lp = length(pix);
qix = [pix(2:end);lx+1]-pix;
for k=1:lp
    n1=x(pix(k));
    n2=x(pix(k)+1); % these two constitute packet number
    mal = min(qix(k),516);
    for n=4:4:mal % there should be at most 520 values in single packet
        try
            if pix(k)+n>lx, break, end
            m1=x(pix(k)+n);
            m2=x(pix(k)+1+n);
        catch
            disp('Warning: crumbled packet structure')
        end
        if n1==m1 & n2==(m2-1)
            qix(k)=n;
            break
        end
    end
end


% if length(pix)>1 | pix(1)~=1,
%     disp('Warning: smaller packets detected')
% end

function print_packet_content(x)
x=reshape(x,4,[])';
for k=1:length(x(:,1))
    disp(sprintf('   %d: %04x %04x %04x %04x',k,x(k,:)))
end

function print_packet_sn(y)
y=reshape(y,2,[])';
for k=1:length(y(:,1))
    disp(sprintf('   %d: %04x %04x',k,y(k,:)))
end


% --- Executes on slider movement.
function averages_control_Callback(hObject, eventdata, handles)
% hObject    handle to averages_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
gas = get(hObject,'Value');
val = round(2^(10*gas));
set(handles.averages_indicator,'String',sprintf('%d',val))

% --- Executes during object creation, after setting all properties.
function averages_control_CreateFcn(hObject, eventdata, handles)
% hObject    handle to averages_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function event_rate_Callback(hObject, eventdata, handles)
% hObject    handle to event_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_rate as text
%        str2double(get(hObject,'String')) returns contents of event_rate as a double


% --- Executes during object creation, after setting all properties.
function event_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drop_rate_Callback(hObject, eventdata, handles)
% hObject    handle to drop_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drop_rate as text
%        str2double(get(hObject,'String')) returns contents of drop_rate as a double


% --- Executes during object creation, after setting all properties.
function drop_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drop_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coinc_rate_Callback(hObject, eventdata, handles)
% hObject    handle to coinc_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coinc_rate as text
%        str2double(get(hObject,'String')) returns contents of coinc_rate as a double


% --- Executes during object creation, after setting all properties.
function coinc_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coinc_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function proc_rate_Callback(hObject, eventdata, handles)
% hObject    handle to proc_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of proc_rate as text
%        str2double(get(hObject,'String')) returns contents of proc_rate as a double


% --- Executes during object creation, after setting all properties.
function proc_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to proc_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in display_mode.
function display_mode_Callback(hObject, eventdata, handles)
% hObject    handle to display_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns display_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display_mode


% --- Executes during object creation, after setting all properties.
function display_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function batch_num_Callback(hObject, eventdata, handles)
% hObject    handle to batch_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batch_num as text
%        str2double(get(hObject,'String')) returns contents of batch_num as a double


% --- Executes during object creation, after setting all properties.
function batch_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batch_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function batch_time_Callback(hObject, eventdata, handles)
% hObject    handle to batch_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batch_time as text
%        str2double(get(hObject,'String')) returns contents of batch_time as a double


% --- Executes during object creation, after setting all properties.
function batch_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batch_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function my_start_Callback(hObject, eventdata, handles)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

btn = findobj(CurFig,'Tag','pushbutton1_start');
switch handles
    case 1
        %         disp(['   ### starting timer',datestr(now,31)])
        set(btn,'String','Start')
        Cfg.start = 0;
        setappdata(CurFig,'Cfg',Cfg);
        drawnow;
    case 2
        %         disp(['   ### triggering action',datestr(now,31)])
        Cfg.start = 0;
        setappdata(CurFig,'Cfg',Cfg);
        drawnow;
    case 3
        %         disp([ '### stopping timer',datestr(now,31)])
end

% --- Executes on button press in pushbutton2_batch_start.
function pushbutton2_batch_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_batch_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

if Cfg.start,
    Cfg.start=0; setappdata(CurFig,'Cfg',Cfg),
    set(hObject,'String','Start'), drawnow;
    return
else,
    set(hObject,'String','Stop'),
end

batchtimestamp = sprintf('%s',datestr(now,30));

if isempty(Cfg.folder), Cfg.folder = pwd; end

batch_root = char(get(handles.batch_root,'String'));
batch_num = round(str2double(get(handles.batch_num,'String')));
batch_time = str2double(get(handles.batch_time,'String'));

if ~isdir(batch_root)
    errordlg('Batch root directory is not valid'), set(hObject,'String','Start'),return
end
if batch_num<1
    errordlg('Batch repetition number is not valid'), set(hObject,'String','Start'),return
end
if batch_time <=0
    errordlg('Batch time is not valid'), set(hObject,'String','Start'),return
end

Cfg.tspm.delay = determine_tspm_delay(Cfg.tspm.delay);
setappdata(CurFig,'Cfg',Cfg), drawnow;

writemode = get(handles.writemode,'value');
if writemode ~=2 & writemode ~=3
%     uiwait(errordlg({'Output mode is not valid','Changing to default'})),
    set(handles.writemode,'value',2), drawnow;
    writemode = get(handles.writemode,'value');
end

t = batch_timer_start(CurFig,batch_time);
for k=1:batch_num
    if strcmp(lower(char(get(hObject,'String'))),'start')
        msgbox('Batch aborted')
        batch_timer_stop(CurFig,t)
        break
    end
    tmp = sprintf('%s executing batch %d',datestr(now,13), k);
    disp(tmp)
    set(handles.batch_status,'string',tmp)
    
    switch writemode,
        case 2
            file2write = sprintf('batch_%s_%04d.evt',batchtimestamp, k);
%             Cfg.eventfile_w = fullfile(Cfg.folder,file2write);
        case 3
            file2write = sprintf('batch_%s_%04d.coinc',batchtimestamp, k);
%             Cfg.coincfile_w = fullfile(Cfg.folder,file2write);
    end
    batch_single_acquisition(CurFig,t,file2write,writemode)
end
batch_timer_stop(CurFig,t)
set(hObject,'String','Start'),


function t = batch_timer_start(CurFig,batch_time)
btn = findobj(CurFig,'Tag','pushbutton1_start');
set(btn,'Enable','inactive')
fh = @my_start_Callback;

t = timer('StartDelay',batch_time);
set(t,'StartFcn',{fh,1});
set(t,'TimerFcn',{fh,2});
set(t,'StopFcn',{fh,3});

function batch_timer_stop(CurFig,t)
btn = findobj(CurFig,'Tag','pushbutton1_start');
set(btn,'Enable','on')
delete(t)

function batch_single_acquisition(CurFig,t,file2write,writemode)
btn = findobj(CurFig,'Tag','pushbutton1_start');
monfile_w = findobj(CurFig,'Tag','monfile_w');
Cfg = getappdata(CurFig,'Cfg');
switch writemode,
    case 2
        Cfg.eventfile_w = fullfile(Cfg.folder,file2write);
        set(monfile_w,'string',file2write)
    case 3
        Cfg.coincfile_w = fullfile(Cfg.folder,file2write);
        set(monfile_w,'string',file2write)
end
setappdata(CurFig,'Cfg',Cfg), drawnow;
start(t)
% disp(['Starting batch into file: ',file2write]);
pushbutton1_start_Callback(btn, [], guidata(CurFig));
stop(t)
% disp(['Stopping batch into file: ',file2write]);
set(monfile_w,'string','')
pause(1)

function batch_root_Callback(hObject, eventdata, handles)
% hObject    handle to batch_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batch_root as text
%        str2double(get(hObject,'String')) returns contents of batch_root as a double
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

if isempty(Cfg.folder), Cfg.folder = pwd; end
set(handles.batch_root,'string',Cfg.folder)

batch_root = uigetdir(Cfg.folder,'Select batch root folder')
if length(batch_root)>1,
    Cfg.folder = batch_root;
end
set(handles.batch_root,'String',Cfg.folder);
setappdata(CurFig,'Cfg',Cfg);

% --- Executes during object creation, after setting all properties.
function batch_root_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batch_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in writemode.
function writemode_Callback(hObject, eventdata, handles)
% hObject    handle to writemode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns writemode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from writemode
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

if Cfg.start,
    tmp = questdlg({'Data acquisition is in progress.','Do you want to stop?'},'Warning');
    if strcmpi(tmp,'yes')
        Cfg.start = 0;
    else
        set(hObject,'Value',Cfg.writemode); drawnow;
        return
    end
end

tmp = get(hObject,'String');
if isempty(tmp), return, end
contents = cellstr(tmp);
selection = get(hObject,'Value');
str = lower(contents{selection});
if isempty(Cfg.folder), Cfg.folder = pwd; end

switch selection
    case 1 % No output
        Cfg.writemode = 1;
        set(handles.monfile_w,'string','')
    case 2 % Write event file
        [fname,paths] = uiputfile(fullfile(Cfg.folder,'*.evt'),'Specify output event file');
        if ischar(fname)
            Cfg.eventfile_w = fullfile(paths,fname);
            Cfg.folder = paths;
        else
            Cfg.eventfile_w=[];
        end
        Cfg.writemode = 2;
        set(handles.monfile_w,'string',fname)
    case 3 % Write coinc file
        [fname,paths] = uiputfile(fullfile(Cfg.folder,'*.coinc'),'Specify output coinc file');
        if ischar(fname)
            Cfg.coincfile_w = fullfile(paths,fname);
            Cfg.folder = paths;
        else
            Cfg.coincfile_w=[];
        end
        Cfg.writemode = 3;
        set(handles.monfile_w,'string',fname)
end
setappdata(CurFig,'Cfg',Cfg);

% --- Executes during object creation, after setting all properties.
function writemode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to writemode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function monfile_r_Callback(hObject, eventdata, handles)
% hObject    handle to monfile_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of monfile_r as text
%        str2double(get(hObject,'String')) returns contents of monfile_r as a double


% --- Executes during object creation, after setting all properties.
function monfile_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monfile_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function monfile_w_Callback(hObject, eventdata, handles)
% hObject    handle to monfile_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of monfile_w as text
%        str2double(get(hObject,'String')) returns contents of monfile_w as a double


% --- Executes during object creation, after setting all properties.
function monfile_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monfile_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5_calibr_start.
function pushbutton5_calibr_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_calibr_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');
tspm_delay = Cfg.tspm.delay;
calibr_add_delay = str2double(get(handles.calibr_add_delay,'String'));

set(handles.calibr_timestamp,'String','');

if strcmp(lower(char(get(hObject,'String'))),'stop'),
    Cfg.start=0; setappdata(CurFig,'Cfg',Cfg),
    set(hObject,'String','Start'), drawnow;
    msgbox('Please wait till calibration sequence stops')
    return
else,
    set(hObject,'String','Stop'),
end

calibr_root = char(get(handles.calibr_root,'String'));

if ~isdir(calibr_root)
    errordlg('Calibration root directory is not valid'),
    set(hObject,'String','Start'), return
end

% derived from HV=400V of Ring 9
Forbiden.Dac = [620 700]; 
Forbiden.Gain = [30 20];
Forbiden.Vector = [diff(Forbiden.Dac), diff(Forbiden.Gain), 0];


g1 = str2double(get(handles.gainstart,'string'));
dg = str2double(get(handles.gainstep,'string'));
g2 = str2double(get(handles.gainstop,'string'));

d1 = str2double(get(handles.dacstart,'string'));
dd = str2double(get(handles.dacstep,'string'));
d2 = str2double(get(handles.dacstop,'string'));

ew = str2double(get(handles.energywindow,'string'));
calibr_duration = str2double(get(handles.calibr_duration,'string'));

gains = g1:dg:g2;
dacs = d1:dd:d2;

DO_NOT_RELOAD_GAINS = 0;
if length(gains)==1
    uiwait(warndlg('Single dac scan. Gains will not be reloaded'))
    DO_NOT_RELOAD_GAINS = 1;
end

if isempty(gains) | isempty(dacs)
    errordlg('Invalid set of gain dac parameters'), return
end

writemode = get(handles.writemode,'value');
if writemode ~=2 & writemode ~=3
    uiwait(errordlg({'Output mode is not valid','Changing to default'})),
    set(handles.writemode,'value',2), drawnow;
    writemode = get(handles.writemode,'value');
end

tmp = get(handles.device_num,'String');
device_type = get(handles.device_num,'value');
if mod(Cfg.defaults.Device.nblk,12), errordlg('Error with detectors'' number. Stop.'), return; end
gaindata = zeros(12*32,4);
% gaindata(:,2)=1; % set vth to 1
gaindata(:,2)=0; % set vth to 1

tspm_timeout = round(500*1e6/25);
set(handles.tspm_timeout,'string','500');
TSPMcontrol('write',36,tspm_timeout); pause(0.1)

%%
% test load to get time right
tic;
if DO_NOT_RELOAD_GAINS
    w=[];
else
    switch Cfg.defaults.Device.nblk
        case 12
            w = prepare_408_words(gaindata);
            for ksure = 1:4,
                calibration_set_gains(w,3,calibr_add_delay);
            end
        case 24
            w = prepare_408_words(gaindata);
            for ksure = 1:4,
                calibration_set_gains(w,3,calibr_add_delay);
                calibration_set_gains(w,4,calibr_add_delay);
            end
        otherwise
            errordlg('Wrong number of detectors'),
            return,
    end
end

%%
% set dacs multiple times for sure
for ksure=1:4,
    calibration_set_dacs(0,min(dacs));
    pause(2*calibr_add_delay),
end
estspm = toc;

estm = length(gains)*length(dacs)*(calibr_duration+estspm);
tmp = sprintf('Estimated total time %d min for %dx%d points(single TSPM delay %4.1f sec)',...
    round(estm/60),length(gains),length(dacs),estspm);
disp(tmp)
set(handles.calibr_status,'string',tmp)
ch = questdlg(sprintf('Asic tester calibration:\nFolder: %s\n%s\nDo you want to continue?',calibr_root,tmp));
if ~strcmp(ch,'Yes')
    set(hObject,'String','Start'),
    set(handles.calibr_status,'string','Calibration sequence canceled')
    return
end

calibr_seq_started = now;
calibr_timestamp = datestr(calibr_seq_started,30);
set(handles.calibr_timestamp,'String',calibr_timestamp);
calibr_seq_estimated = calibr_seq_started + estm/24/60/60;

aborted_sequence = 0;
t = batch_timer_start(CurFig,calibr_duration);

DO_REAL_WORK = 1;

for k=1:length(gains)
    
    if aborted_sequence, break, end
    gaindata(:,4)=gains(k);
    if DO_NOT_RELOAD_GAINS
        w=[];
    else
        w = prepare_408_words(gaindata);
    end
    
    for n=1:length(dacs)
        
        Cfg = getappdata(CurFig,'Cfg');
        
        if strcmp(lower(char(get(hObject,'String'))),'start')
            msgbox('Calibration aborted')
            batch_timer_stop(CurFig,t)
            aborted_sequence = 1;
            break
        end
        
        % check if forbiden dacs and gains
        curdgvector =[dacs(n)-Forbiden.Dac(1), gains(k)-Forbiden.Gain(1), 0];
        crossofvec = cross(curdgvector,Forbiden.Vector);
        
        if crossofvec(3)>0
            if DO_REAL_WORK
                % calibration_set_gains(dacs(n),dacs(n)+ew,w,tspm_delay,calibr_add_delay);
                switch Cfg.defaults.Device.nblk
                    case 12
                        for ksure=1:2,
                            calibration_set_gains(w,3,calibr_add_delay);
                        end
                        for ksure=1:4,
                            calibration_set_dacs(0,dacs(n));
                        end
                    case 24
                        for ksure=1:2,
                            calibration_set_gains(w,3,calibr_add_delay*ksure);
                            calibration_set_gains(w,4,calibr_add_delay*ksure);
                        end
                        for ksure=1:3,
                            calibration_set_dacs(0,dacs(n));
                            pause(calibr_add_delay*ksure)
                        end
                    otherwise
                        errordlg('Wrong number of detectors'),
                        return,
                end
                disp('Threshold scan: disregarding energy window')
            end
        end
        
        tmp = sprintf('%s gain %02d dac %03d %+2d',datestr(calibr_seq_estimated-now,13), gains(k), dacs(n), sign(crossofvec(3)));
        set(handles.calibr_status,'string',tmp)
        disp(tmp)
        
        switch writemode,
            case 2
                file2write = sprintf('calibr %s %02d %03d.evt',calibr_timestamp,gains(k),dacs(n));
            case 3
                file2write = sprintf('calibr %s %02d %03d.coinc',calibr_timestamp,gains(k),dacs(n));
        end
        
        
        if DO_REAL_WORK
            if crossofvec(3)>0
                batch_single_acquisition(CurFig,t,file2write,writemode);
                pause(tspm_delay+calibr_add_delay)
            else
                fclose(fopen(fullfile(Cfg.folder,file2write), 'w'))
            end
        end
    end
end
batch_timer_stop(CurFig,t)
if aborted_sequence
    tmp = sprintf('Sequence aborted at %s',datestr(now-calibr_seq_started,13))
else
    tmp = sprintf('Sequence completed in %s',datestr(now-calibr_seq_started,13))
end
set(handles.calibr_status,'string',tmp)
disp(tmp)

set(hObject,'String','Start'),
set(handles.writemode,'value',1), drawnow;
if ~aborted_sequence
    if DO_REAL_WORK
        pushbutton10_analyze_data_Callback(hObject, eventdata, handles);
    end
end

function calibration_set_dacs(dl,dh)
tdly = 1;
daqval = dl*2^16+dh;
TSPMcontrol('write',11,daqval); pause(tdly)
TSPMcontrol('write',12,daqval); pause(tdly)

function calibration_set_gains(w,spiregister,addfixeddelay)
% nr - number of rings
% reset asics and enable UDP
% tdly = reset_all_tspm(tdly);
tdly = 0.5;
if spiregister<3 | spiregister>4,
    disp('Error: invalid spi register in the call during calibration_set_tspm')
    return
end

%% uncoment this for asic tester
% TSPMcontrol('write',1,hex2dec('00ffffdf')); pause(tdly)
% TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
% TSPMcontrol('write',38,hex2dec('00000020')); pause(tdly)

%     case 1, % single ring
%         TSPMcontrol('write',1,hex2dec('00fff000')); pause(tdly)
%         TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
%         TSPMcontrol('write',38,hex2dec('00000fff')); pause(tdly)
%     case 2, % dual rings
%         TSPMcontrol('write',1,hex2dec('0')); pause(tdly)
%         TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
%         TSPMcontrol('write',38,hex2dec('00ffffff')); pause(tdly)

% load gains
if isempty(w), return, end
if length(w)~=408, return, end
regval = TSPMcontrol('spi write', spiregister, w); pause(tdly+addfixeddelay)





function gainstart_Callback(hObject, eventdata, handles)
% hObject    handle to gainstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gainstart as text
%        str2double(get(hObject,'String')) returns contents of gainstart as a double


% --- Executes during object creation, after setting all properties.
function gainstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gainstep_Callback(hObject, eventdata, handles)
% hObject    handle to gainstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gainstep as text
%        str2double(get(hObject,'String')) returns contents of gainstep as a double


% --- Executes during object creation, after setting all properties.
function gainstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gainstop_Callback(hObject, eventdata, handles)
% hObject    handle to gainstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gainstop as text
%        str2double(get(hObject,'String')) returns contents of gainstop as a double


% --- Executes during object creation, after setting all properties.
function gainstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dacstart_Callback(hObject, eventdata, handles)
% hObject    handle to dacstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dacstart as text
%        str2double(get(hObject,'String')) returns contents of dacstart as a double


% --- Executes during object creation, after setting all properties.
function dacstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dacstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dacstep_Callback(hObject, eventdata, handles)
% hObject    handle to dacstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dacstep as text
%        str2double(get(hObject,'String')) returns contents of dacstep as a double


% --- Executes during object creation, after setting all properties.
function dacstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dacstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dacstop_Callback(hObject, eventdata, handles)
% hObject    handle to dacstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dacstop as text
%        str2double(get(hObject,'String')) returns contents of dacstop as a double


% --- Executes during object creation, after setting all properties.
function dacstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dacstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function energywindow_Callback(hObject, eventdata, handles)
% hObject    handle to energywindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of energywindow as text
%        str2double(get(hObject,'String')) returns contents of energywindow as a double


% --- Executes during object creation, after setting all properties.
function energywindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energywindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function calibr_duration_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calibr_duration as text
%        str2double(get(hObject,'String')) returns contents of calibr_duration as a double


% --- Executes during object creation, after setting all properties.
function calibr_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calibr_root.
function calibr_root_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

if isempty(Cfg.folder), Cfg.folder = pwd; end
set(handles.calibr_root,'string',Cfg.folder)

calibr_root = uigetdir(Cfg.folder,'Select calibration root folder');
if length(calibr_root)>1,
    Cfg.folder = calibr_root;
end
set(handles.calibr_root,'String',Cfg.folder);
setappdata(CurFig,'Cfg',Cfg);


function dac2_set_Callback(hObject, eventdata, handles)
% hObject    handle to dac2_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dac2_set as text
%        str2double(get(hObject,'String')) returns contents of dac2_set as a double


% --- Executes during object creation, after setting all properties.
function dac2_set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dac2_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dac1_set_Callback(hObject, eventdata, handles)
% hObject    handle to dac1_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dac1_set as text
%        str2double(get(hObject,'String')) returns contents of dac1_set as a double


% --- Executes during object creation, after setting all properties.
function dac1_set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dac1_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gains_set_Callback(hObject, eventdata, handles)
% hObject    handle to gains_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gains_set as text
%        str2double(get(hObject,'String')) returns contents of gains_set as a double


% --- Executes during object creation, after setting all properties.
function gains_set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gains_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7_set_tspm.
function pushbutton7_set_tspm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_set_tspm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

% reset TSPM delay if possible
Cfg.tspm.delay = determine_tspm_delay(0.2);
setappdata(CurFig,'Cfg',Cfg), drawnow;

% reset asics and enable UDP
set(handles.tspm_status,'string','Wait until TSPM resets')
tdly = reset_all_tspm(Cfg.tspm.delay);

if tdly>=60
    errordlg({'System reset is too long','Please check network connectivity.'})
    set(handles.tspm_status,'string','TSPM reset failed')
    return
end

tspm_timeout = round(str2double(get(handles.tspm_timeout,'String'))*1e6/25);
% disp(sprintf('TSPM timeout %d %s',tspm_timeout,dec2hex(tspm_timeout)))
TSPMcontrol('write',36,tspm_timeout); pause(tdly)

ResetCase = Cfg.defaults.Device.nring * Cfg.defaults.Device.nblk * Cfg.defaults.Device.nasic/12;

switch ResetCase
    case 1,
        TSPMcontrol('write',1,hex2dec('00fff000')); pause(tdly)
        TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
        TSPMcontrol('write',38,hex2dec('00000fff')); pause(tdly)
    case 2,
        TSPMcontrol('write',1,hex2dec('0')); pause(tdly)
        TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
        TSPMcontrol('write',38,hex2dec('00ffffff')); pause(tdly)
    otherwise,
        TSPMcontrol('write',1,hex2dec('00ffffdf')); pause(tdly)
        TSPMcontrol('write',35,hex2dec('0')); pause(tdly)
        TSPMcontrol('write',38,hex2dec('00000020')); pause(tdly)
end

% set dacs
dl = uint32(str2double(get(handles.dac1_set,'String')));
dh = uint32(str2double(get(handles.dac2_set,'String')));
daqval = dl*2^16+dh;
TSPMcontrol('write',11,daqval); pause(tdly)
TSPMcontrol('write',12,daqval); pause(tdly)
set(handles.tspm_status,'string','TSPM reset completed')

gains_todo = get(handles.gains_todo,'Value');
control_vth_bit = get(handles.control_vth_bit,'Value');

switch gains_todo
    case 1
        gains_load = get(handles.gains_load,'String');
        if isempty(Cfg.gaindata)
            set(handles.tspm_status,'string','Gains were not loaded')
            return
        end
        Cfg.gaindata(:,2)=control_vth_bit;
        w = prepare_408_words(Cfg.gaindata);
    case 2
        Cfg.gaindata = zeros(12*32,4);
        Cfg.gaindata(:,2)=control_vth_bit;
        Cfg.gaindata(:,4)=str2double(get(handles.gains_set,'String'));
        w = prepare_408_words(Cfg.gaindata);
        setappdata(CurFig,'Cfg',Cfg);
    otherwise
        return
end

set(handles.tspm_status,'string','Wait until gains are loaded')
% write gains
switch ResetCase
    case 1,
        regval = TSPMcontrol('spi write', 3, w); pause(tdly)
    case 2,
        TSPMcontrol('spi write', 3, w); pause(tdly)
        TSPMcontrol('spi write', 4, w); pause(tdly)
    otherwise,
        TSPMcontrol('spi write', 3, w); pause(tdly)
end
set(handles.tspm_status,'string','SPI gains loaded')

%%
% reset TSPM to default state
function tdly = reset_all_tspm(inittdly)
tdly = inittdly;
rtn=-1;
while rtn~=256 & tdly < 60
    %     TSPMcontrol('reset all'); pause(tdly)
    TSPMcontrol('write',0,65792); pause(tdly)
    TSPMcontrol('write',0,256); pause(tdly)
    rtn = TSPMcontrol('read',0);
    if isempty(rtn), rtn=-1; end
    disp(sprintf('System reset %d with delay %5.2f',rtn,tdly))
    tdly=tdly+tdly;
end


% --- Executes on selection change in device_num.
function device_num_Callback(hObject, eventdata, handles)
% hObject    handle to device_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns device_num contents as cell array
%        contents{get(hObject,'Value')} returns selected item from device_num
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

device_num = get(hObject,'Value');
contents = cellstr(get(hObject,'String'));
device_name = contents{get(hObject,'Value')};

switch device_num
    case {1,2,3}
        disp(sprintf('Device %s is selected',device_name))
    otherwise
        uiwait(warndlg(sprintf('Unknown device type %d - defaulting to single ring device',device_num)));
        set(hObject,'Value',1)
        device_num = 1;
end

Cfg = initialize_configuration(CurFig, handles, device_num);
setappdata(CurFig,'Cfg',Cfg);

% --- Executes during object creation, after setting all properties.
function device_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to device_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gains_load.
function gains_load_Callback(hObject, eventdata, handles)
% hObject    handle to gains_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

[File2read, Path2read,tmp] = uigetfile([Cfg.folder,'/*.mat'],'Select gain file to read');
if tmp==0, % cancel
    return
else
    try
        q=load([Path2read,File2read]);
        Cfg.gaindata = q.data;
        set(handles.gains_load,'string',File2read)
        setappdata(CurFig,'Cfg',Cfg)
    catch
        errordlg(sprintf('Cannot read file %s',File2read))
    end
    
end

% --- Executes on selection change in gains_todo.
function gains_todo_Callback(hObject, eventdata, handles)
% hObject    handle to gains_todo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gains_todo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gains_todo
chs = get(hObject,'Value');
switch chs
    case 1
        set(handles.gains_set,'enable','off')
    case 2
        set(handles.gains_set,'enable','on')
    case 3
        set(handles.gains_set,'enable','off')
end


% --- Executes during object creation, after setting all properties.
function gains_todo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gains_todo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function w = prepare_408_words(d)
DEBUG_PRINT = 0;

Cfg.NA = 12;
Cfg.NC = 32;

for a=1:Cfg.NA
    
    %     if a==2, DEBUG_PRINT=1, else, DEBUG_PRINT=0; end
    
    offset = (a-1)*(Cfg.NC+2);
    
    ir = (1+(a-1)*Cfg.NC):(a*Cfg.NC);
    
    amon = (d(ir,1)>0); % only zero or one allowed
    vth = (d(ir,2)>0); % only zero or one allowed
    tp = (d(ir,3)>0); %#ok<*NASGU> % only zero or one allowed
    gains = d(ir,4);
    for k=1:32,
        if k<=16,
            w1(2*k-1) = vth(k);
            w1(2*k) = amon(k);
        else
            w2(2*k-1-32) = vth(k);
            w2(2*k-32) = amon(k);
        end
    end
    if DEBUG_PRINT, fprintf('### asic %d\n',a), end
    
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
        w(offset+2+Cfg.NC-c+1) = curw; %
    end
end
tmp = [w(205:408),w(1:204)];
w=tmp;


% --- Executes on button press in pushbutton10_analyze_data.
function pushbutton10_analyze_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10_analyze_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

calibr_timestamp = get(handles.calibr_timestamp,'String');
calibr_root = char(get(handles.calibr_root,'String'));
if ~isdir(calibr_root)
    errordlg('Calibration folder is not valid')
    return
end
NA = Cfg.defaults.Device.nblk * Cfg.defaults.Device.nasic;
NC = prod(Cfg.defaults.Device.bmap);

tmp=[]; for k=1:NA,tmp{k}=sprintf('%2d',k); end
set(handles.calibr_asic,'string',tmp);
tmp=[]; for k=1:NC,tmp{k}=sprintf('%2d',k); end
set(handles.calibr_chan,'string',tmp);
drawnow;

fl = dir(fullfile(calibr_root,sprintf('calibr %s *.evt',calibr_timestamp)));
N= length(fl);
if N==0
    errordlg('Cannot located calibration data')
    return
end
for k=1:N
    tmp = char(regexp(fl(k).name,'\d\d \d\d\d','match'));
    tmp = sscanf(tmp,'%d %d');
    gain(k) = tmp(1);
    dac(k) = tmp(2);
    evr = read_single_evt_file(fullfile(calibr_root,fl(k).name),NC,NA);
    c(k).e = evr;
    total_events(k) = sum(sum(evr));
    saevr = sprintf('%4d ',sum(evr,1));
    disp(sprintf('G %02d D %03d [ %s] with %d events',gain(k),dac(k),saevr,total_events(k)));
    
end
ug = unique(gain);
ud = unique(dac);
NG = length(ug);
ND = N/NG;
gain = reshape(gain,ND,NG);
dac = reshape(dac,ND,NG);
total_events = reshape(total_events,ND,NG);
c =  reshape(c(:),ND,NG);

Calibr.gain = gain(1,:);
Calibr.dac = dac(:,1);
Calibr.tot = total_events;
for nc=1:NC
    for na=1:NA
        for ng=1:NG
            for nd=1:ND
                if isempty(c(nd,ng).e)
                    Calibr.evr(nd,ng,nc,na)=0;
                else
                    Calibr.evr(nd,ng,nc,na) = c(nd,ng).e(nc,na);
                end
            end
        end
    end
end

tmp=[]; for k=1:ND,tmp{k}=sprintf('%3d',Calibr.dac(k)); end
set(handles.calibr_dac,'string',tmp);
tmp=[]; for k=1:NG,tmp{k}=sprintf('%3d',Calibr.gain(k)); end
set(handles.calibr_gain,'string',tmp);
drawnow;    

setappdata(CurFig,'Calibr',Calibr), drawnow;
Calibr.readme = '4D array index order is 1-dacs, 2-gains, 3-channels, 4-asics';
calibrfile = fullfile(calibr_root,sprintf('calibr %s.mat',calibr_timestamp));
Calibr.ded = compute_diff_dac(Calibr);
save(calibrfile,'Calibr');
calibr_show_Callback(hObject, eventdata, handles)

%%
% Read single not too long event file - calibration processing
function evr = read_single_evt_file(filename, NC, NA)
MACHINE_FORMAT = 'ieee-be';
evr=[];

[fid, msg] = fopen(filename,'r',MACHINE_FORMAT);
if ~isempty(msg), disp(sprintf('Error: %s',msg)), return; end

[cx, DataLength] = fread(fid, 'uint16',MACHINE_FORMAT);
if DataLength,
    %     if ~isempty( regexp(filename,'09 300') )
    %         disp('Got it')
    %     end
    
    [pix, qix] = find_packets(cx); % get packet starts and length
    npix = length(pix);
    y=[];x=[];
    for k=1:npix
        y = [y; cx(pix(k):pix(k)+7)]; % get first 8 words for packet info
        x=[x; cx(pix(k)+8:pix(k)+qix(k)-1)];
    end
    
    [etime, easic, echan] = extract_events(x);
    NumberEvents = length(etime);
    
    evr = zeros(NC,NA);
    
    [NC, NA] = size(evr);
    for k=1:NumberEvents
        na = 1+easic(k);
        nc = 1+echan(k);
        if na <= NA & nc <= NC,
            evr(nc,na) = evr(nc,na)+1;
        else
            disp(sprintf('Error read with asic %d and channel %d length of data %d',na, nc, DataLength))
            %             return
        end
    end
    
end
fclose(fid);

% --- Executes on selection change in calibr_show.
function calibr_show_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_show contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_show
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');
Calibr = getappdata(CurFig,'Calibr');
MonitorAxes = handles.monitoraxes;

if isempty(Calibr), return, end
try 
    Calibr.ded;
catch
  Calibr.ded = compute_diff_dac(Calibr);
end

if ishandle(CurFig), set(0, 'CurrentFigure', CurFig);
else CurFig = figure(CurFig); end
set(CurFig,'CurrentAxes',MonitorAxes),

calibr_show = get(handles.calibr_show,'value');
if calibr_show<7, set(handles.calibr_global_max,'enable','off');
else, set(handles.calibr_global_max,'enable','on'); end
calibr_global_max = str2double(get(handles.calibr_global_max,'string'));

tmp = get(handles.calibr_asic,'string');
na = str2num(char(tmp(get(handles.calibr_asic,'value'))));
tmp = get(handles.calibr_chan,'string');
nc = str2num(char(tmp(get(handles.calibr_chan,'value'))));

nd = get(handles.calibr_dac,'value');
ng = get(handles.calibr_gain,'value');

switch calibr_show
    case 1
        imagesc(Calibr.gain,Calibr.dac,log(Calibr.tot))
        set(MonitorAxes, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 2
        imagesc(Calibr.gain,Calibr.dac, sum(Calibr.evr(:,:,:,na),3))
        set(MonitorAxes, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 3
        imagesc(Calibr.gain,Calibr.dac, sum(Calibr.evr(:,:,nc,:),4))
        set(MonitorAxes, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 4
        imagesc(Calibr.gain,Calibr.dac, Calibr.evr(:,:,nc,na))
        set(MonitorAxes, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 5
        plot(Calibr.dac,Calibr.evr(:,ng,nc,na));
        xlabel('Dacs')
        
        if length(Calibr.dac)>5,
%             cevr = smooth(Calibr.evr(:,ng,nc,na),5);
%             dcevr = [0; diff(cevr);]
            [hAx,hLine1,hLine2] = plotyy(Calibr.dac,Calibr.evr(:,ng,nc,na),Calibr.dac,Calibr.ded(:,ng,nc,na));
            xlabel('Dacs')
            ylabel(hAx(1),'Counts') % left y-axis
            ylabel(hAx(2),'Delta') % right y-axis
            grid
        end   
%         figure(11), kk=1; for k=1:4:16,
%             sevr(:,kk) = smooth(Calibr.evr(:,k,nc,na),5);
%             plot(Calibr.dac,Calibr.evr(:,k,nc,na),'linewidth',2); hold on,
%             slgd{kk} = sprintf('gain %4.0f',Calibr.gain(k)); kk=kk+1;
%         end,
%         hold off, grid, xlabel('Dac'), ylabel('Counts'), set(gca,'fontsize',14), legend(slgd)
%         hold on,
%         for k=1:length(slgd)
%             plot(Calibr.dac,sevr(:,k)), hold on
%         end
%         hold off
% 
%         figure(12), kk=1; for k=1:4:16,
%             plot(Calibr.dac(1:end-1),diff(sevr(:,kk)),'linewidth',2); hold on,
%             slgd{kk} = sprintf('gain %4.0f',Calibr.gain(k)); kk=kk+1;
%         end,
%         hold off, grid, xlabel('Dac'), ylabel('delta Counts'), set(gca,'fontsize',14), legend(slgd)

    case 6
        plot(Calibr.gain,Calibr.evr(nd,:,nc,na));
        xlabel('Gains')
%         figure(11), kk=1; for k=1:13:41, 
%             plot(Calibr.gain,Calibr.evr(k,:,nc,na),'linewidth',2); hold on, 
%             slgd{kk} = sprintf('dac %4.0f',Calibr.dac(k)); kk=kk+1;
%         end, 
%         hold off, grid, xlabel('Gain'), ylabel('Counts'), set(gca,'fontsize',14), legend(slgd)
    case 7
        FigureTag = '32 detectors maps'; hfig=findobj('Tag',FigureTag);
        if isempty(hfig), hfig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off'); end;
        sfigure(hfig)
        curevr = Calibr.evr(:,:,:,na); curevr=curevr(:);
        mevr = mean(curevr);
        sevr = std(curevr);
        disp(sprintf('Asic %d has event rate minmax [%6.1f-%6.1f] mean %6.1f(%6.1f)',na, min(curevr),max(curevr),mevr, sevr))
        for k=1:32,
            subplot(4,8,k)
            imagesc(Calibr.gain,Calibr.dac, Calibr.evr(:,:,k,na))
            set(gca, 'ydir', 'reverse','ytick',[],'xtick',[],'clim',[0 calibr_global_max])
        end
    case 8
        FigureTag = '32 detectors diff'; hfig=findobj('Tag',FigureTag);
        if isempty(hfig), hfig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off'); end;
        sfigure(hfig)
        curded = Calibr.ded(:,:,:,na); curded=curded(:);
        mevr = mean(curded);
        sevr = std(curded);
        disp(sprintf('Asic %d has photo peak  minmax [%6.1f-%6.1f] mean %6.1f(%6.1f) ',na, min(curded),max(curded),mevr, sevr))
        for k=1:32,
            subplot(4,8,k)
            imagesc(Calibr.gain,Calibr.dac, Calibr.ded(:,:,k,na))
            set(gca, 'ydir', 'reverse','ytick',[],'xtick',[],'clim',[0 calibr_global_max])
        end
        
    otherwise
        disp('TBD')
end
colormap bone

function devrddac = compute_diff_dac(Calibr)

[ND,NG,NC,NA] = size(Calibr.evr);

if ND>5,
    for na=1:NA,
        for nc=1:NC,
            for ng=1:NG

            cevr = smooth(Calibr.evr(:,ng,nc,na),10);
            devrddac(:,ng,nc,na) = [0; diff(cevr);];
            end
        end
    end
else
    devrddac = [];
end

% --- Executes during object creation, after setting all properties.
function calibr_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in calibr_asic.
function calibr_asic_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_asic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_asic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_asic
calibr_show_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function calibr_asic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_asic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in calibr_chan.
function calibr_chan_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_chan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_chan
calibr_show_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function calibr_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in calibr_dac.
function calibr_dac_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_dac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_dac contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_dac
calibr_show_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function calibr_dac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_dac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in calibr_gain.
function calibr_gain_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_gain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_gain
calibr_show_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function calibr_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in control_vth_bit.
function control_vth_bit_Callback(hObject, eventdata, handles)
% hObject    handle to control_vth_bit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of control_vth_bit



function calibr_add_delay_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_add_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calibr_add_delay as text
%        str2double(get(hObject,'String')) returns contents of calibr_add_delay as a double


% --- Executes during object creation, after setting all properties.
function calibr_add_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_add_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tspm_timeout_Callback(hObject, eventdata, handles)
% hObject    handle to tspm_timeout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tspm_timeout as text
%        str2double(get(hObject,'String')) returns contents of tspm_timeout as a double


% --- Executes during object creation, after setting all properties.
function tspm_timeout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tspm_timeout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% Calculate Singoram
% function [sgr, zgr] = CalculateSinogramOld(Cfg, t, easic, echan, sgr_a, zgr_a, NumAverages)
% % Cfg - various definitions
% % t - event times
% % easic - asic number 0 - 11
% % echan - channel number 0 - 31
% % sgr_a - [48 48 15] old sinogram
% % zgr_a - [17 9] old zinogram
% %
% % sgr - [48 48 15] array singoram
% % zgr - [17 9] array zinogram
%
% cmap = Cfg.defaults.Device.cmap;
% XX = Cfg.defaults.Device.bmap;
%
% detime = diff(double(t));
%
% icoinc = find(detime>Cfg.defaults.Coincidence_time_window(1) & detime<Cfg.defaults.Coincidence_time_window(2));
% detime = detime(icoinc);
% a2 = easic(icoinc+1);
% c2 = echan(icoinc+1);
% a1 = easic(icoinc);
% c1 = echan(icoinc);
%
% sgr=[]; zgr=[];
%
% if isempty(a2) | isempty(a1) | isempty(c1) | isempty(c2),
%     sgr=sgr_a;
%     zgr=zgr_a;
%     return,
% end
%
% if isempty(a2) % all events are included
%     a2=a1(2:end);   c2=c1(2:end);
%     t(end)=[]; a1(end)=[]; c1(end)=[];
% end
%
% % c1 = changem(c1, cmap-1, ((1:length(cmap))-1) );
% % c2 = changem(c2, cmap-1, ((1:length(cmap))-1) );
% c1 = int16(cmap(c1+1)-1)';
% c2 = int16(cmap(c2+1)-1)';
%
% switch Cfg.defaults.Device_type,
%     %%
%     % Single ring 12 elements
%
%     case 1,
%         if max(a1)>11 | max(a2)>11,
%             fprintf('Warining: number of asics exceeds 12\n')
%         end
%
%         try
%             fi1 = mod(c1,XX(2)) + XX(2)*mod(a1,12);
%             %             zi1 = ceil(c1/XX(2));
%             zi1 = ceil(double(c1+1)/XX(2))-1;
%
%             fi2 = mod(c2,XX(2)) + XX(2)*mod(a2,12);
%             %             zi2 = ceil(c2/XX(2));
%             zi2 = ceil(double(c2+1)/XX(2))-1;
%         catch
%             fi1, zi1, fi2, zi2
%             disp('#### sinogram calc work to do ###')
%         end
%
%         % 3D sinogram
%         sgr = double(zeros(48,48,15));
%
%         xfi = 1+ mod(fi1-fi2,48);
%         yfi = 1+ fi1;
%
%         zgr = double(zeros(17,9));
%
%         xzi = 1+ abs(zi1-zi2);
%         yzi = 1+ zi1+zi2;
%
%         for k=1:length(yfi);
%             n1 = yfi(k);
%             n2 = xfi(k);
%             if xzi(k)<=2 % only considering same plane and adjascent planes coincidences
%                 % same plane comes on odd yzi and adjascent planes come on
%                 % even yzi index
%                 n3 = yzi(k);
%                 try
%                     sgr(n1,n2,n3) = sgr(n1,n2,n3)+(1+mod(n3,2))*Cfg.effcorr(n1,n2,n3);
%                 catch
%                     fprintf('### device %d sinogram calc problems
%                     ###\n',Cfg.defaults.Device_type)
%                 end
%                 %                 disp(sprintf('xzi=%d yzi=%d',xzi(k), yzi(k)))
%             end
%             n1 = yzi(k);
%             n2 = xzi(k);
%             zgr(n1,n2) = zgr(n1,n2)+1;
%         end
%     otherwise,
%         disp('Warning: unrecognized device type, cannot display specific info');
% end
%
% % if isempty(sgr_a), sgr_a=0; end
% % if isempty(zgr_a), zgr_a=0; end
%
% if NumAverages>1,
%     sgr = (sgr + sgr_a*(NumAverages-1))/NumAverages;
%     zgr = (zgr + zgr_a*(NumAverages-1))/NumAverages;
% end

%%
% Calculate t0 correction
function [t0, tshift] = Calculate_t0_correction(Cfg, t, easic, echan, t0_a, NumAverages, coinctimevector)
t0=[];

switch Cfg.defaults.Device_type,
    case 1, % Single ring 12 elements
        if numel(t0_a)==1, % this is the first time
            t0=zeros(length(coinctimevector),12,32,12,32); t0_a=t0;
            sb = whos('t0');
            disp(sprintf('Allocated %4.1f Mb for t0 correction array',sb.bytes/1e6))
        end
    otherwise
        disp('Error initializing t0 correction due to unrecognized device type')
        pause(1)
        return
end

t0=t0_a;

DO_NEIGHBORS_PRECOINC =0;

detime = diff(double(t));

icoinc = find(detime>-coinctimevector(end) & detime<coinctimevector(end));

a2 = easic(icoinc+1);
c2 = echan(icoinc+1);
t2 = t(icoinc+1);

a1 = easic(icoinc);
c1 = echan(icoinc);
t1 = t(icoinc);
% whos a1 a2 c1 c2 t1 t2
% d1 = double(a1)*32 + double(c1)+1;
% d2 = double(a2)*32 + double(c2)+1;

if isempty(a2) | isempty(a1) | isempty(c1) | isempty(c2)
    %     disp('... empty ...')
    return,
end

switch Cfg.defaults.Device_type,
    case 1, % Single ring 12 elements
        starttime = tic;
        for ka1=1:12
            for kc1=1:32
                x1 = find(a1==ka1-1 & c1==kc1-1);
                if ~isempty(x1)
                    for ka2=1:12
                        if ka1~=ka2
                            for kc2=1:32
                                x2 = find(a2(x1)==ka2-1 & c2(x1)==kc2-1);
                                if ~isempty(x2)
                                    ptr = x1(x2);
                                    tdiff = double(t2(ptr))-double(t1(ptr));
                                    try
                                        %                                     t0(:,ka1,kc1,ka2,kc2) = (t0_a(:,ka1,kc1,ka2,kc2)*(NumAverages-1) + histc(tdiff(:)',coinctimevector)')/NumAverages;
                                        t0(:,ka1,kc1,ka2,kc2) = t0_a(:,ka1,kc1,ka2,kc2) + histc(tdiff(:)',coinctimevector)';
                                    catch
                                        disp('Error: t0 correction problem')
                                    end
                                end
                                ct0 = t0(:,ka1,kc1,ka2,kc2);
                                mx = max(ct0);
                                ix = find(ct0==mx);
                                tshift(ka1*kc1,ka2*kc2) = mean(coinctimevector(ix));
                            end
                        end
                    end
                end
            end
        end
        eltime = toc(starttime);
        disp(sprintf('t0 hist (%d/%d/%d) for %d evts in %4.3f sec',...
            nnz(t0_a),nnz(t0),numel(t0),length(t1),eltime))
    otherwise;
end







%%
% Calculate Singoram
function [sgr, zgr] = CalculateSinogram(Cfg, t, easic, echan, sgr_a, zgr_a, NumAverages)
% Cfg - various definitions
% t - event times
% easic - asic number 0 - 11
% echan - channel number 0 - 31
% sgr_a - [48 48 15] old sinogram
% zgr_a - [17 9] old zinogram
%
% sgr - [48 48 15] array singoram
% zgr - [17 9] array zinogram

detime = diff(double(t));
LE1 = length(detime);

icoinc = find(detime>Cfg.defaults.In_time_window(1) & detime<Cfg.defaults.In_time_window(2));
detime = detime(icoinc);
a2 = easic(icoinc+1);
c2 = echan(icoinc+1);
a1 = easic(icoinc);
c1 = echan(icoinc);
LE2 = length(detime);

% disp(sprintf('Events: %d -> %d',LE1, LE2))

sgr=[]; zgr=[];

if isempty(a2) | isempty(a1) | isempty(c1) | isempty(c2)
    sgr=sgr_a;
    zgr=zgr_a;
    return,
end
Kye = int32(ones(size(c2)));
SZX = size(Cfg.sindex);

% in this operation indexes should be casted to allow correct math
ns2 = Cfg.sindex(sub2ind(SZX,int32(a1)+1,int32(c1)+1,int32(a2)+1,int32(c2)+1,Kye));
ns1 = Cfg.sindex(sub2ind(SZX,int32(a1)+1,int32(c1)+1,int32(a2)+1,int32(c2)+1,2*Kye));
nz2 = Cfg.sindex(sub2ind(SZX,int32(a1)+1,int32(c1)+1,int32(a2)+1,int32(c2)+1,3*Kye));
nz1 = Cfg.sindex(sub2ind(SZX,int32(a1)+1,int32(c1)+1,int32(a2)+1,int32(c2)+1,4*Kye));

switch Cfg.defaults.Device_type,
    case 1, % Single ring 12 elements
        
        sgr = double(zeros(48,48,15));
        zgr = double(zeros(17,9));
        
        for k=1:length(ns1);
            n1 = ns1(k);
            n2 = ns2(k);
            if nz2(k)<=2 % only considering same plane and adjascent planes coincidences
                % same plane comes on odd yzi and adjascent planes come on
                % even yzi index
                n3 = nz1(k);
                try
                    sgr(n1,n2,n3) = sgr(n1,n2,n3)+(1+mod(n3,2))*Cfg.effcorr(n1,n2,n3);
                catch
                    fprintf('### device %d sinogram calc problems ###\n',Cfg.defaults.Device_type)
                end
                %                 disp(sprintf('xzi=%d yzi=%d',xzi(k), yzi(k)))
            end
            n1 = nz1(k);
            n2 = nz2(k);
            zgr(n1,n2) = zgr(n1,n2)+1;
        end
    case 2, % Double ring 2 x 12 elements
        
        sgr = double(zeros(48,48,31));
        zgr = double(zeros(33,17));
        
        for k=1:length(ns1);
            n1 = ns1(k);
            n2 = ns2(k);
            if nz2(k)<=2 % only considering same plane and adjascent planes coincidences
                % same plane comes on odd yzi and adjascent planes come on
                % even yzi index
                n3 = nz1(k);
                try
                    sgr(n1,n2,n3) = sgr(n1,n2,n3)+(1+mod(n3,2))*Cfg.effcorr(n1,n2,n3);
                catch
                    fprintf('### device %d sinogram calc problems ###\n',Cfg.defaults.Device_type)
                end
                %                 disp(sprintf('xzi=%d yzi=%d',xzi(k), yzi(k)))
            end
            n1 = nz1(k);
            n2 = nz2(k);
            zgr(n1,n2) = zgr(n1,n2)+1;
        end
    case 3, % Large ring composed of 24 elements - plant scanner
        
        sgr = double(zeros(96,96,15));
        zgr = double(zeros(17,9));
        
        for k=1:length(ns1);
            n1 = ns1(k);
            n2 = ns2(k);
            if nz2(k)<=2 % only considering same plane and adjascent planes coincidences
                % same plane comes on odd yzi and adjascent planes come on
                % even yzi index
                n3 = nz1(k);
                try
                    sgr(n1,n2,n3) = sgr(n1,n2,n3)+(1+mod(n3,2))*Cfg.effcorr(n1,n2,n3);
                catch
                    fprintf('### device %d sinogram calc problems ###\n',Cfg.defaults.Device_type)
                end
                %                 disp(sprintf('xzi=%d yzi=%d',xzi(k), yzi(k)))
            end
            n1 = nz1(k);
            n2 = nz2(k);
            zgr(n1,n2) = zgr(n1,n2)+1;
        end    
    otherwise,
        disp('Warning: unrecognized device type, cannot display specific info');
end

% if isempty(sgr_a), sgr_a=0; end
% if isempty(zgr_a), zgr_a=0; end

if NumAverages>1,
    sgr = (sgr + sgr_a*(NumAverages-1))/NumAverages;
    zgr = (zgr + zgr_a*(NumAverages-1))/NumAverages;
end

%%
% Calculate index to transfomr from channel asic space into sinogram space
%
% Sinogram(a,b,c) = sindex(asic,channel)
%
function sindex = CalculateSinogramIndex(Cfg)

cmap = Cfg.defaults.Device.cmap;
XX = Cfg.defaults.Device.bmap;
DZ = XX(1); % detectors along Z axis in each asic
DT = XX(2); % detectors along circumference in each asic

NC = prod(Cfg.defaults.Device.bmap);
NA = Cfg.defaults.Device.nring*Cfg.defaults.Device.nblk*Cfg.defaults.Device.nasic;
NTALL = DT*NA; % total number of detectors in single z virtual ring
NZALL = DZ*Cfg.defaults.Device.nring; % total number of detectors in z directions

switch Cfg.defaults.Device_type,
    case 1, % Single ring 12 elements
        disp('Sinogram for single ring is being computed')
        
        cc1 = cmap(1:NC); % note cc1 start from 1
        cc2 = cc1;
        
        aa1 = 1:NA; % aa1 start from 1
        aa2 = aa1;
        
        sindex = zeros(NA,NC,NA,NC,4,'int8');
        
        for nc1 = 1:NC,
            c1 = cc1(nc1); % remapping
            for na1 = 1:NA,
                a1 = aa1(na1);
                for nc2 = 1:NC
                    c2 = cc2(nc2); % remapping
                    for na2 = 1:NA
                        a2 = aa2(na2);
                        fi1 = mod(c1-1,DT) + DT*mod(a1-1,NA);
                        zi1 = ceil(double(c1)/DT)-1;
                        
                        fi2 = mod(c2-1,DT) + DT*mod(a2-1,NA);
                        zi2 = ceil(double(c2)/DT)-1;
                        
                        xfi = 1+ mod(fi1-fi2,NTALL);
                        yfi = 1+ fi1;
                        
                        xzi = 1+ abs(zi1-zi2);
                        yzi = 1+ zi1+zi2;
                        if numel(a1)>1|numel(c1)>1|numel(a2)>1|numel(c2)>1,
                            disp('Error: number of elements wrong')
                        end
                        if yfi<1 | yfi>NTALL |...
                                xfi<1 | xfi>NTALL |...
                                yzi<1 | yzi>(NZALL*2+1) |...
                                xzi<1 | xzi>(NZALL+1),
                            disp(sprintf('Index problem acac %d %d %d %d -> S %d %d Z %d %d yxyx',...
                                na1,nc1,na2,nc2, yfi, xfi, yzi, xzi));
                        end
                        
                        try
                            sindex(na1,nc1,na2,nc2,1) = int8(xfi);
                            sindex(na1,nc1,na2,nc2,2) = int8(yfi);
                            sindex(na1,nc1,na2,nc2,3) = int8(xzi);
                            sindex(na1,nc1,na2,nc2,4) = int8(yzi);
                        catch
                            disp('Error calculating sinoindex');
                        end
                    end
                end
            end
        end
        %         max(sindex(:))
        %         min(sindex(:))
        sb = whos('sindex'); 
        disp(sprintf('Device %d %d elements sinoindex has %d bytes',...
            Cfg.defaults.Device_type, numel(sindex),sb.bytes));
    case 2, % Double ring 2 x 12 elements
        disp('Sinogram for two rings facing each other is being computed')
        
        cc1 = cmap(1:NC); % note cc1 start from 1
        cc2 = cc1;
        
        aa1 = 1:NA; % aa1 start from 1
        aa2 = aa1;
        
        sindex = zeros(NA,NC,NA,NC,4,'int8');
        
        for nc1 = 1:NC,
            c1 = cc1(nc1); % remapping
            for na1 = 1:NA,
                a1 = aa1(na1);
                for nc2 = 1:NC
                    c2 = cc2(nc2); % remapping
                    for na2 = 1:NA
                        a2 = aa2(na2);
                        fi1 = mod(c1-1,DT) + DT*mod(a1-1,NA);
                        zi1 = ceil(double(c1)/DT)-1;
                        
                        fi2 = mod(c2-1,DT) + DT*mod(a2-1,NA);
                        zi2 = ceil(double(c2)/DT)-1;
                        
                        xfi = 1+ mod(fi1-fi2,NTALL);
                        yfi = 1+ fi1;
                        
                        xzi = 1+ abs(zi1-zi2);
                        yzi = 1+ zi1+zi2;
                        if numel(a1)>1|numel(c1)>1|numel(a2)>1|numel(c2)>1,
                            disp('Error: number of elements wrong')
                        end
                        if yfi<1 | yfi>NTALL |...
                                xfi<1 | xfi>NTALL |...
                                yzi<1 | yzi>(NZALL*2+1) |...
                                xzi<1 | xzi>(NZALL+1),
                            disp(sprintf('Index problem acac %d %d %d %d -> S %d %d Z %d %d yxyx',...
                                na1,nc1,na2,nc2, yfi, xfi, yzi, xzi));
                        end
                        
                        try
                            sindex(na1,nc1,na2,nc2,1) = int8(xfi);
                            sindex(na1,nc1,na2,nc2,2) = int8(yfi);
                            sindex(na1,nc1,na2,nc2,3) = int8(xzi);
                            sindex(na1,nc1,na2,nc2,4) = int8(yzi);
                        catch
                            disp('Error calculating sinoindex');
                        end
                    end
                end
            end
        end
        %         max(sindex(:))
        %         min(sindex(:))
        sb = whos('sindex');
        disp(sprintf('Device %d %d elements sinoindex has %d bytes',...
            Cfg.defaults.Device_type, numel(sindex),sb.bytes));
    case 3, % Large ring 24 asics - plant scanner
        disp('Sinogram for large ring aka plant scanner is being computed')
        
        cc1 = cmap(1:NC); % note cc1 start from 1
        cc2 = cc1;
        
        % asics are positioned sequentially from 1 to 24 over circumference
        aa1 = 1:NA; % aa1 start from 1 to NA
        aa2 = aa1;

        % asics are flipped
%         aa1 = [[1:12],[24:-1:13]]; 
%         aa2 = aa1;
        
        sindex = zeros(NA,NC,NA,NC,4,'int8');
        
        for nc1 = 1:NC,
            c1 = cc1(nc1); % remapping
            for na1 = 1:NA,
                a1 = aa1(na1);
                for nc2 = 1:NC
                    c2 = cc2(nc2); % remapping
                    for na2 = 1:NA
                        a2 = aa2(na2);
                        fi1 = mod(c1-1,DT) + DT*mod(a1-1,NA);
                        zi1 = ceil(double(c1)/DT)-1;
                        
                        fi2 = mod(c2-1,DT) + DT*mod(a2-1,NA);
                        zi2 = ceil(double(c2)/DT)-1;
                        
                        xfi = 1+ mod(fi1-fi2,NTALL);
                        yfi = 1+ fi1;
                        
                        xzi = 1+ abs(zi1-zi2);
                        yzi = 1+ zi1+zi2;
                        if numel(a1)>1|numel(c1)>1|numel(a2)>1|numel(c2)>1,
                            disp('Error: number of elements wrong')
                        end
                        if yfi<1 | yfi>NTALL |...
                                xfi<1 | xfi>NTALL |...
                                yzi<1 | yzi>(NZALL*2+1) |...
                                xzi<1 | xzi>(NZALL+1),
                            disp(sprintf('Index problem acac %d %d %d %d -> S %d %d Z %d %d yxyx',...
                                na1,nc1,na2,nc2, yfi, xfi, yzi, xzi));
                        end
                        
                        try
                            sindex(na1,nc1,na2,nc2,1) = int8(xfi);
                            sindex(na1,nc1,na2,nc2,2) = int8(yfi);
                            sindex(na1,nc1,na2,nc2,3) = int8(xzi);
                            sindex(na1,nc1,na2,nc2,4) = int8(yzi);
                        catch
                            disp('Error calculating sinoindex');
                        end
                    end
                end
            end
        end
        %         max(sindex(:))
        %         min(sindex(:))
        sb = whos('sindex');
        disp(sprintf('Device %d %d elements sinoindex has %d bytes',...
            Cfg.defaults.Device_type, numel(sindex),sb.bytes));

    otherwise,
        disp('Warning: unrecognized device type');
        sindex=[];
end


% --- Executes on selection change in sino_plane_type.
function sino_plane_type_Callback(hObject, eventdata, handles)
% hObject    handle to sino_plane_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sino_plane_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sino_plane_type


% --- Executes during object creation, after setting all properties.
function sino_plane_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sino_plane_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_axial_slice.
function image_axial_slice_Callback(hObject, eventdata, handles)
% hObject    handle to image_axial_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_axial_slice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_axial_slice


% --- Executes during object creation, after setting all properties.
function image_axial_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_axial_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_coronal_slice.
function image_coronal_slice_Callback(hObject, eventdata, handles)
% hObject    handle to image_coronal_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_coronal_slice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_coronal_slice


% --- Executes during object creation, after setting all properties.
function image_coronal_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_coronal_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_sagital_slice.
function image_sagital_slice_Callback(hObject, eventdata, handles)
% hObject    handle to image_sagital_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_sagital_slice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_sagital_slice


% --- Executes during object creation, after setting all properties.
function image_sagital_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_sagital_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_view.
function image_view_Callback(hObject, eventdata, handles)
% hObject    handle to image_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_view
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

N = round(str2double(get(handles.image_zoom,'string')));
set(handles.image_zoom,'string',sprintf('%2d',N));

tmp=[]; for k=1:N, tmp{k}=sprintf('%2d',k); end
set(handles.image_coronal_slice,'string',tmp);
tmp=[]; for k=1:N, tmp{k}=sprintf('%2d',k); end
set(handles.image_sagital_slice,'string',tmp);

image_view = get(handles.image_view,'value');
switch image_view
    case 1 % axial
        set(handles.image_axial_slice,'enable','on');
        set(handles.image_coronal_slice,'enable','off');
        set(handles.image_sagital_slice,'enable','off');
    case 2 % coronal
        set(handles.image_axial_slice,'enable','off');
        set(handles.image_coronal_slice,'enable','on');
        set(handles.image_sagital_slice,'enable','off');
    case 3 % sagital
        set(handles.image_axial_slice,'enable','off');
        set(handles.image_coronal_slice,'enable','off');
        set(handles.image_sagital_slice,'enable','on');
end

Cfg.image.NumberSlices(2)=N;
Cfg.image.NumberSlices(3)=N;

setappdata(CurFig,'Cfg',Cfg);


% --- Executes during object creation, after setting all properties.
function image_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to image_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_zoom as text
%        str2double(get(hObject,'String')) returns contents of image_zoom as a double
image_view_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function image_zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_scaling.
function image_scaling_Callback(hObject, eventdata, handles)
% hObject    handle to image_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_scaling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_scaling
val = get(hObject,'Value');
if val~=3,
    set(handles.image_max,'enable','off');
else
    set(handles.image_max,'enable','on');
end
image_view_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function image_scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in image_all.
function image_all_Callback(hObject, eventdata, handles)
% hObject    handle to image_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of image_all



function calibr_global_max_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_global_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calibr_global_max as text
%        str2double(get(hObject,'String')) returns contents of calibr_global_max as a double
calibr_show_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function calibr_global_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_global_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close CalibratorPET.
function PETmonitor_new_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CalibratorPET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');
Cfg.start=0;
setappdata(CurFig,'Cfg',Cfg);
try
    fclose all;
catch
    disp('Exiting: current operations are aborted by user')
end
delete(hObject);



function image_max_Callback(hObject, eventdata, handles)
% hObject    handle to image_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_max as text
%        str2double(get(hObject,'String')) returns contents of image_max as a double
image_view_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function image_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in effcorr_load.
function effcorr_load_Callback(hObject, eventdata, handles)
% hObject    handle to effcorr_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

[File2read, Path2read,tmp] = uigetfile([Cfg.folder,'/*.mat'],'Select efficiency correction file to read');
if tmp==0, % cancel
    Cfg.effcorr = read_efficiency_corrections(Cfg,[]);
    set(handles.effcorr_load,'string',' ')
    setappdata(CurFig,'Cfg',Cfg)
    return
else
    try
        Cfg.effcorr = read_efficiency_corrections(Cfg,fullfile(Path2read,File2read));
        set(handles.effcorr_load,'string',File2read)
        setappdata(CurFig,'Cfg',Cfg)
    catch
        errordlg(sprintf('Cannot read file %s',File2read))
    end
end



function io_file_start_Callback(hObject, eventdata, handles)
% hObject    handle to io_file_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of io_file_start as text
%        str2double(get(hObject,'String')) returns contents of io_file_start as a double


% --- Executes during object creation, after setting all properties.
function io_file_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to io_file_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function io_file_stop_Callback(hObject, eventdata, handles)
% hObject    handle to io_file_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of io_file_stop as text
%        str2double(get(hObject,'String')) returns contents of io_file_stop as a double


% --- Executes during object creation, after setting all properties.
function io_file_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to io_file_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in control_asic_onoff.
function control_asic_onoff_Callback(hObject, eventdata, handles)
% hObject    handle to control_asic_onoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns control_asic_onoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from control_asic_onoff


% --- Executes during object creation, after setting all properties.
function control_asic_onoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to control_asic_onoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14_asic_enable_disable.
function pushbutton14_asic_enable_disable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14_asic_enable_disable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');

contents = cellstr(get(handles.control_asic_onoff,'String'));
curasic = str2num(contents{get(handles.control_asic_onoff,'Value')});

val = TSPMcontrol('read',1);
if isempty(val), disp('Could not read asic status from TSPM'), return, end
bval = dec2bin(val);

NA = Cfg.defaults.Device.nring * Cfg.defaults.Device.nblk * Cfg.defaults.Device.nasic;
% eval(sprintf('devz = sprintf(''%%0%dd'',0);',NA));
devz = sprintf('%048d',0);

if length(devz)<length(bval)
    disp('Error: asics limit exceeded')
    return
end


for k=length(bval):-1:1
    devz(end-k+1) = bval(end-k+1);
end
bval = devz;

if length(bval)<curasic
    disp('Error: selected asic number is not valid')
    return
end

k = length(bval)-curasic+1;
if bval(k)=='0', bval(k)='1';
elseif bval(k)=='1', bval(k)='0';
else
    disp(['Error: unrecognized asic value: ',bval(k),' '])
end

nval = bin2dec(bval);
TSPMcontrol('write',1,nval);


% --- Executes on selection change in image_reconstruction.
function image_reconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to image_reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_reconstruction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_reconstruction


% --- Executes during object creation, after setting all properties.
function image_reconstruction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% Unfiltered back projection method (schlegel & bille 9.1.2)
% modified after Mark Bangert m.bangert@dkfz.de 2011
function BPI = bpg_unfiltered(sinogram,dtheta)
% dtheta - difference beteeen detecotrs thetas in degrees

% figure out how big our picture is going to be.
numOfParallelProjections = size(sinogram,1);
numOfAngularProjections  = size(sinogram,2); 

% convert thetas to radians
thetas = (0:numOfAngularProjections)*(pi/180)*dtheta;

% set up the backprojected image
BPI = zeros(numOfParallelProjections,numOfParallelProjections);

% find the middle index of the projections
midindex = floor(numOfParallelProjections/2) + 1;

% set up the coords of the image
[xCoords,yCoords] = meshgrid(ceil(-numOfParallelProjections/2):ceil(numOfParallelProjections/2-1));

% loop over each projection
for i = 1:numOfAngularProjections

    % figure out which projections to add to which spots
    rotCoords = round(midindex + xCoords*sin(thetas(i)) + yCoords*cos(thetas(i)));

    % check which coords are in bounds
    indices   = find((rotCoords > 0) & (rotCoords <= numOfParallelProjections));
    newCoords = rotCoords(indices);
    
    % summation
    BPI(indices) = BPI(indices) + sinogram(newCoords,i)./numOfAngularProjections;

end
BPI=rot90(BPI);

%%
% filtered back projection in the spatial domain -> schlegel & bille 9.3.1
% modified after Mark Bangert m.bangert@dkfz.de 2011
%
% note: matlab puts the 0 frequency component of a fourier spectrum _not_
%       in the middle. we need to fumble around with fftshift

function BPI = bpg_filtered_spatial(sinogram,dtheta)

% figure out how big our picture is going to be.
numOfParallelProjections = size(sinogram,1);
numOfAngularProjections  = size(sinogram,2); 

% convert thetas to radians
thetas = (0:numOfAngularProjections)*(pi/180)*dtheta;

% set up the backprojected image
BPI = zeros(numOfParallelProjections,numOfParallelProjections);

% find the middle index of the projections
midindex = floor(numOfParallelProjections/2) + 1;

% set up the coords of the image
[xCoords,yCoords] = meshgrid(ceil(-numOfParallelProjections/2):ceil(numOfParallelProjections/2-1));

% set up filter: now for the spatial domain!!!
filterMode = 'ramLak'; % put either 'sheppLogan' or 'ramLak'

if mod(numOfParallelProjections,2) == 0
    halfFilterSize = floor(1 + numOfParallelProjections/2);
else
    halfFilterSize = floor(numOfParallelProjections/2);
end

if strcmp(filterMode,'ramLak')
    filter = zeros(1,halfFilterSize);
    filter(1:2:halfFilterSize) = -1./([1:2:halfFilterSize].^2 * pi^2);
    filter = [fliplr(filter) 1/4 filter];
elseif strcmp(filterMode,'sheppLogan')
    filter = -2./(pi^2 * (4 * (-halfFilterSize:halfFilterSize).^2 - 1) );
end

% loop over each projection
for i = 1:numOfAngularProjections

    % figure out which projections to add to which spots
    rotCoords = round(midindex + xCoords*sin(thetas(i)) + yCoords*cos(thetas(i)));

    % check which coords are in bounds
    indices   = find((rotCoords > 0) & (rotCoords <= numOfParallelProjections));
    newCoords = rotCoords(indices);
        
    % filter
    filteredProfile = conv(sinogram(:,i),filter,'same');

    % summation
    BPI(indices) = BPI(indices) + filteredProfile(newCoords)./numOfAngularProjections;
   
end
BPI=rot90(BPI);


% --- Executes on button press in pushbutton15_recalc_sino.
function pushbutton15_recalc_sino_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15_recalc_sino (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');
Cfg.sindex = CalculateSinogramIndex(Cfg);
setappdata(CurFig,'Cfg',Cfg);


% --- Executes on button press in pushbutton16_load_calibration.
function pushbutton16_load_calibration_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16_load_calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'hPETmonitor_new');
Cfg = getappdata(CurFig,'Cfg');
try, Calibr = getappdata(CurFig,'Calibr'); catch, Calibr=[]; end

% calibrfile = fullfile(calibr_root,sprintf('calibr %s.mat',calibr_timestamp));
[fname,paths] = uigetfile(fullfile(Cfg.folder,'calibr*.mat'),'Select calibration file');
calibrfile = fullfile(paths,fname);
load(calibrfile,'Calibr');
Calibr.ded = compute_diff_dac(Calibr);

setappdata(CurFig,'Calibr',Calibr), drawnow;
Cfg.folder=paths;
setappdata(CurFig,'Cfg',Cfg);

calibr_status = char(regexp(paths,'[^\\]*\\$','match'))
set(handles.calibr_status,'string',calibr_status(1:end-1))

set(handles.calibr_root,'string',paths);
calibr_timestamp = char(regexp(fname,'\d*T\d*','match'));
set(handles.calibr_timestamp,'String',calibr_timestamp);

NA = Cfg.defaults.Device.nblk * Cfg.defaults.Device.nasic;
NC = prod(Cfg.defaults.Device.bmap);

tmp=[]; for k=1:NA,tmp{k}=sprintf('%2d',k); end
set(handles.calibr_asic,'string',tmp);
tmp=[]; for k=1:NC,tmp{k}=sprintf('%2d',k); end
set(handles.calibr_chan,'string',tmp);
drawnow;

NG = length(Calibr.gain);
ND = length(Calibr.dac);

tmp=[]; for k=1:ND,tmp{k}=sprintf('%3d',Calibr.dac(k)); end
set(handles.calibr_dac,'string',tmp);
tmp=[]; for k=1:NG,tmp{k}=sprintf('%3d',Calibr.gain(k)); end
set(handles.calibr_gain,'string',tmp);
drawnow;

calibr_root = char(get(handles.calibr_root,'String'));

try, dg = Calibr.gain(2) - Calibr.gain(1); catch, dg=0; end
try, dd = Calibr.dac(2) - Calibr.dac(1); catch, dd=0; end

set(handles.gainstart,'string',sprintf('%d',min(Calibr.gain)));
set(handles.gainstep,'string',sprintf('%d',dg));
set(handles.gainstop,'string',sprintf('%d',max(Calibr.gain)));

set(handles.dacstart,'string',sprintf('%d',min(Calibr.dac)));
set(handles.dacstep,'string',sprintf('%d',dd));
set(handles.dacstop,'string',sprintf('%d',max(Calibr.dac)));

calibr_show_Callback(hObject, eventdata, handles)


function calibr_timestamp_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_timestamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calibr_timestamp as text
%        str2double(get(hObject,'String')) returns contents of calibr_timestamp as a double


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over calibr_timestamp.
function calibr_timestamp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to calibr_timestamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function varargout = TSPMcontrol(varargin)
% out = TSPMcontrol(command)
% December 7, 2015 Y.D. Sinelnikov
%
% provides low level control for TSMP system
% out is 1 if command succeed and 0 otherwise
% type TSOMcontrol('help'); for a list of available commands

% Read packet, only the address is meaningful, the data should be zero
pktr = uint8(zeros(12,1));
pktr(1) = hex2dec('DE'); % the “deadbeef header”
pktr(2) = hex2dec('AD');
pktr(3) = hex2dec('BE');
pktr(4) = hex2dec('EF');
pktr(5) = hex2dec('00'); % Register address field is a ushort
pktr(6) = hex2dec('00'); % Register
pktr(7) = hex2dec('00'); % Value field is a uint
pktr(8) = hex2dec('00'); % Value
pktr(9) = hex2dec('00'); % Value
pktr(10) = hex2dec('00'); % Value
pktr(11) = hex2dec('FF'); % last 2 bytes are always 0xff
pktr(12) = hex2dec('FF');

% Once the request is sent, the TSPM will respond with an 18 byte packet consisting of the register address, the value and zeroes. So the response to the previous request would come in on 192.168.120.1 port 32002:
pktw = uint8(zeros(18,1));
pktw(1) = hex2dec('00'); % the register number (ushort)
pktw(2) = hex2dec('00');
pktw(3) = hex2dec('00'); % Value field (uint)
pktw(4) = hex2dec('00');
pktw(5) = hex2dec('00');
pktw(6) = hex2dec('00');
pktw(7) = hex2dec('00');
pktw(8) = hex2dec('00');
pktw(9) = hex2dec('00');
pktw(10) = hex2dec('00');
pktw(11) = hex2dec('00');
pktw(12) = hex2dec('00');
pktw(13) = hex2dec('00');
pktw(14) = hex2dec('00');
pktw(15) = hex2dec('00');
pktw(16) = hex2dec('00');
pktw(17) = hex2dec('00');
pktw(18) = hex2dec('00');

SPILENGTH = 408;
TSPMDELAY = 0.2;


%% Init
out = [];
% return value only if asked
if nargout > 0 ; varargout(1) = {out};  end
% if no input argument, goto help
if nargin==0; command = 'help'; else command = lower(varargin{1}); end


%% Switch on commands
switch command
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% read register
    case 'read'
        if nargin~=2, disp('TSPM read requires register value'); return; end
        regnum = uint8(varargin{2});
        if regnum<0 | regnum>40, disp('TSPM registers are 0 to 40'); return; end
        pktr(6) = regnum;
        
        try
            % Open ports
            u1 = udp('192.168.120.1', 'RemotePort',32001);
            u2 = udp('192.168.120.1', 'RemotePort',32002, 'LocalPort',32002, 'Timeout', 1);
            
            fopen(u1);
            fopen(u2);
            
            % Write to 32001 to indicate which register to read
            fwrite(u1,pktr,'uint8');
            
            % Read from 32002 a value
            regval = fread(u2,18,'uint8');
            
            if ~isempty(regval)
                val = 2.^[24 16 8 0]*regval(3:6);
            else
                val = [];
            end
            
            fclose(u2); delete(u2);
            fclose(u1); delete(u1);
            
            out = val;
            
        catch,
            disp('Read failed')
            TSPMcontrol('clean up');
            out = [];
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TSPM reset excluding Ehternet interface
    case 'reset'
        out = TSPMcontrol('write',0,hex2dec('3e'));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TSPM reset including Ehternet interface
    case 'reset all'
        out = TSPMcontrol('write',0,hex2dec('3f'));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% write value into register
    case 'write'
        
        % Write value into the register
        if nargin~=3, disp('TSPM write requires register number (0-40) and value'); return; end
        regnum = uint8(varargin{2});
        regval = uint32(varargin{3});
        
        if regnum<0 | regnum>40, disp('TSPM registers are 0 to 40'); return; end
        
        pktr(6) = regnum;
        val = typecast(regval,'uint8');
        pktr(7) = val(4); pktr(8) = val(3); pktr(9) = val(2); pktr(10) = val(1);
        
        try
            % Open ports
            u0 = udp('192.168.120.1', 'RemotePort',32000);
            
            fopen(u0);
            
            % Write to 32000
            fwrite(u0,pktr,'uint8');
            
            fclose(u0); delete(u0);

            
        catch,
            disp('Write failed')
            TSPMcontrol('clean up');
            out = [];
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% read spi register
    case 'spi read'
        
        if nargin~=2, disp('TSPM SPI read requires register value'); return; end
        regnum = uint8(varargin{2});
        if regnum<3 | regnum>4, disp('TSPM SPI registers are 3 and 4'); return; end
        pktr(6) = regnum;
        k=0;
        
        try
            % Open ports
            u0 = udp('192.168.120.1', 'RemotePort',32000);
            u1 = udp('192.168.120.1', 'RemotePort',32001);
            u2 = udp('192.168.120.1', 'RemotePort',32002, 'LocalPort',32002, 'Timeout', 2);
            
            fopen(u0);
            fopen(u1);
            fopen(u2);
            
            %%
            % reset FIFO by manipulating with register 2
            %         write_value_into_register(Cfg, 2, 2^31);
            %         write_value_into_register(Cfg, 2, 0);
            pktr(6) = uint8(2);
            val = typecast(uint32(2^31),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            val = typecast(uint32(0),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            pause(TSPMDELAY);
            
            %%
            % initiate read cycle by manipulating with regsiter 2
            %         write_value_into_register(Cfg, 2, 2^(spi_reg+5));
            %         write_value_into_register(Cfg, 2, 0);
            val = typecast(uint32(2^(uint32(regnum)+5)),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            val = typecast(uint32(0),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            pause(TSPMDELAY);
            
            %%
            % push data in fifo: make 408 reads of 32 bits words
            pktr(6) = regnum;
            pktr(7:10) = uint8(0);
            for k=1:SPILENGTH,
                %             regval = read_valuer_from_register(Cfg, spi_reg);
                fwrite(u1,pktr,'uint8');
                regval = fread(u2,18,'uint8');
                
                if ~isempty(regval)
                    val(k) = 2.^[24 16 8 0]*regval(3:6);
                else
                    val(k) = NaN;
                end
                
            end
            
            fclose(u2); delete(u2);
            fclose(u1); delete(u1);
            fclose(u0); delete(u0);
            
            out = val;
            
        catch,
            disp(sprintf('SPI read failed on read %d',k))
            TSPMcontrol('clean up');
            out = [];
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% write 408 values into spi register
    case 'spi write'
        
        if nargin~=3, disp('TSPM SPI writes requires register value and an array of values'); return; end
        regnum = uint8(varargin{2});
        if regnum<3 | regnum>4, disp('TSPM SPI registers are 3 and 4'); return; end
        regval = uint32(varargin{3});
        if length(regval)~=SPILENGTH, disp(sprintf('TSPM SPI write requires %d values',SPILENGTH)); return; end
        
        try
            % Open ports
            u0 = udp('192.168.120.1', 'RemotePort',32000);
%             disp('spi write: openning ports')
            fopen(u0);
            pause(TSPMDELAY);
            
            %%
            % reset FIFO by manipulating with register 2
%             disp('spi write: reseting fifo')
            pktr(6) = uint8(2);
            val = typecast(uint32(2^31),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            val = typecast(uint32(0),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            pause(TSPMDELAY);
            
            
            %%
            % push data in fifo: make 408 writes of 32 bits words
%             disp('spi write: pushing data into fifo')
            pktr(6) = uint8(regnum);
            pktr(7:10) = uint8(0);
            for k=1:SPILENGTH,
                %             regval = read_valuer_from_register(Cfg, spi_reg);
                if ~isnan(regval(k))
                    val = typecast(regval(k),'uint8');
                else
                    val = typecast(0,'uint8');
                end
                pktr(7) = val(4); pktr(8) = val(3); pktr(9) = val(2); pktr(10) = val(1);
                fwrite(u0,pktr,'uint8');
                
            end
            pause(TSPMDELAY)
            
            %%
            % initiate SPI write cycle by manipulating with regsiter 2
%             disp('spi write: initating write cycle')
            pktr(6) = uint8(2);
            val = typecast(uint32(2^(uint32(regnum)-3)),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            val = typecast(uint32(0),'uint8');
            pktr(7) = val(4);pktr(8) = val(3);pktr(9) = val(2);pktr(10) = val(1);
            fwrite(u0,pktr,'uint8');
            pause(TSPMDELAY)
            
%             disp('spi write: closing ports')
            fclose(u0); delete(u0);
            
            out = 1;
            
        catch,
            disp('SPI write failed')
            TSPMcontrol('clean up');
            out = [];
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% clean up udp objects
    case 'clean up'
        
        try
            h=instrfind('type','udp','remoteport',32000); n(1) = length(h);
            if ~isempty(h)
                fclose(h);
                delete(h);
            end
            h=instrfind('type','udp','remoteport',32001);  n(2) = length(h);
            if ~isempty(h)
                fclose(h);
                delete(h);
            end
            h=instrfind('type','udp','remoteport',32002);  n(3) = length(h);
            if ~isempty(h)
                fclose(h);
                delete(h);
            end
            
%             if sum(n), disp(sprintf('%d (%d %d %d) udp objects were removed', sum(n),n(1),n(2),n(3))); end
            
            out = 1;
            
        catch,
            disp('Clean up failed')
            out = [];
        end
    otherwise
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('%% TSPMcontrol.m for control of the TSPM system %%');
        disp('%%     See register map firmware version 415    %%');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp(' ');
        disp('Input arguments should be:');
        disp(' ');
        disp(' ''Help''                               to display this');
        disp(' ''Reset''                              to reset system only');
        disp(' ''Reset all''                          to reset system and network connection');
        disp(' ''Read'' <register>                    to read any register 0-40, one scalar value is returned');
        disp(' ''Write'' <register> <value>           to write value into any register 0-40');
        disp(' ''SPI Read'' <register>                to read SPI FIFO register 3-4, 408 values are returned');
        disp(' ''SPI Write'' <register> <values>      to write 408 values into FIFO register 3-4');
        disp(' ''Clean up''                           to clean up hanging udp objects');
        out = [];
        
end

%% End
if nargout > 0
    varargout(1) = {out};  % return value only if asked
end
