function varargout = Calibration2Gains(varargin)
% CALIBRATION2GAINS M-file for Calibration2Gains.fig
%      CALIBRATION2GAINS, by itself, creates a new CALIBRATION2GAINS or raises the existing
%      singleton*.
%
%      H = CALIBRATION2GAINS returns the handle to a new CALIBRATION2GAINS or the handle to
%      the existing singleton*.
%
%      CALIBRATION2GAINS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATION2GAINS.M with the given input arguments.
%
%      CALIBRATION2GAINS('Property','Value',...) creates a new CALIBRATION2GAINS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Calibration2Gains_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Calibration2Gains_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Calibration2Gains

% Last Modified by GUIDE v2.5 02-Sep-2018 13:59:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Calibration2Gains_OpeningFcn, ...
    'gui_OutputFcn',  @Calibration2Gains_OutputFcn, ...
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


% --- Executes just before Calibration2Gains is made visible.
function Calibration2Gains_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Calibration2Gains (see VARARGIN)

% Choose default command line output for Calibration2Gains
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Calibration2Gains wait for user response (see UIRESUME)
% uiwait(handles.Calibr2Gains);
CurFig = hObject; pause(0.1);
setappdata(0,'SynchroPETCalibr2Gains',CurFig);
Cfg.folder = pwd;
setappdata(CurFig,'Cfg',Cfg);


% --- Outputs from this function are returned to the command line.
function varargout = Calibration2Gains_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_select.
function pushbutton1_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

FilterSpec = 'calibr*.mat';

try  Cfg.folder; catch, Cfg.folder=pwd; end

[InputFile, InputPath, tmp] = uigetfile(FilterSpec,'Select data files',Cfg.folder,'MultiSelect','on');
if tmp==0, % cancel
    set(hObject,'Value',0);
    return
end
Cfg.folder = InputPath;
setappdata(CurFig,'Cfg',Cfg);

if iscell(InputFile)
    BatchNo = length(InputFile);
else
    BatchNo = 1; tmp = InputFile; clear InputFile; InputFile{1} = tmp;
end

set(handles.calibr_files,'String',InputFile);
set(handles.calibr_files,'Value',1);
drawnow

% --- Executes on button press in pushbutton2_gains.
function pushbutton2_gains_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_gains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

calibr_global_max = str2double(get(handles.calibr_plot_max,'string'));

FileNames = get(handles.calibr_files,'string');
NumberOfFiles = length(FileNames);
curk = get(handles.calibr_files,'value');

for k = curk:NumberOfFiles;
    set(handles.calibr_files,'Value',k);
    drawnow
    process_ratcap_calibration_file(Cfg.folder,char(FileNames{k}),calibr_global_max)
end
set(handles.calibr_files,'value',curk)

% --- Executes on selection change in calibr_files.
function calibr_files_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calibr_files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calibr_files


% --- Executes during object creation, after setting all properties.
function calibr_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% top level call to process calibration of ratcap
function process_ratcap_calibration_file(FilePath,FileName,calibr_global_max)
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

Calibr = function_process_calibration_file(FilePath, FileName, calibr_global_max);

if length(Calibr.gain) == 1
    return
end
Calibr = getappdata(CurFig,'Calibr');

figure(11), subplot(221), hist(Calibr.dDdG(:)), grid
xlabel('dDac/dGain'), title(FileName(1:end-4))
subplot(222), hist(Calibr.Doffset(:)), grid
xlabel('Dac offset'), title(datestr(now,31))
subplot(223), imagesc(Calibr.dDdG), colorbar
xlabel('Asic #'), ylabel('Channel #')
subplot(224), imagesc(Calibr.Doffset), colorbar
xlabel('Asic #'),

drawnow, pause(0.1)
saveas(gcf,fullfile(FilePath,[FileName(1:end-4),' hist.jpg']));

%%
% do actual processing and figures
function Calibr = function_process_calibration_file(FilePath, FileName, calibr_global_max)
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

Calibr = process_calibration_threshold_scan(FilePath,FileName);
setappdata(CurFig,'Calibr',Calibr);

pushbutton3_repeat_data_Callback([], [], []);



%%
% detect photo peak line
%%
%%
function xy = line_detect_function(a)

[NX NY] = size(a);
% [y x] = meshgrid(1:NY,1:NX);


[H T R] = hough(a);

% Display the original image.
% subplot(2,1,1);
% imshow(a);
% title('Image');

% Display the Hough matrix.
% subplot(2,1,2);
% imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%       'InitialMagnification','fit');
% title('Hough Transform');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(hot);

P  = houghpeaks(H,1);
% hold on;
% plot(T(P(2)),R(P(1)),'s','color','white');
% hold off

% [T(P(2)),R(P(1))]

for k=2:32
    MinLength = round(NY/k);
    FillGap = round(MinLength/sqrt(2));
    
    lines = houghlines(a,T,R,P,'FillGap',FillGap,'MinLength',MinLength);
    try
        xy = [lines(1).point1; lines(1).point2];
        break;
    catch
        %         [MinLength, FillGap]
        xy=[];
    end
end
% subplot(2,1,1)
% hold on
% plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% set(gca,'ydir','normal')
% axis normal tight
% hold off

%%
% process calibration
function Calibr = process_calibration_threshold_scan(FilePath, FileName)

GainCalcParam = get_gain_calculation_parameters();

calibrfile = fullfile(FilePath, FileName);
load(calibrfile,'Calibr');

Calibr.ded = compute_diff_dac(Calibr);

[ND,NG,NC,NA] = size(Calibr.evr);

Calibr.pfit=[];

for na=1:NA,
    for nc=1:NC,
        for ng=1:NG,
            x = Calibr.dac;
            y = Calibr.ded(:,ng,nc,na);
            
            [p1, r, xcalc,ycalc] = fit_one_detector(x,y,na,nc,Calibr.gain(ng));
            try
                Calibr.pfit(:,ng,nc,na) = [p1, r];
            catch
                p1
                r
                disp(sprintf('Asic %d channel %d: Amplitude %5.2f Sigma %5.2f Offset %5.2f',na,nc,Calibr.pfit(1,ng,nc,na),Calibr.pfit(2,ng,nc,na),Calibr.pfit(3,ng,nc,na)))
            end
        end
    end
end



% datasetname={'Amplitude','FWHM','Shift','Background','Regression fit'};
%
% nsb = ceil(sqrt(NG));
% for n=1:5,
%     figure(n),
%     set(gcf,'name',datasetname{n})
%     for k=1:NG,
%         subplot(nsb,nsb,k)
%         imagesc( shiftdim(Calibr.pfit(k,n,:,:),2) )
%     end
% end

function GainCalcParam = get_gain_calculation_parameters()
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

GainCalcParam.PhotoPeakmV = 511;
GainCalcParam.OperatingLLDmV = 350;

handle=findobj(CurFig,'tag','target_dac');
GainCalcParam.Dac2Go = str2double(get(handle,'string'));

handle=findobj(CurFig,'tag','const_ksi');
GainCalcParam.ksi = str2double(get(handle,'string'));
handle=findobj(CurFig,'tag','const_mv');
GainCalcParam.mv = str2double(get(handle,'string'));
handle=findobj(CurFig,'tag','const_dg');
GainCalcParam.dg = str2double(get(handle,'string'));
handle=findobj(CurFig,'tag','const_B');
GainCalcParam.B = str2double(get(handle,'string'));
handle=findobj(CurFig,'tag','const_T');
GainCalcParam.T = str2double(get(handle,'string'));

function devrddac = compute_diff_dac(Calibr)

[ND,NG,NC,NA] = size(Calibr.evr);

if ND>5,
    for na=1:NA,
        for nc=1:NC,
            for ng=1:NG
                
                cevr = smooth(Calibr.evr(:,ng,nc,na),10);
                cdrv = smooth([0; diff(cevr);],5);
                devrddac(:,ng,nc,na) = cdrv;
            end
        end
    end
else
    devrddac = [];
end

function calibr_show_Callback(Calibr,calibr_show,dgca,calibr_global_max)


if isempty(Calibr), return, end

nd = dgca(1);
ng = dgca(2);
nc = dgca(3);
na = dgca(4);

switch calibr_show
    case 1
        imagesc(Calibr.gain,Calibr.dac,log(Calibr.tot))
        set(gca, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 2
        imagesc(Calibr.gain,Calibr.dac, sum(Calibr.evr(:,:,:,na),3))
        set(gca, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 3
        imagesc(Calibr.gain,Calibr.dac, sum(Calibr.evr(:,:,nc,:),4))
        set(gca, 'ydir', 'reverse')
        axis tight
        xlabel('Gains'), ylabel('Dacs');
    case 4
        imagesc(Calibr.gain,Calibr.dac, Calibr.evr(:,:,nc,na))
        set(gca, 'ydir', 'reverse')
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
        mevr = mean(Calibr.evr(:));
        sevr = std(Calibr.evr(:));
        disp(sprintf('Dataset has mean=%6.1f std=%6.1f',mevr, sevr))
        for k=1:32,
            subplot(4,8,k)
            imagesc(Calibr.gain,Calibr.dac, Calibr.evr(:,:,k,na))
            set(gca, 'ydir', 'reverse','ytick',[],'xtick',[],'clim',[0 calibr_global_max])
        end
    case 8
        FigureTag = '32 detectors diff'; hfig=findobj('Tag',FigureTag);
        if isempty(hfig), hfig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off'); end;
        sfigure(hfig)
        mevr = mean(Calibr.ded(:));
        sevr = std(Calibr.ded(:));
        disp(sprintf('Derivative has mean=%6.1f std=%6.1f',mevr, sevr))
        for k=1:32,
            subplot(4,8,k)
            imagesc(Calibr.gain,Calibr.dac, Calibr.ded(:,:,k,na))
            set(gca, 'ydir', 'reverse','ytick',[],'xtick',[],'clim',[0 calibr_global_max])
        end
    case 9
        FigureTag = 'All detectors diff'; hfig=findobj('Tag',FigureTag);
        if isempty(hfig), hfig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off'); end;
        sfigure(hfig)
        mevr = mean(Calibr.ded(:));
        sevr = std(Calibr.ded(:));
        disp(sprintf('Derivative has mean=%6.1f std=%6.1f',mevr, sevr))
        for na=1:12,
            for k=1:32,
                subplot(12,32,(na-1)*32+k)
                imagesc(Calibr.gain,Calibr.dac, Calibr.ded(:,:,k,na))
                set(gca, 'ydir', 'reverse','ytick',[],'xtick',[],'clim',[0 calibr_global_max])
            end
        end
        
    otherwise
        disp('TBD')
end
colormap bone

function sfigure(h)
if nargin>=1
    if ishandle(h), set(0, 'CurrentFigure', h);
    else h = figure(h); end %#ok<*NASGU>
else h = figure;
end

function [p1, r, xcalc, ycalc] = fit_one_detector(x,y,na,nc, gain)

% coinc_peak_bkg
% p(1) - gaussian amplitude
% p(2) - gaussian sigma
% p(3) - gaussian shift
% p(4) - background time shift
% p(5) - background level to the left
% p(6) - background to the right to the right

% valid region
ivalid = find(x>=100 & x<=400);
x = x(ivalid); xcalc=x;
y = y(ivalid);

y0 = max(y); y0b = mean(y);
x0 = min(x(find(y==y0)));


p0 = [y0   30     x0       y0b       y0b  ];
lb = [y0/2 10     x0-100  -abs(y0b) -abs(y0b) ];
ub = [3*y0 max(x) x0+100   abs(y0b)  abs(y0b)  ];

options = optimset('Display','off','TolX',1e-4,'TolFun',1e-8,'MaxIter',1e10,'MaxFunEvals',1e10);
try
    [p1, r] = lsqnonlin(@coinc_peak_bkg_5arg,p0,lb,ub,options, x, y);
    if isempty(r), r=NaN; end
catch
    p1 = [0 1 0 0 0]; r= NaN;
    disp(sprintf('Warning: Photo-peak for A %d C %d gain %d is not located',na,nc,gain))
end
ycalc = coinc_peak_bkg_5arg(p1,x,0);

% ycalc0 = coinc_peak_bkg_5arg(p0,x,0);
% figure(11), plot(x,y,x,ycalc,x,ycalc0), legend('orig','fit','0'), grid
% drawnow
% 
% [p0;lb;ub;p1] 
% pause;



function yy = coinc_peak(p,x,y)
% p(1) - gaussian amplitude
% p(2) - gaussian sigma
% p(3) - gaussian shift
% p(4) - background level

yy = p(1)*gaussmf(x,[p(2) p(3)]) + p(4) - y;

function yy = coinc_peak_bkg(p,x,y)
% p(1) - gaussian amplitude
% p(2) - gaussian sigma
% p(3) - gaussian shift
% p(4) - background time shift
% p(5) - background level to the left
% p(6) - background to the right to the right

bkg = p(5)*(x<=p(4)) + p(6)*(x>p(4));
yy = p(1)*gaussmf(x,[p(2) p(3)]) + bkg - y;

function yy = coinc_peak_bkg_5arg(p,x,y)
% p(1) - gaussian amplitude
% p(2) - gaussian sigma
% p(3) - gaussian shift
% p(4) - background level to the left
% p(5) - background to the right to the right

yy = p(1)*gaussmf(x,[p(2) p(3)]) + p(4)*(x<=p(3)) + p(5)*(x>p(3)) - y;

function process_verification_scan(Calibr,Gcp,calibr_global_max,curfullfilename)
BaselineDAC = Gcp.B*Gcp.mv;

% single dac scan
ng  = 1;
damax = max(Calibr.ded(:));
aamax = mean(Calibr.evr(:))+3*std(Calibr.evr(:));

damax = calibr_global_max;
% aamax = 3e3;

[ND,NG,NC,NA] = size(Calibr.evr);

SBY = round(NA/12);
SBX = NA/SBY;
for na=1:NA,
    
    aa = shiftdim(Calibr.evr(:,ng,:,na),2);
    da = shiftdim(Calibr.ded(:,ng,:,na),2);

    hf1 = figure(1);
    subplot(SBX,SBY,na)
    plot(Calibr.dac,da), grid, axis tight, 
    set(gca,'ylim',[0 damax],'xtick',[],'ytick',[])
    text(mean(Calibr.dac),damax/2,sprintf('%d',na))
    
    hf2 = figure(2);
    subplot(SBX,SBY,na)
    plot(Calibr.dac,aa), grid, axis tight, set(gca,'ylim',[0 aamax])
    set(gca,'ylim',[0 aamax],'xtick',[],'ytick',[])
    text(mean(Calibr.dac),damax/2,sprintf('%d',na))
    
end
saveas(hf1,[curfullfilename(1:end-4),' photo peaks.jpg'])
saveas(hf2,[curfullfilename(1:end-4),' raw counts.jpg'])

x = Calibr.dac;
y = mean(mean(Calibr.ded,3),4);

[p1, r, xcalc, ycalc] = fit_one_detector(x,y,0,0,0);

EnergyRatio = 2*sqrt(2*log(2))*p1(2)/(BaselineDAC-p1(3));


xkEv = Gcp.PhotoPeakmV*(BaselineDAC-x)/(BaselineDAC-p1(3));
ixcutoff = min(find(abs(xkEv-Gcp.OperatingLLDmV)==min(abs(xkEv-Gcp.OperatingLLDmV))));
xcutoff = x(ixcutoff);
ycutoff = y(ixcutoff);

hf3 = figure(3); plot(x,y,'-b',xcalc,ycalc,'-g',xcutoff,ycutoff,'pr','MarkerSize',14,'LineWidth',2), grid, axis tight
xlabel('Dac (a.u.)','fontsize',14), ylabel('dCounts/dDac','fontsize',14)
legend('Mean','Fit')
title(sprintf('Energy Resolution=%3.0f%%  DAC 511kEv=%4.0f 350kEv=%4.0f',100*EnergyRatio,p1(3),xcutoff),'fontsize',14)
set(gca,'fontsize',14)
saveas(hf3,[curfullfilename(1:end-4),' energy resolution.jpg'])

if NG==1,
    hf4 = figure(4); 
    subplot(221), imagesc(shiftdim(Calibr.pfit(1,1,:,:),2)), colorbar, title('Amplitude')
    subplot(222), imagesc(shiftdim(Calibr.pfit(2,1,:,:),2)), colorbar, title('Sigma')
    subplot(223), imagesc(shiftdim(Calibr.pfit(3,1,:,:),2)), colorbar, title('Offset')
    subplot(224), imagesc(log(shiftdim(Calibr.pfit(6,1,:,:),2))), colorbar, title('log(R)')
    saveas(hf4,[curfullfilename(1:end-4),' statistics.jpg'])
    
    if 0 % temporarily deisabled
        gaincorr = calculate_gain_corrections(shiftdim(Calibr.pfit(3,1,:,:),2), Gcp);
        hf5 = figure(5); imagesc(gaincorr), colorbar, title(sprintf('Gain correction to dac %4.1f', Gcp.Dac2Go))
        saveas(hf5,[curfullfilename(1:end-4),' gain corrections.jpg'])
        
        x=reshape(gaincorr,prod(size(gaincorr)),[]);
        
        [File2read, Path2read,tmp] = uigetfile([pwd,'/*.mat'],'Select gain file to apply corrections');
        if File2read==0, disp('New gains were not generated'), return, end
        
        gfn = fullfile(Path2read, File2read);
        q = load(gfn);
        
        xx=round(q.data(:,4)+x);
        for k=1:length(xx)
            if xx(k)<=0, xx(k)=0; end
            if xx(k)>31, xx(k)=31; end
        end
        data = q.data;
        data(:,4)= xx;
        
        [a b c] = fileparts(gfn);
        newfilename=[b,' CORRECTED ', datestr(now,30), c];
        newfile = fullfile(a,newfilename);
        save(newfile,'data');
        disp(['New gains save in ',newfile])
    end
end

function gaincorr = calculate_gain_corrections(x,Gcp)
% assumes all gains were at zero
% x - 2D matrix of dac offsets
% x0 - desired dac setting

d2 = Gcp.Dac2Go/Gcp.mv;
d1 = x/Gcp.mv;

gaincorr = ((Gcp.B-d2)./(Gcp.B - (d1-Gcp.ksi)) - 1)/Gcp.dg;



function target_dac_Callback(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text1 as text
%        str2double(get(hObject,'String')) returns contents of text1 as a double


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function const_ksi_Callback(hObject, eventdata, handles)
% hObject    handle to const_ksi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of const_ksi as text
%        str2double(get(hObject,'String')) returns contents of const_ksi as a double


% --- Executes during object creation, after setting all properties.
function const_ksi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to const_ksi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function const_mv_Callback(hObject, eventdata, handles)
% hObject    handle to const_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of const_mv as text
%        str2double(get(hObject,'String')) returns contents of const_mv as a double


% --- Executes during object creation, after setting all properties.
function const_mv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to const_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function const_dg_Callback(hObject, eventdata, handles)
% hObject    handle to const_dg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of const_dg as text
%        str2double(get(hObject,'String')) returns contents of const_dg as a double


% --- Executes during object creation, after setting all properties.
function const_dg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to const_dg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function const_B_Callback(hObject, eventdata, handles)
% hObject    handle to const_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of const_B as text
%        str2double(get(hObject,'String')) returns contents of const_B as a double


% --- Executes during object creation, after setting all properties.
function const_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to const_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function const_T_Callback(hObject, eventdata, handles)
% hObject    handle to const_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of const_T as text
%        str2double(get(hObject,'String')) returns contents of const_T as a double


% --- Executes during object creation, after setting all properties.
function const_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to const_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function target_dac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_dac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function calibr_plot_max_Callback(hObject, eventdata, handles)
% hObject    handle to calibr_plot_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calibr_plot_max as text
%        str2double(get(hObject,'String')) returns contents of calibr_plot_max as a double


% --- Executes during object creation, after setting all properties.
function calibr_plot_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibr_plot_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3_repeat_data.
function pushbutton3_repeat_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_repeat_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'SynchroPETCalibr2Gains');
Cfg = getappdata(CurFig,'Cfg');

try, Calibr = getappdata(CurFig,'Calibr');
catch, return; end

calibr_global_max = str2double(get(findobj(0,'tag','calibr_plot_max'),'string'));

FileNames = get(findobj(0,'tag','calibr_files'),'string');
NumberOfFiles = length(FileNames);
curk = get(findobj(0,'tag','calibr_files'),'value');
FileName = char(FileNames{curk});
FilePath = Cfg.folder;

GainCalcParam = get_gain_calculation_parameters();
[ND,NG,NC,NA] = size(Calibr.evr);

if length(Calibr.gain) == 1
    process_verification_scan(Calibr,GainCalcParam,calibr_global_max,fullfile(FilePath,FileName));
    return
end

dDdG_have=[];
Doffset_have=[];

gavg = mean(Calibr.ded(:));
gstd = std(Calibr.ded(:));
% calibr_global_max = gavg + 1.4*gstd;
DeadDetectors = 0;
disp(sprintf('### Analysing gain dac dependence: mean=%6.1f std=%6.1f', gavg, gstd))
for na=1:NA,
    
  
    a = Calibr.ded(:,:,:,na);
    cavg = mean(a(:));
    cstd = std(a(:));
    calibr_global_max = cavg + cstd;
    for k=1:NC,

        handles = guidata(CurFig);
        ALL_MANUAL = get(handles.chk_all_manual,'value');
        drawnow

        a = Calibr.ded(:,:,k,na)';
        cavg = mean(a(:));
        cstd = std(a(:));
        
        figure(1), imagesc(Calibr.dac,Calibr.gain,a), colorbar
        set(gca, 'ydir', 'normal', 'xdir', 'reverse','clim',[0 calibr_global_max])
        title(sprintf('%s asic %d chanel %d',FileName(1:end-4),na,k))
        hold on
        
        
        %         if k==28 & na==7
        %             q=2;
        %         end
        
        if abs(cavg) < abs(gavg)/1e2 | cstd < gstd/1e2
            disp(sprintf('Skip asic %d chanel %d - no data',na, k))
            Calibr.dDdG(k,na) = NaN;
            Calibr.Doffset(k,na) = NaN;
        else
            disp(sprintf('Asic %d chanel %d mean = %6.1f (%6.1f)',na,k,cavg, cstd))
            
            if ALL_MANUAL
                xy=[];
            else
                xy = line_detect_function(a>(cavg+cstd));
            end
            
            if ~isempty(xy)
                plot(Calibr.dac(xy(:,1)),Calibr.gain(xy(:,2)),'LineWidth',2,'Color','yellow');
                
                dDdG = diff(Calibr.dac(xy(:,1)))/diff(Calibr.gain(xy(:,2)));
                Doffset = mean( Calibr.dac(xy(:,1))'- dDdG*Calibr.gain(xy(:,2)) );
                if length(dDdG_have)>3
                    mm_dddg(1) = min(dDdG_have);     mm_dddg(2) = max(dDdG_have);
                    mm_doff(1) = min(Doffset_have);  mm_doff(2) = max(Doffset_have);
                    if dDdG<mm_dddg(1) | dDdG>mm_dddg(2) | Doffset<mm_doff(1) | Doffset>mm_doff(2)
                        disp('### Warning: statistically suspicious data')
                        disp(sprintf('### dDdG=%4.1f [%4.1f %4.1f]',dDdG,mm_dddg(1),mm_dddg(2)))
                        disp(sprintf('### Doffset=%4.0f [%4.0f %4.0f]',Doffset,mm_doff(1),mm_doff(2)))
                        BE_EXTRA_CAREFUL = 1;
                        if BE_EXTRA_CAREFUL
                            dDdG=0; % forse manual processing
                        end
                    end
                else
                    dDdG=0;
                end
            else
                disp(sprintf('About to skip asic %d chanel %d - NO line',na, k))
                Calibr.dDdG(k,na) = NaN;
                Calibr.Doffset(k,na) = NaN;
                dDdG = 0; % forse manual processing anyway
            end
            
            ManualPause = 0;
            IsThisDeadDetector = 0;
            while dDdG==0 | isinf(dDdG) 
                disp('Problem identifying photo-peak line: click two points starting from bottom left')
                dg = ginput(2);
                if isempty(dg)
                    IsThisDeadDetector = 1;
                    break
                else
                    xy(2,1) = min(find(abs(Calibr.dac-dg(1,1))==min(abs(Calibr.dac-dg(1,1)))));
                    xy(1,1) = min(find(abs(Calibr.dac-dg(2,1))==min(abs(Calibr.dac-dg(2,1)))));
                    xy(2,2) = min(find(abs(Calibr.gain-dg(1,2))==min(abs(Calibr.gain-dg(1,2)))));
                    xy(1,2) = min(find(abs(Calibr.gain-dg(2,2))==min(abs(Calibr.gain-dg(2,2)))));
                    dDdG = diff(Calibr.dac(xy(:,1)))/diff(Calibr.gain(xy(:,2)));
                    ManualPause = 1;
                end
            end
            
            if IsThisDeadDetector
                Calibr.dDdG(k,na) = NaN;
                Calibr.Doffset(k,na) = NaN;
            else
                Doffset = mean( Calibr.dac(xy(:,1))'- dDdG*Calibr.gain(xy(:,2)) );
                plot(Calibr.dac(xy(:,1)),Calibr.gain(xy(:,2)),'LineWidth',2,'Color','green');
                xlabel('Dacs'), ylabel('Gains')
                text(500,5,sprintf('D = %5.2f*G + %4.1f',dDdG,Doffset),'color','w','fontsize',18)
                hold off
                
                drawnow
                
                if ManualPause, pause(0.5), end
                
                if dDdG>=0,
                    Calibr.dDdG(k,na) = NaN;
                    Calibr.Doffset(k,na) = NaN;
                else
                    Calibr.dDdG(k,na) = dDdG;
                    Calibr.Doffset(k,na) = Doffset;
                end
                
                dDdG_have = [dDdG_have, dDdG];
                Doffset_have = [Doffset_have, Doffset];
            end
            
            if isnan(Calibr.dDdG(k,na))
                DeadDetectors = DeadDetectors+1;
            end
            
        end
    end
end
disp(sprintf('Number of dead detectors = %d out of %d',DeadDetectors,NA*NC))
computed_gains = [FileName(1:end-4),sprintf(' %s linear equation coeff.mat',datestr(now,30))];
save(fullfile(FilePath,computed_gains),'Calibr');

disp('### Analysing gains')
% now find dac value at which gains are positive and addup to minimum?
cnt=0;
for k=1:length(Calibr.dac)
    gains = (Calibr.dac(k)-Calibr.Doffset)./Calibr.dDdG;
    if isempty(find(gains(:)<0)) & isempty(find(gains(:)>31))
        cnt=cnt+1;
        
        cost_nan(cnt) = sum(sum(isnan(gains)));
        cost_sum(cnt) = sum(gains(find(isnan(gains)==0)));
        
%         gains(isnan(gains)==1) = mean(gains(isnan(gains)==0));
        gains(isnan(gains)==1) = 0; % 10/8/17 change
        
        best_gains(:,:,cnt) = round(gains);
        best_dac(cnt) = Calibr.dac(k);
        
    end
end
if cnt>0
    TOLERATE_NAN=100; % increase from 13 on 10/8/17 when calibrating ring 11
    disp(sprintf('Saving gain files for %s',FileName(1:end-4)))
    for dac = min(best_dac):10:max(best_dac)
        k = find( abs(best_dac-dac)==min(abs(best_dac-dac)) );
        
        if cost_nan(k(1))<=TOLERATE_NAN
            gains = best_gains(:,:,k(1));
            x=reshape(gains,prod(size(gains)),[]);
            data=[x*0,x*0,x*0,x];
            
            gain_file_name = [FileName(1:end-4),sprintf('dac %03d gains %02d-%02d.mat',round(dac),min(gains(:)),max(gains(:)))];
            disp(sprintf('Writing gains for dac %4.1f in a file %s',dac,gain_file_name));
            save(fullfile(FilePath,gain_file_name),'data');
        else
            disp(sprintf('Skipping gains for dac %4.1f: too many missing detectors %d',dac,cost_nan(k(1))));
        end
    end
else
    disp(sprintf('No valid gains exist for %s',FileName(1:end-4)))
end
setappdata(CurFig,'Calibr',Calibr);


% --- Executes on button press in chk_all_manual.
function chk_all_manual_Callback(hObject, eventdata, handles)
% hObject    handle to chk_all_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_all_manual
