function varargout = PendBalistico(varargin)
% PENDBALISTICO MATLAB code for PendBalistico.fig
%      PENDBALISTICO, by itself, creates a new PENDBALISTICO or raises the existing
%      singleton*.
%
%      H = PENDBALISTICO returns the handle to a new PENDBALISTICO or the handle to
%      the existing singleton*.
%
%      PENDBALISTICO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PENDBALISTICO.M with the given input arguments.
%
%      PENDBALISTICO('Property','Value',...) creates a new PENDBALISTICO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PendBalistico_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PendBalistico_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PendBalistico

% Last Modified by GUIDE v2.5 03-Jun-2018 12:23:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PendBalistico_OpeningFcn, ...
                   'gui_OutputFcn',  @PendBalistico_OutputFcn, ...
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


% --- Executes just before PendBalistico is made visible.
function PendBalistico_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PendBalistico (see VARARGIN)
clc
% Choose default command line output for PendBalistico
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PendBalistico wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PendBalistico_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mbala = str2double(get(handles.mBala, 'string'));
mpend = str2double(get(handles.mPend, 'string'));
vinicial = str2double(get(handles.v0Bala, 'string'));

xp = 10;
yp = 5;
yteto = 18; %teto
tetox = [5, 15];
tetoy = [yteto, yteto];
L = yteto-yp; %comprimento do fio
xb = -5; %bala
yb = yp;
g = 9.8;
t= 0;
dt = 0.05;
theta = 0;

while xb < xp
    plot(tetox, tetoy, 'k', 'Linewidth',10); %plot teto
    hold on
    plot([xp xp],[yp yteto], 'k', 'Linewidth', 2); %plot fio
    plot(xp, yp, '.r','markersize', 120); %plot pendulo
    plot(xb, yb,'.','markersize', 20); %plot bala
    hold off
    yb=yp;
    xb=vinicial*t;
    axis ([-10 30 0 30]);
    pause(0.005);
    t = t+0.005;
end

yb = yp;
xb = xp;
v = (mbala/(mbala+mpend))* vinicial %velocidade após a colisão (conservação do momento linear)
theta_f = acos(1-((v*v)/(2*g*L))) %angulo maximo (conservação da enegia mecanica)
vang = sqrt(g/L);
hmax = (v.^2)/(2*g);

set(handles.vApos, 'string', v);
set(handles.h, 'string', hmax);
set(handles.thetamax, 'string', theta_f);

for t = t:dt:50
    theta = theta_f*cos(vang*t-pi/2);
    x = xp + L*sin(theta);
    y = yp + L*(1-cos(theta));
    plot(tetox, tetoy, 'k', 'Linewidth',10); %plot teto
    hold on
    plot([x xp],[y yteto], 'k', 'Linewidth', 2); %plot fio
    plot(x, y, '.r','markersize', 120); %plot pendulo
    plot(x, y,'.','markersize', 20); %plot bala
    hold off
    axis ([-10 30 0 30]);
    pause(dt);
end



function mBala_Callback(hObject, eventdata, handles)
% hObject    handle to mBala (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mBala as text
%        str2double(get(hObject,'String')) returns contents of mBala as a double


% --- Executes during object creation, after setting all properties.
function mBala_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mBala (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mPend_Callback(hObject, eventdata, handles)
% hObject    handle to mPend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mPend as text
%        str2double(get(hObject,'String')) returns contents of mPend as a double


% --- Executes during object creation, after setting all properties.
function mPend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mPend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v0Bala_Callback(hObject, eventdata, handles)
% hObject    handle to v0Bala (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v0Bala as text
%        str2double(get(hObject,'String')) returns contents of v0Bala as a double


% --- Executes during object creation, after setting all properties.
function v0Bala_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v0Bala (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velApos_Callback(hObject, eventdata, handles)
% hObject    handle to velApos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velApos as text
%        str2double(get(hObject,'String')) returns contents of velApos as a double


% --- Executes during object creation, after setting all properties.
function velApos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velApos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function vApos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vApos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uipanel2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
