function varargout = guide(varargin)
% GUIDE MATLAB code for guide.fig
%      GUIDE, by itself, creates a new GUIDE or raises the existing
%      singleton*.
%
%      H = GUIDE returns the handle to a new GUIDE or the handle to
%      the existing singleton*.
%
%      GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE.M with the given input arguments.
%
%      GUIDE('Property','Value',...) creates a new GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guide_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guide_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guide

% Last Modified by GUIDE v2.5 07-Jun-2019 16:02:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guide_OpeningFcn, ...
                   'gui_OutputFcn',  @guide_OutputFcn, ...
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


% --- Executes just before guide is made visible.
function guide_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guide (see VARARGIN)

% Choose default command line output for guide
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guide wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guide_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Tf_Callback(hObject, eventdata, handles)
% hObject    handle to Tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
s1 = "Tf = ";
Tf = get(hObject, 'Value');
s2 = num2str(Tf);
s = s1 + s2;
set(handles.tf_text, 'String', s);

v1 = "2pi^2/Tf^2 = ";
v2 = num2str(2 * pi * pi / (Tf * Tf));
v = v1 + v2;
set(handles.tfb_text, 'String', v);


% --- Executes during object creation, after setting all properties.
function Tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in animacao.
function animacao_Callback(hObject, eventdata, handles)
% hObject    handle to animacao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m = 1;
L = 1;
K = 2; % O user não sabe
g = 10;

n = 10; % Número de períodos
Tf = get(handles.Tf, 'Value');
ti = 0;
dt = 0.1;
tf = Tf * n;

t = ti:dt:tf;
F = cos(2 * pi * t / Tf) / 5;

%% Euler-Cromer

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

for i = 1:length(t) - 1
    if (xb(i) <= 0 && vb(i) < 0) || (xb(i) >= 2 && vb(i) < 0)
        vb(i) = -vb(i);
    end
    
    vb(i + 1) = vb(i) + dt * Fr(K, xb(i), F(i)) / m;
    xb(i + 1) = xb(i) + dt * vb(i + 1);
end

xb_ec = xb;
vb_ec = vb;

%% Verlet

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

xb(2) = xb(1) + vb(1) * dt;

for i = 2:length(t) - 1
    xb(i + 1) = 2 * xb(i) - xb(i - 1) + dt ^2 * Fr(K, xb(i), F(i)) / m;
    vb(i) = (xb(i + 1) - xb(i - 1)) / (2 * dt);
    
    if (xb(i + 1) >= 2 && vb(i) > 0) || (xb(i + 1) <= 0 && vb(i) < 0)
        vb(i) = - vb(i);
        xb(i + 1) = xb(i);
    end
end

xb_v = xb;
vb_v = vb;

%% Plots

axes(handles.p_animacao)
for i = 1:length(t)
    plot(xb_ec(i), 0, 'or', 'MarkerSize', 20)
    hold on
    line([0 0], [-1 1])
    line([2 2], [-1 1])
    line([0 xb_ec(i)], [0 0])
    line([xb_ec(i) 2], [0 0])
    
    plot(xb_v(i), 0, 'ob', 'MarkerSize', 20)
    line([0 xb_v(i)], [0 0])
    line([xb_v(i) 2], [0 0])
    hold off
    drawnow
end
% --- Executes on button press in x_t.
function x_t_Callback(hObject, eventdata, handles)
% hObject    handle to x_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m = 1;
L = 1;
K = 2; % O user não sabe
g = 10;

n = 10; % Número de períodos
Tf = get(handles.Tf, 'Value');
ti = 0;
dt = 0.1;
tf = Tf * n;

t = ti:dt:tf;
F = cos(2 * pi * t / Tf) / 5;

%% Euler-Cromer

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

for i = 1:length(t) - 1
    if (xb(i) <= 0 && vb(i) < 0) || (xb(i) >= 2 && vb(i) < 0)
        vb(i) = -vb(i);
    end
    
    vb(i + 1) = vb(i) + dt * Fr(K, xb(i), F(i)) / m;
    xb(i + 1) = xb(i) + dt * vb(i + 1);
end

xb_ec = xb;
vb_ec = vb;

%% Verlet

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

xb(2) = xb(1) + vb(1) * dt;

for i = 2:length(t) - 1
    xb(i + 1) = 2 * xb(i) - xb(i - 1) + dt ^2 * Fr(K, xb(i), F(i)) / m;
    vb(i) = (xb(i + 1) - xb(i - 1)) / (2 * dt);
    
    if (xb(i + 1) >= 2 && vb(i) > 0) || (xb(i + 1) <= 0 && vb(i) < 0)
        vb(i) = - vb(i);
        xb(i + 1) = xb(i);
    end
end

xb_v = xb;
vb_v = vb;

%% Plots

axes(handles.p_xt)
cla
plot(t, xb_ec)
hold on
plot(t, xb_v)
xlabel("t/s")
ylabel("x/m")
title("Posição em função do tempo")
legend("Euler-Cromer x", "Verlet x")

% --- Executes on button press in em.
function em_Callback(hObject, eventdata, handles)
% hObject    handle to em (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

m = 1;
L = 1;
K = 2; % O user não sabe
g = 10;

n = 10; % Número de períodos
Tf = get(handles.Tf, 'Value');
ti = 0;
dt = 0.1;
tf = Tf * n;

t = ti:dt:tf;
F = cos(2 * pi * t / Tf) / 5;

%% Euler-Cromer

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

for i = 1:length(t) - 1
    if (xb(i) <= 0 && vb(i) < 0) || (xb(i) >= 2 && vb(i) < 0)
        vb(i) = -vb(i);
    end
    
    vb(i + 1) = vb(i) + dt * Fr(K, xb(i), F(i)) / m;
    xb(i + 1) = xb(i) + dt * vb(i + 1);
end

xb_ec = xb;
vb_ec = vb;

%% Verlet

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

xb(2) = xb(1) + vb(1) * dt;

for i = 2:length(t) - 1
    xb(i + 1) = 2 * xb(i) - xb(i - 1) + dt ^2 * Fr(K, xb(i), F(i)) / m;
    vb(i) = (xb(i + 1) - xb(i - 1)) / (2 * dt);
    
    if (xb(i + 1) >= 2 && vb(i) > 0) || (xb(i + 1) <= 0 && vb(i) < 0)
        vb(i) = - vb(i);
        xb(i + 1) = xb(i);
    end
end

xb_v = xb;
vb_v = vb;


%% Energia Mecânica

em_ec = 0.5 * m * vb_ec .^2;
em_v = 0.5 * m * vb_v .^2;

axes(handles.p_xt)
cla
plot(t, em_ec)
hold on
plot(t, em_v)
legend("Euler-Cromer Em", "Verlet Em")


% --- Executes during object creation, after setting all properties.
function tf_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function tfb_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfb_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
