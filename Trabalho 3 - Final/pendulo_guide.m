function varargout = pendulo_guide(varargin)
% PENDULO_GUIDE MATLAB code for pendulo_guide.fig
%      PENDULO_GUIDE, by itself, creates a new PENDULO_GUIDE or raises the existing
%      singleton*.
%
%      H = PENDULO_GUIDE returns the handle to a new PENDULO_GUIDE or the handle to
%      the existing singleton*.
%
%      PENDULO_GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PENDULO_GUIDE.M with the given input arguments.
%
%      PENDULO_GUIDE('Property','Value',...) creates a new PENDULO_GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pendulo_guide_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pendulo_guide_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pendulo_guide

% Last Modified by GUIDE v2.5 03-Jun-2019 22:38:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pendulo_guide_OpeningFcn, ...
                   'gui_OutputFcn',  @pendulo_guide_OutputFcn, ...
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


% --- Executes just before pendulo_guide is made visible.
function pendulo_guide_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pendulo_guide (see VARARGIN)

% Choose default command line output for pendulo_guide
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pendulo_guide wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pendulo_guide_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function m1_Callback(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
s1 = "Massa: ";
s2 = num2str(get(hObject, 'Value'));
s3 = " kg";
s = s1 + s2 + s3;
set(handles.text4, 'String', s);


% --- Executes during object creation, after setting all properties.
function m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function v1_Callback(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
s1 = "v0: ";
s2 = num2str(get(hObject, 'Value'));
s3 = " m/s";
s = s1 + s2 + s3;
set(handles.text3, 'String', s);


% --- Executes during object creation, after setting all properties.
function v1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
s1 = "Massa: ";
s2 = num2str(get(hObject, 'Value'));
s3 = " kg";
s = s1 + s2 + s3;
set(handles.text9, 'String', s);


% --- Executes during object creation, after setting all properties.
function m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function L_Callback(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L as text
%        str2double(get(hObject,'String')) returns contents of L as a double


% --- Executes during object creation, after setting all properties.
function L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function damping_Callback(hObject, eventdata, handles)
% hObject    handle to damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of damping as text
%        str2double(get(hObject,'String')) returns contents of damping as a double


% --- Executes during object creation, after setting all properties.
function damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eulercromer.
function eulercromer_Callback(hObject, eventdata, handles)
% hObject    handle to eulercromer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Constants

% Pendulum
m2 = get(handles.m2, 'Value');
v2 = 0;
L = str2double(get(handles.L, 'String'));
x = 0; y = -L;

% Bullet
m1 = get(handles.m1, 'Value');
v1 = get(handles.v1, 'Value');
xb = -10;
yb = -L;

% Constants
g = 9.81;
theta0 = acos(-((v1 / (1 + (m2 / m1))) ^2 / (2 * g * L)) + 1);

b = str2double(get(handles.damping, 'String'));
M = m1 + m2;

% Plot Constants

tempo = str2double(get(handles.tempo, 'String'));
nPoints = round(tempo * 1001 / 40);
dt = 0.04;
omega = zeros(nPoints, 1);
theta = zeros(nPoints, 1);
t = zeros(nPoints, 1);
theta(1) = theta0;

%% ODE Solver - Euler-Cromer Method
if b == 0
    for i = 1:nPoints - 1
        omega(i + 1) = omega(i) - (g / L) * theta(i) * dt; 
        theta(i + 1) = theta(i) + omega(i + 1) * dt;
        
        t(i + 1) = t(i) + dt;
    end
else
    for i = 1:nPoints - 1
        omega(i + 1) = omega(i) - (g / L) * theta(i) * dt - b * omega(i) * dt; % Last Expression is for the Damping Strength
        theta(i + 1) = theta(i) + omega(i + 1) * dt;
        
        t(i + 1) = t(i) + dt;
    end
end

% 1001 points equals to 40 secunds
% plot(t, theta, 'r')

%% Animation Plots

axes(handles.animation)

% Bullet Impact
while xb < 0
    tf = (x - xb) / v1;
    dt2 = 0.001;
    
    for t1 = 0:dt2:tf
        xb = xb + v1 * dt2;
        
        plot(xb, yb, '.r', 'MarkerSize', 20)
        title("Simulação do Pendulo")
        hold on
        plot(x, y, '.r', 'MarkerSize', 100)
        line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
        axis([-10 10 -L-2 0])
        drawnow update
        hold off
    end
end

% From (0, 0) to initial position
for dTheta = 0:0.02:theta0
    x = L * sin(dTheta);
    y = - L * cos(dTheta);
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    hold on
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

% Simple Pendulum Animation
for i = 1:nPoints
    x = L * sin(theta(i));
    y = - L * cos(theta(i));
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

%% Theta functino of time
axes(handles.theta_time)
plot(t, theta)
title("Variação de theta ao longo do tempo")
xlabel("t / s")
ylabel("theta / rad")

%% Energys

% Kinetic
axes(handles.kinetic_time)
K = 0.5 * M * (omega .* L) .^2;
plot(t, K)
title("Energia Cinética ao longo do tempo")
xlabel("t / s")
ylabel("Ec / J")

% Potential
axes(handles.potential_time)
h = L * (1 - cos(theta));
U = M * g * h;
plot(t, U)
title("Energia Potencial ao longo do tempo")
xlabel("t / s")
ylabel("Ep / J")

% Mechanical
axes(handles.mechanical_time)
E = K + U;
plot(t, E)
title("Energia Mecânica ao longo do tempo")
xlabel("t / s")
ylabel("Em / J")
function tempo_Callback(hObject, eventdata, handles)
% hObject    handle to tempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tempo as text
%        str2double(get(hObject,'String')) returns contents of tempo as a double


% --- Executes during object creation, after setting all properties.
function tempo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Verlet.
function Verlet_Callback(hObject, eventdata, handles)
% hObject    handle to Verlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Constants

% Pendulum
m2 = get(handles.m2, 'Value');
v2 = 0;
L = str2double(get(handles.L, 'String'));
x = 0; y = -L;

% Bullet
m1 = get(handles.m1, 'Value');
v1 = get(handles.v1, 'Value');
xb = -10;
yb = -L;

% Constants
g = 9.81;
theta0 = acos(-((v1 / (1 + (m2 / m1))) ^2 / (2 * g * L)) + 1);

b = str2double(get(handles.damping, 'String'));
M = m1 + m2;

% Plot Constants

tempo = str2double(get(handles.tempo, 'String'));
nPoints = round(tempo * 1001 / 40);
dt = 0.04;
omega = zeros(nPoints, 1);
theta = zeros(nPoints, 1);
t = zeros(nPoints, 1);
theta(1) = theta0;
theta(2) = theta0 - 0.0001;

%% ODE Solver - Verlet Method
if b == 0
    for i = 1:nPoints - 2
        theta(i + 2) = 2 * theta(i + 1) - theta(i) - dt ^2 * (g / L) * theta(i + 1);
        omega(i + 1) = (theta(i + 2) - theta(i)) / (2 * dt);
        
        t(i + 2) = t(i + 1) + dt;
    end
else
    for i = 1:nPoints - 2
        theta(i + 2) = 2 * theta(i + 1) - theta(i) - dt ^2 * (g / L) * theta(i + 1) - dt ^2 * b * omega(i);  % Last Expression is for the Damping Strength
        omega(i + 1) = (theta(i + 2) - theta(i)) / (2 * dt);
        
        t(i + 2) = t(i + 1) + dt;
    end
end

% 1001 points equals to 40 secunds
% plot(t, theta, 'r')

%% Animation Plots

axes(handles.animation)

% Bullet Impact
while xb < 0
    tf = (x - xb) / v1;
    dt2 = 0.001;
    
    for t1 = 0:dt2:tf
        xb = xb + v1 * dt2;
        
        plot(xb, yb, '.r', 'MarkerSize', 20)
        title("Simulação do Pendulo")
        hold on
        plot(x, y, '.r', 'MarkerSize', 100)
        line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
        axis([-10 10 -L-2 0])
        drawnow update
        hold off
    end
end

% From (0, 0) to initial position
for dTheta = 0:0.02:theta0
    x = L * sin(dTheta);
    y = - L * cos(dTheta);
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    hold on
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

% Simple Pendulum Animation
for i = 1:nPoints
    x = L * sin(theta(i));
    y = - L * cos(theta(i));
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

%% Theta functino of time
axes(handles.theta_time)
plot(t, theta)
title("Variação de theta ao longo do tempo")
xlabel("t / s")
ylabel("theta / rad")

%% Energys

% Kinetic
axes(handles.kinetic_time)
K = 0.5 * M * (omega .* L) .^2;
plot(t, K)
title("Energia Cinética ao longo do tempo")
xlabel("t / s")
ylabel("Ec / J")

% Potential
axes(handles.potential_time)
h = L * (1 - cos(theta));
U = M * g * h;
plot(t, U)
title("Energia Potencial ao longo do tempo")
xlabel("t / s")
ylabel("Ep / J")

% Mechanical
axes(handles.mechanical_time)
E = K + U;
plot(t, E)
title("Energia Mecânica ao longo do tempo")
xlabel("t / s")
ylabel("Em / J")
