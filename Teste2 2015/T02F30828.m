function varargout = T02F30828(varargin)
% T02F30828 MATLAB code for T02F30828.fig
%      T02F30828, by itself, creates a new T02F30828 or raises the existing
%      singleton*.
%
%      H = T02F30828 returns the handle to a new T02F30828 or the handle to
%      the existing singleton*.
%
%      T02F30828('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T02F30828.M with the given input arguments.
%
%      T02F30828('Property','Value',...) creates a new T02F30828 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before T02F30828_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to T02F30828_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help T02F30828

% Last Modified by GUIDE v2.5 02-Jun-2015 16:37:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @T02F30828_OpeningFcn, ...
                   'gui_OutputFcn',  @T02F30828_OutputFcn, ...
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


% --- Executes just before T02F30828 is made visible.
function T02F30828_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to T02F30828 (see VARARGIN)

% Choose default command line output for T02F30828
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes T02F30828 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = T02F30828_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radioeuler.
function radioeuler_Callback(hObject, eventdata, handles)
% hObject    handle to radioeuler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioeuler

if get(hObject,'Value')==1
    set(handles.radiocromer,'Value',0);
else
    set(hObject,'Value',1);
    set(handles.radiocromer,'Value',0);
end

% --- Executes on button press in radiocromer.
function radiocromer_Callback(hObject, eventdata, handles)
% hObject    handle to radiocromer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiocromer

if get(hObject,'Value')==1
    set(handles.radioeuler,'Value',0);
else
    set(hObject,'Value',1);
    set(handles.radioeuler,'Value',0);
end

function dtvalue_Callback(hObject, eventdata, handles)
% hObject    handle to dtvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtvalue as text
%        str2double(get(hObject,'String')) returns contents of dtvalue as a double


% --- Executes during object creation, after setting all properties.
function dtvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfvalue_Callback(hObject, eventdata, handles)
% hObject    handle to tfvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfvalue as text
%        str2double(get(hObject,'String')) returns contents of tfvalue as a double


% --- Executes during object creation, after setting all properties.
function tfvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in action.
function action_Callback(hObject, eventdata, handles)
% hObject    handle to action (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x t  % Necessario para a questao 6
 
m=25;
L=5;

dt=str2num(get(handles.dtvalue,'String'));
tf=str2num(get(handles.tfvalue,'String'));

t=0:dt:tf;
x(1)=6;
v(1)=0;

for i=1:(length(t)-1)
    k(i)=10+0.01*(x(i)-L);
    a(i)=-(k(i)/m)*(x(i)-L);
    v(i+1)=v(i)+dt*a(i);
    
    if get(handles.radioeuler,'Value')==1
        x(i+1)=x(i)+dt*v(i);    % Met. Euler
    else
        x(i+1)=x(i)+dt*v(i+1);  % Met. Euler-Cromer
    end
    
    % Animacao
    plot(handles.animacao,[0 0 8],[3 1 1],'-k',[0 x(i+1)],[1.5 1.5],'-b')
    hold(handles.animacao,'on')
    plot(handles.animacao,x(i+1),1.5,'or','MarkerSize',45,'MarkerFaceColor','r')
    hold(handles.animacao,'off')
    axis(handles.animacao,'off')
    
    % Grafico
    E(i+1)=0.5*((10+0.01*(x(i+1)-L))*(x(i+1)-L)^2+m*v(i+1)^2);
    plot(handles.grafico,t(1:i+1),E(1:i+1),'-k')
    xlim(handles.grafico,[0 tf])
    xlabel(handles.grafico,'tempo (s)')
    ylabel(handles.grafico,'energia mecanica (J)')
    pause(0.1)
end

% Optimizacao do Periodo
[md,k]=goldensearch(@tffiting,9,11,1,1e-3);
set(handles.msg,'String',['T = ',num2str(md),' (s)'])
