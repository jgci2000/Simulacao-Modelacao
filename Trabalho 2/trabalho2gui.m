function varargout = trabalho2gui(varargin)
% TRABALHO2GUI MATLAB code for trabalho2gui.fig
%      TRABALHO2GUI, by itself, creates a new TRABALHO2GUI or raises the existing
%      singleton*.
%
%      H = TRABALHO2GUI returns the handle to a new TRABALHO2GUI or the handle to
%      the existing singleton*.
%
%      TRABALHO2GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRABALHO2GUI.M with the given input arguments.
%
%      TRABALHO2GUI('Property','Value',...) creates a new TRABALHO2GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trabalho2gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trabalho2gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trabalho2gui

% Last Modified by GUIDE v2.5 26-May-2019 14:10:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trabalho2gui_OpeningFcn, ...
                   'gui_OutputFcn',  @trabalho2gui_OutputFcn, ...
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


% --- Executes just before trabalho2gui is made visible.
function trabalho2gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trabalho2gui (see VARARGIN)

% Choose default command line output for trabalho2gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trabalho2gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trabalho2gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in video1.
function video1_Callback(hObject, eventdata, handles)
% hObject    handle to video1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mov = VideoReader('bola.mp4');
video=read(mov);
s = size(video);
nFrames = s(4);
framerate = mov.FrameRate;

frame1= video(:,:,:,1) ;
frame1_Vermelho= frame1(:,:,1) ;

G= frame1(:,:,2) ;
B= frame1(:,:,3) ;


F= frame1_Vermelho<85 ;
Imagem=frame1_Vermelho(F);

[yp,xp]=find(F==1);

c1=min(yp);
c2=max(yp);
c=(c2-c1);

bs=zeros(866,2);
b=mean([xp,yp]);

pixeisPorCenti=c/5;
pixeisPorMetro=(c/0.05);

yVideo=680/pixeisPorMetro;
xVideo=300/pixeisPorMetro;

for i=1:nFrames
    
    frame= video(:,:,:,i) ;
    frame_Vermelho= frame(:,:,1) ;
    G= frame(:,:,2) ;
    B= frame(:,:,3) ;
    F=frame_Vermelho<85;
    Imagem=frame_Vermelho(F);
    [yp,xp]=find(F==1);
    b=mean([xp,yp]);
    bs(i,1)=b(1);
    bs(i,2)=b(2);
    x=b(1);
    y=b(2);
    imagesc(handles.axes1,frame)
    hold (handles.axes1,'on')
    plot(handles.axes1,x,y,'r.')
    drawnow update
    
end

bsMetros= bs./pixeisPorMetro;

bsMetros(:,2)= yVideo-bsMetros(:,2);

i=1:866;
dt=1/framerate;

t = i.*dt;

plot(handles.axes5,i*dt,bsMetros(i,2),'r.');
hold on

tempo = linspace(0,6,6*framerate);
interpolacao = spline(t,bsMetros(:,2),tempo);

plot(handles.axes2,tempo,interpolacao,'k-')
pause(1)

[maximos,t_max] = findpeaks(interpolacao,tempo);
[minimos,t_min] = findpeaks(-interpolacao,tempo);

% for i=1:19
%      delta_t=t_min(i+1)- t_min(i);
% end
%      g1 = 8*maximos(1)/(delta_t(1)^2);
%      display(g1)

for j=1:19
        t1 = t_min(j):dt:t_min(j+1);
        t2=round(t1./dt);
        
        y1=interpolacao(t2);
        
        
        poli=polyfit(t2.*dt,y1,2);
        
        g2(j)=-2*poli(1);
end
    g3=sum(g2)/length(g2);
    display(g3)
    
%%---------------------------------PARTE 2---------------------------------
    
g=-9.8;
e=sqrt(maximos(2)/maximos(1));
dt=0.0001;
v(1)=0;
y(1)=1;
t= 0:dt:6;

for i = 1:length(t)-1
    if y(i)<=0 && v(i)<0
        v(i+1)=v(i)*(-e);
    else
        v(i+1) = v(i)+g*dt;
    end
    y(i+1) = y(i)+v(i)*dt;
end


plot(handles.axes3,t,y,'r-')



g=-9.8;
e=sqrt(maximos(2)/maximos(1));
dt=0.0001;
v(1)=0;
y(1)=1;
t= 0:dt:6;

for i = 1:length(t)-1
    if y(i)<=0 && v(i)<0
        v(i+1)=v(i)*(-e);
    else
        v(i+1) = v(i)+g*dt;
    end
    y(i+1) = y(i)+v(i+1)*dt;
end
    
 plot(handles.axes4,t,y,'b--')
