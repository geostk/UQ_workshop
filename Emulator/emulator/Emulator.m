function varargout = Emulator(varargin)
% EMULATOR MATLAB code for Emulator.fig
%      EMULATOR, by itself, creates a new EMULATOR or raises the existing
%      singleton*.
%
%      H = EMULATOR returns the handle to a new EMULATOR or the handle to
%      the existing singleton*.
%
%      EMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMULATOR.M with the given input arguments.
%
%      EMULATOR('Property','Value',...) creates a new EMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Emulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Emulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Emulator

% Last Modified by GUIDE v2.5 06-Mar-2016 13:12:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Emulator_OpeningFcn, ...
                   'gui_OutputFcn',  @Emulator_OutputFcn, ...
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


% --- Executes just before Emulator is made visible.
function Emulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Emulator (see VARARGIN)

% Choose default command line output for Emulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Emulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Initializing parameters
global ypred1;
global ypred2;
global val;
global data;


% --- Outputs from this function are returned to the command line.
function varargout = Emulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ypred1;
global ypred2;
global val;
global data;

corr_method = get(handles.correlation_method,'SelectedObject');
corr_family = get(handles.correlation_family,'SelectedObject');
validation = get(handles.cross_validation,'SelectedObject');

method = get(corr_method,'String');
family = get(corr_family,'String');
val  =get(validation,'String');



%% Reading the data
% x = csvread('lhs_samples.csv');
% y1 = csvread('eruption_height.csv');
% y2 = csvread('particle_flux.csv');
% 
% y = horzcat(y1,y2);



x = data(:,1:2);
y = data(:,3:4);


for i= 1:64
    %data(i,4)
end

%clear data,file

%% generating points for sparse representation
% for first output

for i = 1:16
    xpred(i,1) = min(x(:,1)) + rand(1)*(max(x(:,1))-min(x(:,1)));
    xpred(i,2) = min(x(:,2)) + rand(1)*(max(x(:,2))-min(x(:,2)));
end


F = [ones(64,1) x(:,1) x(:,2) diag(x(:,1)*(x(:,2))')];   
        
Fpred = [ones(16,1) xpred(:,1) xpred(:,2) diag(xpred(:,1)*(xpred(:,2))')];

ypred1 = mperk('X',x, ...
                'Y',y(:,1),...
                'CrossValidation',val,...
                'CorrelationEstimationMethod',method,...
                'CorrelationFamily',family,...
                'Xpred',xpred,...
                'RegressionModel',F,...
                'PredRegressionModel', Fpred);

% for second output            
                

ypred2 = mperk('X',x, ...
                'Y',y(:,2),...
                'CrossValidation','Yes',...
                'CorrelationEstimationMethod','REML',...
                'CorrelationFamily','Gaussian',...
                'Xpred',xpred,...
                'RegressionModel',F,...
                'PredRegressionModel', Fpred);
            

            
Ypred = horzcat(ypred1.preds.ypreds,ypred2.preds.ypreds);

axes(handles.axes1);
cla(handles.axes1,'reset');
tri = delaunay(x(:,1),x(:,2));
trisurf(tri,x(:,1),x(:,2),y(:,1));

% Clean it up
l = light('Position',[-50 -15 29]);
lighting none
shading interp
colorbar EastOutside

hold on
scatter3(xpred(:,1),xpred(:,2),Ypred(:,1),100,'k','filled');
xlabel('Water Fraction');
ylabel('Temperature(K)');
zlabel('Eruption height(m)');
title('Emulator for eruption height');


axes(handles.axes2);
cla(handles.axes2,'reset');
tri = delaunay(x(:,1),x(:,2));
trisurf(tri,x(:,1),x(:,2),y(:,2));

% Clean it up
l = light('Position',[-50 -15 29]);
lighting none
shading interp
colorbar EastOutside

hold on
scatter3(xpred(:,1),xpred(:,2),Ypred(:,2),100,'k','filled');
xlabel('Water Fraction');
ylabel('Temperature(K)');
zlabel('Particle Flux');
title('Emulator for Particle Flux');



% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user dat  a (see GUIDATA)
global ypred1;
global ypred2;
global val;


    
clear x y;
x = ypred1.preds.se;
axes(handles.axes1);
cla(handles.axes1,'reset');

scatter(1:16,x,100,'filled')
xlabel('Data');
ylabel('Standard Error');
title('Eruption Height');




y = ypred2.preds.se;
axes(handles.axes2);
cla(handles.axes2,'reset');

scatter(1:16,y,100,'filled')

xlabel('Data');
ylabel('Standard Error');
title('Particle Flux');










% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ypred1;
global ypred2;

x = ypred1.cv.se;
axes(handles.axes1);
cla(handles.axes1,'reset');

scatter(1:64,x,100,'filled')
xlabel('Data');
ylabel('Standard Error');
title('Eruption Height');

y = ypred2.cv.se;
axes(handles.axes2);
cla(handles.axes2,'reset');

scatter(1:64,y,100,'filled')
xlabel('Data');
ylabel('Standard Error');
title('Particle Flux');






% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ypred1;
global ypred2;

x = ypred1.cv.resids;
axes(handles.axes1);
cla(handles.axes1,'reset');

scatter(1:64,x,100,'filled')
xlabel('Data');
ylabel('Residuals');
title('Eruption Height');


y = ypred2.cv.resids;
axes(handles.axes2);
cla(handles.axes2,'reset');

scatter(1:64,y,100,'filled')
scatter(1:64,y,100,'filled')
xlabel('Data');
ylabel('Residuals');
title('Particle Flux');



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data
file = uigetfile;
data = csvread(file);
