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

% Last Modified by GUIDE v2.5 10-Mar-2016 10:29:24

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
global value;
global prob;


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
set(handles.text20,'String','BUSY');
pause(0.1)


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



%clear data,file

%% generating points for sparse representation
% for first output




for i = 1:16
    xpred(i,1) = min(x(:,1)) + rand(1)*(max(x(:,1))-min(x(:,1)));
    xpred(i,2) = min(x(:,2)) + rand(1)*(max(x(:,2))-min(x(:,2)));
end

length = size(x);

F = [ones(length(1),1) x(:,1) x(:,2) diag(x(:,1)*(x(:,2))')];   
        
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

set(handles.text20,'String','DONE');

var1 = num2str(ypred1.est.sigma2);


set(handles.text25,'String',var1)

res1 = num2str(mean(ypred1.cv.resids));

set(handles.text26,'String',res1)


var2 = num2str(ypred2.est.sigma2);


set(handles.text28,'String',var2)

res2 = num2str(mean(ypred2.cv.resids));

set(handles.text29,'String',res2)







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



function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global value;
value = str2double(get(handles.height,'String'));






% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ypred1;
global ypred2;
global value;
global prob;
param1 = ypred1.est.beta; % Eruption height
param2 = ypred2.est.beta; % particle flux

x1 = 0.01 + 0.04*rand(10000,1);

x2 = 800 + 600*rand(10000,1);

mat = [ones(10000,1) x1 x2 x1.*x2]; 


height = mat*param1;

p=sum(height<value);

prob = p/10000;

probability = num2str(prob)
set(handles.text30,'String',probability);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prob_Callback(hObject, eventdata, handles)
% hObject    handle to prob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prob as text
%        str2double(get(hObject,'String')) returns contents of prob as a double





% --- Executes during object creation, after setting all properties.
function prob_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
