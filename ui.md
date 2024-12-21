function varargout = GUI_pingtai(varargin)
% GUI_PINGTAI MATLAB code for GUI_pingtai.fig
%      GUI_PINGTAI, by itself, creates a new GUI_PINGTAI or raises the existing
%      singleton*.
%
%      H = GUI_PINGTAI returns the handle to a new GUI_PINGTAI or the handle to
%      the existing singleton*.
%   
%      GUI_PINGTAI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PINGTAI.M with the given input arguments.
%
%      GUI_PINGTAI('Property','Value',...) creates a new GUI_PINGTAI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_pingtai_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_pingtai_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_pingtai

% Last Modified by GUIDE v2.5 12-Jul-2024 18:45:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_pingtai_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_pingtai_OutputFcn, ...
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


% --- Executes just before GUI_pingtai is made visible.初始化函数
function GUI_pingtai_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure  当前控件的句柄
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 是一个以GUI中所有控件的Tag属性为字段的结构体，每个字段的取值就是对应控件的句柄.
%类似于C语言中指针，它是某个对象的唯一标识符，通过句柄就可以找到你需要的对象
% varargin   command line arguments to GUI_pingtai (see VARARGIN)

% Choose default command line output for GUI_pingtai
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_pingtai wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_pingtai_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%打开图像
% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path]=uigetfile({'*.jpg';'*.bmp';'*.jpeg';'*.png'}, '打开文件');%uigetfile图像用户界面模块
image=[path file];
handles.file=image;
if (file==0)
    warndlg('请选择一张图片...') ;
end
[fpath, fname, fext]=fileparts(file);
validex=({'.bmp','.jpg','.jpeg','.png'});
found=0;
for (x=1:length(validex))
   if (strcmpi(fext,validex{x}))
       found=1;
     

handles.img=imread(image);
handles.i=imread(image);
h = waitbar(0,'等待...');
steps = 100;

for step = 1:steps
    waitbar(step / steps)
end
close(h) 
axes(handles.g1); 
cla; 
imshow(handles.img);
axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);
break; 
end
end
if (found==0)
     errordlg('文件扩展名不正确，请从可用扩展名[.jpg、.jpeg、.bmp、.png]中选择文件','Image Format Error');
end
%退出
% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
%保存
% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path]= uiputfile('*.jpg','Save Image as');
save=[path file]; imwrite(handles.img,save,'jpg');
%清除
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.img=handles.i;
axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);