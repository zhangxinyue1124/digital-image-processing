%边缘检测prewitt算子
% --- Executes on button press in m3.
function m3_Callback(hObject, eventdata, handles)
% hObject    handle to m3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img=rgb2gray(handles.img);
end
handles.img=edge(handles.img,'prewitt');
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

%边缘检测canny算子
% --- Executes on button press in m5.
function m5_Callback(hObject, eventdata, handles)
% hObject    handle to m5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img=rgb2gray(handles.img);
end
handles.img=edge(handles.img,'canny');
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

%边缘检测roberts算子
% --- Executes on button press in m2.
function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img=rgb2gray(handles.img);
end
handles.img=edge(handles.img,'roberts');
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

%边缘检测log算子
% --- Executes on button press in m4.
function m4_Callback(hObject, eventdata, handles)
% hObject    handle to m4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img=rgb2gray(handles.img);
end
handles.img=edge(handles.img,'log');
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);


%高斯噪声
% --- Executes on button press in n1.
function n1_Callback(hObject, eventdata, handles)
%灰度图变化
img =handles.i;
img = rgb2gray(img);
size(img)
axes(handles.g2);
imshow(img)

% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in n2.
function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img =imadjust(rgb2gray(handles.i));
axes(handles.g2);
imshow(img)
%椒盐噪声
% --- Executes on button press in n3.
function n3_Callback(hObject, eventdata, handles)
%锐化
shadow = handles.i; %待增强图片

shadow_lab = rgb2lab(shadow);
% 亮度值的范围是从 0 到 100。将值缩放到范围 [0 1]，这是数据类型为 double 的图像的预期范围。

max_luminosity = 100;
L = shadow_lab(:,:,1)/max_luminosity;
% 对亮度通道执行三种类型的对比度调整，并保持 a* 和 b* 通道不变。将图像转换回 RGB 颜色空间。


shadow_histeq = shadow_lab;
shadow_histeq(:,:,1) = histeq(L)*max_luminosity;
shadow_histeq = lab2rgb(shadow_histeq);

axes(handles.g2);
imshow(shadow_histeq)
% hObject    handle to n3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%左右翻转
% --- Executes on button press in m8.
function m8_Callback(hObject, eventdata, handles)
% hObject    handle to m8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.img=fliplr(handles.img);
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);
%上下翻转
% --- Executes on button press in m9.
function m9_Callback(hObject, eventdata, handles)
% hObject    handle to m9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.img=flipud(handles.img);
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

%灰度处理
% --- Executes on button press in m6.
function m6_Callback(hObject, eventdata, handles)
% hObject    handle to m6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.img = rgb2gray(handles.img);
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

%二值化
% --- Executes on button press in m7.
function m7_Callback(hObject, eventdata, handles)
% hObject    handle to m7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresh = graythresh(handles.img);     %自动确定二值化阈值
handles.img = im2bw(handles.img,thresh);
axes(handles.g2);
cla;
imshow(handles.img);
guidata(hObject,handles);

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%图形旋转
% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rrv=(get(hObject,'Value'));
handles.rot=handles.img;
handles.rot=imrotate(handles.rot,rrv);%imrotate使图像旋转
axes(handles.g2); cla; imshow(handles.rot);
guidata(hObject,handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%理想低通滤波器
% --- Executes on button press in p2.
function p2_Callback(hObject, eventdata, handles)
% hObject    handle to p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img = rgb2gray(handles.img);
end

I=handles.img;
I=im2double(I);
M=2*size(I,1);
N=2*size(I,2);
u=-M/2:(M/2-1);
v=-N/2:(N/2-1);
[U,V]=meshgrid(u,v);
D=sqrt(U.^2+V.^2);
D0=80;
H=double(D<=D0);
J=fftshift(fft2(I,size(H,1),size(H,2)));
K=J.*H;
L=ifft2(ifftshift(K));
L=L(1:size(I,1),1:size(I,2));
handles.img=L;

axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);

%巴特沃斯高通滤波器
% --- Executes on button press in p3.
function p3_Callback(hObject, eventdata, handles)
% hObject    handle to p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img = rgb2gray(handles.img);
end

I=handles.img;
I=im2double(I);
M=2*size(I,1);%滤波器行数
N=2*size(I,2);%滤波器列数
u=-M/2:(M/2-1);
v=-N/2:(N/2-1);
[U,V]=meshgrid(u,v);
D=sqrt(U.^2+V.^2);
D0=30;%截止频率
n=6;%巴特沃斯滤波器阶数
H=1./(1+(D0./D).^(2*n));
J=fftshift(fft2(I,size(H,1),size(H,2)));
K=J.*H;
L=ifft2(ifftshift(K));
L=L(1:size(I,1),1:size(I,2));
handles.img=L;

axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);

%带阻滤波器
% --- Executes on button press in p4.
function p4_Callback(hObject, eventdata, handles)
% hObject    handle to p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img = rgb2gray(handles.img);
end

I=handles.img;
I=imnoise(I,'gaussian',0,0.01);
I=im2double(I);
M=2*size(I,1);%滤波器行数
N=2*size(I,2);%滤波器列数
u=-M/2:(M/2-1);
v=-N/2:(N/2-1);
[U,V]=meshgrid(u,v);
D=sqrt(U.^2+V.^2);
D0=50;%滤波器D0
W=30;%滤波器带宽
H=double(or(D<(D0-W/2),D>D0+W/2));
J=fftshift(fft2(I,size(H,1),size(H,2)));
K=J.*H;
L=ifft2(ifftshift(K));
L=L(1:size(I,1),1:size(I,2));
handles.img=L;

axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);

%同态滤波器
% --- Executes on button press in p5.
function p5_Callback(hObject, eventdata, handles)
% hObject    handle to p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysize=size(handles.img);
if numel(mysize)>2
    handles.img = rgb2gray(handles.img);
end

I=handles.img;
J=log(im2double(I)+1);
K=fft2(J);
n=5;
D0=0.1*pi;
rh=0.7;
rl=0.4;
[row, column]=size(J);
for i=1:row
    for j=i:column
        D1(i,j)=sqrt(i^2+j^2);
        H(i,j)=rl+(rh/(1+(D0/D1(i,j))^(2*n)));
    end
end
L=K.*H;
M=ifft2(L);
N=exp(M)-1;
handles.img=N;

axes(handles.g2); 
cla; 
imshow(handles.img);
guidata(hObject,handles);
%均值滤波
% --- Executes on button press in f1.
function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=fspecial('average');
handles.img=imfilter(handles.img,h,'replicate');
axes(handles.g2); cla; imshow(handles.img)
guidata(hObject,handles);
%高斯滤波
% --- Executes on button press in f2.
function f2_Callback(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hsize=[8 8]; sigma=1.7;%hsize是滤波器尺寸，sigma是标准差
h=fspecial('gaussian',hsize,sigma);%fspecial生成滤波器
handles.img=imfilter(handles.img,h,'replicate');%imfilter滤波，控制滤波运算的选项类型'replicate'
axes(handles.g2); cla; imshow(handles.img);
guidata(hObject,handles);
%中值滤波
% --- Executes on button press in f3.
function f3_Callback(hObject, eventdata, handles)
% hObject    handle to f3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
r=medfilt2(handles.img(:,:,1));%执行二维中值滤波
g=medfilt2(handles.img(:,:,2));
b=medfilt2(handles.img(:,:,3)); 
handles.img=cat(3,r,g,b);%cat构造三维数组
axes(handles.g2); cla; imshow(handles.img);
guidata(hObject,handles);


%亮度调节
% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=get(hObject,'Value');
img=handles.img;
img=img+x;
axes(handles.g2); cla; imshow(img)
handles.img=img;
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


% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
%裁剪
%裁剪
imgdata = handles.i;
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');
finalRect = rbbox;
point2 = get(gca,'CurrentPoint');
point1 = point1(1,1:2);
point2 = point2(1,1:2);
p1 = min(floor(point1),floor(point2));
p2 = max(floor(point1),floor(point2));
offset = abs(floor(point1)-floor(point2));
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(1) p1(2)];
c = imcrop(imgdata,[p1(1) p1(2) offset(1) offset(2)]); %得到鼠标框选的范围图片
axes(handles.g2)%指定显示区域
imshow(c)
% hObject    handle to pushbutton47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
rrv=(get(hObject,'Value'))
handles.sf=handles.img;
handles.sf=imresize(handles.sf,rrv);%
axes(handles.g2); cla; imshow(handles.sf);
guidata(hObject,handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
A=rgb2gray(handles.img);   %读取到一张图片  
thresh = graythresh(A);     %自动确定二值化阈值

I1 = imbinarize(A,thresh);       %对图像二值化
axes(handles.g2)
imshow(I1)    %显示二值化之后的图片   
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)

I=rgb2gray(handles.img);   %读取到一张图片  
%设置初始阈值   去最大值和最小值的中间值
zmax=max(max(I));
zmin=min(min(I));
T1=(zmax+zmin)/2;
%根据阈值将图像进行分割为前景和背景，分别求出两者的平均灰度  z1和z2
b=1;
[m n]=size(I);
while (b)
        ifg=0;
        ibg=0;
        fnum=0;
        bnum=0;
        for i=1:m
            for j=1:n
                tmp=I(i,j);
                if(tmp>=T1)
                    ifg=ifg+1;
                    fnum=fnum+double(tmp);  %前景像素的个数以及像素值的总和
                else
                    ibg=ibg+1;
                    bnum=bnum+double(tmp);%背景像素的个数以及像素值的总和
                end
            end
        end
        %计算前景和背景的平均值
        z1=fnum/ifg;
        z2=bnum/ibg;
        if(T1==(uint8((z1+z2)/2)))
            b=0;
        else
            T1=uint8((z1+z2)/2);
        end
        %当阈值不变换时，退出迭代
end
thresh = double(T1)/255;
I1=imbinarize(I,thresh);

axes(handles.g2)
imshow(I1)    %显示二值化之后的图片
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
img = handles.i;
%% 低通滤波
for k = 1:3
    image_bin = img(:,:,k);
    [M ,N]=size(image_bin);%图片尺寸
    d0 =80;%滤波半径
    I_f = fft2(double(image_bin));%傅里叶变换得到频谱
    I_f=fftshift(I_f);  %移到中间
    m_mid=floor(M/2);%中心点坐标
    n_mid=floor(N/2);  
    h = zeros(M,N);%低通滤波器构造
    for i = 1:M
        for j = 1:N
            d = ((i-m_mid)^2+(j-n_mid)^2);
            h(i,j) = exp(-(d)/(2*(d0^2)));      
        end
    end
    I_lpf = h.*I_f;
    I_lpf=ifftshift(I_lpf);    %中心平移回原来状态
    I_filter(:,:,k)=uint8(real(ifft2(I_lpf)));  %反傅里叶变换,取实数部分
end
axes(handles.g2)%指定显示区域
imshow(I_filter)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
%聚类


org =handles.i;
k = str2double(get(handles.edit1,'string'));
           
%转化为灰度图
gray=rgb2gray(org);
[m,n]=size(gray); %m,n为所求，p=3为通道数

A = reshape(org(:, :, 1), m*n, 1);   
B = reshape(org(:, :, 2), m*n, 1);
C = reshape(org(:, :, 3), m*n, 1);
data = [A B C]; 
% 第二张图
% kmeans第一个参数: N*P的数据矩阵，N为数据个数，P为单个数据维度
res = kmeans(double(data), k);        
result = reshape(res, m, n);             % 反向转化为图片形式
axes(handles.g2)
% label2rgb功能是转换标记矩阵到RGB图像
imshow(label2rgb(result)),title(strcat('K=',num2str(k),'时RGB通道分割结果')); %显示分割结果


 
 

% 
% segmented_images = cell(1,3);%1*3矩阵
% rgb_label = repmat(label1,[1 1 3]);%把label矩阵编程[1 1 3]的矩阵
% 
% for k = 1:10
%     color = I;
%     color(rgb_label ~= k) = 0;%把label不等于K的点除去
%     segmented_images{k} = color;
% end
% 
% figure(),imshow(segmented_images{1}), title('分割结果——区域1');
% figure(),imshow(segmented_images{2}), title('分割结果——区域2');
% figure(),imshow(segmented_images{3}), title('分割结果——区域3');

