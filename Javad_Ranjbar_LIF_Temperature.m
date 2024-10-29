% Methods of optical measurement in fluid mechanics
% Laser Induced Fluorescence(LIF) Project
% Part2: Temperature Measurment
% Writer: Javad Ranjbar Kermani
%% -----------------------------------------------------------------------
clc;clear;close all; 
tic 
fprintf('Starting teh LIF code(Temperature measurement)...\n')
%File Location of images for temperature measurment 
cd('D:\University\Msc\Term2\Optic\LIF\LIF\T') 
Image1=double(imread('T 280.jpg'));  
[xSize,ySize,zSize]=size(Image1); 
%convert to gray scale 
if zSize>1 %If zSize>1, then the image is RGB, so it shoud convert to grayscale
   Image1=1/3*(Image1(:,:,1)+Image1(:,:,2)+Image1(:,:,3)); 
else % zSize=1, the image is grayscale
   Image1 = double(Image1);
end 
fprintf('reading reference images...\n')
%File location of images at T=T1=21c 
T1_File_Location='D:\University\Msc\Term2\Optic\LIF\LIF\T=21'; 
%File location of images at T=T2=53c 
T2_File_Location='D:\University\Msc\Term2\Optic\LIF\LIF\T=53'; 
%difining matrices and parameters required for computing average of images 
T1Name=dir(fullfile(T1_File_Location,'*.jpg')); 
T2Name=dir(fullfile(T2_File_Location,'*.jpg')); 
NT1=numel(T1Name);%Number of images at T1 in file location
NT2=numel(T2Name);%Number of images at T2 in file location 
sT1I=zeros(xSize,ySize,NT1); %this matrix stores all reference images at T=T1 
sT2I=zeros(xSize,ySize,NT2); %this matrix stores all reference images at T=T1
sum_T1I=zeros(xSize,ySize); %sum of all images at T=T1 
sum_T2I=zeros(xSize,ySize); %sum of all images at T=T2 
%reading reference images at T=T1
fprintf('averaging of reference images at T=21c...\n')
for i=1:1:NT1  
   f1=fullfile(T1_File_Location,T1Name(i).name); 
   T1Image=double(imread(f1)); %read reference images 
   if zSize>1 
     T1Image=1/3*(T1Image(:,:,1)+T1Image(:,:,2)+T1Image(:,:,3)); %convert to gray scale 
   end 
   %Storing grayscale images in a 3d matrix 
   sT1I(:,:,i)=T1Image;
   %Summation of grayscale images to calculate the average at T=T1
   sum_T1I(:,:)=sum_T1I(:,:)+sT1I(:,:,i);
end
ave_T1Image=sum_T1I/NT1; %to calculate the average of reference images at T=T1
fprintf('averaging of reference images at T=53c...\n')
for i=1:1:NT2
  f2=fullfile(T2_File_Location,T2Name(i).name); 
  T2Image=double(imread(f2)); %read reference images 
  if zSize>1 
     T2Image=1/3*(T2Image(:,:,1)+T2Image(:,:,2)+T2Image(:,:,3)); %convert to gray scale 
  end 
  %Storing grayscale images in a 3d matrix 
  sT2I(:,:,i)=T2Image; 
  %Summation of grayscale images to calculate the average at T=T2
  sum_T2I(:,:)=sum_T2I(:,:)+sT2I(:,:,i);
end   
ave_T2Image=sum_T2I/NT2; %to calculate the average of images at T=T2
xPixel=zeros(xSize,1); 
yPixel=zeros(ySize,1); 
for i=1:1:xSize 
    xPixel(i,1)=i; 
end 
for i=1:1:ySize 
    yPixel(i,1)=i;
end 
%Calibration
T1=21+273.15; %value of first temperature for calibration curve 
T2=53+273.15; %value of second temperature for calibration curve 
T=273:1:373; %temperature interval for calibration curve   
%Linear Calibration 
%computing the slop and y-interceot of the linear calibration curve I=aT+b
fprintf('Computing coefficient of linear an exponential calibration...\n')
a=(ave_T2Image-ave_T1Image)/(T2-T1); %slope of linear calibration curve 
b=ave_T1Image+(T1/(T1-T2))*(ave_T2Image-ave_T1Image); %y-intercept of linear calibration curve
%consider four pixel with coordinates of (100,200),(200,100),(200,300),(300,200) 
pixel1=[100,200]; 
pixel2=[200,100]; 
pixel3=[200,300]; 
pixel4=[300,200]; 
a1=a(pixel1(2),pixel1(1)); 
a2=a(pixel2(2),pixel2(1)); 
a3=a(pixel3(2),pixel3(1)); 
a4=a(pixel4(2),pixel4(1)); 
b1=b(pixel1(2),pixel1(1)); 
b2=b(pixel2(2),pixel2(1)); 
b3=b(pixel3(2),pixel3(1)); 
b4=b(pixel4(2),pixel4(1));
%Exponential Calibration
%computing K and beta of exponential calibration curve I=kexp(B/T)
beta=(1/(1/T1-1/T2))*log(ave_T1Image./ave_T2Image); %Computing beta coefficient
K=ave_T1Image./exp(beta/T1);%Computing K coefficient 
beta1=beta(pixel1(2),pixel1(1)); 
beta2=beta(pixel2(2),pixel2(1)); 
beta3=beta(pixel3(2),pixel3(1)); 
beta4=beta(pixel4(2),pixel4(1)); 
K1=K(pixel1(2),pixel1(1)); 
K2=K(pixel2(2),pixel2(1)); 
K3=K(pixel3(2),pixel3(1)); 
K4=K(pixel4(2),pixel4(1));
%omitting margin of the image using parameter OC
OC=0.05; %omitting coefficient for ommit the margin of the image where  intensity=0
Image1_new=Image1(floor(OC*xSize):floor((1-OC)*xSize),floor(OC*ySize):floor((1-OC)*ySize));
xPixel_new=xPixel(floor(OC*xSize):floor((1-OC)*xSize),:); 
yPixel_new=yPixel(floor(OC*ySize):floor((1-OC)*ySize),:); 
a_new=a(floor(OC*xSize):floor((1-OC)*xSize),floor(OC*ySize):floor((1-OC)*ySize)); 
b_new=b(floor(OC*xSize):floor((1-OC)*xSize),floor(OC*ySize):floor((1-OC)*ySize)); 
beta_new=beta(floor(OC*xSize):floor((1-OC)*xSize),floor(OC*ySize):floor((1-OC)*ySize)); 
K_new=K(floor(OC*xSize):floor((1-OC)*xSize),floor(OC*ySize):floor((1-OC)*ySize)); 
LinearT_new=flip((1./a_new).*(Image1_new-b_new),1); 
expT_new=flip(beta_new./log(Image1_new./K_new),1); 
LinearT=flip((1./a).*(Image1-b),1); 
expT=flip(beta./log(Image1./K),1);
dT1 = LinearT(100,200)-273.15;
dT2 = LinearT(200,100)-273.15;
dT3 = LinearT(200,300)-273.15;
dT4 = LinearT(300,200)-273.15;
de1 = expT(100,200)-273.15;
de2 = expT(200,100)-273.15;
de3 = expT(200,300)-273.15;
de4 = expT(300,200)-273.15;
fprintf('plotting the results...\n')
figure (1) 
imshow(mat2gray(ave_T1Image)) 
figure (2) 
imshow(mat2gray(ave_T2Image)) 
figure(3) 
pcolor(yPixel,xPixel,a),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b','YDir','reverse') 
ylabel(colorbar('westoutside','FontSize',10, ...
'fontname','times'),'a','FontSize',13,'fontname','times'); 
figure(4) 
pcolor(yPixel,xPixel,b),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b','YDir','reverse') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'b','FontSize',13,'fontname','times'); 
figure(5) 
pcolor(yPixel_new,xPixel_new,beta_new),shading interp,colormap(jet(200)) 
pbaspect([1 length(xPixel_new)/length(yPixel_new) 1]) 
set(gca,'fontsize',10,'fontweight','b','YDir','reverse')
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'\beta','FontSize',13,'fontname','times'); 
figure(6)
pcolor(yPixel_new,xPixel_new,K_new),shading interp,colormap(jet(200)) 
pbaspect([1 length(xPixel_new)/length(yPixel_new) 1]) 
set(gca,'fontsize',10,'fontweight','b','YDir','reverse')
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'K','FontSize',13,'fontname','times'); 
figure(7) 
plot(T,a1*T+b1,'-k','linewidth',1.1) 
hold on 
plot(T,a2*T+b2,'-m','linewidth',1.1)
hold on 
plot(T,a3*T+b3,'-b','linewidth',1.1)
hold on 
plot(T,a4*T+b4,'-g','linewidth',1.1)  
xlabel('T(K)','fontsize',10,'fontweight','b') 
ylabel('Intensity','fontsize',10,'fontweight','b')  
set(gca,'fontsize',10,'fontweight','b')
legend('pixel(100,200)','pixel(200,100)','pixel(200,300)','pixel(300,200)')
axis([273 373 -40 100])
grid on
figure(8) 
plot(T,K1*exp(beta1./T),'-k','linewidth',1.1) 
hold on 
plot(T,K2*exp(beta2./T),'-m','linewidth',1.1) 
hold on 
plot(T,K3*exp(beta3./T),'-b','linewidth',1.1) 
hold on 
plot(T,K4*exp(beta4./T),'-g','linewidth',1.1)   
xlabel('T(K)','fontsize',10,'fontweight','b') 
ylabel('Intensity','fontsize',10,'fontweight','b')
set(gca,'fontsize',10,'fontweight','b')
legend('pixel(100,200)','pixel(200,100)','pixel(200,300)','pixel(300,200)')
axis([273 373 0 120])
grid on
figure(9) 
LT1=K1*exp(beta1./T); 
eT1=a1*T+b1; 
m=max(max(max(LT1)),max(max(eT1))); 
n=min(min(min(LT1)),min(min(eT1))); 
L=floor(n)-3:1:floor(m)-3; 
firstT=T1*ones(length(L),1); 
secondT=T2*ones(length(L),1); 
plot(T,LT1,'-m','linewidth',1.1) 
hold on 
plot(T,eT1,'-b','linewidth',1.1) 
hold on 
plot(firstT,L,'--k'); 
hold on 
plot(secondT,L,'--k') 
text(T1+1,floor(m)-16,['T_1=',num2str(T1)],'fontname','times','Fontsize',10) 
text(T2+1,floor(m)-16,['T_2=',num2str(T2)],'fontname','times','Fontsize',10) 
set(gca,'fontsize',10,'fontweight','b')    
xlabel('T(K)','fontsize',10,'fontweight','b')
ylabel('Intensity','fontsize',10,'fontweight','b')  
legend('Exponential calibration pixel(100,200)','Linear calibration pixel(100,200)')
axis([273 373 0 100])
grid on
figure (10)
pcolor(yPixel_new,xPixel_new,LinearT_new-273.15),shading interp,colormap(jet(200)) 
pbaspect([1 length(xPixel_new)/length(yPixel_new) 1]) 
set(gca,'fontsize',10,'fontweight','b') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'T(\circC)','FontSize',13,'fontname','times');  
figure (11) 
pcolor(yPixel_new,xPixel_new,expT_new-273.15),shading interp,colormap(jet(200)) 
pbaspect([1 length(xPixel_new)/length(yPixel_new) 1]) 
set(gca,'fontsize',10,'fontweight','b')
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'T(\circC)','FontSize',13,'fontname','times');  
delete_Syring=1; %delete syringe or not (1 correspond to delete)
if (delete_Syring==1) 
   figure(12) 
   linearT_ds=LinearT_new(1:300,:)-273.15; %syringe omission 
   pcolor(yPixel_new,xPixel_new(1:300,:),linearT_ds),shading interp,colormap(jet(200)) 
   pbaspect([1 size(linearT_ds,1)/ySize 1])
   set(gca,'fontsize',10,'fontweight','b') 
   ylabel(colorbar('westoutside','FontSize',10, ... 
   'fontname','times'),'T(\circC)','FontSize',13,'fontname','times'); 
   figure(13) 
   expT_ds=expT_new(1:300,:)-273.15; 
   pcolor(yPixel_new,xPixel_new(1:300,:),expT_ds),shading interp,colormap(jet(200)) 
   pbaspect([1 size(expT_ds,1)/ySize 1]) 
   set(gca,'fontsize',10,'fontweight','b')
   ylabel(colorbar('westoutside','FontSize',10, ... 
   'fontname','times'),'T(\circC)','FontSize',13,'fontname','times');  
end 
fprintf('successfully done!\n')
toc
