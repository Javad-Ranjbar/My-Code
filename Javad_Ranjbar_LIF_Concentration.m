% Methods of optical measurement in fluid mechanics
% Laser Induced Fluorescence(LIF) Project
% Part1: Concentration Measurment
% Writer: Javad Ranjbar Kermani
%% -----------------------------------------------------------------------
clc;clear;close all;
tic
fprintf('Starting teh LIF code(concentration measurement)...\n')
%raeding images for concentration measurment
cd('D:\University\Msc\Term2\Optic\LIF\LIF\C\OutputPics')
%reading images
Image1=double(imread('DSC_4421070.jpg')); 
Image2=double(imread('DSC_4421150.jpg')); 
Image3=double(imread('DSC_4421230.jpg')); 
%convert to gray scale 
[xSize,ySize,zSize]=size(Image1);
if zSize>1 %If zSize>1, then the image is RGB, so it shoud convert to grayscale
  Image1=1/3*(Image1(:,:,1)+Image1(:,:,2)+Image1(:,:,3));
  Image2=1/3*(Image2(:,:,1)+Image2(:,:,2)+Image2(:,:,3)); 
  Image3=1/3*(Image3(:,:,1)+Image3(:,:,2)+Image3(:,:,3));
else % zSize=1, the image is grayscale
  Image1 = double(Image1); 
  Image2 = double(Image2);
  Image3 = double(Image3);
end 
%File location of reference images 
File_Location='D:\University\Msc\Term2\Optic\LIF\LIF\C0= 2.5e-7\OutputPics'; 
RIname=dir(fullfile(File_Location,'*.jpg')); 
NRI=numel(RIname);%Number of reference images in file location
ICI=zeros(NRI,1); %this matrix stores the intensity value of the image center 
IN=zeros(NRI,1); %storing image Number for plot function. 
SRI=zeros(xSize,ySize,NRI); %this matrix stores all reference images 
Sum_RI=zeros(xSize,ySize); %sum of all reference images 
%File location of background images 
background_File_Location='D:\University\Msc\Term2\Optic\LIF\LIF\background for concentration\OutputPics'; 
BIname=dir(fullfile(background_File_Location,'*.jpg')); 
NBI=numel(BIname); %number of background images 
SBI=zeros(xSize,ySize,NBI); %this matrix stores all background images 
Sum_BI=zeros(xSize,ySize); %sum of all background images
fprintf('averaging of reference images...\n')
for i=1:1:NRI  
    IN(i,1)=i; 
    f=fullfile(File_Location,RIname(i).name); 
    RefImage=double(imread(f)); %read reference images 
    if zSize>1 
      RefImage=1/3*(RefImage(:,:,1)+RefImage(:,:,2)+RefImage(:,:,3)); %convert to gray scale 
    end 
    ICI(i,1)=RefImage(floor(xSize/2),floor(ySize/2)); 
    %Storing grayscale images in a 3d matrix  
    SRI(:,:,i)=RefImage;
    %Summation of grayscale images to calculate the average
    Sum_RI(:,:)=Sum_RI(:,:)+SRI(:,:,i);
end 
ave_RImage=Sum_RI/NRI; %to calculate the average of reference images
hIntensity=ave_RImage(xSize/2,:);%horizontal intensity of the averaged reference images
vIntensity=ave_RImage(:,ySize/2);%horizontal intensity of the average reference images 
xPixel=zeros(xSize,1); 
yPixel=zeros(ySize,1); 
for i=1:1:xSize 
   xPixel(i,1)=i;
end 
for i=1:1:ySize 
   yPixel(i,1)=i;
end 
c1=flip(Image1./ave_RImage,1); 
c2=flip(Image2./ave_RImage,1); 
c3=flip(Image3./ave_RImage,1);
fprintf('averaging of background images...\n')
for i=1:1:NBI 
   f1=fullfile(background_File_Location,BIname(i).name); 
   backgroundImage=double(imread(f1)); %read reference images 
   if zSize>1 
      backgroundImage=1/3*(backgroundImage(:,:,1)+backgroundImage(:,:,2)+backgroundImage(:,:,3)); %convert to gray scale 
   end 
   %Storing grayscale images in a 3d matrix  
   SBI(:,:,i)=backgroundImage;
   %Summation of grayscale images to calculate the average
   Sum_BI(:,:)=Sum_BI(:,:)+SBI(:,:,i);
end 
ave_BImage=Sum_BI/NBI; %to calculate the average of background images
%substraction of background from images 
refsubtraction=ave_RImage-ave_BImage; 
imgsubtraction1=Image1-ave_BImage; 
imgsubtraction2=Image2-ave_BImage; 
imgsubtraction3=Image3-ave_BImage; 
%calculating concentration after background subtraction 
c1_new=flip(imgsubtraction1./refsubtraction,1); 
c2_new=flip(imgsubtraction2./refsubtraction,1);
c3_new=flip(imgsubtraction3./refsubtraction,1);  
%plotting the results
fprintf('Plotting the results...\n')
figure(1) 
plot(IN,ICI,'-r*','linewidth',1.5)
axis([1 NRI min(min(ICI))-3 max(max(ICI))+3])
xlabel('Number of reference images.','fontsize',10,'fontweight','b') 
ylabel('Intensity','fontsize',10,'fontweight','b')
set(gca,'fontsize',10,'fontweight','b')
grid on  
figure(2) 
plot(yPixel,hIntensity,'-k','linewidth',1.5) 
xlabel('pixel in x direction','fontsize',10,'fontweight','b') 
ylabel('Intensity','fontsize',10,'fontweight','b') 
grid on
set(gca,'fontsize',10,'fontweight','b') 
figure(3) 
plot(xPixel,vIntensity,'-m','linewidth',1.5) 
xlabel('pixel in y direction ','fontsize',10,'fontweight','b') 
ylabel('intensity','fontsize',10,'fontweight','b')  
set(gca,'fontsize',10,'fontweight','b')  
grid on
figure(4) 
imshow(mat2gray(ave_RImage)) 
figure(5) 
pcolor(yPixel,xPixel,c1),shading interp,colormap(jet(200))
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b')
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b');  
figure(6) 
pcolor(yPixel,xPixel,c2),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b');  
figure(7) 
pcolor(yPixel,xPixel,c3),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b')
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b'); 
figure(8) 
pcolor(yPixel,xPixel,c1_new),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b'); 
figure(9) 
pcolor(yPixel,xPixel,c2_new),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b');  
figure(10) 
pcolor(yPixel,xPixel,c3_new),shading interp,colormap(jet(200)) 
pbaspect([1 xSize/ySize 1]) 
set(gca,'fontsize',10,'fontweight','b') 
ylabel(colorbar('westoutside','FontSize',10, ... 
'fontname','times'),'C/C_0','FontSize',13,'fontweight','b');  
fprintf('successfully done.\n')
toc
