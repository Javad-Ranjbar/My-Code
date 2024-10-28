% Methods of optical measurement in fluid mechanics
% Particle Image Velocimetry (PIV) Project
% Writer: Javad Ranjbar Kermani
% Student ID: 400209795
% Instructor: Dr Shafii
% Daedline : 1401/4/8
%% ----------------------------------------------------------------------------
clc;clear;close all;
fprintf('Starting the PIV code...\n')
tic
%Rmin = 0.7;
Rmin = input('Please enter minimum acceptable value of R corrolation : ');
%Image data loading
%raeding images of jet, mixing or natural convection
%jet
cd('D:\University\Msc\Term2\Optic\PIV\PIV-1401\PIV\1. Jet')
%Mixing
%cd('D:\University\Msc\Term2\Optic\PIV\PIV-1401\PIV\2. Mixing')
%Natural Convection
%cd('D:\University\Msc\Term2\Optic\PIV\PIV-1401\PIV\3. Natural Convection')
Image1 = double(imread('DSC_4032182.jpg')); %reading first image
Image2 = double(imread('DSC_4032182.jpg')); %reading second image
%Image2 = circshift(Image1,[25 25]);
%converting images to grayscale
[xSize,ySize,zSize] = size(Image1); 
if zSize>1  %If zSize>1, then the image is RGB, so it shoud convert to grayscale
    Image1=double((Image1(:,:,1)+Image1(:,:,2)+Image1(:,:,3))/3); 
    Image2=double((Image2(:,:,1)+Image2(:,:,2)+Image2(:,:,3))/3);
else % zSize=1, the image is grayscale
    Image1 =double(Image1); 
    Image2=double(Image2);
end
%Reversing the order of elements in the image 
%matrices to set the velocity vector direction correctly
Image1=flip(Image1,1); 
Image2=flip(Image2,1);
%Grid
% interrogation and search windows size
IWS = 41; % The size of interrogation window
SWS = 71; % The size of search window
%percent of omiting picture(like 0.1)
Omit = 0.1; %This parameter represents percent of omitting the margin of the image
xGrid = floor(Omit*xSize)+floor(IWS/2):IWS:floor((1-Omit)*xSize)-floor(IWS/2);
yGrid = floor(Omit*ySize)+floor(IWS/2):IWS:floor((1-Omit)*ySize)-floor(IWS/2);
xGridSize=size(xGrid,2);% Number of Grids in X direction
yGridSize=size(yGrid,2);% Number of Grids in Y direction
%interrogation and search windows matrix size
IWM = zeros(IWS,IWS); %interrogatio window matrix 
SWM = zeros(SWS,SWS); %search window matrix
AM = zeros(3,3);    %this matrix stores arrays that are around the Rmax
%displacement matrices
xDisp=zeros(xGridSize,yGridSize); %Displacement in X direction to reach Rmax
yDisp=zeros(xGridSize,yGridSize); %Displacement in Y direction to reach Rmax
%the following matrices store the coordinates of where the Rmax occurs
xMax=zeros(xGridSize,yGridSize); %location of Rmax in X direction
yMax=zeros(xGridSize,yGridSize); %location of Rmax in Y direction
%the following matrix stroes the maximum value of R correlation for each interrogation window
Rmax=zeros(xGridSize,yGridSize);
%Velocity field
nvf=10;  %this parameter represents number of velocity filtration.
% velocities with unit of pixel/frame
u_pf=zeros(xGridSize,yGridSize,nvf); %Velocity in X direction with unit of pixel/frame
v_pf=zeros(xGridSize,yGridSize,nvf); %Velocity in Y direction with unit of pixel/frame
V_pf=zeros(xGridSize,yGridSize,nvf); %Velocity magnitude with unit of pixel/frame
% velocities with unit of m/s
u_ms=zeros(xGridSize,yGridSize,nvf); %Velocity in X direction with unit of m/s
v_ms=zeros(xGridSize,yGridSize,nvf); %Velocity in Y direction with unit of m/s
V_ms=zeros(xGridSize,yGridSize,nvf); %Velocity magnitude with unit of m/s
fprintf('\nComputing R correlation matrix...\n')
for i=1:1:xGridSize
    for j=1:1:yGridSize
        xMinI=xGrid(i)-floor(IWS/2); %minimum x coordinate in each interrogation window
        xMaxI=xGrid(i)+floor(IWS/2); %maximum x coordinate in each interrogation window
        yMinI=yGrid(j)-floor(IWS/2); %minimum y coordinate in each interrogation window
        yMaxI=yGrid(j)+floor(IWS/2); %maximum y coordinate in each interrogation window
        xMinS=xGrid(i)-floor(SWS/2); %minimum x coordinate in each search window
        xMaxS=xGrid(i)+floor(SWS/2); %maximum x coordinate in each search window
        yMinS=yGrid(j)-floor(SWS/2); %minimum y coordinate in each search window
        yMaxS=yGrid(j)+floor(SWS/2); %maximum y coordinate in each search window 
        IWM=Image1(xMinI:xMaxI,yMinI:yMaxI);
        SWM=Image2(xMinS:xMaxS,yMinS:yMaxS);
        R = normxcorr2 (IWM,SWM);    %R Correlation matrix
        [xPeak, yPeak] = find (R==max (R(:)));
        RSize=size(R,1);             %R correlation matrix size
        % computing subpixel in which Rmax occurs
        if (xPeak>=2)
            g=1;
        else
            g=0;
        end
        if (yPeak>=2)
            h=1;
        else
            h=0;
        end
        if (xPeak<=RSize-1)
            m=1;
        else
            m=0;
        end
        if (yPeak<=RSize-1)
            n=1;
        else
            n=0;
        end
        p=g+h+m+n;
        if p==4 
           % now we consider a 2dimensional curve to compute subpixel
           % 2D curve is R=ax^2+bxy+cy^2+dx+ey+f
           AM=R(xPeak-1:xPeak+1,yPeak-1:yPeak+1);
           Rc=AM(2,2); %central element of AM matrix
           Rl=AM(2,1); %left element
           Rr=AM(2,3); %right element
           Rt=AM(1,2); %top element
           Rd=AM(3,2); %bottom element
           Rtr=AM(1,3); %top left element
           %computing coefficients of the 2D curve
           f=Rc;
           a=1/2*(Rr+Rl)-f;
           b=Rr-f-a;
           c=1/2*(Rt+Rd)-f;
           d=Rt-f-c;
           e=Rtr-a-b-c-d-f;
           if 4*a*c-e^2<0 %The condition of having a maximum
              xsubpixel=(b*e-2*a*d)/(4*a*c-e^2);
              ysubpixel=(e*d-2*b*c)/(4*a*c-e^2);
           else
              xsubpixel=0;
              ysubpixel=0;
           end
       else
           xsubpixel=0;
           ysubpixel=0;
       end
       %computing the location of xmax & ymax with respect to the origin.
       xMax(i,j)=xPeak+xsubpixel+xGrid(i)-floor(IWS/2)-floor(SWS/2);
       yMax(i,j)=yPeak+ysubpixel+yGrid(j)-floor(IWS/2)-floor(SWS/2);
       %computing displacement
       xDisp(i,j)=xMax(i,j)-xGrid(i)-1;
       yDisp(i,j)=yMax(i,j)-yGrid(j)-1;
       Rmax(i,j)=max(max(R));           
   end
end
fprintf('successfully done.\n')
% filter1
% Deleting R correlation under Rmin and filtering displacement values
filter1=1;
if (filter1==1)
  fprintf('\nFiltering the windows for which R<Rmin...\n') 
  for i=2:xGridSize-1
      for j=2:yGridSize-1
          %central grids
          if Rmax(i,j)<Rmin
              xDisp(i,j)=1/8*(xDisp(i-1,j+1)+xDisp(i-1,j)+xDisp(i-1,j-1)+...
                xDisp(i,j-1)+xDisp(i+1,j-1)+xDisp(i+1,j)+xDisp(i+1,j+1)+...
                xDisp(i,j+1));
            
              yDisp(i,j)=1/8*(yDisp(i-1,j+1)+yDisp(i-1,j)+yDisp(i-1,j-1)+...
                yDisp(i,j-1)+yDisp(i+1,j-1)+yDisp(i+1,j)+yDisp(i+1,j+1)+...
                yDisp(i,j+1));
          end
          %top boundary grids
          if Rmax(1,j)<Rmin
              xDisp(1,j)=1/5*(xDisp(1,j-1)+xDisp(1,j+1)+xDisp(2,j-1)+...
                xDisp(2,j)+xDisp(2,j+1));
            
              yDisp(1,j)=1/5*(yDisp(1,j-1)+yDisp(1,j+1)+yDisp(2,j-1)+...
                yDisp(2,j)+yDisp(2,j+1));
          end
          %bottom boundary grids
          if Rmax(xGridSize,j)<Rmin
              xDisp(xGridSize,j)=1/5*(xDisp(xGridSize,j-1)+...
               xDisp(xGridSize,j+1)+xDisp(xGridSize-1,j-1)+...
               xDisp(xGridSize-1,j)+xDisp(xGridSize-1,j+1));
           
             yDisp(xGridSize,j)=1/5*(yDisp(xGridSize,j-1)+...
              yDisp(xGridSize,j+1)+yDisp(xGridSize-1,j-1)+...
              yDisp(xGridSize-1,j)+yDisp(xGridSize-1,j+1));
          end
          %right boundary grids
          if Rmax(i,yGridSize)<Rmin
              xDisp(i,yGridSize)=1/5*(xDisp(i+1,yGridSize)+...
                xDisp(i-1,yGridSize)+xDisp(i+1,yGridSize-1)+...
                xDisp(i,yGridSize-1)+xDisp(i-1,yGridSize-1));
            
              yDisp(i,yGridSize)=1/5*(yDisp(i+1,yGridSize)+...
                yDisp(i-1,yGridSize)+yDisp(i+1,yGridSize-1)+...
                yDisp(i,yGridSize-1)+yDisp(i-1,yGridSize-1));
          end
          %left boundary grids
          if Rmax(i,1)<Rmin
              xDisp(i,1)=1/5*(xDisp(i+1,1)+xDisp(i-1,1)+...
               xDisp(i+1,2)+xDisp(i,2)+xDisp(i-1,2));  
             
              yDisp(i,1)=1/5*(yDisp(i+1,1)+yDisp(i-1,1)+...
               yDisp(i+1,2)+yDisp(i,2)+yDisp(i-1,2));
          end
          %corners
          %upper left corner
          if Rmax(1,1)<Rmin
              xDisp(1,1)=1/3*(xDisp(1,2)+xDisp(2,1)+xDisp(2,2));
              yDisp(1,1)=1/3*(yDisp(1,2)+yDisp(2,1)+yDisp(2,2));
          end
          %lower left corner
          if Rmax(xGridSize,1)<Rmin
              xDisp(xGridSize,1)=1/3*(xDisp(xGridSize,2)...
               +xDisp(xGridSize-1,2)+xDisp(xGridSize-1,1));
            
              yDisp(xGridSize,1)=1/3*(yDisp(xGridSize,2)...
                +yDisp(xGridSize-1,2)+yDisp(xGridSize-1,1));
          end
          %upper right corner
          if Rmax(1,yGridSize)<Rmin
             xDisp(1,yGridSize)=1/3*(xDisp(1,yGridSize-1)...
                +xDisp(2,yGridSize-1)+xDisp(2,yGridSize));
            
             yDisp(1,yGridSize)=1/3*(yDisp(1,yGridSize-1)...
                 +yDisp(2,yGridSize-1)+yDisp(2,yGridSize));
          end
          %lower right
          if Rmax(xGridSize,yGridSize)<Rmin
              xDisp(xGridSize,yGridSize)=1/3*(xDisp(xGridSize-1,yGridSize)+...
               xDisp(xGridSize,yGridSize-1)+xDisp(xGridSize-1,yGridSize-1));
           
              yDisp(xGridSize,yGridSize)=1/3*(yDisp(xGridSize-1,yGridSize)+...
               yDisp(xGridSize,yGridSize-1)+yDisp(xGridSize-1,yGridSize-1));
          end       
      end
  end
end
fprintf('successfully done.\n')
%computing velocities(pixel/frame)
deltat = 1;
fprintf('\nComputing velocities...\n')
for i=1:1:xGridSize
    for j=1:1:yGridSize
        u_pf(i,j,1)=yDisp(i,j)/deltat;
        v_pf(i,j,1)=xDisp(i,j)/deltat;
        V_pf(i,j,1)=sqrt( (u_pf(i,j,1))^2+(v_pf(i,j,1))^2);
    end
end
fprintf('successfully done.\n')
%filter2
VF=1.20; % parameter for velocity filteration
filter2=1;
if (filter2==1)
    fprintf('\nstarting velocity filteration...\n')
    fprintf('number of velocity flitration:')
    for k=1:1:nvf
        fprintf('\n')
        fprintf(num2str(k))
        for i=2:xGridSize-1
            for j=2:yGridSize-1
                %central grids
                if (V_pf(i,j,k)>VF*1/8*(V_pf(i+1,j+1,k)...
                       +V_pf(i+1,j,k)+V_pf(i+1,j-1,k)+V_pf(i,j-1,k)+...
                        V_pf(i-1,j-1,k)+V_pf(i-1,j,k)+V_pf(i-1,j+1,k)+...
                        V_pf(i,j+1,k)))
                    
                     u_pf(i,j,k+1)=1/8*(u_pf(i+1,j+1,k)+u_pf(i+1,j,k)+...
                        u_pf(i+1,j-1,k)+u_pf(i,j-1,k)+u_pf(i-1,j-1,k)+...
                        u_pf(i-1,j,k)+u_pf(i-1,j+1,k)+u_pf(i,j+1,k));
                    
                     v_pf (i,j,k+1)=1/8*(v_pf(i+1,j+1,k)+v_pf(i+1,j,k)+...
                        v_pf(i+1,j-1,k)+v_pf(i,j-1,k)+v_pf(i-1,j-1,k)+...
                        v_pf(i-1,j,k)+v_pf(i-1,j+1,k)+v_pf(i,j+1,k));
                    
                     V_pf (i,j,k+1)=sqrt(u_pf(i,j,k+1)^2+v_pf(i,j,k+1)^2);
                else
                    u_pf(i,j,k+1)=u_pf(i,j,k);
                    v_pf(i,j,k+1)=v_pf(i,j,k);
                    V_pf(i,j,k+1)=V_pf(i,j,k);
                end
                %top boundary grids
                if (V_pf(1,j,k)>VF*1/5*(V_pf(1,j+1,k)+V_pf(1,j-1,k)+...
                        V_pf(2,j-1,k)+V_pf(2,j,k)+V_pf(2,j+1,k)))
                    
                    u_pf(1,j,k+1)=1/5*(u_pf(1,j+1,k)+u_pf(1,j-1,k)+...
                        u_pf(2,j-1,k)+u_pf(2,j,k)+u_pf(2,j+1,k));
                    
                    v_pf(1,j,k+1)=1/5*(v_pf(1,j+1,k)+v_pf(1,j-1,k)+...
                        v_pf(2,j-1,k)+v_pf(2,j,k)+v_pf(2,j+1,k));
                    
                    V_pf(1,j,k+1)=sqrt(u_pf(1,j,k+1)^2+v_pf(1,j,k+1)^2);
                else
                    u_pf(1,j,k+1)=u_pf(1,j,k); 
                    v_pf(1,j,k+1)=v_pf(1,j,k);
                    V_pf(1,j,k+1)=V_pf(1,j,k);
                end
                %bottom boundary grids
                if (V_pf(xGridSize,j,k)>VF*1/5*(V_pf(xGridSize,j-1,k)...
                        +V_pf(xGridSize,j+1,k)+V_pf(xGridSize-1,j-1,k)+...
                         V_pf(xGridSize-1,j,k)+V_pf(xGridSize-1,j+1,k)))
                    
                    u_pf(xGridSize,j,k+1)=1/5*(u_pf(xGridSize,j-1,k)...
                        +u_pf(xGridSize,j+1,k+1)+u_pf(xGridSize-1,j-1,k)+...
                         u_pf(xGridSize-1,j,k+1)+u_pf(xGridSize-1,j+1,k));
                    
                    v_pf(xGridSize,j,k+1)=1/5*(v_pf(xGridSize,j-1,k)...
                      +v_pf(xGridSize,j+1,k+1)+v_pf(xGridSize-1,j-1,k)+... 
                       v_pf(xGridSize-1,j,k+1)+v_pf(xGridSize-1,j+1,k));
                  
                    V_pf(xGridSize,j,k+1)=sqrt(u_pf(xGridSize,j,k+1)^2 ...
                      +v_pf(xGridSize,j,k+1)^2); 
                else
                    u_pf(xGridSize,j,k+1)=u_pf(xGridSize,j,k); 
                    v_pf(xGridSize,j,k+1)=v_pf(xGridSize,j,k);
                    V_pf(xGridSize,j,k+1)=V_pf(xGridSize,j,k);
                end
                %right boundary grids
                if (V_pf(i,yGridSize,k)>VF*1/5*(V_pf(i-1,yGridSize,k)...
                    +V_pf(i+1,yGridSize,k)+V_pf(i+1,yGridSize-1,k)+...
                     V_pf(i,yGridSize-1,k)+V_pf(i-1,yGridSize-1,k)))
                
                    u_pf(i,yGridSize,k+1)=1/5*(u_pf(i-1,yGridSize,k)+...
                     +u_pf(i+1,yGridSize,k)+u_pf(i+1,yGridSize-1,k)+...
                       u_pf(i,yGridSize-1,k)+u_pf(i-1,yGridSize-1,k));
             
                    v_pf(i,yGridSize,k+1)=1/5*(v_pf(i-1,yGridSize,k)+...
                     +v_pf(i+1,yGridSize,k)+v_pf(i+1,yGridSize-1,k)+... 
                      v_pf(i,yGridSize-1,k)+v_pf(i-1,yGridSize-1,k));
            
                    V_pf(i,yGridSize,k+1)=sqrt(u_pf(i,yGridSize,k+1)^2+...
                     v_pf(i,yGridSize,k+1)^2);
                else
                    u_pf(i,yGridSize,k+1)=u_pf(i,yGridSize,k);
                    v_pf(i,yGridSize,k+1)=v_pf(i,yGridSize,k);
                    V_pf(i,yGridSize,k+1)=V_pf(i,yGridSize,k);
                end
                %left boundary grids
                if (V_pf(i,1,k)>VF*1/5*(V_pf(i+1,1,k)+V_pf(i-1,1,k)+...
                     V_pf(i+1,2,k)+V_pf(i,2,k)+V_pf(i-1,2,k)))
                 
                    u_pf(i,1,k+1)=1/5*(u_pf(i+1,1,k)+u_pf(i-1,1,k)+...
                     u_pf(i+1,2,k)+u_pf(i,2,k)+u_pf(i-1,2,k));
               
                    v_pf(i,1,k+1)=1/5*(v_pf(i+1,1,k)+v_pf(i-1,1,k)+...
                     v_pf(i+1,2,k)+v_pf(i,2,k)+v_pf(i-1,2,k));
               
                    V_pf(i,1,k+1)=sqrt(u_pf(i,1,k+1)^2+v_pf(i,1,k+1)^2);
                else
                    u_pf(i,1,k+1)=u_pf(i,1,k);
                    v_pf(i,1,k+1)=v_pf(i,1,k);
                    V_pf(i,1,k+1)=V_pf(i,1,k);
                end
                %corners
                %upper left grid
                if (V_pf(1,1,k)>VF*1/3*(V_pf(1,2,k)+V_pf(2,1,k)+...
                        V_pf(2,2,k)))
                    
                    u_pf(1,1,k+1)=1/3*(u_pf(1,2,k)+u_pf(2,1,k)+u_pf(2,2,k));
                    
                    v_pf(1,1,k+1)=1/3*(v_pf(1,2,k)+v_pf(2,1,k)+v_pf(2,2,k));
                    
                    V_pf(1,1,k+1)=sqrt(u_pf(1,1,k+1)^2+v_pf(1,1,k+1)^2);
                else
                    u_pf(1,1,k+1)=u_pf(1,1,k);
                    v_pf(1,1,k+1)=v_pf(1,1,k);
                    V_pf(1,1,k+1)=V_pf(1,1,k); 
                end
                %lower left grid
                if (V_pf(xGridSize,1,k)>VF*1/3*(V_pf(xGridSize,2,k)+...
                        V_pf(xGridSize-1,1,k)+V_pf(xGridSize-1,2,k)))
                    
                    u_pf(xGridSize,1,k+1)=1/3*(u_pf(xGridSize,2,k)...
                      +u_pf(xGridSize-1,1,k)+u_pf(xGridSize-1,2,k));
                  
                    v_pf(xGridSize,1,k+1)=1/3*(v_pf(xGridSize,2,k)...
                      +v_pf(xGridSize-1,1,k)+v_pf(xGridSize-1,2,k));
                  
                    V_pf(xGridSize,1,k+1)=sqrt(u_pf(xGridSize,1,k+1)^2 ...
                      +v_pf(xGridSize,1,k+1)^2);
                else
                    u_pf(xGridSize,1,k+1)=u_pf(xGridSize,1,k);
                    v_pf(xGridSize,1,k+1)=v_pf(xGridSize,1,k);
                    V_pf(xGridSize,1,k+1)=V_pf(xGridSize,1,k); 
                end
                %lower right grid
                if (V_pf(xGridSize,yGridSize,k)>VF*1/3*...
                        (V_pf(xGridSize,yGridSize-1,k) ...
                        +V_pf(xGridSize-1,yGridSize-1,k)...
                        +V_pf(xGridSize-1,yGridSize,k)))
                    
                    u_pf(xGridSize,yGridSize,k+1)=1/3*...
                     (u_pf(xGridSize,yGridSize-1,k)...
                     +u_pf(xGridSize-1,yGridSize-1,k)...
                     +u_pf(xGridSize-1,yGridSize,k));
                 
                    v_pf(xGridSize,yGridSize,k+1)=1/3*...
                     (v_pf(xGridSize,yGridSize-1,k)...
                     +v_pf(xGridSize-1,yGridSize-1,k)...
                     +v_pf(xGridSize-1,yGridSize,k));
              
                    V_pf(xGridSize,yGridSize,k+1)=sqrt(...
                     u_pf(xGridSize,yGridSize,k+1)^2 ...
                     +v_pf(xGridSize,yGridSize,k+1)^2);
                else
                    u_pf(xGridSize,yGridSize,k+1)=u_pf(xGridSize,yGridSize,k);
                    v_pf(xGridSize,yGridSize,k+1)=v_pf(xGridSize,yGridSize,k);
                    V_pf(xGridSize,yGridSize,k+1)=V_pf(xGridSize,yGridSize,k);
                end
                %upper right grid
                if (V_pf(1,yGridSize,k)>VF*1/3*(V_pf(2,yGridSize,k)+...
                   +V_pf(2,yGridSize-1,k)+V_pf(1,yGridSize-1,k)))
               
                    u_pf(1,yGridSize,k+1)=1/3*(u_pf(2,yGridSize,k)+...
                    +u_pf(2,yGridSize-1,k)+u_pf(1,yGridSize-1,k));
                
                    v_pf(1,yGridSize,k+1)=1/3*(v_pf(2,yGridSize,k)+...
                    +v_pf(2,yGridSize-1,k)+v_pf(1,yGridSize-1,k));
               
                    V_pf(1,yGridSize,k+1)=sqrt(u_pf(1,yGridSize,k+1)^2+...
                     v_pf(1,yGridSize,k+1)^2);
                else
                    u_pf(1,yGridSize,k+1)= u_pf(1,yGridSize,k);
                    v_pf(1,yGridSize,k+1)= v_pf(1,yGridSize,k);
                    V_pf(1,yGridSize,k+1)= V_pf(1,yGridSize,k);
                end
            end
        end
    end
end
fprintf('\n')
fprintf('\nMaximum R correlation is=')
fprintf(num2str(max(max(Rmax))))
fprintf('\nMinimum R correlation=')
fprintf(num2str(min(min(Rmax))))
fprintf('\nAverage of R correlation is =')
fprintf(num2str(mean2(Rmax)))
fprintf('\n')
% converting unit of velocity from pixel/frame to m/s
% it is known that 1 pixel = 2.64583e-4
% and time interval between the images is 1/30s(because the camera is 30fps)
dX=2.64583e-4;%converting pixel to meter
dT=deltat/30; %time interval between the images
Convert_unit=0;%convert data to SI units or not(m/s)
if (Convert_unit==1)
    for k=1:1:nvf+1
        for i=1:1:xGridSize
            for j=1:1:yGridSize
                u_ms(i,j,k)=u_pf(i,j,k)*dX/dT;
                v_ms(i,j,k)=v_pf(i,j,k)*dX/dT;
                V_ms(i,j,k)=V_pf(i,j,k)*dX/dT;
            end
        end
    end
    
    figure (1)
    quiver(yGrid,xGrid,u_ms(:,:,1),v_ms(:,:,1))
    axis([0 ySize 0 xSize])
    pbaspect([1 xSize/ySize 1])
    xlabel('Pixel in x direction','fontsize',10,'fontweight','b')
    ylabel('Pixel in y direction','fontsize',10,'fontweight','b')
    title('Velocity field without filtration','fontsize',6,'fontweight','b')
    set(gca,'fontsize',10,'fontweight','b')
    
    if (filter2==1)
        figure (2)
        quiver(yGrid,xGrid,u_ms(:,:,nvf),v_ms(:,:,nvf))
        axis([0 ySize 0 xSize])
        pbaspect([1 xSize/ySize 1])
        xlabel('Pixel in x direction','fontsize',10,'fontweight','b')
        ylabel('Pixel in y direction','fontsize',10,'fontweight','b')
        title('Velocity field after filtration','fontsize',10,'fontweight','b')
        set(gca,'fontsize',10,'fontweight','b')

        figure (3)
        pcolor(yGrid,xGrid,V_ms(:,:,nvf)),shading interp,colormap(parula(200))
        pbaspect([1 xSize/ySize 1])
        set(gca,'fontsize',10,'fontweight','b') 
        caxis([0 max(max(V_ms(:,:,nvf)))])
        ylabel(colorbar('westoutside'),'Velocity Magnitude (m/s)'...
        ,'FontSize',10,'fontweight','b','fontname','times','Rotation',90);
    end
    
else    
    figure (1)
    quiver(yGrid,xGrid,u_pf(:,:,1),v_pf(:,:,1))
    axis([0 ySize 0 xSize])
    pbaspect([1 xSize/ySize 1])
    xlabel('Pixel in x direction','fontsize',10,'fontweight','b')
    ylabel('Pixel in y direction','fontsize',10,'fontweight','b')
    title('Velocity field without filtration','fontsize',10,'fontweight','b')
    set(gca,'fontsize',10,'fontweight','b') 
    
    if (filter2==1)
        
        figure (2)
        quiver(yGrid,xGrid,u_pf(:,:,nvf),v_pf(:,:,nvf))
        axis([0 ySize 0 xSize])
        pbaspect([1 xSize/ySize 1])
        xlabel('Pixel in x direction','fontsize',10,'fontweight','b')
        ylabel('Pixel in y direction','fontsize',10,'fontweight','b')
        title('Velocity field after filtration','fontsize',10,'fontweight','b')
        set(gca,'fontsize',10,'fontweight','b') 
    
        figure (3)
        pcolor(yGrid,xGrid,V_pf(:,:,nvf)),shading interp,colormap(parula(200))
        pbaspect([1 xSize/ySize 1])
        set(gca,'fontsize',10,'fontweight','b') 
        caxis([0 max(max(V_pf(:,:,nvf)))])
        ylabel(colorbar('westoutside'),'Velocity Magnitude'...
        ,'FontSize',11,'fontname','times','Rotation',90);
    end
end
toc

