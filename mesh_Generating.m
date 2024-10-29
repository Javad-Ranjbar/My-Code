clc;
clear;
%mesh generating
delta_zeta=1;
delta_eta=1;
R=1.5;
IM=601;
JM=61;
k=0;
s1=0;
s2=0;
x(1:IM,1:JM)=0;
y(1:IM,1:JM)=0;
a(1:IM,1:JM)=1;
b(1:IM,1:JM)=1;
c(1:IM,1:JM)=1;

%Initialize type_1
for i=2:IM-1
    for j=2:JM-1
        x(i,j)=0;%+s1;
        y(i,j)=0;%+s2;
        s1=x(i,j);
        s2=y(i,j);
    end
end

%BC_1 to 2
for i=1:IM
    x(i,1)=2.5+R*cos(2*pi*(i-1)/IM);
    y(i,1)=2.5+R*sin(2*pi*(i-1)/IM);
end
x(1,1)=4;
y(1,1)=2.5;
x(IM,1)=4;
y(IM,1)=2.5;

%BC_2 to 3 ,8 to 1
%for j=1:JM
    %x(1,j)=4+(j-1)*6/60;
    %y(1,j)=2.5;
    %x(IM,j)=4+(j-1)*6/60;
    %y(IM,j)=2.5;
%end

%BC_8 to 7
for i=1:51
    x(i,JM)=10;
    y(i,JM)=2.5+(i-1)*2.5/50;
end

%BC_7 to 6
for i=51:251
    x(i,JM)=10-(i-51)*10/200;
    y(i,JM)=5;
end

%BC_6 to 5
for i=251:351
    x(i,JM)=0;
    y(i,JM)=5-(i-251)*5/100;
end

%BC_5 to 4
for i=351:551
    x(i,JM)=(i-351)*10/200;
    y(i,JM)=0;
end

%BC_4 to 3
for i=551:601
    x(i,JM)=10;
    y(i,JM)=(i-551)*2.5/50;
end

%Inner Nodes
while k<20000
    xx=x;
    yy=y;
    for i=2:IM-1
        for j=2:JM-1
            x(i,j)=(a(i,j)/(delta_zeta)^2.*(x(i+1,j)+x(i-1,j))+c(i,j)./(delta_eta)^2.*(x(i,j+1)+x(i,j-1))-b(i,j)/(2*delta_zeta*delta_eta).*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1)))./(2*(a(i,j)./delta_zeta^2+c(i,j)./delta_eta^2));
            y(i,j)=(a(i,j)/(delta_zeta)^2.*(y(i+1,j)+y(i-1,j))+c(i,j)./(delta_eta)^2.*(y(i,j+1)+y(i,j-1))-b(i,j)/(2*delta_zeta*delta_eta).*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1)))./(2*(a(i,j)./delta_zeta^2+c(i,j)./delta_eta^2));
            a(i,j)=((x(i,j+1)-x(i,j-1))/(2*delta_eta)).^2+((y(i,j+1)-y(i,j-1))/(2*delta_eta)).^2;
            b(i,j)=(x(i+1,j)-x(i-1,j))/(2*delta_zeta).*(x(i,j+1)-x(i,j-1))/(2*delta_eta)+(y(i+1,j)-y(i-1,j))/(2*delta_zeta).*(y(i,j+1)-y(i,j-1))/(2*delta_eta);
            c(i,j)=((x(i+1,j)-x(i-1,j))/(2*delta_zeta)).^2+((y(i+1,j)-y(i-1,j))/(2*delta_zeta)).^2;
        end
    end
    
    %Boundary of 1 to 8 and 2 to 3
    for j=2:JM-1
        i=1;
        x(1,j)=(a(i,j)/(delta_zeta)^2.*(x(2,j)+x(IM-1,j))+c(i,j)./(delta_eta)^2.*(x(1,j+1)+x(1,j-1))-b(i,j)/(2*delta_zeta*delta_eta).*(x(2,j+1)-x(2,j-1)+x(IM-1,j-1)-x(IM-1,j+1)))./(2*(a(i,j)./delta_zeta^2+c(i,j)./delta_eta^2));
        %y(1,j)=(a(i,j)/(delta_zeta)^2.*(y(2,j)+y(IM-1,j))+c(i,j)./(delta_eta)^2.*(y(1,j+1)+y(1,j-1))-b(i,j)/(2*delta_zeta*delta_eta).*(y(2,j+1)-x(2,j-1)+y(IM-1,j-1)-y(IM-1,j+1)))./(2*(a(i,j)./delta_zeta^2+c(i,j)./delta_eta^2));
        y(1,j)=2.5;
        x(IM,j)=x(1,j);
        y(IM,j)=y(1,j);
    end
    
    Error_X=sum(sum(abs(x-xx)));
    Error_Y=sum(sum(abs(y-yy)));
    Error_T=Error_X+Error_Y;
    if Error_X<0.05 && Error_Y<0.05
        break
    else
        k=k+1;
    end
end

t=1;
for j=1:JM
    for i=1:IM
        xxx(t,1)=x(i,j);
        yyy(t,1)=y(i,j);
        t=t+1;
    end
end
p=IM*JM;
d=(IM-1)*(JM-1);
Grid=fopen('2D_Meshing.dat','wt');
fprintf(Grid,'TITLE="2D_Meshing.dat"\nVARIABLES="X","Y"\n');
fprintf(Grid,'ZONE NODES=3721, ELEMENTS=3600, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n')

for s=1:p
    fprintf(Grid,'%.8f\t%.8f\n',xxx(s,1),yyy(s,1)');
end
for s=1:d
    fprintf(Grid,'%d\t%d\t%d\t%d\n',s,s+1,s+1+IM,s+IM);
end        
fclose(Grid);
k
figure(1)
plot(x,y,'k',x',y','k','linewidth',0.1)
%figure(2)
%plot(xxx,yyy,'k','linewidth',0.1)
%col_header={'X','Y'};
%xlswrite('Mesh_Report_2.xls',[xxx,yyy]);
            