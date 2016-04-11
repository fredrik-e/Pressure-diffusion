function [S,times]=Crank_Nich_2D(data,dP,D,mmres,mmperpx,time,framestep,P_initial)

% This function solves the 2-D diffusion equation by the Crank-Nicholson
% method. The left boundary and fracture pattern in data is held at
% constant overpressure dP, while the right boundary is held at 0. The top
% and bottom boundaries reflects the closest interior points.
% Inputs: data=binary image of fracture pattern
%         dP=overpressure
%         D=diffusion constant mm^2/s
%         mmres=the wanted space steps in mm
%         mmperpx=the original mm per pixel value
%         time=duration time in seconds for simulation
%         framestep=time step between state records
%         P_initial=initial state of the pressure field (optional)
% Outputs: S=cell of stored states, times=list of time at stored states

scale=mmperpx/mmres;
P=double(fliplr(imresize(data,scale,'nearest')));
if ~exist('P_initial')
    P_initial=zeros(size(P));
else
    P_initial=imresize(P_initial,size(P),'nearest');
end
dx=mmres;       %x-step
dt=((dx^2)/(2*D)); %define max timestep (alpha_min=0.5)
a=(D*dt)/(dx^2); %factor in iteration scheme

%Boundary functions
P=padarray(P,[0 1],0,'both');
Frac=(P==1);    %fracture
P=padarray(P_initial,[0 1],0,'both');
[~,J]=size(P);
P(:,1)=dP;       %left
P(:,J)=0;       %right
P(Frac)=dP;     %set pressure

I=51;J=51;
%Left side tridiagonal matrix with fringes
d=[-I,-1,0,1,I]';
D0=ones(I*J,1)*((4*a)+2);
D1=ones(I*J,1)*(-a);
Diags=cat(2,D1,D1,D0,D1,D1);
A=spdiags(Diags,d,I*J,I*J);
%Right side tridiagonal matrix with fringes
d=[-I,-1,0,1,I]';
D0=ones(I*J,1)*(2-(4*a));
D1=ones(I*J,1)*(a);
Diags=cat(2,D1,D1,D0,D1,D1);
B=spdiags(Diags,d,I*J,I*J);
clear D0 D1 Diags d
%Weight matrix for Crank-Nicholson vector equation
F=A\B;
clear A B

%filter kernel
T=reshape(full(F(((I*J)+1)/2,:)),I,J);
clear F
[y,x]=find(T>1e-10);
T=T(min(y):max(y),min(x):max(x));
[~,J]=size(P);


t=0;
ft=0;%State save iteration
S=cell(1,floor(time/framestep)); %States storage cell
times=zeros(1,floor(time/framestep)); %Times storage list
h=waitbar(0,['Numerical time: ' num2str(t*dt) ' s, of total ' num2str(time) ' s']);
while t*dt<=time
    if t*dt>=ft*framestep
        S{1,ft+1}=P(:,2:J-1); %Save state
        times(ft+1)=t*dt; %Save time
        ft=ft+1;
    end
    P=imfilter(P,T,'replicate');
    P(:,1)=dP;      %left BC
    P(:,J)=0;       %right BC
    P(Frac)=dP;     %keep pressure in fracture
    t=t+1;
    waitbar((t*dt)/time,h,['Numerical time: ' num2str(t*dt) ' s, of total ' num2str(time) ' s']);
end
S{1,ft+1}=P; %Save an extra state just in case
times(ft+1)=t*dt; %Save the extra time 
close(h)
end