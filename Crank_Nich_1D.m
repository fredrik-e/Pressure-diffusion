function [S,times,pos]=Crank_Nich_1D(L,points,dt,D,time,framestep)

% This function solves the 1-D Diffusion equation by the Crank-Nicholson 
% method and returns the evolution in a movie. The left boundary has 
% constant concentration of 1 and the right boundary is always 0.
% inputs: L=length of grid, points=number of grid points, dt=time step, 
%         D=diffusion constant,
%         time=duration, framestep=timestep between state records
% output: S=saved concentration states over time (cell) 
%         times=times of saved states in S
%         pos=position vector for grid points in S

if round(points/2)*2==points
    points=points+1; %Ensure odd number of points
end
cpos=(points+1)/2; %Center position of grid

U=zeros(points,1);%Column vector for grid, the grid length is 1000*dx
U(1)=1;%Left BC
U(length(U))=0;%Right BC
dx=L/(points-1);
pos=(0:length(U)-1)*dx;%Position vector for gridpoints
if dt>((dx^2)/(2*D))
    dt=((dx^2)/(2*D)); %Defines a maximum timestep
end
a=(D*dt)/(dx^2);

%Tridiagonal matrix for left side of equation
Aa=ones(1,3*length(U)-1)*(-a);%sub diagonal entries for interior points
Ab=ones(1,3*length(U))*((2*a)+2);%diagonal entries for interior points
Ac=Aa;%super diagonal entries for interior points
A=gallery('tridiag',Aa,Ab,Ac); %Tridiagonal matrix

%Tridiagonal matrix for right side of equation
Ba=ones(1,3*length(U)-1)*(a);%sub diagonal entries for interior points
Bb=ones(1,3*length(U))*(2-(2*a));%diagonal entries for interior points
Bc=Ba;%super diagonal entries for interior points
B=gallery('tridiag',Ba,Bb,Bc); %Tridiagonal matrix

%Tridiagonal matrix for vector equation to find concentrations for dt+1
F=A\B;
F=F(length(U)+1:2*length(U),length(U)+1:2*length(U));

%Implement boundary conditions in weighting matrix
F(1,:)=0;
F(1,1)=1;
F(size(F,1),:)=0;
F(size(F,1),size(F,2))=1;
R=F(cpos,cpos+1:points);
for i=2:length(R)
    F(i,1)=sum(R(i-1:length(R))); %Weights left of matrix into col. 1
end
for i=size(F,2)-length(R)-1:size(F,2)-1
    F(i,size(F,2))=sum(R(size(F,2)-i:length(R)));%Weights right of F
end

t=0; %Time iteration
ft=0;%State save iteration
h=waitbar(0,['Numerical time: ' num2str(t*dt) ' s, of total ' num2str(time) ' s']);
S=cell(1,floor(time/framestep));
times=zeros(1,floor(time/framestep));
while t*dt<time
    if t*dt>=ft*framestep
        S{1,ft+1}=U; %Save state
        times(ft+1)=t*dt; %Save time
        ft=ft+1;
    end
    t=t+1;
    U=F*U;
    waitbar((t*dt)/time,h,['Numerical time: ' num2str(t*dt) ' s, of total ' num2str(time) ' s']);
end
S{1,ft+1}=U; %Save state
times(ft+1)=t*dt; %Save time
close(h)
end