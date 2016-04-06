
% This function solves the Diffusion equation by the Crank-Nicholson method

U=zeros(1000,1);%Column vector of values
LBC=1;%Left BC
RBC=0;%Right BC
U(1)=LBC;
U(length(U))=RBC;
dx=1;%x-step
dt=0.01;%t-step
D=5;%Diffusion constant
if dt>((dx^2)/(2*D))
    dt=((dx^2)/(2*D));
end
a=(D*dt)/(dx^2);

%Tridiagonal matrix for left side of equation
Aa=ones(1,length(U)-1)*(-a);%sub diagonal entries for interior points
Ab=ones(1,length(U))*((2*a)+2);%diagonal entries for interior points
Ac=Aa;%super diagonal entries for interior points
A=gallery('tridiag',Aa,Ab,Ac); %Tridiagonal matrix

%Tridiagonal matrix for right side of equation
Ba=ones(1,length(U)-1)*(a);%sub diagonal entries for interior points
Bb=ones(1,length(U))*(2-(2*a));%diagonal entries for interior points
Bc=Ba;%super diagonal entries for interior points
B=gallery('tridiag',Ba,Bb,Bc); %Tridiagonal matrix

%Tridiagonal matrix for vector equation
F=A\B;
%Make boundary point equations trivial
F(1,:)=0;
F(1,1)=1;
F(size(F,1),:)=0;
F(size(F,1),size(F,2))=1;

for t=1:1000
    plot(U)
    M(t)=getframe;
    U=F*U;
    %U(1)=LBC;%Left BC
    %U(length(U))=RBC;%Right BC
end