function M=Crank_Nich_1D(dx,dt,D,time)

% This function solves the Diffusion equation by the Crank-Nicholson method
% and returns the evolution in a movie. The left boundary has constant
% concentration of 1 and the right boundary is always 0.
% inputs: dx=position step, dt=time step, D=diffusion constant,
%         time=duration
% output: M=Matlab movie

U=zeros(1000,1);%Column vector for positions, the grid length is 1000*dx
U(1)=1;%Left BC
U(length(U))=0;%Right BC
if dt>((dx^2)/(2*D))
    dt=((dx^2)/(2*D)); %Defines a maximum timestep
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

%Tridiagonal matrix for vector equation to find concentrations for dt+1
F=A\B;

%Fix smeared out diagonals
test=F(500,499:501);
test=test+(1-sum(test))/3;
test=repmat(test,size(F,1),1);
sub=test(2:size(test,1),1); %sub diagonal
main=test(1:size(test,1),2);%main diagonal
sup=test(2:size(test,1),3); %super diagonal
F=gallery('tridiag',sub,main,sup);

%Make boundary point equations trivial (keep left/right boundaries static)
F(1,:)=0;
F(1,1)=1;
F(size(F,1),:)=0;
F(size(F,1),size(F,2))=1;

t=0;
h=waitbar(0,['Numerical time: ' num2str(t*dt) ' s']);
h.Position=[320.25,36.75,270,57];
while t*dt<time
    plot(U)
    M(t+1)=getframe;
    t=t+1;
    U=F*U;
    waitbar((t*dt)/time,h,['Numerical time: ' num2str(t*dt) ' s']);
end
close(h)
end