function [S,iters]=Laplace_2D(data,dP,mmres,mmperpx,framestep,P_initial)

scale=mmperpx/mmres;
P=double(fliplr(imresize(data,scale,'nearest')));
if ~exist('P_initial')
    P_initial=zeros(size(P));
else
    P_initial=imresize(P_initial,size(P),'nearest');
end

%Boundary functions
P=padarray(P,[0 1],0,'both');
Frac=(P==1);    %fracture
[~,J]=size(P);
P=padarray(P_initial,[0 1],0,'both');
P(:,1)=dP;       %left
P(:,J)=0;        %right
P(Frac)=dP;      %set pressure

%Filter kernel
T=[0 1 0;1 0 1;0 1 0]/4;

iter=0;
conv=1e-8;%Convergence tolerance
rmse=1; %Initial error
ft=0;%State save iteration
S={};%States storage cell
iters=[]; %Iterations storage list
N=size(P,1)*size(P,2);%Number of grid points
h=waitbar(0,['Solving Laplace equation, Iteration: ' num2str(iter) '. RMSE: ' num2str(rmse)]);
while rmse>=conv % convergence criterion
    Ptmp=P; %Store current state
    if iter>=ft*framestep
        S{1,ft+1}=P(:,2:J-1); %Save state
        iters(ft+1)=iter; %Save time
        ft=ft+1;
    end
    P=imfilter(P,T,'replicate');
    P(:,1)=dP;      %left BC
    P(:,J)=0;       %right BC
    P(Frac)=dP;     %keep pressure in fracture
    rmse=sum(sqrt((P(:)-Ptmp(:)).^2)/(N*dP)); %Normalized RMSE between frames
    iter=iter+1;
    waitbar(conv/rmse,h,['Solving Laplace equation, Iteration: ' num2str(iter) '. RMSE: ' num2str(rmse)]);
end
S{1,ft+1}=P(:,2:J-1); %Save an extra state just in case
iters(ft+1)=iter; %Save the extra time
close(h)
end