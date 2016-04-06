P=double(imresize(fliplr(data),0.1,'nearest'));
P=padarray(P,[1 1],0,'both');
[I,J]=size(P);
dP=150;%kPa

%Boundary functions
P(1,:)=P(2,:);
P(size(P,1),:)=P(size(P,1)-1,:);
P(:,1)=1;
P(:,size(P,2))=0;
Frac=(P==1);
P(Frac)=dP;
dx=7.48;%x-step
dt=0.001;%t-step
D=3.2e5;%Diffusion constant mm^2/s

if dt>((dx^2)/(2*D))
    dt=((dx^2)/(2*D));
end
a=(D*dt)/(dx^2);

d=[-I,-1,0,1,I]';
D0=ones(I*J,1)*((4*a)+2);
D1=ones(I*J,1)*(-a);
Diags=cat(2,D1,D1,D0,D1,D1);
A=spdiags(Diags,d,I*J,I*J);
d=[-I,-1,0,1,I]';
D0=ones(I*J,1)*(2-(4*a));
D1=ones(I*J,1)*(a);
Diags=cat(2,D1,D1,D0,D1,D1);
B=spdiags(Diags,d,I*J,I*J);
clear D0 D1 Diags d
F=A\B;
clear A B

pvec=P(:);
t=0;
while t*dt<0.096% t=1:10000
    imagesc(P)
    M(t+1)=getframe;
    pvec=F*pvec;
    P=(reshape(pvec,I,J));
    %Boundary functions
    P(1,:)=P(2,:);
    P(size(P,1),:)=P(size(P,1)-1,:);
    P(:,1)=1;
    P(:,size(P,2))=0;
    P(Frac)=dP;
    pvec=P(:);
    t=t+1;
end

%%%% SUB BLOCKS
A_diagblock=A(1:I,1:I);
A_subblock=A(I+1:2*I,1:I);
A_subblock(1,I)=0;
A_superblock=A(1:I,I+1:2*I);
A_superblock(I,1)=0;
%imagesc(A_diagblock),figure,imagesc(A_subblock),figure,imagesc(A_superblock)
%imagesc(blcktridiag(A_diagblock,A_subblock,A_superblock,2))
B_diagblock=B(1:I,1:I);
B_subblock=B(I+1:2*I,1:I);
B_subblock(1,I)=0;
B_superblock=B(1:I,I+1:2*I);
B_superblock(I,1)=0;
A=blktridiag(A_diagblock,A_subblock,A_superblock,2);
B=blktridiag(B_diagblock,B_subblock,B_superblock,2);
F=A\B;
F_diagblock=F(1:I,1:I);
F_subblock=F(I+1:2*I,1:I);
F_superblock=F(1:I,I+1:2*I);
F=blktridiag(F_diagblock,F_subblock,F_superblock,J);
imagesc(F)