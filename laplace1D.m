U=zeros(1000,1);%Column vector for positions, the grid length is 1000*dx
U(1)=1;%Left BC
U(length(U))=0;%Right BC

for i=1:1000
    plot(U);
    M(i)=getframe;
    L=U(1:length(U)-2);
    R=U(3:length(U));
    C=(L+R)/2;
    U=[1;C;0];
end
plot(U)
M(i+1)=getframe;
