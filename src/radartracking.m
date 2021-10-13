dt=0.05;
N=20/dt;
i3=eye(3);
TXk=zeros(3,N);
EXk=zeros(3,N);
tempEX=zeros(3,N);
TX0=[0;100;1000];
EX0=[-100;200;2000];
Wr=3*randn(N,1);
Mk=zeros(N,1);
Pk=zeros(3,3,N+1);
Pk(:,:,1)=[49 0 0;0 49 0;0 0 49];
fid=fopen('radardata2.txt','w');
rkR=25;
rkQ=[0.1*0.1*dt*dt 0 0;
    0 0.1*0.1*dt*dt 0;
    0  0  0];
Fk=[1 dt 0;
    0 1 0;
    0 0 1];
for k=1:N
    TXk(:,k)=[k*dt*100;100;TX0(3)];
    Mk(k)=sqrt(TXk(1,k)^2+TXk(3,k)^2)+Wr(k);
end
for k=1:N
    if k==1
        tempEX(:,k)=Fk*EX0;
    else
        tempEX(:,k)=Fk*EXk(:,k-1);
    end
    tempP=Fk*Pk(:,:,k)*Fk'+rkQ;
    denorm=sqrt(tempEX(1,k)^2+tempEX(3,k)^2);
    Hk=[tempEX(1,k)/denorm 0 tempEX(3,k)/denorm];
    K=tempP*Hk'*inv(Hk*tempP*Hk'+rkR);
    Pk(:,:,k+1)=(i3-K*Hk)*tempP;
    %EXk(:,k)=tempEX(:,k)+K*(Mk(k)-Hk*tempEX(:,k));
    EXk(:,k)=tempEX(:,k)+K*(Mk(k)-denorm);
end
 for k=1:N
     fprintf(fid,'%f %f %f %f %f %f %f %f\n',k*dt,TXk(1,k),TXk(2,k),Mk(k),EXk(1,k),EXk(2,k),EXk(3,k),sqrt(EXk(3,k)));
 end
