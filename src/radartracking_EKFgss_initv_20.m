clear all;
dt=0.05;
N=20/dt;
fid=fopen('radardata_gss_initv_20.txt','w');
%Generate state truth and measurement.
TXk=zeros(3,N);
Wr=3*randn(N,1);
Mk=zeros(N,1);
TX0=[0;100;1000];
TimeK=zeros(N,1);
for k=1:N
    TXk(:,k)=[k*dt*100;100;TX0(3)];
    Mk(k)=sqrt(TXk(1,k)^2+TXk(3,k)^2)+Wr(k);
    TimeK(k)=k;
end
%EKF process
rkR=9;
rkQ=[0.1*0.1*dt*dt 0 0;
    0 0.1*0.1*dt*dt 0;
    0  0  0];
Fk=[1 dt 0;
    0 1 0;
    0 0 1];
i3=eye(3);
EXk=zeros(3,N);
tempEX=zeros(3,N);
EX0=[-100;20;2000];
Pk=zeros(3,3,N+1);
Pk(:,:,1)=[49 0 0;0 49 0;0 0 49];
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
    EXk(:,k)=tempEX(:,k)+K*(Mk(k)-denorm);
end
%EKF in GSS Process
rkRgss=[9,0;0,9];
rkQgss=[0.1*0.1*dt*dt 0 0 0 ;
    0    0.1*0.1*dt*dt 0 0 ;
    0    0 0.1*0.1*dt*dt 0;
     0  0  0 0 ];
Fkgss=[1 0 dt 0;
    0 1 dt 0;
    0 0 1 0;
    0 0 0 1];
i4=eye(4);
EXkgss=zeros(4,N);
tempEXgss=zeros(4,N);
EXkgss(:,1)=[EXk(1,1)-EXk(2,1)*dt;EXk(1,1);EXk(2,1);EXk(3,1)];
Pkgss=zeros(4,4,N);
Pkgss(:,:,1)=[49 0 0 0;0 Pk(1,1,2) 0 0;0 0 Pk(2,2,2) 0;0 0 0 Pk(3,3,2)];
for k=2:N
    tempEXgss(:,k)=Fkgss*EXkgss(:,k-1);
    tempPgss=Fkgss*Pkgss(:,:,k-1)*Fkgss'+rkQgss;
    denorm1=sqrt(tempEXgss(1,k)^2+tempEXgss(4,k)^2);
    denorm2=sqrt(tempEXgss(2,k)^2+tempEXgss(4,k)^2);
    Hkgss=[tempEXgss(1,k)/denorm1 0 0 tempEXgss(4,k)/denorm1;
    0 tempEXgss(2,k)/denorm2 0  tempEXgss(4,k)/denorm2];
    Kgss=tempPgss*Hkgss'*inv(Hkgss*tempPgss*Hkgss'+rkRgss);
    Pkgss(:,:,k)=(i4-Kgss*Hkgss)*tempPgss;
    EXkgss(:,k)=tempEXgss(:,k)+Kgss*[Mk(k)-denorm1;Mk(k-1)-denorm2];
end
figure(1);
plot(TimeK*dt,EXk(1,:)-TXk(1,:),TimeK*dt,EXkgss(1,:)-TXk(1,:));
ylabel('The Estimation Error of X');xlabel('Time(s)');
legend('EKF','GSSM');%
figure(2);
plot(TimeK,EXk(3,:),TimeK,EXkgss(4,:));
ylabel('The Estimation of Height');xlabel('Time(s)');
legend('EKF','GSSM');%
figure(3);
plot(TimeK,EXk(2,:),TimeK,EXkgss(3,:));
ylabel('The Estimation of Velocity');xlabel('Time(s)');
legend('EKF','GSSM');%
 for k=1:N
     fprintf(fid,'%f %f %f %f %f %f %f %f %f %f  %f %f  %f\n',k*dt,TXk(1,k),TXk(2,k),Mk(k),EXk(1,k),EXk(1,k)-TXk(1,k),EXk(2,k),EXk(3,k),EXkgss(1,k),EXkgss(1,k)-TXk(1,k),EXkgss(2,k),EXkgss(3,k),EXkgss(4,k));
 end
 fclose(fid);
