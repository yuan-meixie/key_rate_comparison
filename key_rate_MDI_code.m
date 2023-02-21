clear all
clc
tic
%The calculation method is the same as PRA.130.012402, without the random selection process
lb=[800,  1,  1,400,  1, 70,700, 500, 500,450,  1, 70];
ub=[999,400,150,999,250,450,999,3000,3000,999,250,420];
%Upper and lower boundary correspond mu_a, nu_a, omega_a, Pmu_a, Pnu_a, Pomega_a, mu_b, nu_a/nu_b, omega_a/omega_b, Pmu_b, Pnu_b, Pomega_b, 
IntCon=[1,2,3,4,5,6,7,8,9,10,11,12];
nonlcon=[];
Ldata=[0:10:450,450:1:490];
para_result=zeros(12,length(Ldata));
para_resultbest=zeros(12,length(Ldata));
Rfin=zeros(1,length(Ldata));

Fre=4e9;
hour_set=22;
eta_fiber=0.16;
eta_int=0;

d=zeros(10,1);%Experimental parameter assignment
d(1)=Fre*hour_set*3600/2;%N
d(2)=1.1;%Error correction factor
d(3)=0.1/Fre;%Dark count rate
d(4)=0;%Separate xz, z basis, misligment error rate
d(5)=0.04;%Separate xz, x basis, misligment error rate
d(6)=10^(-10);
d(7)=10^(-10);
d(8)=10^(-10);
d(9)=10^(-10);
d(10)=10^(-10);
A=[0,0,0,1,1,1,0,0,0,0,0,0;
    0,-1,1,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,1,1,1;
    0,0,0,0,0,0,0,1,-1,0,0,0];%decoy state condition, va/vb<wa/wb
b=[999;0;999;0];

p=parpool(96);
parfor j=1:length(Ldata)
    l0=zeros(2,1);
    l0(1)=0.8*10^(-eta_fiber*(Ldata(j))/20);
    l0(2)=0.8*10^(-eta_fiber*(Ldata(j))/20);
    tot=3;
    fun=@(x)mdi1(x,d,l0);
    for i=1:tot
        [para_result(:,j)]=ga(fun,12,A,b,[],[],lb,ub,nonlcon,IntCon);
        Mp=-fun(para_result(:,j))/d(1);
        if Mp>Rfin(j)
            Rfin(j)=Mp;
            para_resultbest(:,j)=para_result(:,j);
        end
    end
end
delete(p);

semilogy(Ldata,Rfin,'r')
save('Rfin_10_8_9_s.mat','Rfin')%Save code rate calculation parameters
save('Ldata_10_8_9_s.mat','Ldata')%Save distance taking parameters
save('para_resultbest_10_8_9_s.mat','para_resultbest')%Save the results of the optimal parameters of the search
toc

function R=mdi1(para,d,l)

d1(1)=d(1);%N
d1(2)=d(2);%f
d1(3)=d(3);%pd
d1(4)=d(4);%edz
d1(5)=d(5);%edx
d1(6)=l(1);%la (distance between Alice and Charlie)
d1(7)=l(2);%lb (distance between Bob and Charlie)
d1(8)=d(6);%epsilion PA
d1(9)=d(7);%epsilion cor
d1(10)=d(8);
d1(11)=d(9);
d1(12)=d(10);
R1=mdi(para,d1);
R=R1(1);

end
function R=mdi(para,d)

para1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
for i=1:3
    para1(i)=para(i)/1000;%mu_a, nu_a, omega_a
end
para1(4)=0;%Vacuum state from Alice
for j=1:3
    para1(j+4)=para(j+3)/1000;%probility Pmu_a, Pnu_a, Pomega_a
end
para1(8)=1-para1(5)-para1(6)-para1(7);%probility of vacuum state from Alice

para1(9)=para(7)/1000;
para1(10)=(para(2)/(para(8)/1000))/1000;
para1(11)=(para(3)/(para(9)/1000))/1000;
%mu_b, nu_b, omega_b

para1(12)=0;%Vacuum state from Bob
for n=1:3
    para1(n+12)=para(n+9)/1000;%probility Pmu_b, Pnu_b, Pomega_b
end
para1(16)=1-para1(13)-para1(14)-para1(15);%probility of vacuum state from Bob

R=-l(para1,d(1),d(2),d(3),d(4),d(5),d(6),...
    d(7),d(8),d(9),d(10),d(11),d(12));

end
function l=l(para,N,f,pd,edz,edx,etaa,etab,ecor,esec,e1,ec,ePA)
%para has 16 elements, light intensity u,v,w,0 and corresponding probabilities.
%16 elements of vector "para" correspond to mu_a, nu_a, omega_a, o_a Pmu_a,
%Pnu_a, Pomega_a, Po_a mu_b, nu_b, omega_b, o_b, Pmu_b, Pnu_b, Pomega_b, Po_b

l=[0,0,0,0];
%l(1) code length; l(2) phase error; l(3) Sz11 observation; l(4) Sx11 observation;
epsilon=esec/26;%Failure probability related counts

MQxww=Qxkakb(pd,edx,etaa,etab,para(3),para(11));
mxww=N*para(7)*para(15)*MQxww(2);
ML=chernoff(mxww,1,0,epsilon);
MU=chernoff(mxww,0,0,epsilon);
%The double scan method involves upper and lower bounds on the expected value of the first variable M

MQxw0=Qxkakb(pd,edx,etaa,etab,para(3),0);
MQx0w=Qxkakb(pd,edx,etaa,etab,0,para(11));
MQx00=Qxkakb(pd,edx,etaa,etab,0,0);
nxw0=N*para(7)*para(16)*MQxw0(1);
nx0w=N*para(8)*para(15)*MQx0w(1);
nx00=N*para(8)*para(16)*MQx00(1);
a1=exp(para(3))/(N*para(7)*para(16));
a2=exp(para(11))/(N*para(8)*para(15));
a3=0;
if a1<a2
    b1=a1;
    c1=nxw0;
    b2=a2;
    c2=nx0w;
else
    b1=a2;
    c1=nx0w;
    b2=a1;
    c2=nxw0;
end
if b2<a3
    b3=a3;
    c3=0;
else
    if b1<a3
        b3=b2;
        c3=c2;
        b2=a3;
        c2=0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a3;
        c1=0;
    end
end
FL=b1*chernoff(c1+c2+c3,1,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,1,0,epsilon)+...
    (b3-b2)*chernoff(c3,1,0,epsilon);
HL=FL-chernoff(nx00,0,0,epsilon)/(N*para(8)*para(16));

a4=exp(para(3))/(N*para(7)*para(16));
a5=exp(para(11))/(N*para(8)*para(15));
a6=0;
if a4<a5
    b1=a4;
    c1=nxw0;
    b2=a5;
    c2=nx0w;
else
    b1=a5;
    c1=nx0w;
    b2=a4;
    c2=nxw0;
end
if b2<a6
    b3=a6;
    c3=0;
else
    if b1<a6
        b3=b2;
        c3=c2;
        b2=a6;
        c2=0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a6;
        c1=0;
    end
end
FU=b1*chernoff(c1+c2+c3,0,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,0,0,epsilon)+...
    (b3-b2)*chernoff(c3,0,0,epsilon);
HU=FU-chernoff(nx00,1,0,epsilon)/(N*para(8)*para(16));
%The double scan method involves upper and lower bounds for the second variable H. The upper and lower bounds are obtained from the joint constraint
if MU<ML
    l(1)=0;
    l(2)=0;
    l(3)=0;
    l(4)=0;
    return
end
if HU<HL
    l(1)=0;
    l(2)=0;
    l(3)=0;
    l(4)=0;
    return
end

% %original method begin
% para_result1=zeros(2,1);
% rest(1:16)=para;
% rest(9+8)=N;
% rest(10+8)=f;
% rest(11+8)=pd;
% rest(12+8)=edz;
% rest(13+8)=edx;
% rest(14+8)=etaa;
% rest(15+8)=etab;
% rest(16+8)=ecor;
% rest(17+8)=esec;
% rest(18+8)=e1;
% rest(19+8)=ec;
% rest(20+8)=ePA;
% rest(21+8)=ML;
% rest(22+8)=MU;
% rest(23+8)=HL;
% rest(24+8)=HU;
% 
% 
% myfun=@(x)l1(x,rest);
% Rfin=myfun(para_result1(:,1));
% 
% ESx11=fESx11(para,pd,edx,etaa,etab,N,ML,...
%     HU,epsilon);
% 
% phiz11z=fphiz11(para,ESx11,N,MU,...
%     HL,epsilon);
% %original method end

%begin double-scanning
ML1=0;
HL1=0;
MU1=500;
HU1=500;
lb=[ML1,HL1];
ub=[MU1,HU1];
IntCon=[1,2];
nonlcon=[];
para_result=zeros(2,1);
para_result1=zeros(2,1);
para_resultbest1=zeros(2,1);
rest(1:16)=para;
rest(9+8)=N;
rest(10+8)=f;
rest(11+8)=pd;
rest(12+8)=edz;
rest(13+8)=edx;
rest(14+8)=etaa;
rest(15+8)=etab;
rest(16+8)=ecor;
rest(17+8)=esec;
rest(18+8)=e1;
rest(19+8)=ec;
rest(20+8)=ePA;
rest(21+8)=ML;
rest(22+8)=MU;
rest(23+8)=HL;
rest(24+8)=HU;

min=N;
tot=3;
for i=1:tot
    myfun=@(x)l1(x,rest);
    [para_result(:,1)]=ga(myfun,2,[],[],[],[],lb,ub,...
        nonlcon,IntCon);
    if min>=myfun(para_result(:,1))
        min=myfun(para_result(:,1));
        para_result1(:,1)=para_result(:,1);
    end
end
Rfin=myfun(para_result1(:,1));

para_resultbest1(1,1)=para_result1(1,1)*(MU-ML)/500+ML;
para_resultbest1(2,1)=para_result1(2,1)*(HU-HL)/500+HL;

ESx11=fESx11(para,pd,edx,etaa,etab,N,para_resultbest1(1,1),...
    para_resultbest1(2,1),epsilon);

phiz11z=fphiz11(para,ESx11,N,para_resultbest1(1,1),...
    para_resultbest1(2,1),epsilon);
%end double-scanning
phiz11=phiz11z(1);
Sz11=phiz11z(2);
Sx11=phiz11z(3);

l(1)=Rfin;
l(2)=phiz11;
l(3)=Sz11;
l(4)=Sx11;
%searching condition
% if l(1)==N
%     l(1)=0;
% end
if l(1)<0
    l(1)=0;
end
if imag(l(1))~=0
    l(1)=0;
end
if phiz11>0.5
    l(1)=0;
end

end
function l1=l1(x,rest)
%x is the vector with M, H as elements
%vector "rest" is vector "para" (for light intensity with probability parameters), N,f,pd,edz,edx,
%etaa,etab,ecor,esec,e1,ec,ePA。


para=rest(1:16);
N=rest(9+8);
f=rest(10+8);
pd=rest(11+8);
edz=rest(12+8);
edx=rest(13+8);
etaa=rest(14+8);
etab=rest(15+8);
ecor=rest(16+8);
esec=rest(17+8);
e1=rest(18+8);
ec=rest(19+8);
ePA=rest(20+8);
% %original method begin
% ML=rest(21+8);
% MU=rest(22+8);
% HL=rest(23+8);
% HU=rest(24+8);
% 
% epsilon=esec/26;%Failure probability related counts
% QEz=Qzkakb(pd,edz,etaa,etab,para(1),para(9));
% Qzuu=QEz(1);
% Ezuu=QEz(2)/QEz(1);
% nzuu=N*para(5)*para(13)*Qzuu;
% 
% Sz0u=fSz0mu(para,pd,edz,etaa,etab,N,epsilon);
% ESx11=fESx11(para,pd,edx,etaa,etab,N,ML,HU,epsilon);
% 
% phiz11z=fphiz11(para,ESx11,N,MU,HL,epsilon);
% phiz11=phiz11z(1);
% Sz11=phiz11z(2);
% C=log2(8/ecor)+2*log2(2/(e1*ec))+2*log2(1/(2*ePA));
% H1=-phiz11*log2(phiz11)-(1-phiz11)*log2(1-phiz11);
% H2=-Ezuu*log2(Ezuu)-(1-Ezuu)*log2(1-Ezuu);
% %original method end

%double scanning method begin
%16 elements of vector "para" correspond to mu_a, nu_a, omega_a, o_a Pmu_a,
%Pnu_a, Pomega_a, Po_a mu_b, nu_b, omega_b, o_b, Pmu_b, Pnu_b, Pomega_b, Po_b
M=x(1)*(rest(22+8)-rest(21+8))/500+rest(21+8);
H=x(2)*(rest(24+8)-rest(23+8))/500+rest(23+8);

epsilon=esec/26;%Failure probability related counts
QEz=Qzkakb(pd,edz,etaa,etab,para(1),para(9));
Qzuu=QEz(1);
Ezuu=QEz(2)/QEz(1);
nzuu=N*para(5)*para(13)*Qzuu;

Sz0u=fSz0mu(para,pd,edz,etaa,etab,N,epsilon);
ESx11=fESx11(para,pd,edx,etaa,etab,N,M,H,epsilon);

phiz11z=fphiz11(para,ESx11,N,M,H,epsilon);
phiz11=phiz11z(1);
Sz11=phiz11z(2);
C=log2(8/ecor)+2*log2(2/(e1*ec))+2*log2(1/(2*ePA));
H1=-phiz11*log2(phiz11)-(1-phiz11)*log2(1-phiz11);
H2=-Ezuu*log2(Ezuu)-(1-Ezuu)*log2(1-Ezuu);
%double scanning method end

l1=(Sz0u+Sz11*(1-H1)-nzuu*f*H2-C);
if Ezuu>1/2
    l1=0;
end
%searching condition
% if l1<0
%     l1=N;
% end
% if imag(l1)~=0
%     l1=N;
% end
% if phiz11>0.5
%     l1=N;
% end

end
function phiz11=fphiz11(para,ESx11,N,M,H,epsilon)
%16 elements of vector "para" correspond to mu_a, nu_a, omega_a, o_a Pmu_a,
%Pnu_a, Pomega_a, Po_a mu_b, nu_b, omega_b, o_b, Pmu_b, Pnu_b, Pomega_b, Po_b

phiz11=[0,0,0];

Etx11=M-exp(-para(3)-para(11))*N*para(7)*para(15)*H/2;
ESz11=para(1)*para(9)*para(5)*para(13)/(para(3)*para(11)*para(7)*para(15))...
    *exp((para(3)-para(1))+(para(11)-para(9)))*ESx11;
Sz11=chernoff(ESz11,1,1,epsilon);
% tx11=chernoff(Etx11*ESz11/ESx11,0,1,epsilon);
% Sx11=chernoff(ESx11,1,1,epsilon);
% 
% phiz11(1)=tx11/Sz11;
%tx11=chernoff(Etx11,0,1,epsilon);
Sx11=chernoff(ESx11,1,1,epsilon);

phiz11(1)=chernoff(Sz11*(Etx11/ESx11),0,1,epsilon)/Sz11;
phiz11(2)=Sz11;
phiz11(3)=Sx11;
end
function ESx11=fESx11(para,pd,edx,etaa,etab,N,M,H,epsilon)
%16 elements of vector "para" correspond to mu_a, nu_a, omega_a, o_a Pmu_a,
%Pnu_a, Pomega_a, Po_a mu_b, nu_b, omega_b, o_b, Pmu_b, Pnu_b, Pomega_b, Po_b

A=N*para(7)*para(15)*para(3)*para(11)*exp(-para(3)-para(11))/...
    (para(2)*para(10)*para(3)*para(11)*(para(2)-para(3)));
a=para(2)^2*para(10)*exp(para(3)+para(11))/(N*para(7)*para(15));
b=para(3)^2*para(11)*exp(para(10))/(N*para(14)*para(8));
c=para(3)^2*para(11)*exp(para(2))/(N*para(6)*para(16));
d=para(3)^2*para(11)*exp(para(2)+para(10))/(N*para(6)*para(14));
f=para(3)^2*para(11)/(N*para(8)*para(16));
m=para(2)^2*para(10)*exp(para(3)+para(11))/(N*para(7)*para(15));
h=para(2)^2*para(10);
%A is the overall factor, a is the coefficient of EUn~xww, b is the coefficient of EUnxv0, and c is the coefficient of EUnx0v.
%d is the coefficient of ELnxvv, f is the coefficient of ELnx00, m is the coefficient of M, and h is the coefficient of H.

MQxww=Qxkakb(pd,edx,etaa,etab,para(3),para(11));
nxww=N*para(7)*para(15)*MQxww(1);
mxww=N*para(7)*para(15)*MQxww(2);
nCxww=nxww-mxww;
%First variable

MQxv0=Qxkakb(pd,edx,etaa,etab,para(2),0);
nxv0=N*para(6)*para(16)*MQxv0(1);
%The 3rd variable

MQx0v=Qxkakb(pd,edx,etaa,etab,0,para(10));
nx0v=N*para(8)*para(14)*MQx0v(1);
%The 2ed variable

MQxvv=Qxkakb(pd,edx,etaa,etab,para(2),para(10));
nxvv=N*para(6)*para(14)*MQxvv(1);
%The 4th variable

MQx00=Qxkakb(pd,edx,etaa,etab,0,0);
nx00=N*para(8)*para(16)*MQx00(1);
%The 5th variable

if a<b
    b1=a;
    c1=nCxww;
    b2=b;
    c2=nx0v;
else
    b1=b;
    c1=nx0v;
    b2=a;
    c2=nCxww;
end
if b2<c
    b3=c;
    c3=nxv0;
else
    if b1<c
        b3=b2;
        c3=c2;
        b2=c;
        c2=nxv0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=c;
        c1=nxv0;
    end
end
FL=b1*chernoff(c1+c2+c3,1,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,1,0,epsilon)+...
    (b3-b2)*chernoff(c3,1,0,epsilon);
%Joint Binding
a6=0;
if d<f
    b1=d;
    c1=nxvv;
    b2=f;
    c2=nx00;
else
    b1=f;
    c1=nx00;
    b2=d;
    c2=nxvv;
end
if b2<a6
    b3=a6;
    c3=0;
else
    if b1<a6
        b3=b2;
        c3=c2;
        b2=a6;
        c2=0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a6;
        c1=0;
    end
end
FU=b1*chernoff(c1+c2+c3,0,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,0,0,epsilon)+...
    (b3-b2)*chernoff(c3,0,0,epsilon);
%The first three variables use the joint constraint, 
%and the last two variables are complemented by a zero coefficient multiplied by zero variables to continue using the joint constraint

ESx11=A*(FL-FU+m*M-h*H);

end
function Sz0mu=fSz0mu(para,pd,edz,etaa,etab,N,epsilon)
%16 elements of vector "para" correspond to mu_a, nu_a, omega_a, o_a Pmu_a,
%Pnu_a, Pomega_a, Po_a mu_b, nu_b, omega_b, o_b, Pmu_b, Pnu_b, Pomega_b, Po_b

Qz0u=Qzkakb(pd,edz,etaa,etab,0,para(9));
nz0mu=N*para(8)*para(13)*Qz0u(1);
ELnz0mu=chernoff(nz0mu,1,0,epsilon);

ELnz=exp(-para(1))*para(5)/para(8)*ELnz0mu;
nz=chernoff(ELnz,1,1,epsilon);

Sz0mu=nz;

end
function Qxkakb=Qxkakb(pd,edx,etaa,etab,ka,kb)

Qxkakb=[0,0];

if (ka~=0)&&(kb~=0)
    I0 = besseli(0,sqrt(ka*etaa*kb*etab));
    I01 = besseli(0,sqrt(ka*etaa*kb*etab)/2);
else
    I0 =1;
    I01=1;
end

y=(1-pd)*exp(-(ka*etaa+kb*etab)/4);
Qc=2*y^2*(y^2-2*y*I01+I0);
Qe=2*y^2*(1+y^2-2*y*I01);
Qxkakb(1)=(Qc+Qe)/2;
Qxkakb(2)=(edx*Qc+(1-edx)*Qe)/2;
%Qxkakb(2) is equal to Exkakb multipled by Qxkakb

end
function Qzkakb=Qzkakb(pd,edz,etaa,etab,ka,kb)

Qzkakb=[0,0];
if (ka~=0)&&(kb~=0)
    I0 = besseli(0,sqrt(ka*etaa*kb*etab));
else
    I0 =1;
end
A=I0-(1-pd)*exp(-(ka*etaa+kb*etab)/2);
B=(1-(1-pd)*exp(-(ka*etaa)/2))*(1-(1-pd)*exp(-(kb*etab)/2));
Qc=2*4*(1-pd)^2*exp(-(ka*etaa+kb*etab)/2)*B;
Qe=2*4*pd*(1-pd)^2*exp(-(ka*etaa+kb*etab)/2)*A;
Qzkakb(1)=(Qc+Qe)/4/2;
Qzkakb(2)=(edz*Qc+(1-edz)*Qe)/4/2;
%Qzkakb(1)equal to Qzkakb
%Qzkakb(2) is equal to Ezkakb multiplied by Qzkakb

end
function chernoff=chernoff(x,U,E,epsilon)
%U=0 means to find the upper limit and U=1 means to find the lower limit. 
%e=0 means from the observed value to the expected value and E=1 vice versa.

beta=log(1/epsilon);
triangleL=(beta+sqrt(8*beta*x+beta^2))/2;
triangleU=beta+sqrt(2*beta*x+beta^2);
deltaL=sqrt(2*beta/x);
deltaU=(beta+sqrt(8*beta*x+beta^2))/(2*x);
if E==0
    if U==0
        chernoff=x+triangleU;
    else
        chernoff=x-triangleL;
        if chernoff<0
            chernoff=0;
        end
    end
else
    if x==0
        chernoff=0;
        return
    end
    if U==0
        chernoff=x+deltaU*x;
    else
        chernoff=x-deltaL*x;
        if chernoff<0
            chernoff=0;
        end
    end
end
end
