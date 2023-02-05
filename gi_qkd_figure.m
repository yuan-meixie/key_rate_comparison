clear all
clc
tic
% Search functions for decoy-state QKD

Ldata=[0:10:470,471:510];
lb=[1,1,1,400,400,1,1];
ub=[400,500,300,800,1000,500,500];   
 
A=[0,0,1,-1,0,0,0;
    0,0,0,0,1,1,1];
b=[-1,999];

IntCon=[1,2,3,4,5,6,7];
nonlcon=[]; 

leng_data=length(Ldata);
para_result=zeros(leng_data,length(lb));
para_resultbest=zeros(leng_data,length(lb));
Rqkd_figure_paper=zeros(leng_data,1);

Fre=4e9;
hour_set=22;
alfa=0.16;
eta_int =2;
N=Fre*hour_set*3600;
pd=0.1/(Fre);
etad=0.8;
datapara=[alfa,eta_int,N,pd,pd,etad,etad];

for j=1:leng_data
    d=Ldata(j);
    tot=5;
    fun=@(x)qkd_figure(x,d,datapara);
    for i=1:tot
        [para_result(j,:)]=ga(fun,length(lb),A,b,[],[],lb,ub,nonlcon,IntCon);
        mmp=-fun(para_result(j,:));
        if mmp>Rqkd_figure_paper(j)
            Rqkd_figure_paper(j)=mmp;
            para_resultbest(j,:)=para_result(j,:);
        end
    end
end
plob=-log2(1-10.^(-alfa*(Ldata)/10));
semilogy(Ldata,plob,Ldata,Rqkd_figure_paper,'r')
save Rqkd_figure_paper Rqkd_figure_paper
toc



