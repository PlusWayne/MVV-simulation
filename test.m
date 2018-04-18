clear all

iterationForh=50;
iterationForuser=10;
beta=3;         %beta参数
K=[1 2 3 4 5 6 7 8 9 10];           %用户集合

y1=zeros(6,iterationForh,iterationForuser);
y2=zeros(6,iterationForh,iterationForuser);
y3=zeros(6,iterationForh,iterationForuser);
y4=zeros(6,iterationForh,iterationForuser);
% y5=[];
n_0=1;          %白噪声参数
E_b=5;          %基站合成view能量
E_u=5*ones(1,length(K));       %用户端合成view需要的能量，1*K个不同的值
precision=0.1;  %离散情况下两个view之间的间隔
V=[1:0.1:5];           %一共有多少个view (离散情况) 1:precision:N
delta=ones(1,length(K));       %每个用户接受的精度范围
T=0.1;            %所有任务要在T内完成
D=1;            %数据块大小
h=sqrt(1/2)*(randn(iterationForh,length(K))+sqrt(-1)*randn(iterationForh,length(K)));         %信道参数 1*K个
for B=1:6
    
for tt=1:iterationForh
for j=1:iterationForuser
     r=V(randi(length(V),1,length(K)));         %用户请求的view 1*K个
     y1(B,tt,j)=mainWithPenalty(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
     y2(B,tt,j)=main2(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
     y3(B,tt,j)=mainAllsynthesized(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
     y4(B,tt,j)=mainRound(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
%    y5=[y5,optsearch(i,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)];
end
end
save(['B=',num2str(B)]);
end
save('All')
% plot(B,y1,'*-')
% hold on
% plot(B,y2,'^-')
% hold on
% plot(B,y3,'o-')
% hold on
% plot(B,y4,'x-')
% hold on
% plot(B,y5,'^-')
% a=mean(mean(y2,3),2);
% plot(1:3,a(1:3))
% hold on
% a=mean(mean(y3,3),2);
% plot(1:3,a(1:3))
% a=mean(mean(y4,3),2);
% plot(1:3,a(1:3))