clear all

iterationForh=10;
iterationForuser=10;
beta=3;         %beta参数
y1=zeros(6,iterationForh,iterationForuser);
% y2=zeros(6,iterationForh,iterationForuser);
% y3=zeros(6,iterationForh,iterationForuser);
y4=zeros(6,iterationForh,iterationForuser);
y5=zeros();
n_0=1;          %白噪声参数
precision=0.1;  %离散情况下两个view之间的间隔
T=0.1;            %所有任务要在T内完成
D=1;            %数据块大小
V=[1:0.1:5];           %一共有多少个view (离散情况) 1:precision:N
E_b=5;          %基站合成view能量
t1=[];
t4=[];
t5=[];
for k=2:4
K=1:k;           %用户集合
h=sqrt(1/2)*(randn(iterationForh,length(K))+sqrt(-1)*randn(iterationForh,length(K)));         %信道参数 1*K个
E_u=5*ones(1,length(K));       %用户端合成view需要的能量，1*K个不同的值
delta=ones(1,length(K));       %每个用户接受的精度范围
for tt=1:iterationForh
for j=1:iterationForuser
     r=V(randi(length(V),1,length(K)));         %用户请求的view 1*K个
	 tic
     y1(k,tt,j)=mainWithPenalty(10,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
	 t1=[t1 toc];
%      y2(k,tt,j)=main2(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
%      y3(k,tt,j)=mainAllsynthesized(B+9,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
	 tic
     y4(k,tt,j)=mainRound(10,beta,h(tt,:),r,n_0,E_b,E_u,precision,V,K,delta,T,D);
	 t4=[t4 toc];
	 tic
     y5(k,tt,j)=optsearch(10,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D);
	 t5=[t5 toc];
end
save(['K=',num2str(k)]);
end
end
save('All_k')