clear all

iterationForh=10;
iterationForuser=10;
beta=3;         %beta����
y1=zeros(6,iterationForh,iterationForuser);
% y2=zeros(6,iterationForh,iterationForuser);
% y3=zeros(6,iterationForh,iterationForuser);
y4=zeros(6,iterationForh,iterationForuser);
y5=zeros();
n_0=1;          %����������
precision=0.1;  %��ɢ���������view֮��ļ��
T=0.1;            %��������Ҫ��T�����
D=1;            %���ݿ��С
V=[1:0.1:5];           %һ���ж��ٸ�view (��ɢ���) 1:precision:N
E_b=5;          %��վ�ϳ�view����
t1=[];
t4=[];
t5=[];
for k=2:4
K=1:k;           %�û�����
h=sqrt(1/2)*(randn(iterationForh,length(K))+sqrt(-1)*randn(iterationForh,length(K)));         %�ŵ����� 1*K��
E_u=5*ones(1,length(K));       %�û��˺ϳ�view��Ҫ��������1*K����ͬ��ֵ
delta=ones(1,length(K));       %ÿ���û����ܵľ��ȷ�Χ
for tt=1:iterationForh
for j=1:iterationForuser
     r=V(randi(length(V),1,length(K)));         %�û������view 1*K��
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