function [r1,E_u1,delta1,h1,K1]=sortAllvector(r,E_u,delta,h,K)
% r=[1 1.5 2 1.5];
% E_u=[1 2 3 4];
% delta=[0.5 0.3 0.2 0.6];
% h=[2 1 3 4];
% K=[1 2 3 4];

r=r';
E_u=E_u';
delta=delta';
h=h';
K=K';
h_abs=abs(h);

A=sortrows([r,E_u,delta,K,h,h_abs],5);

r1=A(:,1)';
E_u1=A(:,2)';
delta1=A(:,3)';
K1=A(:,4)';
h1=A(:,5)';