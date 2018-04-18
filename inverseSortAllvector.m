function [r1,E_u1,delta1,h1,y0]=inverseSortAllvector(r,E_u,delta,h,K,y)
% r=[1 1.5 2 1.5];
% E_u=[1 2 3 4];
% delta=[0.5 0.3 0.2 0.6];
% h=[2 1 3 4];

r=r';
E_u=E_u';
delta=delta';
h=h';
K=K';
A=sortrows([r,E_u,delta,h,K,y],5);

r1=A(:,1)';
E_u1=A(:,2)';
delta1=A(:,3)';
h1=A(:,4)';
y0=A(:,6:end);