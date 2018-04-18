function P=mainAllsynthesized(B,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)
n=length(V);
m=length(K);
y=zeros(m,n);
for i=1:length(K)
    y(i,int8((floor(r(i))-1)/precision+1))=1;
    y(i,int8((ceil(r(i))-1)/precision+1))=1;
end


%% cvx 
cvx_begin quiet
    variables t(n) z(m,n)
    %优化
    for v=1:length(V)
        for k=1:length(K)
        temp1(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
        end
    end
    term1 = sum(max(temp1));
    minimize (term1)
    subject to
       % constraint 1
        sum(t) == T;
       % constraint 2
        for k=1:length(K)
            for v=1:length(V)
                (log(2)*y(k,v)*D*n_0)/(B*abs(h(k))^2)+rel_entr(n_0*t(v)/abs(h(k))^2,z(k,v)) <=0;
            end
        end
        % t_v>0  <=> t_v>=0
        tempSum=0;
        for v=1:length(V)
            t(v)>=0;
        end
cvx_end
term2=0; %基站端不会合成
term3=0; 
for j=1:length(r)
    if(floor(r(j))~=r(j))
        term3=term3+beta*E_u(j);
    end
end
P=cvx_optval+term3;