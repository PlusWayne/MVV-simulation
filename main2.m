function P=main2(B,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)
%% 考虑用户全部不合成,请求什么发什么
n=length(V);
m=length(K);
y=zeros(m,n);
for i=1:length(K)
    y(i,int8((r(i)-1)/precision+1))=1;
end

%% cvx 
cvx_begin quiet
    variables t(n) z(m,n)
    for v=1:length(V)
        for k=1:length(K)
        temp1(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
        end
    end
    term1 = sum(max(temp1));
    %优化
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
        for v=1:length(V)
            t(v)>=0;
        end
cvx_end
term2=0;
term3=0; %用户端不会合成
for j=1:length(r)
    if(floor(r(j))~=r(j))
        term2=term2+E_b;
    end
end
P=cvx_optval+term2;