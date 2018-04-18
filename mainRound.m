function P=mainRound(B,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)
n=length(V);
m=length(K);
%% cvx 求解松弛问题
cvx_begin quiet
%     cvx_solver mosek
    variables t(n) z(m,n)
    variable y(m,n)
    %目标函数第一项
    for v=1:length(V)
        for k=1:length(K)
        temp1(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
        end
    end
    term1 = sum(max(temp1));
    %目标函数第二项
    term2=0;
    for i=1:length(V)
        if V(i)~=floor(V(i))
         term2=term2+max(y(:,int8((V(i)-1)/precision+1)))*E_b;
        end
    end
    %目标函数第三项
      term3=0;
    for i=1:length(K)
        temp1=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
        for j=1:length(temp1)
            if(temp1(j)~=r(i))
                term3=term3+y(i,int8((temp1(j)-1)/precision+1))*E_u(i);
            end
        end
    end
    term3=beta*term3;
    %优化
    minimize (term1+term2+term3)
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
        % y_{k,v} 左侧
        for i=1:length(r)
        left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
        for j=1:length(left)
           tempSum=tempSum+y(i,int8((left(j)-1)/precision+1));
        end
           tempSum == 1;
           tempSum=0;
            left=[];
        end
        % y_{k,v} 右侧
        for i=1:length(r)
        right=unique([r(i):precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
        for j=1:length(right)
           tempSum=tempSum+y(i,int8((right(j)-1)/precision+1));
        end
           tempSum == 1;
           tempSum=0;
            right=[];
        end
         for i=1:length(K)
             for j=1:length(V)
                 0<= y(i,j) <=1;
             end
         end
cvx_end

%%
for i=1:length(K)
  left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
  leftNum=int8((left-1)./precision+1);
  right=unique([r(i):precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
  rightNum=int8((right-1)./precision+1);
    for j=1:length(V)
     if ismember(j,leftNum)==0 && ismember(j,rightNum)==0
         y(i,j)=0;
     end
    end
    left=[];
    leftNum=[];
    right=[];
    rightNum=[];
end

for i=1:length(K)
  left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)-precision]);
  leftNum=int8((left-1)./precision+1);
  right=unique([r(i)+precision:precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
  rightNum=int8((right-1)./precision+1);
  l1=find(y(i,leftNum)==max(y(i,leftNum)));
  r1=find(y(i,rightNum)==max(y(i,rightNum)));

       if y(i,int8((r(i)-1)/precision+1))>=y(i,leftNum(l1(1))) && y(i,int8((r(i)-1)/precision+1))>=y(i,rightNum(r1(1)))
          for j=1:length(V)
               y(i,j)=0;
          end
           y(i,int8((r(i)-1)/precision+1))=1;
       else
          for j=1:length(V)
               y(i,j)=0;
          end
          y(i,leftNum(l1(1)))=1;
          y(i,rightNum(r1(1)))=1;
       end
%        if j==leftNum(l1(1)) || j==rightNum(r1(1))
%         y(i,j)=1;
%        else
%         y(i,j)=0;
%        end
%       end
    left=[];
    leftNum=[];
    right=[];
    rightNum=[];
end

%% 重新求解
cvx_begin quiet
    clear temp1 temp2
    variables t(n) z(m,n)
    %优化
    %目标函数第一项
    for v=1:length(V)
        for k=1:length(K)
        temp2(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
        end
    end
    term1 = sum(max(temp2));
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
for v=1:length(V)
    for k=1:length(K)
    temp2(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
    end
end
term1 = sum(max(temp2));
term2=0; %基站端
term3=0; %用户端
x=max(y);
sendingNum=find(x==1);
sending=(sendingNum-1)*precision+1;
for j=1:length(r)
    if(ismember(int8((r(j)-1)/precision+1),sendingNum)==0)
        term3=term3+beta*E_u(j);
    end
end

for j=1:length(sending)
    if(floor(sending(j))~=sending(j))
        term2=term2+E_b;
    end
end
P=term1+term2+term3;
