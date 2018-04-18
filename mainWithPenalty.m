function P=mainWithPenalty(B,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)
%% 参数(离散情况下)
rho=10;  %惩罚函数系数
n=length(V);
m=length(K);
optimal_y=[];
optimal_value=inf;
%% cvx开始迭代优化 随机取10个可行点，分别进行迭代，取结果最好的
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
%     term2=0;
    for i=1:length(V)
        if V(i)~=floor(V(i))
         temp2(i)=max(y(:,int8((V(i)-1)/precision+1)))*E_b;
        end
    end
    term2=sum(temp2);
    %目标函数第三项
      term3=0;
    for i=1:length(K)
        temp3=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
        for j=1:length(temp3)
            if(temp3(j)~=r(i))
                term3=term3+y(i,int8((temp3(j)-1)/precision+1))*E_u(i);
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
y0=y;
indicator =0;
%% 开始迭代优化，理论上一定向减小方向进行
while 1
clear temp1;
cvx_begin quiet
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
        temp3=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
        for j=1:length(temp3)
            if(temp3(j)~=r(i))
                term3=term3+y(i,int8((temp3(j)-1)/precision+1))*E_u(i);
            end
        end
    end
    term3=beta*term3;
    %penalty 可以考虑之前的系数
    penalty=0;
    for i=1:length(K)
       for j=1:length(V)
           penalty=penalty+(1-2*y0(i,j))*y(i,j)+y0(i,j)^2;
       end
    end
    %优化
    minimize (term1+term2+term3+rho*penalty)
    subject to
       % constraint 1
        sum(t) == T;
       % constraint 2
        for k=1:length(K)
            for v=1:length(V)
                (log(2)*y(k,v)*D*n_0)/(B*abs(h(k))^2)+rel_entr(n_0*t(v)/abs(h(k))^2,z(k,v)) <=0; % 利用转换可以写成这样
            end
        end
        % t_v>0  <=> t_v>=0
        tempSum=0;
        for v=1:length(V)
             t(v)>=0;
        end
        % y_{k,v} 左侧
        for i=1:length(r)
        left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);% 考虑左侧可行的区间，利用unique函数去掉重复的值
            for j=1:length(left)
               tempSum=tempSum+y(i,int8((left(j)-1)/precision+1));%转换到实际矩阵中的位置，转换公式(r-1)/precision+1
            end
            tempSum == 1;% 对于每个左侧加和要等于1
            tempSum=0;   %重新初始化这两个变量以用来下次变量使用，用完之后一定是0
            left=[];
        end
        % y_{k,v} 右侧
        for i=1:length(r)
        right=unique([r(i):precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]); % 考虑右侧可行的区间，利用unique函数去掉重复的值
            for j=1:length(right)
               tempSum=tempSum+y(i,int8((right(j)-1)/precision+1)); %转换到实际矩阵中的位置，转换公式(r-1)/precision+1
            end
            tempSum == 1;  % 对于每个右侧加和要等于1
            tempSum=0;    %重新初始化这两个变量以用来下次变量使用，用完之后一定是0
            right=[];
        end
         for i=1:length(K)
             for j=1:length(V)
                 0<= y(i,j) <=1;
             end
         end
cvx_end
 y0=y;
%% 重新计算目标函数值
%     %目标函数第一项
%     for v=1:length(V)
%         for k=1:length(K)
%         temp1(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
%         end
%     end
%     term1 = sum(max(temp1));
%     %目标函数第二项
%     term2=0;
%     for i=1:length(V)
%         if V(i)~=floor(V(i))
%          term2=term2+max(y(:,int8((V(i)-1)/precision+1)))*E_b;
%         end
%     end
%     %目标函数第三项
%     term3=0;
%     for i=1:length(K)
%         temp3=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
%         for j=1:length(temp3)
%             if(temp3(j)~=r(i))
%                 term3=term3+y(i,int8((temp3(j)-1)/precision+1))*E_u(i);
%             end
%         end
%     end
%     term3=beta*term3;
%     %penalty 可以考虑之前的系数
%     penalty=0;
%     for i=1:length(K)
%        for j=1:length(V)
%            penalty=penalty+(1-2*y0(i,j))*y(i,j)+y0(i,j)^2;
%        end
%     end
%     cvx_optval=term1+term2+term3+rho*penalty;
    indicator = indicator +1;
    %% 比较
    if indicator ==5
        break
    end
    if abs(cvx_optval-optimal_value) > 0.1
    
        if cvx_optval<=optimal_value
        optimal_y=y;
        optimal_t=t;
        optimal_value=cvx_optval;
        end
    else
        if cvx_optval<=optimal_value
        optimal_y=y;
        optimal_t=t;
        optimal_value=cvx_optval;
        end
        break;
    end
end
P=optimal_value-rho*penalty;



% left=[];
% right=[];
%% 对最后的解进行适当的修改
% for i=1:length(K)
%   left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
%   leftNum=int8((left-1)./precision+1);
%   right=unique([r(i):precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
%   rightNum=int8((right-1)./precision+1);
%     for j=1:length(V)
%      if ismember(j,leftNum)==0 && ismember(j,rightNum)==0
%          y(i,j)=0;
%      end
%     end
%     left=[];
%     leftNum=[];
%     right=[];
%     rightNum=[];
% end

% for i=1:length(K)
%   left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)]);
%   leftNum=int8((left-1)./precision+1);
%   right=unique([r(i):precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
%   rightNum=int8((right-1)./precision+1);
%     for j=1:length(V)
%         l1=find(y(i,leftNum)==max(y(i,leftNum)));
%         r1=find(y(i,rightNum)==max(y(i,rightNum)));
%        if j==leftNum(l1(1)) || j==rightNum(r1(1))
%         y(i,j)=1;
%        else
%         y(i,j)=0;
%        end
%     end
%     left=[];
%     leftNum=[];
%     right=[];
%     rightNum=[];
% end

%% 重新求解
% cvx_begin quiet
%     clear temp1 temp2
%     variables t(n) z(m,n)
%     %优化
%     %目标函数第一项
%     for v=1:length(V)
%         for k=1:length(K)
%         temp2(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
%         end
%     end
%     term1 = sum(max(temp2));
%     minimize (term1)
%     subject to
%        % constraint 1
%         sum(t) == T;
%        % constraint 2
%         for k=1:length(K)
%             for v=1:length(V)
%                 (log(2)*y(k,v)*D*n_0)/(B*abs(h(k))^2)+rel_entr(n_0*t(v)/abs(h(k))^2,z(k,v)) <=0;
%             end
%         end
%         % t_v>0  <=> t_v>=0
%         for v=1:length(V)
%             t(v)>=0;
%         end
% cvx_end
% for v=1:length(V)
%     for k=1:length(K)
%     temp2(k,v)=z(k,v)-n_0/abs(h(k))^2*t(v);
%     end
% end
% term1 = sum(max(temp2));
% term2=0; %基站端
% term3=0; %用户端
% x=max(y);
% sendingNum=find(x==1);
% sending=(sendingNum-1)*precision+1;
% for j=1:length(r)
%     if(ismember(int8((r(j)-1)/precision+1),sendingNum)==0)
%         term3=term3+beta*E_u(j);
%     end
% end
% 
% for j=1:length(sending)
%     if(floor(sending(j))~=sending(j))
%         term2=term2+E_b;
%     end
% end
% P=term1+term2+term3;

