function optimalE=optsearch(B,beta,h,r,n_0,E_b,E_u,precision,V,K,delta,T,D)
% beta=3;         %beta参数
% h=[1 1 1 1 1];        %信道参数 1*K个
% r=[1.3 1.5 1.6 1.8 2.5];         %用户请求的view 1*K个
% B;            %带宽，TDMA只用一个带宽
% n_0=1;          %白噪声参数
% E_b=1;          %基站合成view能量
% E_u=[1 1 1 1 1];       %用户端合成view需要的能量，1*K个不同的值
% precision=0.1;  %离散情况下两个view之间的间隔
% V=[1:0.1:3];           %一共有多少个view (离散情况) 1:precision:N
% K=[1 2 3 4 5];           %用户集合
% delta=[0.2 0.2 0.2 0.2 0.2]+0.4;       %每个用户接受的精度范围
% T=1;            %所有任务要在T内完成
% D=1;            %数据块大小
n=length(V);
m=length(K);
[r,E_u,delta,h,K]=sortAllvector(r,E_u,delta,h,K);

y = zeros(m,n);
node_left = zeros(m,m);
node_right = zeros(m,m);
for i=1:length(r)
    temNum = 2;
    node_left(i,1) = r(i);
    left=unique([floor(r(i)),max(r(i)-delta(i),V(1)):precision:r(i)-precision]);
    for j=1:length(left)
       if((any(abs(r - left(j))<0.01))||(left(j) == round(left(j))))
           node_left(i,temNum) = left(j);
           temNum = temNum + 1;
       end
    end
end

for i=1:length(r)
    temNum = 1;
    right=unique([r(i)+precision:precision:min(r(i)+delta(i),V(length(V))),ceil(r(i))]);
    for j=1:length(right)
        if((any(abs(r - right(j))<0.01)) || (right(j) == round(right(j))))
            node_right(i,temNum) = right(j);
            temNum = temNum + 1;
       end
    end
end
left_length = sum(node_left ~= 0,2);
right_length = sum(node_right ~= 0,2);
left_index = ones(m,1);
right_index = ones(m,1);
temy = zeros(n,1);
optimaly = zeros(m,n);
current_user = 1;
optimalE = 10000;
while(current_user ~= 0)
    if(left_index(current_user) == 1)
        y(current_user,:)=0;
        y(current_user,int8((node_left(current_user,left_index(current_user))-V(1))/precision + 1)) = 1;
        left_index(current_user) = left_index(current_user) + 1;
        current_user = current_user + 1;
    else
            if(left_index(current_user) <= left_length(current_user) || right_index(current_user) <= right_length(current_user))
                y(current_user,:)=0;
                
                if(current_user == 1)
                    temy = zeros(n,1);
                elseif(current_user == 2)
                    temy = y(1,:);
                else
                    temy = sum(y(1:(current_user-1),:));
                end
                
                for i = 1:length(node_left(current_user,:))
                    if((node_left(current_user,i) ~= 0) && (temy(int8((node_left(current_user,i)-V(1))/precision + 1)) ~= 0))
                        left_index(current_user) = left_length(current_user) + 2;
                        y(current_user,int8((node_left(current_user,i)-V(1))/precision + 1)) = 1;
                        if(i==1)
                        right_index(current_user) = right_length(current_user) + 2;
                        end
                        break;
                    end
                end


                for i = 1:length(node_right(current_user,:))
                    if((node_right(current_user,i) ~= 0) && (temy(int8((node_right(current_user,i)-V(1))/precision + 1)) ~= 0))
                         y(current_user,int8((node_right(current_user,i)-V(1))/precision + 1)) = 1;
                         right_index(current_user) = right_length(current_user) + 2;
                        break;
                    end
                end
                
                if(left_index(current_user) <= left_length(current_user))
                    y(current_user,int8((node_left(current_user,left_index(current_user))-V(1))/precision + 1)) = 1;
                end
                
                
                if(right_index(current_user) <= right_length(current_user))
                    y(current_user,int8((node_right(current_user,right_index(current_user))-V(1))/precision + 1)) = 1;
                end
                
                if(right_index(current_user) <= right_length(current_user))
                    right_index(current_user) = right_index(current_user) + 1;
                    current_user = current_user + 1;
                elseif(right_index(current_user) == right_length(current_user)+1)
                   left_index(current_user) = left_index(current_user) + 1;
                   if(left_index(current_user) <= left_length(current_user))
                        right_index(current_user) =  1;
                   end
                else
                        left_index(current_user) = left_index(current_user) + 1;
                        current_user = current_user + 1;
                end
            else
                left_index(current_user) = 1;
                right_index(current_user) = 1;
                current_user = current_user - 1;
            end
    end
    if(current_user == m+1)
        S_v = [];
        h_min = [];
        current_user = current_user - 1;
        
        for i = 1:n
            for j = 1:m
                if(y(j,i) == 1)
                    S_v = [S_v V(i)];
                    h_min = [h_min h(j)];
                    break;
                end
            end
        end
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
        x=max(y);
        sendingNum=find(x==1);
        sending=(sendingNum-1)*precision+1;
        term2=0;
        term3=0;
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
        E_all=cvx_optval+term2+term3;
        if(E_all <optimalE)
            optimalE = E_all;
            optimaly = y;
        end
    end
end
