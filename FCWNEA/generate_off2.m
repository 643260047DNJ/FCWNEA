%% 跟新子代个体
% Off 子代种群   Ind_off 子代变量索引  
function [Population,Ind_dec,Off,Off_ind,Gra,act_net]=generate_off2(Population,Ind_dec,val_net,w_net,act_net,lif_net,Problem,alph)
    %% 
    
    %种群去重，补齐去重后缺失的解
    N=Problem.N;
    [Population,Ind_dec] = ridoff(Population,Ind_dec,N,val_net);
    
    %计算拥挤度距离，主要是服务与结构拓展生成解时，选择拥挤度小得的解
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
%     Crowdis = CrowdingDistance(Population.decs,FrontNo);
    Crowdis = CrowdingDistance(Population.objs,FrontNo);
%     S = randperm(N,N);

%     tau = sin(8*pi.*(alph-1))/10+exp((alph-1));
    tau = exp((alph-1));
    F2 = 2*floor(0.3*N*tau/2);
    F1 = round(tau * (N-F2));
%     if rand() <0.1
%         F2 = 2*floor(N/4);
%         F1 = round( (N-F2)./2);
%     end
%     alph=alph/4+0.5;     
%     F2 = 2*ceil(50*alph/2);
%     F1 = ceil( (N-F2)/2 );  % 0->F1 个解为开发， F1->N-F2 个解为为探索 N-F2+1 ->N个解用于结构拓展
    
    % 重新调整种群中个体的先后次序
    [ss,ss_ind] = sort(Crowdis);  % 
%     [ss,ss_ind] = sort(Crowdis,'descend');  %
    Population = [Population(ss_ind(F2+1:N)),Population(ss_ind(1:F2))];
    Ind_dec = [Ind_dec(ss_ind(F2+1:N),:);Ind_dec(ss_ind(1:F2),:)];
    

    S = [randperm(N-F2,N-F2),N-F2+[1:F2]];
%         S = [randperm(N-F2,N-F2),N-F2+randperm(F2,F2)]; %生成类似matepool序列
    P_dec = Population.decs;
    Off_ind = Ind_dec;
    Off_dec = P_dec;
    vitality = abs(lif_net);
    
    %% 计算变量值域访问方差
    for i =1:length(act_net)
        if i==1
            act_var(i)=var(act_net{i});
        else
            act_var(i)=var(var(act_net{i}));
        end
    end
    
   %% 开发搜素， 搜索高权重附近的值域
     sel_var = select(vitality);

   for i =1:F1
       parent = Ind_dec(S(i),:);
       % 选择要调整的变量
       
       sel_var = select(vitality);
%        sel_var = randperm(length(vitality),1);
       for k=1:length(sel_var)
           vitality(sel_var(k)) = vitality(sel_var(k))-vitality(sel_var(k))./2;

           if sel_var(k)==1
                val=val_net{1};
                Pro = w_net{1};
           else
                val=val_net{sel_var(k)};
                if rand()<1
                    if rand() <0
                        nn=randperm(length(Off_ind(:,1)),1);
                        Pro = w_net{sel_var(k)}(Off_ind(nn,sel_var(k)-1),:);
                    else
                        Pro = w_net{sel_var(k)}(parent(sel_var(k)-1),:);
                    end             
                else
                    Pro = sum(w_net{sel_var(k)},1);
                end
           end

            Pro(parent(sel_var(k)))=0;


            mark = routlette(Pro);

            Off_ind(S(i),sel_var(k)) = mark;
            Off_dec(S(i),sel_var(k)) = val(mark);
        
       end
        
       %更新变量值域访问量 
        SS = Off_ind(S(i),:);
        for j=1:length(SS)
            if j==1
                act_net{j}(SS(j) )= act_net{1}(SS(j) ) +1;
            else
                act_net{j}(SS(j-1),SS(j)) = act_net{j}(SS(j-1),SS(j))+1;
            end
        end
        
   end
% %    
   % 探索搜索， 搜索活跃度底的值域
   for i=  F1: N % N-F2
       parent = Ind_dec(S(i),:);
       % 选择要调整的变量
%        if rand() <0
%            sel_var = select(vitality);
%        else
%            ss =abs(vitality);
%            sel_var = select(1+ max(ss)-ss);   
%        end

     
%        sel_var = select(act_var);
       sel_var = randperm(length(vitality),1);
        
    
       for k=1:length(sel_var)
           vitality(sel_var(k)) = vitality(sel_var(k))-vitality(sel_var(k))./2;
           if sel_var(k)==1
                val=val_net{1};
                Pro = w_net{1};
                Act = act_net{1};
           else
                val=val_net{sel_var(k)};
                Pro = w_net{sel_var(k)}(parent(sel_var(k)-1),:);
                Act = act_net{sel_var(k)}(parent(sel_var(k)-1),:);
           end

           % 生成变量时临时根据不同值域的访问次数  选择访问次数较少的值域
           Act = max(Act)+1 -Act;
           Act(sel_var(k))=0;
           mark = routlette(Act);
           while mark== Off_ind(S(i),sel_var(k))
               mark = routlette(Act);
           end

           if  Off_ind(S(i),sel_var(k)) == mark
                aa=1;
           end
           Off_ind(S(i),sel_var(k)) = mark;
           Off_dec(S(i),sel_var(k)) = val(mark);
       end
       %更新变量值域访问量 
       SS = Off_ind(S(i),:);
        for j=1:length(SS)
            if j==1
                act_net{j}(SS(j) )= act_net{1}(SS(j) ) +1;
            else
                act_net{j}(SS(j-1),SS(j)) = act_net{j}(SS(j-1),SS(j))+1;
            end
        end
       
   end
   
   Gra = Off_ind - Ind_dec;
    Off = SOLUTION(Off_dec);
end

%% 种群去重： 把相同的解去掉
function [Pop, Ind] = ridoff(Population,Ind_dec,N,val_net)
    obj = Population.objs;
    [obj,S_ind] = unique(obj,'row');
    Pop = Population(S_ind);
    Ind = Ind_dec(S_ind,:);
    
    ext_pop=[];
    for i =1: N-length(Pop)
        par = randperm(length(Pop),2);
        off_ind = round( (Ind_dec(par(1),:) + Ind_dec(par(2),:))/2 ) ;
        Ind = [Ind;off_ind];
        off_dec = indtodec(off_ind,val_net);
        off = SOLUTION(off_dec);

        Pop = [Pop,off];
    end

end

%% 选取活跃度最高的变量
function Ind = select(lif_net)
    m=max(lif_net);
    index = find(lif_net==m);
    a = randperm(length(index),1);
    Ind = index(a);
end

 %% 轮盘选择
function Ind = routlette(P)
%       pro = P.^2;
    pro = P;
    miu = rand()*sum(pro);
    i=1;
    miu = miu - pro(i);
    while miu>0
        i=i+1;
        miu = miu - pro(i);
    end
    Ind = i;
end

%% 索引转数值
function dec = indtodec(ind,val_net)
    dec=[];
    for i=1:length(ind)
        dec = [dec,val_net{i}(ind(i))];
    end
end
