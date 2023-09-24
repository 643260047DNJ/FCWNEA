%% 种群初始化
function [pop,Ind_dec] = initial_pop(val_net,w_net,N)
    Dec=[];
    Ind_dec=[];
    for i=1:N
       dec=[];
       mark = routlette(w_net{1});
       ind_dec=[mark];
       dec = [dec, val_net{1}(mark) ];
       for j=2:length(val_net)
          w = w_net{j}(mark,:);
          mark = routlette(w);
          ind_dec = [ind_dec,mark];
          dec = [dec, val_net{j}(mark) ];
       end
       Dec = [Dec;dec];
       Ind_dec = [Ind_dec;ind_dec];
    end
    pop=SOLUTION(Dec);
end

    
%% 轮盘选择
function Ind = routlette(P)
    pro = P.^2;
%     pro = P;
    miu = rand()*sum(pro);
    i=1;
    miu = miu - pro(i);
    while miu>0
        i=i+1;
        miu = miu - pro(i);
    end
    Ind = i;
end
