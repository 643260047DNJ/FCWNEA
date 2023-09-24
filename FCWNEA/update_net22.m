function [w_net] = update_net22(w_net, Pop,Ind_dec,theta)
%UPDATE_NET2 此处显示有关此函数的摘要
%   此处显示详细说明
    obj=Pop.objs;
    dec=Pop.decs;
    dis_O = sum(obj.^2,2);
    dis_O = max(dis_O) - dis_O;
    w1 = dis_O./sum(dis_O);
%     w2 = Crowdis./sum(Crowdis); 
    [~,ind]= sort(w1);
    w1=w1;
    ll=round(length(ind)/2);
        normv = 4;   % 权重更新时的方差参数
    for i=1:round(length(ind)/2)
       mm = Ind_dec(ind(i),:);
       for j=1:length(mm)
           if j==1
               w_net{j} = w_net{j}+ normpdf([1:length(w_net{j})],mm(j),normv)'.^4;
           else
               w_net{j}(mm(j-1),:) = w_net{j}(mm(j-1),:)+ normpdf([1:length( w_net{j}(mm(j-1),:)  )],mm(j),normv).^4;
           end       
       end
        w_net = utility_net(w_net);
    end
end

function [nw_net]=utility_net(w_net)
    nw_net = w_net;
    for i=1:length(w_net)
        if i==1
            nw_net{i} = nw_net{i}./sum(nw_net{1});
        else
            nw_net{i} = nw_net{i}./ repmat(sum(nw_net{i},2),1,length(nw_net{i}(1,:))  )  ;
        end
    end
        
end
