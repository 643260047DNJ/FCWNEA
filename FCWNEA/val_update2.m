function [val_net,w_net] = val_update(Ind_dec,Population,Problem,act_net,val_net,w_net)
%UPDATE_NET2  对连续变量的值域取值进行更新
%   此处显示详细说明
    ID = Problem.var_id;
    mm=find(ID(1:Problem.D)==1);
    del = length(val_net{1});
    mark=0;
    for i =1:length(mm)
        ind = Ind_dec(:,mm(i));
        vnet = val_net{mm(i)};
        new_vnet =[];

        
        ind_state = tabulate(ind);
        tl = length(ind_state(:,1));
        if tl~=del
            ind_state = [ind_state;zeros(del-tl,3)];
        end
        rate = ind_state(:,3);  % 当前种群中变量的各个值域的访问频率
     
        rr=[];
        kk=1;
        % 找出没有被访问过的值域点， rr(i) 表示该节点为连续没被访问过的数量
        for j=1:length(rate)
           if rate(j)==0
              rr(j)=kk;
              kk=kk+1
           else
              kk=1;
              rr(j) = 0;
           end
        end 
        
        % 标记出连续未被访问节点数量超过一定量 的片段区域
        temp_r = zeros(1,length(rr));   
        for j=1:length(rr)-1
            if rr(j)~=1 & rr(j+1)==0
                if rr(j)>11
                    temp_r(j-rr(j)+1:j)=1;
                end 
            end     
        end
        if rr(end)>10
            temp_r(j-rr(end)+1:del)=1;
        end   
        R = rr.*temp_r;
        
        % 重新标记， 标记连续访问的片段并计数， 未被访问节点标记未0
        kk=1;
        for j=1:length(R)
            if R(j)==0
              R(j)=kk;
              kk=kk+1;
           else
              kk=1;
              R(j) = 0;
           end
        end
        
%         if length(R)==1
%             new_vnet = repmat(vnet(R),1,51);
%         end
        
        for j =1:length(R)-1
            if (R(j)~=0 & R(j+1)==0 ) || ( R(j+1)~=0 && j==length(R)-1)
                    mark =1;
                    if R(j+1)~=0
                        nn=round(del*R(j+1)./sum(R~=0));
                    else
                        nn=round(del*R(j)./sum(R~=0));
                    end
                    
                    if j-R(j)>1 
                        v1 = (vnet(j-R(j))+vnet(j-R(j)+1))/2; 
                    else
                        v1 = vnet(1);
                    end
                    if j+2>length(R)
                        v2 = vnet(end);
                    else 
                        v2 = (vnet(j)+vnet(j+1))/2;
                    end
                    %更新可取值域
                    if v1==v2
                        new_vnet = [v1,v2];
                    else
                        new_vnet = [new_vnet,v1:(v2-v1)/(nn-1):v2];
                    end

                    %更新响应的权重                  
            end
        end
        
         if length(new_vnet) ~=del
                        if length(new_vnet)<del
                            while length(new_vnet)~=del
                                mark1 = randperm(length(new_vnet)-1,1);
                                new_vnet = [new_vnet(1:mark1),(new_vnet(mark1)+new_vnet(mark1+1))/2,new_vnet(mark1+1:end)];
                            end
                        else
                            while length(new_vnet)~=del
                                mark1 = 1+randperm(length(new_vnet)-2,1);
                                new_vnet = [new_vnet(1:mark1-1),new_vnet(mark1+1:end) ];
                            end
                        end
                    end
        
        if mark ==1 
            val_net{i}= new_vnet;
            w_net{i}=0.02.* ones(size(w_net{i}));
            mark=0;
        else
           w_ent{i}=w_net{i};
        end
    end

end
