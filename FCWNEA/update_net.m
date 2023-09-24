function [w_net, act_net] = update_net(w_net, Gra, ind_dec,theta)
%  
%   此处显示详细说明
    N_Gra = Gra;
    N_Gra(N_Gra~=0) =   N_Gra(N_Gra~=0)./abs(N_Gra(N_Gra~=0));    
    aa = sum(abs(N_Gra),2);
    
    ind_dec = ind_dec(aa~=0,:);
    N_Gra = N_Gra(aa~=0,:);
    
    ind_dec = ind_dec + N_Gra;
    ind_dec(ind_dec<1) = 1;
    ind_dec(ind_dec>51)=51;
    theta = 0.02;
    for i=1:length(ind_dec(:,1))
        gra = Gra(i,:);
        ind = ind_dec(i,:);
        for j =1 :length(ind)
           if  gra(j)~=0
               if j==1
%                    
                    m = max(w_net{j});
                    w_net{j} = w_net{j}./2;
                    w_net{j}(ind(j)) = m;
               else
%                    
                    m = max(w_net{j}(ind(j-1),:));
                    w_net{j}(ind(j-1),:) = w_net{j}(ind(j-1),:)./2;
                    w_net{j}(ind(j-1),ind(j) ) = m ;
               end
           end
        end
    end
end

