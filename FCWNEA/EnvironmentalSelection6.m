function [Pop,Crowdis,Ind_Dec,FFrontNo,N_Gra,w_net,lif_net] = EnvironmentalSelection6(Population,N,Z,Zmin,Ind_dec,w_net,Gra,lif_net)
% The environmental selection of NSGA-III

%------------------------------- Copyright --------------------------------
% 增加了拥挤距离的计算
%--------------------------------------------------------------------------

    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end
    aa=lif_net;
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);
    Next(Last(Choose)) = true;
    % Population for next generation
%     Population = Population(Next);

    %% 
    Pop =[];
    N_Gra = [];
    FFrontNo =[];

    D_state = zeros(2,length(Ind_dec(1,:)));    %记录变量扰动对个体质量影响情况
    Ind_Dec = [];
    theta = 0.02;
    L = length(Ind_dec(1,:));
    normv = 5;   % 权重更新时的方差参数
    for i=1:N
        gra = Gra(i,:);
        m1 = Ind_dec(i,:); 
        m2 = Ind_dec(i+N,:);
        
        [j]= find(gra~=0);
        
        m_dom = sum(Population(i).obj>Population(i+N).obj);
        
        if Next(i) == 1 & Next(N+i) ==1   % 父个体、子个体都被选入下一代，对权重进行调整
            for k =1:length(j)
                if j(k)==1
%                     w_net{j(k)} = w_net{j(k)}./2+ normpdf([1:51],m1(j(k)),normv)';               
                    w_net{j(k)} = w_net{j(k)}./10+ normpdf([1:length(w_net{j(k)})],m1(j(k)),normv)';     
                else
%                     w_net{j(k)}(m1(j(k)-1),:) = w_net{j(k)}(m1(j(k)-1),:)./2+ normpdf([1:51],m1(j(k)),normv);          
                    w_net{j(k)}(m1(j(k)-1),:) = w_net{j(k)}(m1(j(k)-1),:)./10+ normpdf([1:length( w_net{j(k)}(m1(j(k)-1),:))],m1(j(k)),normv); 
%                      w_net{j(k)}(m1(j(k)-1),:) =  normpdf([1:length( w_net{j(k)}(m1(j(k)-1),:))],m1(j(k)),normv).^2;  
                end
                %更新变量活跃度
                if m_dom==3
                    D_state(1,j(k)) = D_state(1,j(k))+20;
                else
                    if lif_net(j(k)) >100
                        D_state(2,j(k)) = D_state(2,j(k))+1;
                    else 
                        D_state(1,j(k)) = D_state(1,j(k))+1;
                    end
                end
            end
                Pop = [Pop,Population(i),Population(N+i)];
                N_Gra = [N_Gra; zeros(1,L);gra];
                Ind_Dec = [Ind_Dec; Ind_dec(i,:);Ind_dec(i+N,:)];
                FFrontNo = [FFrontNo,FrontNo(i),FrontNo(N+i)];

        end      
        if Next(i) == 1 & Next(N+i) == 0    %父个体被选入、子个体没有
            for k=1:length(j)
                if j(k)==1
                   w_net{j(k)}(m1(j(k))) =  w_net{j(k)}(j(k))*0.5;
                else
                   w_net{j(k)}(m1(j(k)-1),m1(j(k))) = w_net{j(k)}(m1(j(k)-1),m1(j(k)))*1 ;               
                end
                Gra(i,j(k))=0;
                %更新变量活跃度
                if rand() <0.4
                    D_state(1,j(k)) = D_state(1,j(k))+4;
                else
                    if lif_net(j(k)) >1
                        D_state(2,j(k)) = D_state(2,j(k))+1;
                    end
                end

            end
                Pop = [Pop,Population(i)];
                N_Gra = [N_Gra;-gra];
                Ind_Dec = [Ind_Dec; Ind_dec(i,:)];
                FFrontNo = [FFrontNo,FrontNo(i)];
                
        end
        if Next(i) == 0 & Next(N+i) ==1     % 父个体未被选、 子个体被选
            for k=1:length(j)
                if j(k)==1
%                     w_net{(k)} = w_net{(k)}./2+ normpdf([1:51],m1(j(k)),normv)';
                    w_net{(k)} = w_net{(k)}./1+ normpdf([1:length(w_net{j(k)})],m1(j(k)),normv)';
                else
%                     w_net{j(k)}(m1(j(k)-1),:) = w_net{j(k)}(m1(j(k)-1),:)./2+ normpdf([1:51],m1(j(k)),normv);
                    w_net{j(k)}(m1(j(k)-1),:) = w_net{j(k)}(m1(j(k)-1),:)./10+ normpdf([1:length( w_net{j(k)}(m1(j(k)-1),:))],m1(j(k)),normv);
%                      w_net{j(k)}(m1(j(k)-1),:) =  normpdf([1:length( w_net{j(k)}(m1(j(k)-1),:))],m1(j(k)),normv).^2;

                end
                %更新变量活跃度
                if m_dom==3
                    D_state(1,j(k)) = D_state(1,j(k))+20;
                else
                    if lif_net(j(k)) >100
                        D_state(2,j(k)) = D_state(2,j(k))+1;
                    else 
                        D_state(1,j(k)) = D_state(1,j(k))+1;
                    end
                end
                
            end 
                Pop = [Pop,Population(i+N)];
                N_Gra =[N_Gra; gra];
                Ind_Dec = [Ind_Dec; Ind_dec(N+i,:)];
                FFrontNo = [FFrontNo,FrontNo(N+i)];
                
        end
        if Next(i) == 0 & Next(N+i) == 0   % 都未被选
%             for k=1:length(j)                 
%                 Gra(i,j(k))=0;
%                 %更新变量活跃度
%                 if rand() <0.4
%                     D_state(1,j(k)) = D_state(1,j(k))+1;
%                     lif_net(j(k))= lif_net(j(k))+50;
%                 else
%                     if lif_net(j(k)) >1
%                         D_state(2,j(k)) = D_state(2,j(k))+1;
%                         lif_net(j(k))= lif_net(j(k))-1;
%                     end
%                 end
% 
%             end
        end
    end
%     w_net = utility_net(w_net);
    [FrontNo,MaxFNo] = NDSort(Pop.objs,Pop.cons,N);
    Crowdis = CrowdingDistance(Population.objs,FrontNo);
    %更新变量活跃度
    a= exp(1)+D_state;
    a= 1+D_state;
    lif_net = aa+(a(1,:).^3)./a(2,:);
%     lif_net = log(a(1,:))./a(2,:);
    lif_net = lif_net./sum(lif_net);
    lif_net = aa+lif_net;
    
    
end


function [nw_net]=utility_net(w_net)
    nw_net = w_net;
    for i=1:length(w_net)
        if i==1
            nw_net{i} = nw_net{i}./sum(nw_net{1});
        else
            nw_net{i} = nw_net{i}./ repmat(sum(nw_net{i},2),1,length(nw_net{i}(1,: ))  )  ;
        end
    end
        
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end
