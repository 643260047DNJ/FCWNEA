%% 初始化概率网   
%del:变量切分粒度； N_x 决策变量个数； L 变量下届； U 变量上届。
% del=100;
% N_x=10;
% [,bb]=initial_net(del,N_x,zeros(N_x,1),ones(N_x,1))
% clc;

function [val_net,w_net,act_net,lif_net]=initial_net(Problem)
    del = 150;
    N_x = Problem.D;
    L = Problem.lower;
    U = Problem.upper;
    var_id = Problem.var_id;
%     o_val = Problem.o_val;
    val_net = {};
    mm=1;
    for i=1:N_x
        if var_id(i)==1
            val_net = [val_net;[L(i): (U(i)-L(i))/del: U(i)] ];
        else 
%           
            val_net = [val_net;[L(i): U(i)] ];  
            mm=mm+1;
        end
    end
    
    w_net= {};
    act_net ={};
    w = 1/length(val_net{1}).*ones(length(val_net{1}),1);
    act = ones(length(val_net{1}),1);
    w_net = [w_net;w];
    act_net = [act_net;act];
    for i=1:N_x-1
        n=length(val_net{i});
        m=length(val_net{i+1});
        w_net = [w_net; 1/m .* ones(n,m) ];
        act = ones(n,m);
        act_net = [act_net;act];
    end
    lif_net =ones(1,N_x);
    
end