classdef A4_FCWNEDA < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate the reference points and random population
            
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            [val_net, w_net,act_net,lif_net ] = initial_net(Problem);
           [Population, Ind_dec] = initial_pop(val_net,w_net,Problem.N);
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            record_lif = [];
            %% Optimization
            kk=1;
            NN = [];
            while Algorithm.NotTerminated(Population)
                kk=kk+1;
                 theta =0.08; 
                 phi = Problem.FE/Problem.maxFE;
                 
              
                [Population,Ind_dec,Offspring,Ind_decoff,Gra,act_net]  = generate_off2(Population,Ind_dec,val_net,w_net,act_net,lif_net,Problem,phi);
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                [Population,Crowdis,Ind_dec,FrontNo,N_Gra,w_net,lif_net] = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin,[Ind_dec;Ind_decoff],w_net,Gra,lif_net,phi);
               
               
%                [w_net] = update_net(w_net,N_Gra, Ind_dec,theta);
                [w_net] = update_net2(w_net,Population,Ind_dec,phi);
                
                if Problem.FE/Problem.maxFE>0.6 & mod(kk,10)==0
                    [val_net,w_net] = val_update2(Ind_dec,Population,Problem,act_net,val_net,w_net);
                    NN = [NN;val_net{2}];
                end
                record_lif = [record_lif;lif_net];
                
            end
            a=0;
        end
    end
end