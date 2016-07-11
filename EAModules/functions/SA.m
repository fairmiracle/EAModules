%% FUNCTION SA
%   simulated annealing for active modules identification in molecular interaction networks
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   modulesize: maximal module size
%   subnetwork of size i
%   T: starting temperature
%   T_end: ending temperature
%   iteration: maximal iteration
%% OUTPUT
%   Module membership
%
%% LICENSE
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2015 - 2016 Dong Li and Shan He
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via dxl466@cs.bham.ac.uk
%
%   Last modified on June 1, 2016.
%
%% Related papers
%  [1] Discovering regulatory and signalling circuits in molecular interaction networks. Trey Ideker et al, Bioinformatics 2002
%  [2] Memetic algorithm for finding active connected subnetworks in intracellular networks. Dong Li et al, 2016

function [corrected_subnet_score, fsubset,func] = SA(G, array_basic_z, randomscore,modulesize,T,T_end,iteration)

if nargin < 6
    error('\n Inputs: G, array_basic_z, randomscore,T,T_end,iteration should be specified!\n');
end

prob_active_inactive = 0.1;

N = length(G);
v = randperm(N);
index_subnet = v(1:modulesize);

temp_step = 1 - (T_end/T)^(1.0/iteration);

[corrected_subnet_score,fsubset] = topscore(G,array_basic_z,randomscore,index_subnet);
k = length(index_subnet);
func = zeros(iteration,1);
for i = 1:iteration
    T = T*(1 - temp_step);
    new_node =ceil(rand()*N);
    new_index_subnet = index_subnet;
    if ismember(new_node,index_subnet)
        new_index_subnet(new_index_subnet==new_node)=[];
        new_k = k-1;
    else
        if (k < modulesize)
            new_index_subnet = [new_index_subnet, new_node];
            new_k = k+1;
        end
    end
    [new_corrected_subnet_score,fsubset] = topscore(G,array_basic_z,randomscore,new_index_subnet);
    
    if (new_corrected_subnet_score > corrected_subnet_score)
        corrected_subnet_score = new_corrected_subnet_score;
        index_subnet = new_index_subnet;
        k = new_k;
    else
        prob_sa = exp((new_corrected_subnet_score - corrected_subnet_score )/T);
        ran_num =rand(1);
        if (ran_num < prob_sa)
            corrected_subnet_score = new_corrected_subnet_score;
            index_subnet = new_index_subnet;
            k = new_k;
        end
    end
    func(i) = corrected_subnet_score;
    disp(['completed ' num2str(i/iteration)]);
end
[corrected_subnet_score, fsubset] = topscore(G,array_basic_z,randomscore,index_subnet);