%% FUNCTION TOPSCOREMULTILAYER
%   calculated the highest-scoring connected component with given nodes
%   set, in a multilayer network
%
%% INPUT
%   TG: a k*n*n tensor, each slice is a n*n adjacency matrix
%   nodesocre_z: k*n matrix, each column is a vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   nodeset: given nodes set

%% OUTPUT
%   module score and component nodes set
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
%   Copyright (C) 2017 Dong Li
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via donggeat@gmail.com
%
%   Last modified on April 18, 2017.
%
%% Related papers
%  [1] Discovering regulatory and signalling circuits in molecular interaction networks. Trey Ideker et al, Bioinformatics 2002
%  [2] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.

function [s,subset] = topscore(TG,nodesocre_z,randomscore,nodeset)

if nargin < 4
    error('\n Inputs: G, array_basic_z, randomscore, nodeset should be specified!\n');
end
k = size(TG,1);
G = TG(k,:,:);
array_basic_z = nodesocre_z(k,:);

s = -inf;
[Lnew, Cnew] = graph_conn_comp(G(nodeset,nodeset));
labels = unique(Lnew);

for i = 1:length(labels)
    nodeList = nodeset(find(Lnew==labels(i)));
    k=length(nodeList);
    compscore = sum(array_basic_z(nodeList))/sqrt(k);
    %aggregate_score = sum(array_basic_z(nodeList))/sqrt(k);
    %compscore = (aggregate_score - randomscore(k,1))/randomscore(k,2);
    if compscore > s
        s = compscore;
        subset = nodeList;
    end
    
end