%% FUNCTION TOPSCORE
%   calculated the highest-scoring connected component with given nodes set
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
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
%   Copyright (C) 2015 - 2016 Dong Li
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via donggeat@gmail.com
%
%   Last modified on June 1, 2016.
%
%% Related papers
%  [1] Discovering regulatory and signalling circuits in molecular interaction networks. Trey Ideker et al, Bioinformatics 2002
%  [2] Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. Dong Li et al, BMC Genomics 2017

function [s,subset] = topscore(G,array_basic_z,randomscore,nodeset)

if nargin < 4
    error('\n Inputs: G, array_basic_z, randomscore, nodeset should be specified!\n');
end

s = -inf;
[Lnew, Cnew] = graph_conn_comp(G(nodeset,nodeset));
labels = unique(Lnew);

% for i = 1:length(labels)
%     nodeList = nodeset(find(Lnew==labels(i)));
%     k=length(nodeList);
%     compscore = sum(array_basic_z(nodeList))/sqrt(k);
%     %aggregate_score = sum(array_basic_z(nodeList))/sqrt(k);
%     %compscore = (aggregate_score - randomscore(k,1))/randomscore(k,2);
%     if compscore > s
%         s = compscore;
%         subset = nodeList;
%     end
%     
% end

for i = 1:length(labels)
    nodeList = nodeset(find(Lnew==labels(i)));
    k=length(nodeList);
    if k > s
        s = k;
        subset = nodeList;
    end
    
end