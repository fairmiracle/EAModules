%% FUNCTION TOPSCOREMULTILAYER
%   calculated the highest-scoring connected component with given nodes
%   set, in a multilayer network. The current strategy is to find the 
%   highest-scoring connected component of each layer, and find the maximal
%   consensus one as the final result
%
%% INPUT
%   TG: a k*n*n tensor, each slice is a n*n adjacency matrix G
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
%  [2] Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. Dong Li et al, BMC Genomics 2017
%  [3] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.

function [s,subset] = topscore(TG,nodesocre_z,randomscore,nodeset)

if nargin < 4
    error('\n Inputs: G, array_basic_z, randomscore, nodeset should be specified!\n');
end
k = size(TG,1);
C = cell(k,2)
for layers = 1:k
    G = TG(k,:,:);
    array_basic_z = nodesocre_z(k,:);
    s = -inf;
    [Lnew, Cnew] = conncomp(G(nodeset,nodeset));
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
    C{k,1}=subset;
    C{k,2}=s;
end

% give a set of component, find the maximal consensus one

%http://stackoverflow.com/questions/16883367/how-to-find-connected-components-in-matlab
function [S,C] = conncomp(G)
  [p,q,r] = dmperm(G'+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end

%http://stackoverflow.com/questions/21078445/find-connected-components-in-a-graph
global vistable;
function dfs(v,vistable,G)
    idxu = G(v,:);
    for i=1:length(idxu)
        if vistable(idxu(i))==0
            vistable(idxu(i))=1
            dfs(idxu(i),vistable,G)
        end
    end
end

function [S,C] = dfsconncomp(G)
    S = 0;
    C = zeros(n,1);
    n = size(G,1);
    vistable = zeros(n,1);
    for i = 1:n
        if vistable(i)==0
            vistable(i)=1;
            S = S+1;
            %C(i) = S;
            dfs(i,vistable,G)
        end
    end
end