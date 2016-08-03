%% FUNCTION topscore
%   find the maximal module score for given node set, based on connected
%   component finding on current node set.
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   subnetwork of size i
%   nodeset: numberic vector, current node id set
%
%% OUTPUT
%   maximal module score
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
%     % put module size constraint
%     % if resulted module is too large then partition. But this is not best
%     if k > maxmodulesize
%         p = randperm(k);
%         for j = 1:ceil(k/maxmodulesize)
%             tmpnodeset = nodeList(p((j-1)*maxmodulesize+1:min(k,j*maxmodulesize)));
%             [tmpLnew, tmpCnew] = graph_conn_comp(G(tmpnodeset,tmpnodeset));
%             tmplabels = unique(tmpLnew);
%             
%             for insidek = 1:length(tmplabels)
%                 tmpnodeList = nodeList(find(tmpLnew==tmplabels(insidek)));
%                 tmpk = length(tmpnodeList);
%                 aggregate_score = sum(array_basic_z(tmpnodeList))/sqrt(tmpk);
%                 compscore = (aggregate_score - randomscore(tmpk,1))/randomscore(tmpk,2);
%                 if compscore > s
%                     s = compscore;
%                     subset = tmpnodeList;
%                 end
%             end
%         end
%         
%     else
%         aggregate_score = sum(array_basic_z(nodeList))/sqrt(k);
%         compscore = (aggregate_score - randomscore(k,1))/randomscore(k,2);
%         if compscore > s
%             s = compscore;
%             subset = nodeList;
%         end
%     end
%     
% end
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