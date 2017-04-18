%% FUNCTION MA
%   memtic algorithm for active modules identification in molecular interaction networks
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   minmodulesize: minimal module size
%   maxmodulesize: maximal module size
%   subnetwork of size i
%   popsize: population size
%   crossrate: cross rate of population
%   lsrate: local search rate of population
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

function [corrected_subnet_score, fsubset,func] = MA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,popsize,crossrate,lsrate,lsiter,iteration)

if nargin < 8
    error('\n Inputs: G, array_basic_z, randomscore,popsize,crossrate,lsrate,lsiter,iteration should be specified!\n');
end
N = length(G);
Pop = cell(popsize,1);
for i = 1:popsize
    v = randperm(N);
    nodeSet = v(1:maxmodulesize);
    nodelist = zeros(1,N);
    nodelist(nodeSet) = 1;
    indiv.nodes = nodelist;
    [s1,topset] = topscore(G,array_basic_z,randomscore,nodeSet);
    indiv.score = s1;
    Pop{i} = indiv;
end
tic;
times=zeros(iteration,1);
func = zeros(iteration,1);
for T = 1:iteration
    % selection
    allscore = zeros(popsize,1);
    scoresum = 0;
    for i = 1:popsize
        allscore(i) = Pop{i}.score;
        scoresum = scoresum+Pop{i}.score;
    end
    [B,I] = sort(allscore,'descend');
    bestindiv = Pop{I(1)};
    B = B/scoresum;
    
    % accumulated score
    accumscore = zeros(popsize,1);
    accumscore(1) = B(1);
    for i = 2:popsize
        accumscore(i) = B(i)+accumscore(i-1);
    end
    
    Popnew = cell(popsize,1);
    for i = 1:popsize
        idx = locatep(accumscore,rand());
        Popnew{i} = Pop{idx};
    end
    
    %cross
    for i = 1:ceil(popsize*crossrate)
        p1 = ceil(popsize*rand());
        p2 = ceil(popsize*rand());
        nodes1 = Popnew{p1}.nodes;
        nodes2 = Popnew{p2}.nodes;
        crosspoint =  ceil(N*rand());
        for j = 1:crosspoint
            tmp = nodes1(j);
            nodes1(j) = nodes2(j);
            nodes2(j) = tmp;
        end
        Popnew{p1}.nodes = nodes1;
        Popnew{p2}.nodes = nodes2;
    end
    
    % local search
    for i = 1:ceil(popsize*lsrate)
        %p1 = ceil(popsize*rand());
        nodes1 = Popnew{p1}.nodes;
        %nodes1 = bestindiv.nodes;
        
        bestnodes1 = nodes1;
        nodeSet = find(bestnodes1==1);
        [bstscore,bstset] = topscore(G,array_basic_z,randomscore,nodeSet);
        for j = 1:lsiter
            mutatepoint =  ceil(N*rand());
            if nodes1(mutatepoint) == 1
                if (sum(nodes1) > minmodulesize)
                    nodes1(mutatepoint) = 0;
                end
            else
                if (sum(nodes1) < maxmodulesize)
                    nodes1(mutatepoint) = 1;
                end
            end
            nodeSet = find(nodes1==1);
            [topkscore,topkset] = topscore(G,array_basic_z,randomscore,nodeSet);
            if bstscore < topkscore
                bstscore = topkscore;
                bestnodes1 = nodes1;
            end
        end
        Popnew{p1}.nodes =  bestnodes1;
    end
    
    Pop = Popnew;
    Pop{1} = bestindiv;
    
    %compute fitness score
    for i = 1:popsize
        nodeSet = find(Pop{i}.nodes==1);
        [s1,topset] = topscore(G,array_basic_z,randomscore,nodeSet);
        Pop{i}.score = s1;
    end
    times(T) = toc;
    func(T) = Pop{1}.score;
    disp(['completed ' num2str(T/iteration)]);
end

for i = 1:popsize
    allscore(i) = Pop{i}.score;
    scoresum = scoresum+Pop{i}.score;
end
totaltime = toc;

[B,I] = sort(allscore,'descend');
nodeSet = find(Pop{I(1)}.nodes==1);
[corrected_subnet_score, fsubset] = topscore(G,array_basic_z,randomscore,nodeSet);