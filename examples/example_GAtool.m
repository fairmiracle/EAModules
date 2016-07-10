%% Examples
%   genetic algorithm for active modules identification in molecular interaction networks,
% 	using Matlab GA toolbox
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   subnetwork of size i
%
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

clear all
addpath('../EAModules/utils/')
addpath('../EAModules/functions/')
addpath('../data/')
global G;
global array_basic_z;
global randomscore;
load('galFiltered.mat')
popsize = 1000;% import in ga toolbox
load('galdata.mat')
popsize = 100;
crossrate = 0.9;
mutrate = 0.5;
iteration = 500;
Nvar = length(array_basic_z);
lb = zeros(Nvar,1);
ub = ones(Nvar,1);
%[corrected_subnet_score, fsubset,func] = GA(G, array_basic_z, randomscore,popsize,crossrate,mutrate,iteration);
IntCon = 1:Nvar;
options = gaoptimset('Generations',1000,'PopulationSize',1000,'PopulationType','bitstring',...
    'CreationFcn',@gacreationuniform ,'CrossoverFcn',@crossoversinglepoint,...
    'MutationFcn',{@mutationuniform,0.5},'PlotFcns',@gaplotbestf,'TolFun',1e-19);
[x,fval,exitflag] = ga(@fitness,Nvar,[],[],[],[],[],[],[],[],options);
nodeSet = find(x==1);
[corrected_subnet_score, fsubset] = topscore(G,array_basic_z,randomscore,nodeSet);