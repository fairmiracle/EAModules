%% Examples
%   genetic algorithm for active modules identification in molecular interaction networks
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
load('galdata.mat')
popsize = 100;
crossrate = 0.9;
lsrate = 0.5;
iteration = 10000;
lsiter = 10;
minmodulesize = 10;
maxmodulesize = 100;

%[corrected_subnet_score_MA, fsubset_MA,func_MA] = MA(G, array_basic_z, randomscore,popsize,crossrate,lsrate,lsiter,iteration);

for i = 1:10
[corrected_subnet_score_MA, fsubset_MA,func_MA] = MA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,popsize,crossrate,lsrate,lsiter,iteration);

fid = fopen('result.txt', 'a+');
fprintf(fid, '%f \n', corrected_subnet_score_MA);
fprintf(fid, '%s \n\n', num2str(fsubset_MA));
fclose(fid);

disp(GeneNames(fsubset_MA));


G(fsubset_MA,:)=[];
G(:,fsubset_MA)=[];
array_basic_z(fsubset_MA)=[];
randomscore(fsubset_MA,:)=[];
GeneNames(fsubset_MA) = [];
end
