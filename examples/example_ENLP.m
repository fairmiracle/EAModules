%% Examples
%   genetic algorithm for active modules identification in molecular interaction networks
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
%   Last modified on May 22, 2017.
%
%% Related papers
%  [1] Discovering regulatory and signalling circuits in molecular interaction networks. Trey Ideker et al, Bioinformatics 2002
%  [2] Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. Dong Li et al, BMC Genomics 2017

clear all
addpath('../EAModules/utils/')
addpath('../EAModules/functions/')
addpath('../data/')
load('galdata.mat')
minmodulesize = 10;
maxmodulesize = 100;
for i = 1:10
[corrected_subnet_scorelp, fsubsetlp] = ENLP(G, array_basic_z,randomscore, minmodulesize, maxmodulesize,20);
fid = fopen('result.txt', 'a+');
fprintf(fid, '%f \n', corrected_subnet_scorelp);
fprintf(fid, '%s \n\n', num2str(fsubsetlp'));
fclose(fid);

G(fsubsetlp,:)=[];
G(:,fsubsetlp)=[];
array_basic_z(fsubsetlp)=[];
end
