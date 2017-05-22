%% Examples
%   Elastic net constrained linear programming for active modules identification in molecular interaction networks
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: 1*n vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   minmodulesize: minimal module size
%   maxmodulesize: maximal module size
%   iteration: maximal iteration for searching proper a
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
%  [3] Active module identification in multi-layer intracellular networks using constrained linear programming. Dong Li et al, in prearing
%  [4] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.

function [corrected_subnet_scorelp, fsubsetlp] = ENLP(G, z,randomscore, minmodulesize, maxmodulesize,iteration)
abegin = 0;
aend = 1;
n = size(G,1);

for iter=1:iteration
    a=(abegin+aend)/2
    cvx_begin
        variable x(n)
        maximize( z*x )
        subject to
        a*norm( x, 1 )+(1-a)*norm( x, 2 ) <= 1;
        x >= 0;
    cvx_end
    index_subnetlp=find(x > 1e-6);
    [corrected_subnet_scorelp, fsubsetlp] = topscore(G,z,randomscore,index_subnetlp);
    %length(fsubsetlp)
    if length(fsubsetlp) > maxmodulesize
        abegin = (abegin+aend)/2;
    elseif length(fsubsetlp) < minmodulesize
        aend = (abegin+aend)/2;
    else
        break;
    end
end