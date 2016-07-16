function aggregate_score = aggregate (array_basic_z, index_subnet);  
% input : array of the z score values in subnetwork A
%output  : one value of score that represents subnetwork A 
k = length(index_subnet);
aggregate_score =sum(array_basic_z(index_subnet))/sqrt(k);  %sum elements in array z only index in subnet