function [mu,sigma] = monte(G,k,array_basic_z,iter_monte);       
% monte carlo approach 
%input : G = graph, k = number of node in subgraph, array of p value,
%            array_basic_z(all genes are already computed basic z score), 
%            iteration of monte
%output: mu, sigma

%method: for i = 1: iter_monte 
%           generate subgraph with k nodes
%           compute and keep z_a score (p_value-> z score -> aggregate score z_a)
%        end
% compute mu and sigma from all z_a score


z_monte = zeros(1,iter_monte);          %temp array for keep value
aggregate_score_monte = zeros(1,iter_monte);  
for i = 1:iter_monte
    %generate subnetwork size k
    index_subnet_monte = gen_subnetwork(G,k);
    aggregate_score_monte(i) = aggregate(array_basic_z,index_subnet_monte);
end
mu = mean(aggregate_score_monte);
sigma = std(aggregate_score_monte);

%s = (z - mu)/sigma;