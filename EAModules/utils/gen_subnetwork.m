function index_subnet = gen_subnetwork(G,k);
%input G = adjacent matrix representing network;
%k = size of subnetwork ( number of nodes in subnetwork)
%G = [1,0,0,0,0,1,1,0;0,0,0,0,0,1,1,0;0,0,0,0,0,1,1,0;0,0,0,0,0,1,1,0;0,0,0,0,0,1,1,0;1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1;0,0,0,0,0,1,1,0;]

%G = [0,1,1,0,0,0; 1,0,0,0,0,1;1,0,0,1,1,0;0,0,1,0,0,0;0,0,1,0,0,0;0,1,0,0,0,0]
%k = 3;
temp  = size(G);
number_of_node = temp(1);
v = randperm(number_of_node);  %total
v = v(1:k);          % node of subgraph
index_subnet = sort(v);
%    index_subnet = [1,2,3]
%{
    G_w = zeros(number_of_node);
    for i = 1:length(index_subnet)
        for j = 1: length(index_subnet)
           G_w(index_subnet(i),index_subnet(j)) = G(index_subnet(i),index_subnet(j)); 
           G_w(index_subnet(j),index_subnet(i)) = G(index_subnet(j),index_subnet(i));
        end
    end
%}

