function [y] = fitness( x )
global G;
global array_basic_z;
global randomscore;
s = -inf;
nodeSet = find(x==1);
[Lnew, Cnew] = graph_conn_comp(G(nodeSet,nodeSet));
labels = unique(Lnew);
for i = 1:length(labels)
    nodeList = nodeSet(find(Lnew==labels(i)));
    k=length(nodeList);
    aggregate_score =sum(array_basic_z(nodeList))/sqrt(k);
    compscore = (aggregate_score - randomscore(k,1))/randomscore(k,2);
    if compscore > s
        s = compscore;
        subset = nodeList;
    end
end
 y = -s;
end