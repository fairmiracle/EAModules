%The algebraic connectivity of a graph G is the second-smallest eigenvalue of the Laplacian matrix of G
%This eigenvalue is greater than 0 if and only if G is a connected graph. 
%Laplacian matrix, L = D-A
function w = algebraic(adj)
if length(adj)==1
    w = 1;
    return;
end
L = diag(sum(adj)) - adj;
lambda = eig(L);
lambda=sort(lambda,'ascend');
if lambda(2) > 0
    w = 1;          %connected
else
    w = 0;          %not connected
end


