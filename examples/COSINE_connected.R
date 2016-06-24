library(COSINE)
data(scaled_node_score)
data(scaled_edge_score)
data(PPI)
GA_result2<-GA_search_PPI(lambda=0.5,scaled_node_score,scaled_edge_score,PPI,
num_iter=100, muCh=0.05, zToR=10, minsize=50)

GA_result <- GA_search_connected(lambda=0.5,scaled_node_score,scaled_edge_score,PPI,
num_iter=100, muCh=0.05, zToR=10, minsize=50)

## visualized by igraph
selected = c()
for (i in 1:dim(PPI)[1]){
    if (PPI[i,1] %in% GA_result2$Subnet || PPI[i,2] %in% GA_result2$Subnet)
        selected = c(selected,i)
}
library(igraph)
g <- graph.data.frame(PPI[selected,], directed=FALSE)
V(g)$color <- "blue"
V(g)$color[match(GA_result2$Subnet,V(g)$name)] <- "red"
layout <- layout.reingold.tilford(g, circular=T)
plot(g,layout=layout, vertex.size=5, vertex.label.cex=5,vertex.label.dist=0.5,edge.width=5)

# visualized by qgraph
library(qgraph)
x <- PPI[selected,]
nodes <- unique(union(x[,1],x[,2]))
adj <- matrix(0,nrow = length(nodes), ncol = length(nodes))
rownames(adj) <- as.character(nodes)
colnames(adj) <- as.character(nodes)
