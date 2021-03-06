## Example to illustrate the differences between module identified by COSINE 
## and the connected ensured method
library(COSINE)
data(scaled_node_score)
data(scaled_edge_score)
data(PPI)

## GA from COSINE
GA_result<-GA_search_PPI(lambda=0.5,scaled_node_score,scaled_edge_score,PPI,
num_iter=100, muCh=0.05, zToR=10, minsize=50)

## GA that ensures connectedness of resulted module
GA_result2 <- GA_search_connected(lambda=0.5,scaled_node_score,scaled_edge_score,PPI,
num_iter=100, muCh=0.05, zToR=10, minsize=50)

## visualized by igraph
visualizedModule <- function(GA_result2,PPI,savefilename){
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
    #layout <- layout.fruchterman.reingold(g)
    png(file = savefilename,width = 4096, height = 3652)
    plot(g,layout=layout, vertex.size=5, vertex.label.cex=5,vertex.label.dist=0.5,edge.width=5)
    dev.off()
    plot(g,layout=layout, vertex.size=5, vertex.label.dist=0.5)
}

visualizedModule(GA_result2,PPI,'ppi.png')

## venn digram
grid.newpage()
draw.triple.venn(area1 = 34, area2 = 55, area3 = 55, n12 = length(intersect(saset,gaset)), n23 = length(intersect(maset,gaset)), n13 = length(intersect(saset,maset)), 
    n123 = length(intersect(maset,intersect(saset,gaset))), category = c("SA", "GA", "MA"), lty = "blank", 
    fill = c("skyblue", "pink1", "mediumorchid"))