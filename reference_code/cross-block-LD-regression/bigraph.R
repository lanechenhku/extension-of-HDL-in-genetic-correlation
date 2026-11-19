library(igraph)

#### create example input data
S1=c("SEC16B", "WLS", "WLS", "WLS", "TGFB2", "NUAK2", "TGFB2", "NUAK2", "C1orf141", "WLS", "METTL13", "JMJD4", "CD34", "GCSAML", "ZNF678", "CTPS1", "TGFB2", "IL23R", "COLGALT2", "ZNF678", "JMJD4", "SLFNL1", "COLGALT2", "COLGALT2", "WLS", "TGFB2", "JMJD4", "NFASC", "SCMH1", "TGFB2", "ZNF678", "ZNF678", "SNAP47", "PTGFR", "SNAP47", "DCST2", "PCSK9", "ZBTB7B", "TGFB2", "COLGALT2", "CTPS1", "TGFB2", "SCMH1", "DSTYK", "TMCC2", "NUAK2", "C1orf127")
nodesSet1 <- unique(S1)
S2=c("TMEM18", "PKDCC", "SPTBN1", "EML6", "EFEMP1", "ANKMY1", "DIS3L2", "RNPEPL1", "ATG16L1", "EN1", "EHD3", "DIS3L2", "CCDC141", "EHD3", "EFEMP1", "ZNF638", "NPPC", "ATG16L1", "LTBP1", "DIS3L2", "LTBP1", "ZNF638", "EFEMP1", "DIS3L2", "C2orf73", "LTBP1", "EFEMP1", "ANKMY1", "ZNF638", "BOK", "NPPC", "LTBP1", "DIS3L2", "TMEM18", "LTBP1", "CYBRD1", "APOB", "CYBRD1", "ZNF638", "NPPC", "BOK", "EIF2AK3", "BOK", "RNPEPL1", "RNPEPL1", "GPC1", "WNT10A")
nodesSet2 <- unique(S2)
edgeList <- data.frame(S1=S1,
                       S2=S2)
####

### PREMISE :
### graph.bipartite function simply create a graph and add a boolean 'type' attribute
### to each vertex that indicate the partition of the vertex (e.g. TRUE first partition, 
### FALSE second partition).
### So, it's not strictly necessary to use that function to get a bipartite graph, but 
### you can use any method you like (or feel easier) as long as you add a 'type' 
### attribute. 
### Hence, in the following code I won't use graph.bipartite since I don't like it :)


# first we give prefixes to the nodes to discern the two partition
g <- graph.empty()
g <- add.vertices(g,nv=length(nodesSet1),attr=list(name=paste0('',nodesSet1),
                                                   type=rep(TRUE,length(nodesSet1))))
g <- add.vertices(g,nv=length(nodesSet2),attr=list(name=paste0('',nodesSet2),
                                                   type=rep(FALSE,length(nodesSet2))))

# we need to turn edgeList into a vector (and using names instead of indexes)
edgeListVec <- as.vector(t(as.matrix(data.frame(S1=paste0('',edgeList$S1),
                                                S2=paste0('',edgeList$S2)))))
g <- add.edges(g,edgeListVec)

# check if is recognized as bipartite
is.bipartite(g)

bipartite.mapping(g)
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
E(g)$color <- "lightgray"
par(mar=c(0,0,0,0))
plot(g, vertex.label.cex = 0.7,vertex.size=12,edge.color='black',edge.width=2,edge.arrow.size=0.2, vertex.label.color = "black",margin=-c(-0.1,-0.1,-0.1,-0.1))

# # let's plot it !
# png('bigraph.png',width = 5000,height=600)
# #par(mar=c(0,0,0,0))
# plot.igraph(g, layout=layout.bipartite,mark.col = 'black',vertex.label.cex=1,vertex.size=20,edge.color='black',edge.width=2,edge.arrow.size=0.3,
#             margin=c(-0.1,-0.1,-0.1,-0.1), 
#             vertex.color=c("orange","green")[V(g)$type+1])
# dev.off()

#### create example input data
S1=c("SEC16B", "WLS", "WLS", "WLS", "TGFB2", "NUAK2", "TGFB2", "NUAK2", "C1orf141", "WLS", "METTL13", "JMJD4", "CD34", "GCSAML", "ZNF678", "CTPS1", "TGFB2", "IL23R", "COLGALT2", "ZNF678", "JMJD4", "SLFNL1", "COLGALT2", "COLGALT2", "WLS", "TGFB2", "JMJD4", "NFASC", "SCMH1", "TGFB2", "ZNF678", "ZNF678", "SNAP47", "PTGFR", "SNAP47", "DCST2", "PCSK9", "ZBTB7B", "TGFB2", "COLGALT2", "CTPS1", "TGFB2", "SCMH1", "DSTYK", "TMCC2", "NUAK2", "C1orf127")
nodesSet1 <- unique(S1)
S2=c("TMEM18", "PKDCC", "SPTBN1", "EML6", "EFEMP1", "ANKMY1", "DIS3L2", "RNPEPL1", "ATG16L1", "EN1", "EHD3", "DIS3L2", "CCDC141", "EHD3", "EFEMP1", "ZNF638", "NPPC", "ATG16L1", "LTBP1", "DIS3L2", "LTBP1", "ZNF638", "EFEMP1", "DIS3L2", "C2orf73", "LTBP1", "EFEMP1", "ANKMY1", "ZNF638", "BOK", "NPPC", "LTBP1", "DIS3L2", "TMEM18", "LTBP1", "CYBRD1", "APOB", "CYBRD1", "ZNF638", "NPPC", "BOK", "EIF2AK3", "BOK", "RNPEPL1", "RNPEPL1", "GPC1", "WNT10A")
nodesSet2 <- unique(S2)
SS=data.frame(
  stringsAsFactors = FALSE,
                S1 = c("ANKMY1","APOB","APOB",
                          "ATG16L1","C1orf141","C2orf73","COLGALT2",
                          "DIS3L2","DSTYK","EFEMP1","EIF2AK3","GPC1","JMJD4",
                          "LTBP1","NFASC","SEC16B","WLS"),
                S2 = c("RNPEPL1","GPC1",
                          "PCSK9","IL23R","IL23R","SPTBN1","METTL13","EIF2AK3",
                          "TMCC2","TGFB2","PCSK9","WNT10A","SNAP47","TGFB2",
                          "SPTBN1","TMEM18","WNT10A")
   )
edgeList <- SS
  
####

### PREMISE :
### graph.bipartite function simply create a graph and add a boolean 'type' attribute
### to each vertex that indicate the partition of the vertex (e.g. TRUE first partition, 
### FALSE second partition).
### So, it's not strictly necessary to use that function to get a bipartite graph, but 
### you can use any method you like (or feel easier) as long as you add a 'type' 
### attribute. 
### Hence, in the following code I won't use graph.bipartite since I don't like it :)


# first we give prefixes to the nodes to discern the two partition
g <- graph.empty()
g <- add.vertices(g,nv=length(nodesSet1),attr=list(name=paste0('',nodesSet1),
                                                   type=rep(TRUE,length(nodesSet1))))
g <- add.vertices(g,nv=length(nodesSet2),attr=list(name=paste0('',nodesSet2),
                                                   type=rep(FALSE,length(nodesSet2))))

# we need to turn edgeList into a vector (and using names instead of indexes)
edgeListVec <- as.vector(t(as.matrix(data.frame(S1=paste0('',edgeList$S1),
                                                S2=paste0('',edgeList$S2)))))
g <- add.edges(g,edgeListVec)

# check if is recognized as bipartite
is.bipartite(g)

bipartite.mapping(g)
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
E(g)$color <- "lightgray"
par(mar=c(0,0,0,0))
plot(g, vertex.label.cex = 0.7,vertex.size=12,edge.color='black',edge.width=2,edge.arrow.size=0.2, vertex.label.color = "black",margin=-c(-0.1,-0.1,-0.1,-0.1))

