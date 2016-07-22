
# this is the work flow for single cell RNAseq Transcriptome heterogeneity analysis using SCDE R package

# load required R library
library("scde")
library("DBI")

# set the working directory
setwd("./")

# load raw data and filter low expressed genes
count_tab <- read.table("data/sample_data.txt",head=T,row.names =1)
cd <- clean.counts(count_tab, min.lib.size = 1800, min.reads = 30, min.detected = 10)

# find cell group inforamtion from cell name
x <- gsub("(.*)\\.mouse.*", "\\1",colnames(cd))
two_groups <- c("athero", "healthy")[as.integer(factor(x, levels = c("athero", "healthy")))]

# construct error models for individual cells in two groups (athero and healthy)
knn <- knn.error.models(cd, groups = two_groups, k = 15, n.cores = 20, min.count.threshold = 3, min.nonfailed = 10, max.model.plots = 10, save.model.plots = TRUE, verbose = 1)  # turn on verbosity

# normalize out expected levels of technical and intrinsic biological noise
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 10, n.cores = 20, plot = T)
# Controlling for sequencing depth, how many gene is present in each cell
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

# save results
write.table(varinfo$mat,file="tab_varinfo_mat",quote=F,sep="\t")
write.table(varinfo$matw,file="tab_varinfo_matw",quote=F,sep="\t")
write.table(varinfo$arv,file="tab_varinfo_arv",quote=F,sep="\t")
write.table(varinfo$avmodes,file="tab_varinfo_avmodes",quote=F,sep="\t")
write.table(sort(varinfo$arv, decreasing = TRUE)[1:1000],file="tab_top_1000_overdis_genes",sep="\t",quote=F)

# load gene ontology data base
library(org.Mm.eg.db)
# translate gene names to Entrez id
ids <- unlist(lapply(mget(rownames(cd), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids

go.env <- eapply(org.Mm.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env <- go.env[unlist(lapply(go.env, length))>5]

# adding Go term description
library(GO.db)
desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
names(go.env) <- paste(names(go.env), desc)  # append description to the names
# convert to an environment
go.env <- clean.gos(go.env)
go.env <- list2env(go.env)

# calculate weighted first principal component magnitudes for each GO gene set in the provided environment
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 20, verbose = 1)

# save top overdispersed gene sets
df <- pagoda.top.aspects(pwpca,return.table = TRUE, plot = TRUE, z.score = 1.96)
write.table(df,file="tab_sig_aspects",sep="\t",quote=F)

# save top overdispersed genes
df <- pagoda.top.aspects(pwpca,return.genes = TRUE, plot = TRUE, z.score = 1.96)
write.table(df,file="tab_genes_driving_significant_aspects",sep="\t",quote=F)

# get significant aspects of heterogeneity
tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))

# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo, include.aspects = T)
pdf("plot_cluster.pdf",width = 9, height = 5)
plot(hc)
dev.off()

# get the details about cell clustering tree
hc_details <- pagoda.cluster.cells(tam, varinfo, return.details =1)
# hc_details$clustering$order
# hc_details$clustering$merge
# hc_details$distance
write.table(hc_details$clustering$order,file="tab_cell_order",sep="\t",quote=F)

# combine pathways that are driven by the same sets of genes:
pdf("plot_tamr.pdf",width = 9, height = 5)
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca,  plot=T)
dev.off()

# combine aspects that show similar patterns
pdf("plot_tamr2.pdf",width = 9, height = 5)
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.85, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
dev.off()

# translate group and sample source data into color codes.
x <- gsub("(.*)_Plate.*", "\\1",colnames(cd))
l2cols <- c("coral4", "olivedrab3", "skyblue2", "darkorchid4")[as.integer(factor(x, levels = c("athero.mouse1.large", "athero.mouse2.large", "athero.mouse2.medium", "healthy.mouse.medium")))]

# review top aspect
pdf("plot_top_aspect.pdf",width = 9, height = 5)
pagoda.view.aspects(tam, cell.clustering = hc, row.clustering = hclust(dist(tam$xv)), box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))
pagoda.view.aspects(tamr, cell.clustering = hc, row.clustering = hclust(dist(tamr$xv)), box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))
pagoda.view.aspects(tamr2, cell.clustering = hc, row.clustering = hclust(dist(tamr2$xv)), box = TRUE, labCol = NA, margins = c(0.5, 20),  col.cols = rbind(l2cols))
dev.off()

# specific genes heat map
pdf("important_genes.pdf",width =9, height =5)
genes = c("Abca1","Ly6a", "Pdgfrb","Smc1a")
d <-pagoda.show.pathways(genes, varinfo, go.env, cell.clustering = hc, two.sided =T, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE, n.genes=200,cexRow=1,return.details=1)
dev.off()

# specific Goterms
pdf("important_genes.pdf",width =9, height =5)
genes = c("Myh11", "Tagln", "Cnn1", "Myocd", "Acta2","Lgals3" , "Lamp2" , "Cd68", "Abca1","Ly6a", "Pdgfrb","Smc1a", "Smc1b","Smc2","Smc3","Smc4", "Smc5","Smc6")
d <-pagoda.show.pathways(genes, varinfo, go.env, cell.clustering = hc, two.sided =T, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE, n.genes=200,cexRow=1,return.details=1)
dev.off()


