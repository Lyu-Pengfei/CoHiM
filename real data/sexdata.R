library(tidyverse)
library(STAREG)
library(scales)
library(patchwork)
library(adaFilter)
library("snpGeneSets")
library(qqman)

save_file = FALSE

data.a <- read.table("./1. Data/DIAGRAM.Morris2012.males.txt", sep = "", header = TRUE)
data.b <- read.table("./1. Data/DIAGRAM.Morris2012.females.txt", sep = "", header = TRUE) # q = 1e-5

# data.a <- read.table("./1. Data/DIAGRAM.website.GWAS.metabochip.txt", sep = "", header = TRUE)
# data.b <- read.table("./1. Data/DIAGRAMv3.2012DEC17.txt", sep = "", header = TRUE)


data.example <- read.table("./1. Data/example.txt", sep = "", header = TRUE)

# data.a <- data.a[data.a$Freq1 >= 0.05, ]
# data.b <- data.b[data.b$Freq1 >= 0.05, ]

rownames(data.a) <- data.a[, "ID"]
rownames(data.b) <- data.b[, "ID"]
snpID <- intersect(data.a$ID, data.b$ID)


pa <- data.a[snpID, "P_VALUE"]
pb <- data.b[snpID, "P_VALUE"]

source("./Yan's codes/Simulations/RepLIS_0321.R")
q = 1e-5
J = length(snpID)

# EB_HMM
tic <- Sys.time()
res.hmm <- RepLIS(pa, pb)
toc <- Sys.time()
toc - tic
padj.hmm <- res.hmm$fdr
snp.hmm <- snpID[which(padj.hmm <= q)]
rLIS.hmm <- res.hmm$repLIS
thres.hmm = max(rLIS.hmm[padj.hmm <= q])



# EB (STAREG)
# source("./R/RepLfdr.R")
# res.eb <- RepLfdr(pa, pb)
tic <- Sys.time()
res.eb <- stareg(pa, pb)
toc <- Sys.time()
toc - tic
padj.eb <- res.eb$fdr
snp.stareg <- snpID[which(padj.eb<=q)]
LIS.eb = res.eb$Lfdr
thres.eb = max(LIS.eb[padj.eb <= q])

# MaxP
maxp <- pmax(pa, pb)
padj.maxp <- p.adjust(maxp, method = "BH")
snp.maxp <- snpID[which(padj.maxp<=q)]
thres.maxp = max(maxp[padj.maxp <= q])

# ad hoc BH
padj1.bh <- p.adjust(pa, method = "BH")
padj2.bh <- p.adjust(pb, method = "BH")
snp.bh <- snpID[which(padj1.bh<=q & padj2.bh<=q)]

# MaRR
x <- rank(pa)
y <- rank(pb)
max.rank = apply(cbind(x,y),1,max)
res.marr <- MaRR(max.rank, alpha = q)
marr.rej <- rep(0, J)
# marr.rej[res.marr$which.sig] = 1
snp.marr = snpID[res.marr$which.sig]

# Bogomolov & Heller 2018 (Biometrika)
library(radjust)
res.rv18 <- radjust_sym(pa, pb, input_type = "all", directional_rep_claim = FALSE, 
                        variant = "adaptive", alpha = q)
rv18 <- rep(1, J)
rv18[as.numeric(res.rv18$results_table$name)] <- res.rv18$results_table$r_value
snp.radjust <- snpID[which(rv18<=q)]


# JUMP
library(JUMP)
jump.obj <- JUMP(pa, pb, q)
jump.thr <- jump.obj$jump.thr
p.max <- jump.obj$p.max
snp.jump <- snpID[which(p.max<=jump.thr)]

# adaFilter
p_matrix = cbind(pa, pb)
ada_result <- adaFilter(p_matrix, r = 2, type.I.err = "FDR", alpha = q)
ada_result %>% head()
snp.ada <- snpID[which(ada_result$decision == 1)]
# snp.discovery = snp.discovery %>% append(list(ada = snp.ada))

# summary
snp.discovery <- list(hmm = snp.hmm, stareg = snp.stareg, radjust = snp.radjust, 
                      jump = snp.jump, bh = snp.bh, maxp = snp.maxp, marr = snp.marr)
file.wd = "C:/Users/pl19f/OneDrive - Florida State University/5. dissertation/11. Hidder Markov Model/3. Codes/2. results/"
if(save_file){
  save(snp.discovery,
     file = paste0(file.wd, "snp.discovery_20240919"))
}# end if save_file
load(paste0(file.wd, "snp.discovery_20240919"))

num_rep <- c(length(snp.hmm), length(snp.stareg), length(snp.maxp), length(snp.bh), 
             length(snp.radjust), length(snp.jump), length(snp.marr))
names(num_rep) <- c("CoHiM", "STAREG", "MaxP", "ad hoc BH", "radjust", "JUMP", "MaRR")
q
num_rep

{
  snp.bh = snp.discovery$bh; snp.maxp = snp.discovery$maxp; snp.jump = snp.discovery$jump
  snp.radjust = snp.discovery$radjust; snp.stareg = snp.discovery$stareg; snp.hmm = snp.discovery$hmm
  snp.marr = snp.discovery$marr; snp.ada = snp.discovery$ada
}

library(UpSetR)
listInput <- list(Ad_hoc_BH = snp.discovery$bh, MaxP = snp.discovery$maxp, JUMP = snp.discovery$jump, 
                  radjust = snp.discovery$radjust, STAREG = snp.discovery$stareg, HMM = snp.discovery$hmm,
                  MaRR = snp.discovery$marr, AdaFilter = snp.discovery$ada)
names(listInput)[1] <- "ad hoc BH"
names(listInput)[6] <- "CoHiM"
scale.max <- listInput %>% lapply(length) %>% unlist() %>% max()
# pdf(file="Data analysis/Mouse olfactory bulb (ST Rep1 + ST Rep8)/result.pdf",onefile = FALSE,width=7.5,height=5)
upset(fromList(listInput), nsets=8, sets.x.label = "# replicable SNPs", 
      matrix.color = "darkblue", main.bar.color = "darkblue", 
      sets.bar.color = "darkblue", order.by = "freq", text.scale = 1.8,
      set_size.show = TRUE, set_size.scale_max = scale.max*1.4,
      point.size = 2.5)

# library(VennDiagram)
# venn.diagram(
#   x = list(snp.hmm, snp.bh, snp.maxp, snp.stareg, snp.rv18),
#   category.names = c("rLIS", "BH", "MaxP", "STAREG", "RV"),
#   fill = c("#8FBC8F", "#F4B183", "#56B4E9", "#9370DB", "gray"),
#   filename = 'venn.tiff',
#   output = TRUE,
#   margin = 0.05,
#   cex = 1.5,
#   cat.cex = 1.5
# )

f.dframe <- data.frame(pa = sort(pa), f1 = sort(res.hmm$f1, decreasing = TRUE),
                      pb = sort(pb), f2 = sort(res.hmm$f2, decreasing = TRUE))
p1 = ggplot(f.dframe, aes(x = pa, y = f1)) + 
  geom_step(linewidth = 2, color = 'blue') +
  scale_x_continuous(limits = c(1e-21, 1),
                     trans = log2_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_y_continuous(limits = c(0, 4e17))+
  xlab("p-value") + 
  ylab("Non-null density") +
  ggtitle("Male") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
print(p1)
  
p2 = ggplot(f.dframe, aes(x = pb, y = f2)) + 
  geom_step(linewidth = 2, color = 'blue') +
  scale_x_continuous(limits = c(2e-20, 1),
                     trans = log2_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_y_continuous(limits = c(0, 1e15))+
  xlab("p-value") + 
  ylab("Non-null density") +
  ggtitle("Female") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
print(p2)
print(p1 + p2)

hmm.sex = res.hmm
snp.others = snp.stareg %>% union(snp.maxp) %>% union(snp.bh) %>% union(snp.radjust) %>% union(snp.jump) %>% union(snp.marr)
snp.unique = snp.hmm %>% setdiff(snp.others)
file_a <- data.a[snp.unique, ]
file_b <- data.b[snp.unique, ]
file.wd = "C:/Users/pl19f/OneDrive - Florida State University/5. dissertation/11. Hidder Markov Model/3. Codes/2. results/"
if(save_file){
  write.table(file_a, file = paste0(file.wd, "sex_a_0607.txt"), row.names = TRUE)
  write.table(file_b, file = paste0(file.wd, "sex_b_0607.txt"), row.names = TRUE)
} # end if save_file

# validation -----------------------------------
library(readr)
file.wd = "C:/Users/pl19f/OneDrive - Florida State University/5. dissertation/11. Hidder Markov Model/3. Codes/1. Data/"

# validation1 = readr::read_tsv(paste0(file.wd, "gwas-association-downloaded_2023-12-09-MONDO_0005148.tsv"))
validation1 = readr::read_tsv(paste0(file.wd, "gwas-association-downloaded_2024-02-27-MONDO_0005148.tsv"))
snp.validated <- snp.unique %>% intersect(validation1$SNPS)
length(snp.validated)

loci.validated <- as.integer(data.a[snp.unique, "POSITION"]) %>% intersect(strtoi(validation1$CHR_POS))
length(loci.validated)

snp.2gene <- snp.unique %>% setdiff(snp.validated)
gene.map = snp2Gene(snp.2gene)
gene.id = gene.map$map$gene_id %>% unique()
gene.name = getGeneMap(gene.id)$gene_map$gene_name %>% sort()
gene.name

eg.name = c("JAZF1", "ADAMTS9", "NOTCH2", "CDC123", "THADA")
eg.id = c(221895, 100507098,      4853,      8872,     63892)
for (i in 1:5) {
  print(paste("SNPs mapped to gene", eg.name[i], ':'))
  print(gene.map$map[gene.map$map$gene_id==eg.id[i], 1])
}


pi = res.hmm$A %>% t() %>% eigen()
pi = pi$vectors[,1]
pi = pi/sum(pi)
pi %>% round(3)
res.hmm$A %>% round(4)

# Manhattan plot --------------------------------------------
data.common = cbind(data.a[snpID, c(1:3, 9)], data.b[snpID, 9], maxp, padj.eb, padj.hmm)
colnames(data.common)[4:8] <- c("pa", "pb", "maxp", "LIS", "rLIS")
par(mfrow = c(1, 3))
manhattan(data.common, chr = "CHROMOSOME", bp = "POSITION", snp = "ID", genomewideline = F,
          col = c('grey', 'orange'), ylim = c(0, 30), cex.lab = 1.5, cex.axis = 1.5, mgp = c(2.5, 1, 0),
          suggestiveline = -log10(thres.maxp),
          p = 'maxp', ylab = expression(-log[10](p[max])))
manhattan(data.common, chr = "CHROMOSOME", bp = "POSITION", snp = "ID", genomewideline = F,
          col = c('grey', 'orange'), ylim = c(0, 30), cex.lab = 1.5, cex.axis = 1.5, mgp = c(2.5, 1, 0),
          suggestiveline = -log10(thres.eb),
          p = 'LIS', ylab = expression(-log[10](Lfdr)))
manhattan(data.common, chr = "CHROMOSOME", bp = "POSITION", snp = "ID", genomewideline = F,
          col = c('grey', 'orange'), ylim = c(0, 30),, cex.lab = 1.5, cex.axis = 1.5, mgp = c(2.5, 1, 0),
          suggestiveline = -log10(thres.hmm),
          p = 'rLIS', ylab = expression(-log[10](rLIS)))

