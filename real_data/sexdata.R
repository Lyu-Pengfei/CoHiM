library(tidyverse)
# processing the data ----------------------------------------------------
data.a <- read.table("./real_data/DIAGRAM.Morris2012.males.txt", sep = "", header = TRUE)
data.b <- read.table("./real_data/DIAGRAM.Morris2012.females.txt", sep = "", header = TRUE) # q = 1e-5

rownames(data.a) <- data.a[, "ID"]
rownames(data.b) <- data.b[, "ID"]
snpID <- intersect(data.a$ID, data.b$ID)

pa <- data.a[snpID, "P_VALUE"]
pb <- data.b[snpID, "P_VALUE"]

q = 1e-5
J = length(snpID)

# CoHIM-----------------------------------------------------------------------
res.hmm <- CoHiM(rbind(pa, pb), q)
snp.hmm <- snpID[res.hmm$rej_index]
rej_num.hmm = sum(res.hmm$rej_index)
rLIS.hmm <- res.hmm$rLIS_mat
thres.hmm = sort(rLIS.hmm)[rej_num.hmm]

# Estimations----------------------------------------
# density function estimations
f.dframe <- data.frame(pa = sort(pa), f1 = sort(res.hmm$f1_mat, decreasing = TRUE),
                      pb = sort(pb), f2 = sort(res.hmm$f2_mat, decreasing = TRUE))
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

# estimated stationary and transition probabilities
pi = res.hmm$pi_mat
pi %>% round(3)
res.hmm$A_list[[1]] %>% round(4)
                                          
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

