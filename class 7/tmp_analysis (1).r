
# One time only install of bio3d package
#install.packages(bio3d)

# loads bio3d functions
library(bio3d)

## Read your alignment file
raw <- read.fasta("alignment.fa", rm.dup=FALSE)

new <- raw
new$id <- paste0(1:length(raw$id), "_", raw$id)
write.fasta(new, file="alignment_cleanIDs.fasta")

## Sequence identity
ide <- seqidentity(new)

hc.all <- hclust(as.dist(1-ide))

## Trim to your region 49 to 157
trim.aln <- new
trim.aln$ali[,49:157]

## Sequence identity of trimed alignment section.
trim.ide <- seqidentity(trim.aln)

hc.trim <- hclust(as.dist(1-trim.ide))

plot(hc.trim)

pdf("plot_trim_tree.pdf", width=14)
plot(hc.trim)
dev.off()

pdf("plot_full_tree.pdf", width=14)
plot(hc.all)
dev.off()



