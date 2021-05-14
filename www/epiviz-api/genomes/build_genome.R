library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

library(ensembldb)
library(Homo.sapiens)
library(OrganismDbi)
library(Mus.musculus)

# use a GENOME, txdb or orgdb or ensembl
# org <- makeOrganismDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
#                                        dataset="hsapiens_gene_ensembl",
#                                        transcript_ids=NULL,
#                                        circ_seqs=DEFAULT_CIRC_SEQS,
#                                        filter="",
#                                        id_prefix="ensembl_",
#                                        host="www.ensembl.org",
#                                        port=80,
#                                        miRBaseBuild=NA,
#                                        keytype = "ENSEMBL",
#                                        orgdb = NA)

# Extract genes and exon positions
org <- Mus.musculus


library(BSgenome.Cjacchus.UCSC.calJac3)
# marmoset
org <- makeOrganismDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                       dataset="cjacchus_gene_ensembl",
                                       transcript_ids=NULL,
                                       circ_seqs=DEFAULT_CIRC_SEQS,
                                       filter="",
                                       id_prefix="ensembl_",
                                       host="www.ensembl.org",
                                       port=80,
                                       miRBaseBuild=NA,
                                       keytype = "ENSEMBL",
                                       orgdb = NA)

org <- TxDb.Cjacchus.BioMart.ENSEMBLMARTENSEMBL.ASM275486v1

org <- TxDb.Hsapiens.UCSC.hg38.knownGene
# for homo.sapiens, use GENEID, SYMBOL
# for txdb objects, use GENEID
genome <- GenomicFeatures::genes(org, columns=c("GENEID"))

exons <- GenomicFeatures::exonsBy(org, by="gene")

ids <- as.character(genome$GENEID)
exons <- reduce(ranges(exons)[ids])

exons_start <- unname(lapply(start(exons), paste, collapse=","))
exons_end <- unname(lapply(end(exons), paste, collapse=","))

# create data frame for the genome
genome$exons_start <- exons_start
genome$exons_end <- exons_end
genome$Gene <- genome$SYMBOL

genome <- sortSeqlevels(genome)
genome <- sort(genome, ignore.strand=TRUE)


genome_df <- as.data.frame(genome)
genome_df$gene <- genome_df$GENEID
# genome_df = subset(genome_df, select = -c(SYMBOL))

# save genome as tsv -> index by tabix!
write.table(genome_df, file="marmoset.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  