suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(VariantAnnotation))
fl <- ("hg18_refGene.gtf")
txdb <- makeTxDbFromGFF(file=fl,format="gtf")
fl <- system.file("extdata", "SRR1518358_default.vcf",package="VariantAnnotation")
vcf <- readVcf('SRR1518358_default.vcf', "hg18")

#head(seqlevels(vcf))
#head(seqlevels(txdb))
#intersect(seqlevels(vcf), seqlevels(txdb))

#vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))

intersect(seqlevels(vcf), seqlevels(txdb))

region <- IntergenicVariants(upstream=70000, downstream=70000)
loc_int <- locateVariants(vcf, txdb, region)
mcols(loc_int)[c("LOCATION", "PRECEDEID", "FOLLOWID")]

p_ids <- unlist(loc_int$PRECEDEID, use.names=FALSE)
exp_ranges <- rep(loc_int,  elementNROWS(loc_int$PRECEDEID))

distance(exp_ranges, txdb, id=p_ids, type="gene")

cdsbytx <- cdsBy(txdb)
var_cds <-locateVariants(vcf, cdsbytx, CodingVariants())
df <- data.frame()

intbytx <- intronsByTranscript(txdb)
locateVariants(vcf, intbytx, IntronVariants())

fa <- open(FaFile('all.fa'))
coding1 <- predictCoding(vcf,txdb,fa)
final <- data.frame(seqnames=seqnames(coding1),start=start(coding1),ref=coding1$REF,alt=coding1$ALT,type=coding1$CONSEQUENCE,gene_id=coding1$GENEID)
test <- final
final <- test[,c("seqnames","start","ref","alt.value","type","gene_id")]
ans <- final %>% group_by(seqnames,start,ref,alt.value,type) %>% summarise(gene_id = toString(gene_id))
write.table(ans, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)
