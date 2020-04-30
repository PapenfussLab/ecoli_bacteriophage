library(StructuralVariantAnnotation)
library(tidyverse)

# argv = parse_args(argp, argv=c("--input", "D:/dev/ecoli_phage_resistence/gridss/bl21_de3_.ASM956v1_gridss.vcf", "--normalordinal", "1", "--tumourordinal", "2", "3", "--output", "D:/dev/ecoli_phage_resistence/gridss/123_bl21_de3_.ASM956v1_gridss.somatic.vcf", "-f", "D:/dev/ecoli_phage_resistence/gridss/123_bl21_de3_.ASM956v1_gridss.somatic.full.vcf", "--scriptdir", "D:/dev/gridss/scripts", "--gc"))
# argv = parse_args(argp, argv=c("--input", "D:/dev/ecoli_phage_resistence/gridss/bl21_de3_.ASM956v1_gridss.vcf", "--normalordinal", "4", "--tumourordinal", "5", "6", "--output", "D:/dev/ecoli_phage_resistence/gridss/456_bl21_de3_.ASM956v1_gridss.somatic.vcf", "-f", "D:/dev/ecoli_phage_resistence/gridss/456_bl21_de3_.ASM956v1_gridss.somatic.full.vcf", "--scriptdir", "D:/dev/gridss/scripts", "--gc"))

for (file in c(
	"D:/dev/ecoli_phage_resistence/gridss/123_bl21_de3_.ASM956v1_gridss.somatic.full.vcf",
	"D:/dev/ecoli_phage_resistence/gridss/456_bl21_de3_.ASM956v1_gridss.somatic.full.vcf")) {
	svvcf = readVcf(file)
	af = as.matrix(info(svvcf)$TAF)
	af[is.na(af)] = 0
	# at least 10% AF (20% for single breakends)
	svvcf = svvcf[rowSums(af) > 0.2]
	writeVcf(svvcf, file=paste0(file, ".af10.vcf"))
}
svvcf = readVcf("D:/dev/ecoli_phage_resistence/gridss/123_bl21_de3_.ASM956v1_gridss.vcf")
bpgr = breakpointRanges(svvcf)
begr = breakendRanges(svvcf)
bpgr$vaf1 = gridss_bp_af(bpgr, svvcf, 1)
bpgr$vaf2 = gridss_bp_af(bpgr, svvcf, 2)
bpgr$vaf3 = gridss_bp_af(bpgr, svvcf, 3)
bpgr$vaf4 = gridss_bp_af(bpgr, svvcf, 4)
bpgr$vaf5 = gridss_bp_af(bpgr, svvcf, 5)
bpgr$vaf6 = gridss_bp_af(bpgr, svvcf, 6)
begr = breakendRanges(vcf)
begr$vaf1 = gridss_be_af(begr, svvcf, 1)
begr$vaf2 = gridss_be_af(begr, svvcf, 2)
begr$vaf3 = gridss_be_af(begr, svvcf, 3)
begr$vaf4 = gridss_be_af(begr, svvcf, 4)
begr$vaf5 = gridss_be_af(begr, svvcf, 5)
begr$vaf6 = gridss_be_af(begr, svvcf, 6)
gr = c(bpgr, begr)
table(gr$vaf1 > 0.3)
table(gr$vaf2 > 0.3)
table(gr$vaf3 > 0.3)
table(gr$vaf4 > 0.3)
table(gr$vaf5 > 0.3)
table(gr$vaf6 > 0.3)
gr[gr$vaf6 + gr$vaf5 > 0.2 & gr$vaf4 < 0.1]
gr[gr$vaf3 + gr$vaf4 > 0.2 & gr$vaf1 < 0.1]



snvvcf = readVcf("D:/dev/ecoli_phage_resistence/variants/bl21_de3_.ASM956v1.bb.vcf")
for (ordinal in c(2,3,5,6)) {
	normal_ordinal = c(1,1,1,4,4,4)[ordinal]
	filtered_vcf = snvvcf[geno(snvvcf)$PF[,ordinal] == "PASS" & geno(snvvcf)$PF[,normal_ordinal] == "FAIL"]
	writeVcf(filtered_vcf, paste0("D:/dev/ecoli_phage_resistence/variants/bl21_de3_.ASM956v1.bb.vcf_somatic", ordinal, ".vcf"))
}

