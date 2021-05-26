#Script for implementing additional filters to vcf files. All variants are written to a text file. 
#The script also produces plots of variant allele frequency against read depth in each sample

library(VariantAnnotation)
library(parallel)
library(ggplot2)
library(jsonlite)


# A function to create sensible mutation types:
mutType <- function(ref, alt, collapse=TRUE){
  if(length(ref) != length(alt)) stop('ref and alt length differ')
  ref <- as.character(ref)
  alt <- as.character(alt)
  res <- sprintf('%s>%s', ref, alt)
  res[ref==alt] <- NA
  if(identical(collapse, TRUE)){
    mut.types <- c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
    res[res=='G>T'] <- 'C>A'
    res[res=='G>C'] <- 'C>G'
    res[res=='G>A'] <- 'C>T'
    res[res=='A>T'] <- 'T>A'
    res[res=='A>G'] <- 'T>C'
    res[res=='A>C'] <- 'T>G'
  } else {
    mut.types <- sort(c(c('G>T', 'G>C', 'G>A', 'A>T', 'A>G', 'A>C'), c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')))
  }
  return(factor(res, levels=mut.types))
}

#Read in sample id mappings
sample.map <- fromJSON("../metadata/sample-map.json") 
sample.map.df <- read.csv("../metadata/sample-map.csv", sep = "")



#calculate the 'new' VAF with the reads passing Mutect2 filtering - this new VAF will then be used for the additional filtering 

#function for calculating AF
calc_AF <- function(vcf){
sample.id <- colnames(geno(vcf)$AD)
AD <- sapply(geno(vcf)$AD[,sample.id], '[[',2 )
total.depth <- geno(vcf)$DP[,sample.id]
return(AD/total.depth)
}

#will test if MBQ meets filtering requirements
MBQ_filter <- function(vcf){
  MBQ <- info(vcf)$MBQ
  MBQ.ok <- as.matrix(MBQ >= params$min.MBQ)
  return(apply(MBQ.ok, 1, all)) #need to make sure that all samples meet this quality 
}  

#will test depth meets filter
depth_filter <- function(vcf){
  depth <- geno(vcf)$DP
  depth.ok <- depth >= params$min.depth
  return(apply(depth.ok, 1, all)) #make sure all samples meet this depth
} 

depth_filter_2 <- function(vcf){
  depth.mal <- geno(vcf)$DP[, malignant.id] #extract the depth for the malignant sample
  depth.normals <- apply(as.matrix(geno(vcf)$DP[, normal.id]), 1, mean) #extract the mean depth from the two normals 
  depth <- cbind(depth.mal, depth.normals) #combine the two so they can be tested against the threshold
  depth.ok <- depth >= params$min.depth
  return(apply(depth.ok, 1, all))
} 


#will test VAF in normal sample meets filters
VAF_normal_filter <- function(vcf){
  VAF <- geno(vcf)$AF[, normal.id]
  return(VAF <= params$max.normAF)
} 

#will test VAF in normal sample meets filters - need to average the two normals 
VAF_normal_filter_2 <- function(vcf){
  VAF <- as.data.frame(geno(vcf)$AF[, normal.id])
  VAF <- as.data.frame(lapply(VAF, unlist))
  VAF.avg <- apply(VAF, 1, mean) #need to extract the mean and then test this against the threshold 
  return(VAF.avg <= params$max.normAF)
} 

#Will test VAF in malignant sample meets filters
VAF_malignant_filter <- function(vcf){
  depth <- geno(vcf)$DP[, malignant.id]
  VAF <- geno(vcf)$AF[, malignant.id]
  res <- rep(NA, length(depth))
  res[depth<30] <-  VAF[depth<30] >= params$min.malignantAF.DPthirty #need to do this before <20 
  res[depth<20] <-  VAF[depth<20] >= params$min.malignantAF.DPtwenty #this will overwrite the <30 filter for those whose depth is below 20  
  res[depth>=30] <- VAF[depth>=30] >= params$min.malignantAF
  return(res)
}


#create empty data frames for recording all fitlered variants  
d_all_variants <- data.frame(stringsAsFactors = FALSE)


#for loop for writing tables and creating graphs 
  for (i in names(sample.map)){

#read in input vcf files
  io <- list('input_NT.vcf' = paste('../', 'output-vcf/',i, 'NT-annotated.vcf.gz', sep = ""),
  'input_raw_NT.vcf' = paste('../', 'output-vcf/', i, 'NT-mutect.vcf.gz', sep = ""),
  'input_ND.vcf' = paste('../', 'output-vcf/', i, 'ND-annotated.vcf.gz', sep = ""),
  'input_raw_ND.vcf' = paste('../', 'output-vcf/', i, 'ND-mutect.vcf.gz', sep = "")
  )

  params <- list(
  "sample" = i,
  'malignant.id' = c(sample.map[i][[1]]$tumour$malignant, sample.map[i][[1]]$dediff$malignant),
  'normal.id' = as.vector(sample.map.df[sample.map.df$sample == i & sample.map.df$type == "normal",1]),
  'min.depth' = 10,
  'max.normAF' = 0.05, # 5%
  'min.malignantAF' = 0.05, # 5% - this is the base filter
  'min.malignantAF.DPtwenty' = 0.2, #20% - this is for reads whose read depth is <20 
  'min.malignantAF.DPthirty' = 0.1, #10% - this is for reads whose read depth is <30
  'min.MBQ' = 30) #set minimum median base quality to 30 
 

#No conventional melanoma from 'PD42799' individual
if(i != "PD42799"){
  # Load the input VCF files - normal vs tumour:
message(sprintf('reading annotated input VCF from "%s"...', io$input_NT.vcf))
suppressWarnings(input_NT.v <- readVcf(io$input_NT.vcf))
message(sprintf('%s variants read from file.', prettyNum(nrow(input_NT.v), big.mark=',')))

message(sprintf('reading raw input VCF from "%s"...', io$input_raw_NT.vcf))
suppressWarnings(input_raw_NT.v <- readVcf(io$input_raw_NT.vcf))
message(sprintf('%s variants read from file.', prettyNum(nrow(input_raw_NT.v), big.mark=',')))

# Check that the normal and malignant ids match one of the ids described in the parameter list:
#create a variable for the malignant and normal id depending on which are present in both the vcf header and the parameter list 

if(length(intersect(params$malignant.id, samples(header(input_NT.v)))) == 0) stop(sprintf('ERROR: malignant sample %s not in VCF file', params$normal.id))
malignant.id <- intersect(params$malignant.id, samples(header(input_NT.v)))

if(length(intersect(params$normal.id, samples(header(input_NT.v)))) == 0) stop(sprintf('ERROR: normal sample %s not in VCF file', params$normal.id))
normal.id <- intersect(params$normal.id, samples(header(input_NT.v)))

message(sprintf('normal sample: %s ', normal.id))
message(sprintf('malignant sample: %s ', malignant.id))

# Pull out the normal and malignant IDs (in case there are multiple other samples):
v_annotated_NT <- input_NT.v[,c(normal.id, malignant.id)]
v_raw_NT <- input_raw_NT.v[,c(normal.id, malignant.id)]


#Require that there are only SNVs 
v_annotated_NT <- v_annotated_NT[isSNV(v_annotated_NT)]
v_raw_NT <-  v_raw_NT[isSNV(v_raw_NT)]

# For the annotated VCF - Require that the filter tag is "PASS":
v_annotated_NT <- v_annotated_NT[filt(v_annotated_NT) == 'PASS',]

#Create a copy of the annotated VCF which will have our extra filters applied 
v_filtered_NT <- v_annotated_NT  


#apply function and get out the results
new_VAF <- calc_AF(v_filtered_NT)

#change the vcf allele frequencies for the new calculated frequencies 
geno(v_filtered_NT)$AF <- new_VAF

#Apply the filter functions to the vcf
MBQ.ok <- MBQ_filter(v_filtered_NT)

#for depth - apply depth_filter if there is only 1 normal, or depth_filter_2 if there are two normals 
depth.ok <- as.array(if(length(normal.id) == 1){
  depth_filter(v_filtered_NT)} else {
  depth_filter_2(v_filtered_NT)
  })

#for VAF - apply VAF_normal_filter if there is only 1 normal, or VAF_normal_filter_2 if there are two normals 
VAF_normal.ok <- as.array(if(length(normal.id) ==1){
  VAF_normal_filter(v_filtered_NT)} else {
     VAF_normal_filter_2(v_filtered_NT)
  })

VAF_malignant.ok <- VAF_malignant_filter(v_filtered_NT)

#compile the results from all the filters 
all_filters <- rbind(MBQ.ok, depth.ok, VAF_normal.ok, VAF_malignant.ok)

#for each variant check all filters have been passed 
filters.ok <- apply(all_filters, 2, all)

#only include variants for which all filters have been passed 
v_filtered_NT <- v_filtered_NT[filters.ok, ]
message(sprintf('%s variants after filtering.', prettyNum(nrow(v_filtered_NT), big.mark=',')))  

#data frame of variants from the normal-tumour vcfs
v.data_NT <- do.call(rbind, lapply(colnames(v_raw_NT), function(sample.id){
  message(sprintf('  reformatting sample %s...', sample.id))
  res <- data.frame(
    'id' = c(rownames(rowData(v_raw_NT)), rownames(rowData(v_annotated_NT)), rownames(v_filtered_NT)),
   'sample' = factor(sample.map.df[sample.map.df$id == sample.id, 'sample'], levels=sort(unique(sample.map.df$sample))),
   'sample_id' = sample.id,
    'type' = factor(sample.map.df[sample.map.df$id == sample.id, 'type'], levels=c('normal', 'tumour', 'dediff')),
    'stage' = factor(c(rep("raw", nrow(rowData(v_raw_NT))), rep("annotated", nrow(rowData(v_annotated_NT))), rep("filtered", nrow(rowData(v_filtered_NT)))), levels = c('raw', "annotated", "filtered")), #stage of filtering - here all will be raw 
    'ref' = factor(c(as.character(ref(v_raw_NT)), as.character(ref(v_annotated_NT)), as.character(ref(v_filtered_NT))), levels=DNA_BASES),
    'alt' = factor(c(as.character(unlist(alt(v_raw_NT))), as.character(unlist(alt(v_annotated_NT))), as.character(unlist(alt(v_filtered_NT)))), levels=DNA_BASES),
   'GT' = c(geno(v_raw_NT)$GT[,sample.id],geno(v_annotated_NT)$GT[,sample.id], geno(v_filtered_NT)$GT[,sample.id]), 
    'mutation' = NA,
    'depth' = c(geno(v_raw_NT)$DP[,sample.id],geno(v_annotated_NT)$DP[,sample.id], geno(v_filtered_NT)$DP[,sample.id]), 
    'ref.depth' = c(sapply(geno(v_raw_NT)$AD[,sample.id], '[[', 1), sapply(geno(v_annotated_NT)$AD[,sample.id], '[[', 1), sapply(geno(v_filtered_NT)$AD[,sample.id], '[[', 1)),
    'alt.depth' = c(sapply(geno(v_raw_NT)$AD[,sample.id], '[[', 2),sapply(geno(v_annotated_NT)$AD[,sample.id], '[[', 2), sapply(geno(v_filtered_NT)$AD[,sample.id], '[[', 2)),
    'VAF' = c(unlist(geno(v_raw_NT)$AF[,sample.id]), unlist(geno(v_annotated_NT)$AF[,sample.id]), unlist(geno(v_filtered_NT)$AF[,sample.id])),
    "chr" = c(v_raw_NT@rowRanges@seqnames, v_annotated_NT@rowRanges@seqnames, v_filtered_NT@rowRanges@seqnames),
    "pos" = c(v_raw_NT@rowRanges@ranges@start, v_annotated_NT@rowRanges@ranges@start, v_filtered_NT@rowRanges@ranges@start)
  )
  res <- res[!is.na(res$VAF),]
  return(res)
}))
 
v.data_NT$mutation <- mutType(v.data_NT$ref, v.data_NT$alt, collapse=TRUE) #add mutation type using the mutType function
v.data_NT$VAF <- v.data_NT$alt.depth / v.data_NT$depth #change the VAF so it is now the 'new' true VAF 
v.data_NT[v.data_NT$alt.depth == 0,]$VAF <- NA #if VAF is 0 then change it to NA (plot uses log scale so it cannot be 0)
v.data_NT <- v.data_NT[!is.na(v.data_NT$VAF),]



#add data to the df containing all variants 
d_all_variants <- rbind(d_all_variants, v.data_NT)

}
  
#Filter ND file second 
# Load the input VCF files:
message(sprintf('reading annotated input VCF from "%s"...', io$input_ND.vcf))
suppressWarnings(input_ND.v <- readVcf(io$input_ND.vcf))
message(sprintf('%s variants read from file.', prettyNum(nrow(input_ND.v), big.mark=',')))

message(sprintf('reading raw input VCF from "%s"...', io$input_raw_ND.vcf))
suppressWarnings(input_raw_ND.v <- readVcf(io$input_raw_ND.vcf))
message(sprintf('%s variants read from file.', prettyNum(nrow(input_raw_ND.v), big.mark=',')))


# Check that the normal and malignant ids match one of the ids described in the parameter list:
#create a varible for the malignant and normal id depending on which are present in both the vcf header and the parameter list 

if(length(intersect(params$malignant.id, samples(header(input_ND.v)))) == 0) stop(sprintf('ERROR: malignant sample %s not in VCF file', params$normal.id))
malignant.id <- intersect(params$malignant.id, samples(header(input_ND.v)))

if(length(intersect(params$normal.id, samples(header(input_ND.v)))) == 0) stop(sprintf('ERROR: normal sample %s not in VCF file', params$normal.id))
normal.id <- intersect(params$normal.id, samples(header(input_ND.v)))

message(sprintf('normal sample: %s ', normal.id))
message(sprintf('malignant sample: %s ', malignant.id))

    
# Pull out the normal and malignant IDs (in case there are multiple other samples):
v_annotated_ND <- input_ND.v[,c(normal.id, malignant.id)]
v_raw_ND <- input_raw_ND.v[,c(normal.id, malignant.id)]

#Require that there are only SNVs 
v_annotated_ND <- v_annotated_ND[isSNV(v_annotated_ND)]
v_raw_ND <-  v_raw_ND[isSNV(v_raw_ND)]

# For the annotated VCF - Require that the filter tag is "PASS":
v_annotated_ND <- v_annotated_ND[filt(v_annotated_ND) == 'PASS',]

#Create a copy of the annotated VCF which will have our extra filters applied 
v_filtered_ND <- v_annotated_ND

#apply function and get out the results
new_VAF <- calc_AF(v_filtered_ND)

#change the vcf allele frequencies for the new calculated frequencies 
geno(v_filtered_ND)$AF <- new_VAF
  
#Apply the filter functions to the vcf
MBQ.ok <- MBQ_filter(v_filtered_ND)

#filter depends on number of normal samples in vcf 
depth.ok <- as.array(if(length(normal.id) == 1){
  depth_filter(v_filtered_ND)} else {
  depth_filter_2(v_filtered_ND)
  })

#filter depends on number of normal samples in vcf 
VAF_normal.ok <- as.array(if(length(normal.id) ==1){
  VAF_normal_filter(v_filtered_ND)} else {
     VAF_normal_filter_2(v_filtered_ND)
  })

VAF_malignant.ok <- VAF_malignant_filter(v_filtered_ND)

#compile the results from all the filters 
all_filters <- rbind(MBQ.ok, depth.ok, VAF_normal.ok, VAF_malignant.ok)

#for each variant check all filters have been passed 
filters.ok <- apply(all_filters, 2, all)

#only include variants for which all filters have been passed 
v_filtered_ND <- v_filtered_ND[filters.ok, ]
message(sprintf('%s variants after filtering.', prettyNum(nrow(v_filtered_ND), big.mark=',')))


#create dataframe of all variants from the normal-dediff vcfs
v.data_ND <- do.call(rbind, lapply(colnames(v_raw_ND), function(sample.id){
  message(sprintf('  reformatting sample %s...', sample.id))
  res <- data.frame(
    'id' = c(rownames(rowData(v_raw_ND)), rownames(rowData(v_annotated_ND)), rownames(v_filtered_ND)),
    'sample' = factor(sample.map.df[sample.map.df$id == sample.id, 'sample'], levels=sort(unique(sample.map.df$sample))),
    'sample_id' = sample.id,
    'type' = factor(sample.map.df[sample.map.df$id == sample.id, 'type'], levels=c('normal', 'tumour', 'dediff')),
    'stage' = factor(c(rep("raw", nrow(rowData(v_raw_ND))), rep("annotated", nrow(rowData(v_annotated_ND))), rep("filtered", nrow(rowData(v_filtered_ND)))), levels = c('raw', "annotated", "filtered")), #stage of filtering - here all will be raw 
    'ref' = factor(c(as.character(ref(v_raw_ND)), as.character(ref(v_annotated_ND)), as.character(ref(v_filtered_ND))), levels=DNA_BASES),
    'alt' = factor(c(as.character(unlist(alt(v_raw_ND))), as.character(unlist(alt(v_annotated_ND))), as.character(unlist(alt(v_filtered_ND)))), levels=DNA_BASES),
    'GT' = c(geno(v_raw_ND)$GT[,sample.id],geno(v_annotated_ND)$GT[,sample.id], geno(v_filtered_ND)$GT[,sample.id]), 
    'mutation' = NA,
    'depth' = c(geno(v_raw_ND)$DP[,sample.id],geno(v_annotated_ND)$DP[,sample.id], geno(v_filtered_ND)$DP[,sample.id]), 
    'ref.depth' = c(sapply(geno(v_raw_ND)$AD[,sample.id], '[[', 1), sapply(geno(v_annotated_ND)$AD[,sample.id], '[[', 1), sapply(geno(v_filtered_ND)$AD[,sample.id], '[[', 1)),
    'alt.depth' = c(sapply(geno(v_raw_ND)$AD[,sample.id], '[[', 2),sapply(geno(v_annotated_ND)$AD[,sample.id], '[[', 2), sapply(geno(v_filtered_ND)$AD[,sample.id], '[[', 2)),
    'VAF' = c(unlist(geno(v_raw_ND)$AF[,sample.id]), unlist(geno(v_annotated_ND)$AF[,sample.id]), unlist(geno(v_filtered_ND)$AF[,sample.id])),
    "chr" = c(v_raw_ND@rowRanges@seqnames, v_annotated_ND@rowRanges@seqnames, v_filtered_ND@rowRanges@seqnames),
    "pos" = c(v_raw_ND@rowRanges@ranges@start, v_annotated_ND@rowRanges@ranges@start, v_filtered_ND@rowRanges@ranges@start)
  )
  res <- res[!is.na(res$VAF),]
  return(res)
}))



 
v.data_ND$mutation <- mutType(v.data_ND$ref, v.data_ND$alt, collapse=TRUE) #add mutation type using the mutType function
v.data_ND$VAF <- v.data_ND$alt.depth / v.data_ND$depth #change the VAF so it is now the 'new' true VAF 
v.data_ND[v.data_ND$alt.depth == 0,]$VAF <- NA #if VAF is 0 then change it no NA (plot uses log scale so does not like zeros)
v.data_ND <- v.data_ND[!is.na(v.data_ND$VAF),]

#add data to the df containing all variants 
d_all_variants <- rbind(d_all_variants, v.data_ND)
return(d_all_variants)

}

#change colnames and write the shared variant table and C>T proportion table to file 

d_all_variants <- d_all_variants[!is.na(d_all_variants$VAF),]

#alter d_all_variants so that the sample annotations are correct
d_all_variants[d_all_variants$sample_id == "PD42798a", 'type'] <- "dediff"
d_all_variants[d_all_variants$sample_id == "PD45781c", 'type'] <- "dediff"

#write all variants which have been filtered to a csv file 
write.table(d_all_variants[d_all_variants$stage == "filtered",], file = "all_variants_filtered.csv", sep = " ", row.names = F)






