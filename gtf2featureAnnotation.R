#!/usr/bin/env Rscript

# This script parses the GTF file to create a feature-wise annotation file with
# mitochondrial features flagged, to assist in annotation and QC of single-cell
# expression data analysis.

suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(optparse))

ucfirst <- function (str) {
  paste(toupper(substring(str, 1, 1)), tolower(substring(str, 2)), sep = "")
}

die <- function(message){
  write(message, stderr())
  q(status = 1)
}

cleanlist <- function(str){
  tolower(unlist(strsplit(str, ',')))
}

cl <- commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option(
    c("-g", "--gtf-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to a valid GTF file"
  ),
  make_option(
    c("-t", "--feature-type"),
    action = "store",
    default = 'gene',
    type = 'character',
    help = 'Feature type to use (default: gene)'
  ),
  make_option(
    c("-f", "--first-field"),
    action = "store",
    default = 'gene_id',
    type = 'character',
    help = 'Field to place first in output table (default: gene_id)'
  ),
  make_option(
    c("-r", "--no-header"),
    action = "store_false",
    default = TRUE,
    type = 'logical',
    help = 'Suppress header on output'
  ),
  make_option(
    c("-l", "--fields"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Comma-separated list of output fields to retain (default: all)'
  ),
  make_option(
    c("-m", "--mito"),
    action = "store_true",
    default = FALSE,
    type = 'character',
    help = 'Mark mitochondrial elements with reference to chromsomes and biotypes'
  ),
  make_option(
    c("-n", "--mito-chr"),
    action = "store",
    default = 'mt,mitochondrion_genome,mito,m,chrM,chrMt',
    type = 'character',
    help = 'If specified, marks in a column called "mito" features on the specified chromosomes (case insensitive)'
  ),
  make_option(
    c("-p", "--mito-biotypes"),
    action = "store",
    default = 'mt_trna,mt_rrna,mt_trna_pseudogene',
    type = 'character',
    help = 'If specified,  marks in a column called "mito" features with the specified biotypes (case insensitve)'
  ),
  make_option(
    c("-c", "--parse-cdnas"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Provide a cDNA file for extracting meta info and/or filtering.'
  ),
  make_option(
    c("-y", "--parse-cdna-names"),
    action = "store_true",
    default = NULL,
    type = 'character',
    help = 'Where --parse-cdnas is specified, parse out info from the Fasta name. Will likely only work for Ensembl GTFs'
  ),
  make_option(
    c("-d", "--parse-cdna-field"),
    action = "store",
    default = 'transcript_id',
    type = 'character',
    help = 'Where --parse-cdnas is specified, what field should be used to compare to identfiers from the FASTA?'
  ),
  make_option(
    c("-i", "--fill-empty"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Where --fields is specified, fill empty specified columns with the content of the specified field. Useful when you need to guarantee a value, for example a gene ID for a transcript/gene mapping. '
  ),
  make_option(
    c("-e", "--filter-cdnas-output"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Where --parse-cdnas is specified, filter sequences and output to the specified file. No file will be output if this is not specified (for example for use of --dummy-from-cdnas only).'
  ),
  make_option(
    c("-u", "--version-transcripts"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = 'Where the GTF contains transcript versions, should these be appended to transcript identifiers? Useful when generating transcript/gene mappings for use with transcriptomes. NOTE: if --filter-cdnas-field is set, the setting of this field is set to best match transcript identifiers in the cDNAs FASTA.'
  ),
  make_option(
    c("-o", "--output-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file path'
  )
)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

if (is.na(opt$gtf_file)){
  die('ERROR: No input GTF file specified')
}

if (is.na(opt$output_file)){
  die('ERROR: No output file specified')
}

# Try to get some annotation fields from the the transcript names themselves-
# this will likely not work for non-Ensembl GTFs

parse_ensembl_fasta_transcript_info <- function(tname){
  description <- sub(".*description:(.*)", "\\1", tname)
  tname <- sub(" description:.*", '', tname)
  tsplit <- unlist(strsplit(tname, ' '))
  names(tsplit) <- lapply(tsplit, function(x){
    if (grepl(':', x[1])){
      cnames <- unique(unlist(lapply(strsplit(x, ':'), function(y) y[1])))
      if (length(cnames) == 1){
        return(cnames)
      }else{
        return(NULL)
      }
    }
  })
  names(tsplit)[1] <- 'transcript_id'
  tsplit <- tsplit[names(tsplit) != 'NULL']
  tsplit <- sub('.*:(.*)', "\\1", tsplit)
  tsplit['description'] <- description
  names(tsplit) <- sub('^gene$', 'gene_id', names(tsplit))
  names(tsplit) <- sub('^gene_symbol$', 'gene_name', names(tsplit))
  tsplit['gene_id'] <- sub('\\.[0-9]+', '', tsplit['gene_id'])
  
  data.frame(t(tsplit), stringsAsFactors = FALSE)
}

# Import the GTF

print(paste('Reading', opt$gtf_file, 'elements of type', opt$feature_type))
gtf <- import(opt$gtf_file, feature.type = opt$feature_type )

# Combine basic info (chromosomes, coordinates) with annotation found in GTF attributes

anno <- cbind(chromosome = seqnames(gtf), as.data.frame(ranges(gtf)), elementMetadata(gtf))
print(paste('Found', nrow(anno), 'features'))

# Mark mitochondrial features

if (opt$mito){
  anno$mito <- ucfirst(as.character(tolower(anno$gene_biotype) %in% cleanlist(opt$mito_biotypes) | tolower(anno$chromosome) %in% cleanlist(opt$mito_chr)))
}

# If specified, put the desired field first

if (! is.na(opt$first_field)){
  if (! opt$first_field %in% colnames(anno)){
    die(paste(first_field, 'is not a valid field'))
  }
  anno <- anno[,c(opt$first_field, colnames(anno)[colnames(anno) != opt$first_field])]
}

# If we're filtering a list of cDNAs, find their IDs while considering the
# transcript ID versioning

cdna_transcript_names <- NULL
if (! is.null(opt$parse_cdnas)){
  
  suppressPackageStartupMessages(require(Biostrings))
  cdna <- readDNAStringSet(opt$parse_cdnas)
  cdna_transcript_names <- unlist(lapply(names(cdna), function(x) unlist(strsplit(x, ' '))[1]  ))
}  

# Version transcripts

if ( opt$feature_type == 'transcript' &&  all(c('transcript_id', 'transcript_version') %in% colnames(anno) )){
  has_transcript_version <- ! is.na(anno$transcript_version)
  
  versioned_transcripts <- anno$transcript_id

  if ( any(has_transcript_version) ){
    versioned_transcripts[has_transcript_version] <- paste(anno$transcript_id[has_transcript_version], anno$transcript_version[has_transcript_version], sep='.')
  }

  # Force transcript versioning if that's what's required to match the cDNAs

  if ( (! is.null(cdna_transcript_names)) && length(intersect(cdna_transcript_names, versioned_transcripts)) > length(intersect(cdna_transcript_names, anno$transcript_id)) ){
    opt$version_transcripts <- TRUE
    print('Versioning transcripts')
  } else{
    opt$version_transcripts <- FALSE
    print('Not versioning transcripts (they do not look versioned)')
  }

  if ( opt$version_transcripts ){
    anno$transcript_id <- versioned_transcripts
  } 
}

# If specified, filter down a provided cDNA FASTA file

if (! is.null(opt$parse_cdnas)){
  
  filtered_cdna <- cdna[which(cdna_transcript_names %in% anno[[opt$parse_cdna_field]])]
  cdna_only =  cdna[which(! cdna_transcript_names %in% anno[[opt$parse_cdna_field]])]
  
  if (! is.null(opt$filter_cdnas_output)){
    print(paste("Filtering", opt$parse_cdnas, "to match the GTF"))
    
    # Filter out cDNAs without matching transcript entries in the GTF
    
    if (! any(cdna_transcript_names %in% anno[[opt$parse_cdna_field]])){
      die(paste("ERROR: None of the input sequences have matching", opt$parse_cdna_field, 'values in the GTF file'))
    }
    
    print(paste('Storing filtered sequences to', opt$filter_cdnas_output))
    writeXStringSet(x = filtered_cdna, filepath = opt$filter_cdnas_output, compress = 'gzip')
  }
  
  # If we have transcripts with no annotation, see if we can get it from the fasta names

  if (length(cdna_only) > 0 && ! is.null(opt$parse_cdna_names)){
    print(paste("Info missing from GTF for", length(cdna_only), "supplied cDNAs, trying to extract it from the cDNA FASTA headers"))
    tinfo <- plyr::rbind.fill(lapply(names(cdna_only), parse_ensembl_fasta_transcript_info))
    anno <- plyr::rbind.fill(as.data.frame(anno), tinfo[,colnames(tinfo) %in% colnames(anno)])
  }
}

# If specified, subset to desired fields

if (! is.null(opt$fields) && opt$fields != ''){
  fields <- unlist(strsplit(opt$fields, ','))
  if (any(! fields %in% colnames(anno))){
    die(paste('ERROR:', fields, 'contains invalid field(s)'))
  }
  anno <- anno[,fields, drop = FALSE]
  if (! is.null(opt$fill_empty)){
    for (f in fields){
      empty_vals <- is.na(anno[[f]])
      anno[[f]][empty_vals] <- anno[[opt$fill_empty]][empty_vals]
    } 
  }
  anno <- anno[apply(anno, 1, function(x) all(! is.na(x))), ]
}

print(paste('Storing output to', opt$output_file))
write.table(anno, file = opt$output_file, sep = "\t", quote=FALSE, row.names = FALSE, col.names = opt$no_header)
