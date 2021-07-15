# Gene annotation manipulation
[![install with Conda](https://img.shields.io/conda/v/ebi-gene-expression-group/atlas-gene-annotation-manipulation)](https://anaconda.org/ebi-gene-expression-group/atlas-gene-annotation-manipulation)

This package currently holds one script: gtf2featureAnnotation.R (though it may contain more in future). gtf2featureAnnotation.R takes a GTF annotation file and creates a table of annotation by feature, optionally filtering a cDNA file supplied at the same time. 

## Install

The recommended method for script installation is via the Conda package:

```
conda install -c ebi-gene-expression-group atlas-gene-annotation-manipulation
```

## Usage

For full usage instructions run:


```
gtf2featureAnnotation.R --help
```

```
Usage: ./gtf2featureAnnotation.R [options]


Options:
    -g GTF-FILE, --gtf-file=GTF-FILE
        Path to a valid GTF file

    -t FEATURE-TYPE, --feature-type=FEATURE-TYPE
        Feature type to use (default: gene)

    -f FIRST-FIELD, --first-field=FIRST-FIELD
        Field to place first in output table (default: gene_id)

    -r, --no-header
        Suppress header on output

    -l FIELDS, --fields=FIELDS
        Comma-separated list of output fields to retain (default: all)

    -m, --mito
        Mark mitochondrial elements with reference to chromsomes and biotypes

    -n MITO-CHR, --mito-chr=MITO-CHR
        If specified, marks in a column called "mito" features on the specified chromosomes (case insensitive)

    -p MITO-BIOTYPES, --mito-biotypes=MITO-BIOTYPES
        If specified,  marks in a column called "mito" features with the specified biotypes (case insensitve)

    -c PARSE-CDNAS, --parse-cdnas=PARSE-CDNAS
        Provide a cDNA file for extracting meta info and/or filtering.

    -y, --parse-cdna-names
        Where --parse-cdnas is specified, parse out info from the Fasta name. Will likely only work for Ensembl GTFs

    -d PARSE-CDNA-FIELD, --parse-cdna-field=PARSE-CDNA-FIELD
        Where --parse-cdnas is specified, what field should be used to compare to identfiers from the FASTA?

    -i FILL-EMPTY, --fill-empty=FILL-EMPTY
        Where --fields is specified, fill empty specified columns with the content of the specified field. Useful when you need to guarantee a value, for example a gene ID for a transcript/gene mapping. 

    -e FILTER-CDNAS-OUTPUT, --filter-cdnas-output=FILTER-CDNAS-OUTPUT
        Where --parse-cdnas is specified, filter sequences and output to the specified file. No file will be output if this is not specified (for example for use of --dummy-from-cdnas only).

    -u, --version-transcripts
        Where the GTF contains transcript versions, should these be appended to transcript identifiers? Useful when generating transcript/gene mappings for use with transcriptomes. NOTE: if --filter-cdnas-field is set, the setting of this field is set to best match transcript identifiers in the cDNAs FASTA.

    -o OUTPUT-FILE, --output-file=OUTPUT-FILE
        Output file path

    -h, --help
        Show this help message and exit
```

## Example commands

This script should be able to produce tabular text files compiling any of the information present in a GTF file, e.g. for passing to downstream analysis.

### Make a table of all gene annotations in a GTF file

```
gtf2featureAnnotation.R \
    --gtf-file <INPUT GTF> \
    --feature-type "gene" \
    --first-field "gene_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE>
```

### Make a gene ID/ gene symbol mapping

```
gtf2featureAnnotation.R \
    --gtf-file <INPUT GTF> \
    --feature-type "gene" \
    --fields "gene_id,gene_name"
    --first-field "gene_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE> \
```

### Make a transcript-to-gene ID mapping for use with e.g. Kallisto, Salmon etc:

This is a common thing to need: get a mapping of transcript ID to gene ID. So we tell the script we want to extract annotations for 'transcript' feature types, and just have transcript and gene IDs in the output table.

```
gtf2featureAnnotation.R \
    --gtf-file <INPUT GTF> \
    --version-transcripts \
    --feature-type "transcript" \
    --first-field "transcript_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE> \
    --fields "transcript_id,gene_id" \
    --no-header
```

### Synchronise a cDNA file with the annotation

Sometimes it's important that there are no transcripts in a FASTA-format transcriptome that cannot be matched to a transcript/gene mapping. Alevin, for example, produces (or at least used to produce) errors when this mismatch was present. We can synchronise the cDNA file, removing mismatches:

```
gtf2featureAnnotation.R \
    --gtf-file <INPUT GTF> \
    --version-transcripts \
    --parse-cdnas <REFERENCE CDNA FASTA> \
    --parse-cdna-field "transcript_id" \
    --feature-type "transcript" \
    --first-field "transcript_id" 
    --output-file <OUTPUT TEXT FILE FOR TABLE> \
    --fields "transcript_id,gene_id" \
    --no-header \
    --filter-cdnas-output <FILTERED CDNA OUTPUT FILE>
```

Output will be a an annotation table, and a FASTA-format cDNAs file with unannotated transcripts removed.

When running in this way, the annotation talbe that's output needs to match the cDNAs file. So the `--version-transcripts` setting is overriden internally using whichever of versioned or unversioned identifiers match the transcript identifiers. We have found Ensembl to be a little inconsistent in this.

### Synchronise an annotation file with the cDNA

Alternatively, you may want to do the converse to the above, and make the annotation file match the cDNA:

```
gtf2featureAnnotation.R \
    --gtf-file <INPUT GTF> \
    --version-transcripts \
    --parse-cdnas <REFERENCE CDNA FASTA> \
    --parse-cdna-field "transcript_id" \
    --feature-type "transcript" \
    --parse-cdna-names \
    --fill-empty transcript_id \
    --first-field "transcript_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE> \
    --fields "transcript_id,gene_id" \
    --no-header
```

There's no `--filter-cdnas-output` here so the Fasta won't be filtered. The --parse-cdna-names options will allow the Fasta headers to be searched for the annotation info- which is often present with Ensembl data. Failing that, the `--fill-empty` parameter will just copy the `transcript_id` to empty `gene_id` values. 
