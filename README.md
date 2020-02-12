# Gene annotation manipulation

This package currently holds one script: gtf2featureAnnotation.R (though it may contain more in future). gtf2featureAnnotation.R takes a GTF annotation file and creates a table of annotation by feature, optionally filtering a cDNA file supplied at the same time. For full usage instructions run:

```
gtf2featureAnnotation.R --help
```

```
Usage: gtf2featureAnnotation.R [options]


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

    -c FILTER-CDNAS, --filter-cdnas=FILTER-CDNAS
        If specified, sequences in the provided FASTA-format cDNAs file will be filtered to remove entries not present in the annotation

    -d FILTER-CDNAS-FIELD, --filter-cdnas-field=FILTER-CDNAS-FIELD
        Where --filter-cdnas is specified, what field should be used to compare to identfiers from the FASTA?

    -e FILTER-CDNAS-OUTPUT, --filter-cdnas-output=FILTER-CDNAS-OUTPUT
        Where --filter-cdnas is specified, what file should the filtered sequences be output to?

    -u, --version-transcripts
        Where the GTF contains transcript versions, should these be appended to transcript identifiers? Useful when generating transcript/gene mappings for use with transcriptomes. NOTE: if --filter-cdnas-field is set, the setting of this field is set to best match transcript identifiers in the cDNAs FASTA.

    -o OUTPUT-FILE, --output-file=OUTPUT-FILE
        Output file path

    -h, --help
        Show this help message and exit
```

## Example commands

This script should be able to produce tabular text files compiling any of the information present in a GTF file, e.g. for passing to downstream analysis.

### Make a trasnscript-to-gene ID mapping for use with e.g. Kallisto, Salmon etc:

This is a common thing to need: get a mapping of transcript ID to gene ID. So we tell the script we want to extract annotations for 'transcript' feature types, and just have transcript and gene IDs in the output table.

```
gtf2featureAnnotation.R \
    --gtf-file <input GTF> 
    --no-header \
    --version-transcripts \ 
    --feature-type "transcript" \
    --fields "transcript_id,gene_id" \
    --first-field "transcript_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE>
```

Parameter meanings:

 - `--no-header`: Do not apply a header to the output file
 - `--version-transcripts`: A flag which specifies that the transcript versions commonly found in Ensembl GTF files, if present, should be appended to transcript identifiers.
 - `--fields`: Comma-separated list of attributes to extract
 - `--first-field`: Which field to place first in output? 


### Synchronise a cDNA file with the annotation

Sometimes it's important that there are no transcripts in a FASTA-format transcriptome that cannot be matched to a transcript/gene mapping. Alevin, for example, produces (or at least used to produce) errors when this mismatch was present. We can synchronise the cDNA file, removing mismatches:

```
gtf2featureAnnotation.R \
    --gtf-file <input GTF> 
    --no-header \
    --version-transcripts \ 
    --feature-type "transcript" \
    --fields "transcript_id,gene_id" \
    --first-field "transcript_id" \
    --output-file <OUTPUT TEXT FILE FOR TABLE> \
    --filter-cdnas reference.fa.gz \
    --filter-cdnas-field "transcript_id" \
    --filter-cdnas-output <OUTPUT FASTA FILE> 
```

Output will be a an annotation table, and a FASTA-format cDNAs file with unannotated transcripts removed.

The additional fields here specify:

 - `--filer-cdnas`: FASTA-format cDNAs file to filter for matches with the annotation
 - `--filter-cdnas-field`: which field in the annotation should match the transript identifiers

When running in this way, the annotation talbe that's output needs to match the cDNAs file. So the `--version-transcripts` setting is overriden internally using whichever of versioned or unversioned identifiers match the transcript identifiers. We have found Ensembl to be a little inconsistent in this.
