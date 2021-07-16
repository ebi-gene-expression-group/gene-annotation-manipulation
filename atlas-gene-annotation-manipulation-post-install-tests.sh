#!/usr/bin/env bats

# Extract the test data
setup() {
    test_cdna_uri='http://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz'
    test_gtf_uri='http://ftp.ensembl.org/pub/release-104/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.104.gtf.gz'

    test_dir="post_install_tests"
    data_dir="${test_dir}/data"
    output_dir="${test_dir}/outputs"

    test_cdna=${data_dir}/$(basename $test_cdna_uri)
    test_gtf=${data_dir}/$(basename $test_gtf_uri)
    test_gene_to_remove=WBGene00000001
    part_test_gtf=${data_dir}/part.gtf

    gene_anno=${output_dir}/gene_anno.txt
    enriched_gene_anno=${output_dir}/enriched_gene_anno.txt
    gene_id_to_symbol=${output_dir}/gene_id_to_symbol.txt
    t2gene=${output_dir}/t2gene.txt
    t2gene_part=${output_dir}/t2gene_part.txt
    t2gene_part_matched=${output_dir}/t2gene_part_matched.txt
    filtered_cdnas=${output_dir}/filtered.fa.gz

    if [ ! -d "$data_dir" ]; then
        mkdir -p $data_dir
    fi

    if [ ! -d "$output_dir" ]; then
        mkdir -p $output_dir
    fi
}

@test "Download test data from FTP" {
    if  [ "$resume" = 'true' ] && [ -f "$test_gtf" ]; then
        skip "$test_gtf exists"
    fi

    run rm -f $test_gtf && rm -f $test_cdna && wget -P $test_dir/data $test_cdna_uri && wget -P $test_dir/data $test_gtf_uri && gunzip -c $test_gtf | grep -v $test_gene_to_remove > $part_test_gtf

    [ "$status" -eq 0 ]
    [ -f "$test_gtf" ]
}

@test "Make a table of all gene annotation in a GTF file" {
    if  [ "$resume" = 'true' ] && [ -f "$gene_anno" ]; then
        skip "$gene_anno exists"
    fi

    run rm -rf $gene_anno && gtf2featureAnnotation.R --gtf-file $test_gtf --feature-type "gene" --first-field "gene_id" --output-file $gene_anno

    [ "$status" -eq 0 ]
    [ -f "$gene_anno" ]
}

@test "Make a gene ID/ gene symbol mapping from a GTF file" {
    if  [ "$resume" = 'true' ] && [ -f "$gene_id_to_symbol" ]; then
        skip "$gene_id_to_symbol exists"
    fi

    run rm -rf $gene_id_to_symbol && gtf2featureAnnotation.R --gtf-file $test_gtf --feature-type "gene" --first-field "gene_id" --output-file $gene_id_to_symbol --first-field "gene_id" --fields "gene_id,gene_name"

    [ "$status" -eq 0 ]
    [ -f "$gene_id_to_symbol" ]
}

@test "Make a transcript to gene file from GTF only" {
    if  [ "$resume" = 'true' ] && [ -f "$t2gene" ]; then
        skip "$t2gene exists"
    fi

    run rm -rf $t2gene && gtf2featureAnnotation.R --gtf-file $test_gtf --version-transcripts --feature-type "transcript" --first-field "transcript_id" --output-file $t2gene --fields "transcript_id,gene_id" --no-header

    [ "$status" -eq 0 ]
    [ -f "$t2gene" ]
}

@test "Make a transcript to gene file (using transcriptome, some missing GTF lines)" {
    if  [ "$resume" = 'true' ] && [ -f "$t2gene_part" ]; then
        skip "$t2gene_part exists"
    fi

    run rm -rf $t2gene_part && gtf2featureAnnotation.R --gtf-file $part_test_gtf --version-transcripts --parse-cdnas $test_cdna  --parse-cdna-field "transcript_id" --feature-type "transcript" --parse-cdna-names --fill-empty transcript_id --first-field "transcript_id" --output-file $t2gene_part --fields "transcript_id,gene_id" --no-header

    [ "$status" -eq 0 ]
    [ -f "$t2gene_part" ]
}

@test "Check agreement of GTF and Fasta-derived gene IDs" {
    run diff <(cat $t2gene | sort) <(cat $t2gene_part | sort)

    [ "$status" -eq 0 ]
}

@test "Make a transcript to gene file and filter cDNAs to match" {
    if  [ "$resume" = 'true' ] && [ -f "$filtered_cdnas" ]; then
        skip "$filtered_cdnas exists"
    fi

    run rm -rf $filtered_cdnas && gtf2featureAnnotation.R --gtf-file $part_test_gtf --version-transcripts --parse-cdnas $test_cdna  --parse-cdna-field "transcript_id" --feature-type "transcript" --first-field "transcript_id" --output-file $t2gene_part_matched --fields "transcript_id,gene_id" --no-header --filter-cdnas-output $filtered_cdnas

    [ "$status" -eq 0 ]
    [ -f "$filtered_cdnas" ]
}

@test "Make sure transcripts were successfully filtered" {
    run eval "zcat $filtered_cdnas | grep $test_gene_to_remove"  
    [ "$status" -eq 1 ]
}

@test "Make a gene-level annotation file and include missing gene info from cDNA FASTA headers" {
    if  [ "$resume" = 'true' ] && [ -f "$enriched_gene_anno" ]; then
        skip "$enriched_gene_anno  exists"
    fi

    run rm -rf $enriched_gene_anno && gtf2featureAnnotation.R --gtf-file $part_test_gtf --version-transcripts --parse-cdnas $test_cdna  --parse-cdna-field "gene_id" --feature-type "gene" --parse-cdna-names --first-field "gene_id" --output-file $enriched_gene_anno

    [ "$status" -eq 0 ]
    [ -f "$enriched_gene_anno" ]
}
