# hiv-haystack

## About
hiv-haystack is a Python tool derived and adapted from [epiVIA](https://github.com/wangwl/epiVIA) (by Wenliang Wang from the Vahedi Lab at the University of Pennsylvania; [published in PNAS](https://www.pnas.org/content/117/10/5442)), sharing many of the same pipeline steps. The purpose of hiv-haystack is to:
1. parse scATACseq BAM files (from cellranger-atac which have been aligned to both human and HIV genomes)
2. identify proviral reads (both mates are in the provirus)
3. identify host chimeric reads (both mates are in human, but one mate has a soft clip containing viral LTR)
4. identify viral chimeric reads (both mates are in virus, but one mate has a soft clip containing human DNA)
5. identify umapped reads which are a dual chimera (one mate in human and one mate in virus)

I highly recommend visiting epiVIA's repo and paper to see their diagram regarding this pipeline.

While hiv-haystack is built upon epiVIA, the main differences between hiv-haystack and epiVIA are as follows (items below indicate functionalites added to hiv-haystack):
- Allows for BAM files that were aligned to more than one provirus reference genome, which is useful for samples derived from persons living with HIV (espeically with antiretroviral treatment). This is due to the high mutational rate which misses out on viral reads if using just the default HXB2 HIV reference genome in combination with cellranger-atac's bwa mem aligner.
- Has more stringent checking of host reads with potential viral chimeras to reduce chances of false positives.
  - In the case of a potential viral soft clip not mapping continously to the end of LTR (either 5' or 3' in either direction), hiv-haystack will check the adjacent host base pairs to see if those reads could align to both virus and host. If there is no alignment to virus, this read is discarded and not considered a read with a valid chimera since the viral sequence did not extent to the end of the LTR as would be needed for a correct HIV-host genome junction site.


## Getting started

There are two ways to run hiv-haystack depending on the number of viral sequences that were included in the chimeric reference genome used during cellranger-atac alignment.

### (A) If there is only one viral sequence:
```bash
python main.py \
  --bamfile=namesorted.bam \
  --outputDir=run0Output \
  --viralFasta=suma_TF1.fasta \
  --LTRpositions="1,633,9094,9726" \
  --hostGenomeIndex=refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa
```

### (B) If there are more than one viral sequence:
```bash
# build BLAST db
makeblastdb -in viralSequences.fasta -dbtype nucl -out viralSeqs

# run blastn. hxb2_ltr5.fa is the 5' LTR sequence that we'll use to find any matches with our viral sequences
blastn -db viralSeqs -query hxb2_ltr5.fa -out viralSeqs_possible.table -outfmt 6

# optional. If your viral sequences don't include HXB2 but 
# they were used in the reference genome during cellranger-atac
# alignment, you will need to manually add HXB2 LTR sequences to
# the output table from the above command.
#
# hxb2    chrHXB2    -1    -1    -1    -1    1    634    1    634    0.0    -1
# hxb2    chrHXB2    -1    -1    -1    -1    1    634    9086    9719    0.0    -1

# run hiv-haystack
python main.py \
  --bamfile=namesorted.bam \
  --outputDir=run1Output \
  --viralFasta=viralSequences_withHXB2.fa \
  --LTRmatches=viralSeqs_possible.table \
  --hostGenomeIndex=refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa

```

## Parameters

- `--bamfile` *(required)* Namesorted BAM file from cellranger-atac. Note that the default output from cellranger-atac is position sorted. You will need to run name sorting via samtools.
- `--outputDir` *(required)* Directory for output files.
- `--viralFasta` *(required)* Viral fasta file of all (and only) viral sequences that were part of the reference chimeric genome used for the initial alignment with cellranger-atac. Can have multiple sequences in the file.
- `--hostGenomeIndex` *(required)* Prefix of bwa indexed host reference genome (NO provirus sequences included). Can use the 10X Genomics cellranger-atac reference genome which should be bwa indexed.
- `--topNReads` Integer value for the first *n* reads from the BAM file. Default value is all reads (n = -1).
- `--LTRmatches` blastn table output format for LTR matches to HXB2 LTR. This is required when running with multiple autologous sequences (i.e. if there are multiple fasta sequences in the file associated with the `--viralFasta` argument.
- `--LTRpositions` LTR positions when running with only one viral sequence (i.e. only one fasta sequence in the file associated with the `--viralFasta` arugment). LTR positions should be provided as 1-indexed positions: 5' start, 5' end, 3' start, 3'end (example: 1,634,9086,9719)
- `--LTRClipLen` Number of basepairs to extend into LTR from a chimeric fragment. The default value is 11 as used by epiVIA.
- `--hostClipLen` Number of basepairs to extend into the host genome from a chimeric fragment. The default value is 17 as used by epiVIA.

## Outputs
- `proviralReads.bam`: all reads from namesorted bam with both mates aligning to the viral sequences.
- `hostWithPotentialChimera.bam`: all reads from namesorted bam where soft clip present in host genome read that passes the requirements set in the arguments.
- `unmappedWithPotentialChimera.bam`: all reads from unmapped reads where soft clip is present.
- `hostWithValidChimera.bam`: all reads from `hostWithPotentialChimera.bam` where the soft clip has a confirmed alignemnt to a LTR region.
- `viralReadHostClipFasta.fa`: soft clip sequences from viral aligned reads that need to be chcked for alignemnt to host genome.
- `unmappedHostClipFasta.fa`: soft clip sequences from unmapped reads that need to be chcked for alignemnt to host genome.
- `integrationSites.tsv`: tsv file of valid integration sites.

  | cbc | chr | orient | pos |
  |---|---|---|---|
  | cellbarcode | chromosome | orientation (+ or -) | basepair position |

- `integrationSites_viralFrags.tsv`: tsv file of viral fragments associated with the valid integration sites. See description for output format (same as viralFrags.tsv excpet for the the `alreadyRecordedInIntegration` column)
- `viralFrags.tsv`: tsv file of all viral fragments. See below for output format.

  | cbc | seqname | startBp | endBp | readname | usingAlt | confirmedAlt | alreadyRecordedInIntegration |
  |---|---|---|---|---|---|---|---|
  | cell barcode | name of viral sequence | start basepair (0-index) | end basepair (0-index; inclusive) | readname of BAM record | alternative alignments | confirmed alternative | this read is already recorded in `integrationSites_viralFrags.tsv` |