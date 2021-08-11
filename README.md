# hiv-haystack

## About
hiv-haystack is a Python tool derived and adapted from [epiVIA](https://github.com/wangwl/epiVIA) (by Wenliang Wang from the Vahedi lab at the University of Pennsylvania; [published in PNAS](https://www.pnas.org/content/117/10/5442)), sharing many of the same pipeline steps. The purpose of hiv-haystack is to:
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
More to come about this...
