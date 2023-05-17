# vcfByGene

## Table of contents

   * [About](#About)
      * [Filtering pipeline overview](#filtering-pipeline-overview)
      * [Resources](#resources)
   * [Running the vcfByGene pipeline](#running-the-vcfbygene-pipeline)
   * [Running VEP on the genetics cluster](#running-vep-on-the-genetics-cluster)
   * [Creating patient-variant mappings from annotated files](#creating-patient-variant-mappings-from-annotated-files)
      * [Logic overview](#logic-overview)
      * [Running the script](#running-the-script)
   * [Contributors](#contributors)
   * [Copyrights](#copyrights)
   * [How to provide feedback](#how-to-provide-feedback)
   * [Questions](#questions)

## About

This repository contains code necessary to filter WES data from the UK BioBank for further analysis. It requires access to the UKB Research Analysis Platform (RAP) on DNAnexus through a valid UK BioBank Application ID and linked RAP project.

### Filtering pipeline overview

Exon coordinates were determined for genes of interest using MANE transcripts, with an additional 5nt retained up- and downstream of each coding region to capture splice acceptor and donor region variants. Gene-level VCF files were extracted from the UK Biobank WES joint-called pVCFs using bcftools. The VCF files were then normalized to flatten multiallelic sites and align variants to the GRCh38 reference genome. Variants located in NIST Genome in a Bottle “difficult regions” were removed from analysis, as were variants with an alternate allele frequency greater than 0.1% in the UK Biobank cohort. Further filtering removed variants where more than 10% of samples were missing genotype calls (Szustakowski et al. 2021) and variants that did not appear in the UK Biobank cohort. To mitigate differences in sequencing coverage between individuals who were sampled at different phases of the UK Biobank project, variants were only retained in the final set if at least 90% of their called genotypes had a read depth of at least 10 (UK Biobank Whole Exome Sequencing 300k Release: Analysis Best Practices). All filtering took place by running bcftools through the Swiss Army Knife app on the UK Biobank Research Analysis Platform.

### Resources

Szustakowski, J. D. et al. [Advancing human genetics research and drug discovery through exome sequencing of the UK Biobank](https://www.nature.com/articles/s41588-021-00885-0). Nat Genet 53, 942–948 (2021). 

[UK Biobank Whole Exome Sequencing 300k Release: Analysis Best Practices](https://biobank.ndph.ox.ac.uk/crystal/ukb/docs/UKB_WES_AnalysisBestPractices.pdf)

## Running the vcfByGene Pipeline

1. To begin, ensure you are in a conda environment with Jupyter lab installed (see instructions in "Laptop Setup" file).
2. Clone this repository from GitHub
3. From your conda environment, launch a Jupyter lab session by running `jupyter lab`
    1. This should launch in your browser. If not, just navigate directly to the session URL. It should usually just be http://localhost:8888/lab (but sometimes the port number will be different).
4. Locate the notebook file `vcfByGene.ipynb`. Once you're here the process is fairly straightforward. There's just a few places you'll have to make changes:
    1. First, try loading the libraries. You might need to install a few within your conda environment the first time you do this.
    2. In the code cell below the heading “Set Parameters,” you'll need to set parameters. (haha, who would've guessed?) The ones that change most often are:
        1. `tag_str`: this is a tag that will be attached to all the UKB RAP jobs you run in this notebook. This allows you to search for this set of jobs after the fact, to diagnose problems or to go into the logs and remember the parameters you used for them. Also super useful if you mess up somehow and need to batch terminate the jobs! You can find them by the tag instead of terminating one-by-one or terminating all running jobs and accidentally killing one that wasn't part of the batch.
        2. `dx_vcf_out_path`: the folder on UKB in which you want to store the outputs. If creating both VCF and carrier files, you'll end up with two different out paths, just give them different names.
    3. Next, specify the genes you want to extract and filter. Use the code cell below “List of gene symbols as input.” You have a couple options here:
        1. If your list is short, just type it out as a list of strings.
        2. If you have a lot of genes, it's probably easier to read them in from a text file.
    4. The bulk of the work occurs in the section “Run swiss-army-knife on DNAnexus.”
        1. The `known_large_genes` list is something I've created to preempt memory issues with genes that I know from experience are large and require more job memory. Feel free to add to it as needed.
        2. Inspect the command strings in the cell and adjust them as needed. This cell as I've formatted it actually executes two commands per gene. Both commands follow the same filtering steps. In one, the individual genotypes are dropped, yielding a vcf of UKB variants. In the other, the variants are retained along with a pipe-separated list of their UKB carriers. In the filtering step, you can change things like max-allele filters or allele frequency filters. In the formatting step, you can change how the carrier output looks or add or remove columns (eg allele frequency) as needed.
        3. After running this cell (it may take a while if you have lots of genes), double-check and note which genes were not found in MANE. Follow up with whoever asked you to generate these files to figure out if finding data for these genes is critical or not.
5. Once all of the jobs have finished, I like to be safe and search through all their logs to ensure that non of them ran into any memory issues.
    1. Run bash script `search_job_logs.sh`. This may take a while depending on how many jobs you have to search through.
    2. When the script is finished running, inspect the generated text file. Re-run the original jobs for any genes that appear in this file.
        1. A nice way to do this would be to load this file in and create that genes list in the `vcfByGene` notebook. Then add these genes to the `known_large_genes` list and re-run the swiss-army-knife cell. This will re-run these jobs on higher-memory servers. Once they're done, you can be extra safe by re-running the search_job_logs.sh script one more time to make sure the jobs were all successful.

## Running VEP on the genetics cluster

Once all your filtered vcf files are ready, you can download them to the cluster to annotate with VEP. **This is definitely a spot for future improvement — ideally, we want to do all of these steps directly on the RAP.** For now, you can follow the following steps:

1. Bulk-download the filtered vcf files (**with patient genotypes removed**) from the RAP. Assuming you are logged in to the dx toolkit, run the following on a local terminal: `dx download -r vcf_out_directory`
2. Use CyberDuck to upload this downloaded folder to a destination of your choosing on the BWH genetics cluster. In general, we've been storing files within `/net/ukbb`.
3. Log into the BWH genetics cluster through the command line. You can then modify and run the wrapper script stored at `/net/ukbb/500k_processing_lara/parse_vcf_parallel.sh` (I suggest doing any modifications needed by opening the files through CyberDuck). To run this script and submit batch jobs, run `qsub -t 1-N /net/ukbb/500k_processing_lara/parse_vcf_parallel.sh`, where N is the number of gene vcf files to annotate. This will run VEP over all the files in your specified directory. It will yield a csv file of VEP-annotated variants for each gene. The two main components of this script are:
    1. Running VEP using the following command:
    ```
    vep --cache --force_overwrite --offline --dir /net/data/vep --assembly GRCh38 --plugin dbNSFP,/net/ukbb/data/dbnsfp_content/dbNSFP4.0a.gz,ALL --plugin LoF,loftee_path:/net/data/vep/loftee-grch38,human_ancestor_fa:/net/data/vep/loftee-grch38/human_ancestor.fa.gz --vcf --max_af --af_gnomade --af_gnomadg --no_check_variants_order --canonical --no_stats --quiet -I $VCF_IN_FILE -o $VEP_FILE
    ```
    2. Parsing the VEP output by running `/net/ukbb/500k_processing_lara/vep_parser_af_lof.py`. Note that this parser requires that the VEP annotations include (a) max allele frequencies, both the default from dbNSFP (`allele_frequency`), as well as `gnomADe_MAX_AF` and `gnomADg_MAX_AF`, and (b) LoF confidence scores, obtained from the `loftee` VEP plugin.
4. Download the directory of parsed VEP annotations to your local machine. Then upload them to the RAP: `dx upload -r path/to/directory --destination RAP/path/to/store`

##  Creating patient-variant mappings from annotated files

### Logic overview

For many of our analyses, we require csv/ssv files oriented by patient, as opposed to by variant. Further, we normally want only one variant per patient (per gene), and therefore, we wish to find the most severe variant that each UKB patient carries in each gene. We choose the most severe according to the following decision tree:

1. Choose “deleterious” variant if present (any of frameshift, splice acceptor, splice donor, or stop gained)
    1. If multiple deleterious, choose one in the following priority order, moving to next option if all deleterious variants present are equal at the given level:
        1. Variant with CADD score (most deleterious variants, however, don’t have CADD scores annotated)
        2. Variant with a high confidence (HC) LoF_confidence flag, from the loftee annotations
            1. *Note: the updated script now applies a filter to only retain HC variants prior to the ranked choice occurs. This step is therefore effectively skipped; only HC deleterious variants are retained from the outset.*
        3. Variant with lower allele frequency
            1. Selected for each variant as the maximum of max gnomAD exomes AF, max gnomAD genomes AF, and dbNSFP allele_frequency
        4. Random choice
2. If no deleterious variants, choose the missense variant with the highest CADD score
    1. If two or more variants have the same CADD score, choose the variant with the lower allele frequency (where AF is determined as described above).
    2. If all missense variants have the same allele frequency, choose one randomly.
3. If no missense variants, choose the synonymous variant with the lowest allele frequency.
    1. If all synonymous variants have the same allele frequency, choose one randomly.

### Running the script

To generate the patient-variant mapping, run the script stored on the RAP at `patient_variant_extracts/scripts/var-pt-map.py`. To do so, you will need to create a Docker image based on Python 3.6 that properly installs Pandas and Numpy. A Docker image with those libraries can be found at `patient_variant_extracts/scripts/pandas_numpy.tar.gz`. Then, run the following commands on a local terminal where Python SDK/dx-toolkit has been properly configured:

```
cmd="python3 var-pt-map.py gene_name -v blend_ukb/resources/parsed_bc_genes -s -e VEST4_score,CADD,REVEL_score,ukb_af"
dx run swiss-army-knife -iin="file-GQVXx8jJqBj0Vg55zP03kXj5" -icmd="$cmd" -iimage_file="patient_variant_extracts/scripts/pandas_numpy.tar.gz" --destination "patient_variant_extracts/your_targest_directory" -y
```

Ensure your target directory already exists. The file ID specified as the Swiss Army Knife input corresponds to the aforementioned `var-pt-map.py` script. Replace LDLR with your gene of choice, or run in a for-loop over multiple genes. The Python script is set up to accept the following flags:

| Flag | Full name | Description | Data type |
| :- | :- | :- | :-: |
| gene | n/a | Required; symbol of gene to access | str 
| -v | --var_path | Name of directory in which annotated variant file is stored. Assumes file name is {gene}_parsed.csv | str
| -c | --car_path | Name of directory in which variant-patient `.ssv` mapping is stored. Assumes file name is {gene}.ssv | str
| -s | --most_severe | Includes only most severe variant per patient when specified | store_true
| -e | --extra_cols | Comma delimited list of additional VEP columns to include in output. Use to include, eg, computational scores like VEST4 or CADD9 in the output file. ukb_af can be included if it was retained at initial generation of var-pt SSVs. | comma separated str (converted to list) |

The generated csv file(s) will have, at minimum, the following columns:

| Column name | Description | Data type |
| :- | :- | :-: |
| eid | ID of carrier in the UK Biobank, for your Application ID | object
| variant_id/most_severe_variant | ID of the variant described in the remainder of the row. Format Chrom-Pos-Ref-Alt | object
| Synonymous | 1 if variant is synonymous; 0 otherwise| int64
| Missense | 1 if variant is missense; 0 otherwise| int64
| Deleterious | 1 if variant is one of: `stop_gained, start_lost, splice_acceptor_variant, splice_donor_variant, frameshift_variant`; 0 otherwise| int64
| coding_position | relative position of base pair in coding sequence | float64
| AA_name | for non-synonymous variants, lists the amino acid substitution as RefLocationAlt | object

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START -->
| Contributions | Name |
| ----: | :---- |
| [💻](# "Code") [🤔](# "Ideas and Planning") | [Tian Yu](https://github.com/TYTYBU) |
| [💻](# "Code") [🤔](# "Ideas and Planning") [📖](# "Documentation")| [Lara Brown](https://github.com/larabbrown) |
| [💻](# "Code") [🤔](# "Ideas and Planning") | Vineel Bhat |
| [📆](# "Project Management") [🤔](# "Ideas and Planning") | [Chris Cassa](https://github.com/cassalab) |

<!-- ALL-CONTRIBUTORS-LIST:END -->

(For a key to the contribution emoji or more info on this format, check out [“All Contributors.”](https://allcontributors.org/docs/en/emoji-key))

## Copyrights

All code in this repo is licensed with a BSD 3-Clause License. Please see the [license file](LICENSE) for details. 

All written materials are licensed with a Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0). Please see [this license](https://creativecommons.org/licenses/by-sa/3.0/) for details.

## How to Provide Feedback

Questions, bug reports, and feature requests can be submitted to this repo's [issue queue](https://github.com/TYTYBU/vcfByGene/issues).

## Questions

Please reach out to any of the contributors above with any questions you may have.