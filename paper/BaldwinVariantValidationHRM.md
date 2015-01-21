Workflow for PoolSeq Variant Validation  by Optimised High Resolution Melting Assay Design
===================================================================================================

S. Baldwin<sup>1</sup>, S Thomson<sup>1</sup>, M Pither-Joyce<sup>1</sup>, K Wright<sup>1</sup>, B
Warren<sup>1</sup>, S Liu<sup>1</sup>, Z Dwight<sup>2</sup> , R. Lee<sup>3</sup>, R Macknight<sup>3</sup>, J
McCallum<sup>1 \#</sup>

<sup>1</sup> The New Zealand Institute for Plant & Food Research Limited, Private
Bag 4704, Christchurch, NEW ZEALAND

<sup>2</sup> Department of Pathology, School of Medicine, University of Utah,
Salt Lake City, UT 84132, USA.

<sup>3</sup> Biochemistry Department, University of Otago, Dunedin, NEW ZEALAND

\# Corresponding Author <john.mccallum@plantandfood.co.nz>

Abstract
--------

### Background

Improvements in sequencing technologies and services enable large-scale variant
discovery even for small genetics projects in non-model
systems. Scaling validation of selected variants is now the principal barrier to
applying sequence data in downstream genetic studies. It is desirable to identify SNP validation
strategies that can be performed using PCR instrumentation and computational resources 
readily in most laboratories.

### Results

We report a workflow for automated design of small-amplicon
high-resolution melting (HRM) assays for validation of SNPs revealed by
standard variant detection software used in next-generation sequencing
analysis. Publically-available tools for PCR primer pair design in Unix
command-line or Galaxy workflow environment were enhanced by addition of
melt prediction through the uMelt webservice. Melt prediction enabled
large scale screens of candidate primer sets to select those predicted
to reveal larger melt differences between reference and variant
homozygotes. This workflow was evaluated by validating variants from
bulked segregant RNA PoolSeq (BSR-seq) of F<sub>2</sub> progeny pools of a wide
onion cross segregating for bolting (precocious flowering). The
validation rates for HRM markers selected using optimised assay design
were 69%, yielding 77 resolvable polymorphic HRM markers. Using a
selective bin mapping strategy these markers validated the BSR-seq
approach with enrichment for chromosome 1 and *AcBlt1* associated
markers.

### Conclusions

This method provides an efficient approach for medium-scale validation
of variants from Pool-seq or other genetics approaches employing
second-generation sequencing. The large SNP set and markers developed in
this reproducible onion stock will be useful for marker-assisted
selection and comparative mapping for this important adaptive trait.

Keywords
---------

PCR, diagnostics, bioinformatics, onion, PoolSeq

Background
----------

Global and context-specific detection of genetic variants by sequencing
is now the standard approach even for non-model organisms. PoolSeq is
the term given to sequencing pools of individuals by next-generation
sequencing NGS to identify SNPs or estimate population genetic
parameters [1, 2]. Although most widely applied to genomic DNA pools,
PoolSeq can also be applied to RNA pools, which may be more practical
for large genomes or species lacking reference genomes. PoolSeq has been
used for bulked segregant analysis (BSA) for mapping quantitative trait
loci (QTL) and simultaneously revealing candidate and linked variants.
Bulked segregant RNA-Seq (BSR-seq) or bulked segregant transcriptome
analysis (BSTA) methods have been described for plants [3-5] and catfish
[6].

Liu et al. [3] used bulked segregant RNA-Seq (BSR-seq) to identify the
*gl3* gene in an F<sub>2</sub> maize population where the recessive mutant is
involved in the accumulation of epicuticular wax. Trick et al. [4] used
BSA and RNA-seq to identify and map SNPs to an interval containing the
*GPC-B1* locus in wheat. In catfish the BSR-seq technique was used to
identify genes involved in disease resistance and mapped using the
reference genome [6]. In sunflower, Livaja *et al*. [5] used BSTA by
pooling F<sub>2</sub> lines that were homozygous resistant or susceptible and
homozygous for markers flanking a known resistance associated interval
to generate bulks. These were then used for pooled transcriptome
sequencing and markers designed to resistance gene candidates mapped
back to the interval.

Although deep sequencing and bioinformatics can readily reveal and rank
candidate variants, SNP validation remains an expensive step which is a
major barrier to smaller laboratories. Commercial SNP platforms such as
iPlex (Sequenom ) and Kaspar (LGC) have been widely employed for this
purpose. Although they are less universal and scalable than such
technologies, PCR-based methods are appealing because all genetics
laboratories have access to PCR machines and most have real-time PCR
instruments with capability for high-resolution melting (HRM; [7]).
However, there is a lack of usable modular, scalable, open-source
componentry that can be used to enable bulk design of PCR-based assays
from standard NGS variant formats.

HRM of small amplicons is one of the simplest closed-tube genotyping
methodologies [8]. HRM markers are rapid to screen, high-throughput and
can be used for scanning different sequence polymorphisms, such as SNP,
SSRs and indels. HRM analysis has been used for many applications from
research-based gene scanning [9] to clinical diagnostics [10]. It is an
attractive method for genotyping in small labs because it does not
require specialised reagents and can be performed on most real-time PCR
instruments [7]. Although HRM is very efficient for detecting
heterozygotes, resolving different homozygous classes can be a
limitation. A series of web-based ‘uMelt’ tools ([Dwight *et al.*,
2011](#_ENREF_4)) have been developed to aid in the development of HRM
diagnostics, which provide prediction of melting behaviour for sets of
amplicons, as well as assay design for dbSNP variants. However these
tools have not been previously employed in larger-scale assay design
pipelines.

Onion and shallot (*Allium cepa* L.) are internationally important
vegetable crops that have been domesticated for over 5000 years [11].
Trancriptome-based genetic maps [12-14] are available but the huge
nuclear genome size (1.8 x 10<sup>10</sup>  bp) has so far prevented development
of a genome reference. Onion is an outcrossing diploid, maintained
traditionally as open-pollinated populations that are selected for
adaptive traits but maintain high levels of heterozygosity [15, 16].
Bienniality, where the plants form a bulb in the first year and flower
in the following year, is a key adaptive trait selected during the
breeding and domestication of onion. Key genes involved in regulation of
bulb and flower formation have been described [17-19] and we recently
reported evidence for QTL in a wide onion cross conditioning ‘bolting’
(precocious flowering; [18]). This work was conducted with a low density
map and did not reveal any candidate genes linked with the major QTL
*AcBlt1.*

In the current study, we conducted BSR-seq on pools of RNA from F<sub>2</sub>
progeny used in validation of the *AcBlt1* QTL in the previous study. To
validate leads we developed an improved HRM-based SNP validation
workflow that incorporates amplicon melt prediction [20] for optimal
primer pair selection. This strategy was found to be a practical and
efficient means to screen large numbers of variant leads from NGS with
minimal benchwork and has provided gene and genome regions leads for
future physiological and genetic studies.

Results and Discussion
----------------------

### RNA-Seq

The combined assembly of RNA-Seq from leaves of the doubled-haploid line
‘CUDH2150’ before and after transfer to long days with earlier leaf plus
shoot RNA-Seq [14] provided a new reference assembly (GenBank
GBGJ00000000.1) comprising 140764 contigs of mean length 762 bp and N50
of 1247 bp. ‘CUDH2150’ was not only a parent of the population but being
a doubled haploid provided an ideal homozygous reference transcriptome
assembly as paralogous sequences could be identified.

Sequencing of pooled RNA from bulbs of the 12 bulks of F<sub>2</sub> progeny
yielded in the range of 19-27 million 100bp paired-end reads per pool
(Table 1). Variant detection using Freebayes revealed 281,385 variants
after filtering for those variants with minor allele frequencies between
0.3 and 0.7. Replication, normalisation and stringent filtering (MAF and
FDR) was used to minimise the effects of unequal expression of
individuals within pools, allele specific expression and sequencing
errors on the allele frequency results and comparisons [21]. Replicated
comparisons of normalised allele counts between the three pool types by
a Cochran–Mantel–Haenszel (CMH) test using Popoolation2 revealed that
1443 contigs (1% of total) exhibited significant allele frequency
differences in one or more comparisons (Figure 1). These contigs
contained 13,096 filtered variants (Supplementary File
filtered\_sign\_contigs.vcf.gz). Preliminary screening focussed on the
263 contigs which were significant in all three comparison types (Figure
1).

### Tool and Workflow Development

We previously reported that use of small amplicon HRM for SNP marker
development in this population was relatively inefficient, due to
difficulties in reliably discriminating homozygotes [14]. To improve
outcomes we modified the scripts used for designing primer pairs to
flank features of interest to provide automated estimation of melting
temperature in reference and alternate alleles through the uMelt [20]
web service provided at University of Utah. We employed a criterion of a
predicted melting temperature difference of 0.6 °C or higher during
selection of primer sets for validation assays. These updated design
scripts, as well as additional components for marker design from
standard vcf files, expand the utility of our modular PCR marker design
toolset for usage at the command-line and in the Galaxy workflow
environment.

Although Galaxy provides an excellent graphical environment for
organising large-scale variant analysis and assay design, we have
utilised the iPython Notebook (<http://ipython.org/notebook.html>) as a
means for desktop prototyping, documentation and sharing of design
workflows. Although widely-used in the physical sciences, such notebooks
have yet to gain wide usage for bioinformatics analysis. Our experience
shows they are a highly visual, portable and expressive means to combine
multiple tools in a reproducible manner. Since marker design from modern
sequencing data will typically require multiple steps employing a
variety of tools, the ability to use diverse tools and document this in
a single, web-ready JSON file is a major aid to reproducibility and
communication. While configuration of Galaxy tools requires server and
systems administration resources, iPython notebooks may be configured
for desktop or server usage. We have provided explanations of format
requirements and command-line usage in the wiki at
<https://github.com/cfljam/galaxy-pcr-markers/wiki/An-Introduction-to-the-PCR-Markers-Tools>.
The wiki includes links to worked examples in iPython notebooks. In
addition, we provide a simple Vagrant (https://www.vagrantup.com/)
virtual machine definition that can be used to create a local iPython
notebook server configured to run the notebook examples.

### HRM Marker Screening

By selecting primer sets with the criterion of large predicted melting
difference between homozygotes, we observed increased efficiency of SNP
validation by HRM compared to our previous report [14]. We designed 138
HRM markers to 112 contigs of the 263 identified as being significantly
associated with bolting (Table 2). For 23 contigs multiple markers were
designed, including five contigs where two HRM assays were designed
across the same SNP. Of the 115 assays that showed simple melt profiles,
77 (or 81%) could be used for mapping because the homozygous parental
alleles could be resolved. The median melting temperature differences
between homozygotes for these was 0.5 °C. This outcome represents an
overall validation rate of 56% across all assays and 67% for those
amplifying a single locus. These values compared to 56% polymorphic and
9% resolved in the previous study [14]. Therefore the modifications made
to the HRM design have significantly improved the validation rates.
Conversion rates of SNPs identified from RNA-Seq data to markers in
other studies have included 13% for Salmon using HRMA [22], 48% in onion
of SNP to Kaspar (LGC) [13] and 28% to Sequenom assays for rainbow trout
[23].

**HRM markers to *AcBlt1***

While the primary focus of this publication is the ability of the
pipeline to provide discriminating HRM markers, the data to date
confirms that it has provided meaningful leads for further functional
and genetic analysis. Bin mapping assignment and inspection of allele
frequencies reveals two markers closely linked to the *AcBlt1* regions
of interest (Table 3). Based on allele frequencies, ACPA40 is tightly
linked to the ACP267 marker, while ACPA59 may be closer to the QTL. The
ACPA59 marker matches 18 contigs in the assembly with similarity to
serine-arginine rich RNA binding proteins including *Arabidopsis*
ATSRP30. The regulation of this gene has been shown to be
temperature-dependent in *Arabidopsis* [24]. Since genetic variation in
onion bolting is conditioned by variation in vernalization response
[11], the gene or genes that these transcripts are derived from are
plausible candidates for the *AcBlt1* locus. The HRM markers that were
more closely linked to the *AcBlt1* locus than ACP267 could be tested in
wider germplasm for evaluation for marker-assisted selection (MAS) in
breeding programmes. This trait is binary but controlled by multiple QTL
and it would be interesting to test the strategy for a quantitative
trait.

Conclusions
-----------

HRM is a simple yet powerful technology for genetic analysis in
wet-bench genetics that can be readily scaled for screening and
validating leads from NGS. Our results confirm that the same
methodologies supported by the current generation of uMelt web
applications for clinical diagnostic assay design can be scaled to
enable SNP validation on a scale typically conducted on proprietary
platforms. We see considerable scope for improving the current toolkit
of scripts and services in areas such as development of multi-allelic
markers, probe design and secondary structure prediction.

Although initial development of our marker design toolkit focused on the
Galaxy bioinformatics environment, the iPython notebook has since
emerged as a powerful means to develop and share scientific code, output
and narrative as ‘executable explanations’. We suggest that iPython
notebook provides an ideal ‘electronic lab book’ for management and
exchange of complex code-based workflows such as marker design and is
worthy of wider usage by the bioinformatics community. Although the
filtering and manipulation of variant files and design of assays does
not require large computational resources, configuration and
customization of the software tools can be challenging for scientists.
Automated provision of software appliances custom configured for assay
design using Vagrant or containerization technologies such as Docker
(<https://www.docker.com/>) will be the most practical means in future
to deliver up-to-date Galaxy or iPython environments to scientists
desktops.

Evaluation of the leads identified to date in our BSR-seq data suggest
that the catalog of transcriptome and variant data we have developed
will provide a rich resource for genetic and physiological studies of
major adaptive traits in onion.

Methods
-------

### Tissue Sampling and RNA Sequencing

####  ‘CUDH2150’ Parental Reference Line

Plants of the homozygous doubled haploid onion line ‘CUDH2150’ (Cornell
University;[25, 26] were grown in short days (8h day) at 20 C to *circa*
11 leaves and then transferred to long days (16 hours, 20 degrees).
Three biological replicate samples were prepared from 1-2mm sections of
the innermost fully expanded leaf at transfer and after three days in
long days. RNA was isolated using Qiagen RNeasy Plant Mini Kit as per
the manufacturer’s instructions. Sequencing libraries were prepared
using TruSeq RNA Sample Preparation v2.0 (Illumina) and sequenced on an
Illumina HiSeq 2000 by NZ Genomics Ltd. Reads were submitted to Genbank
SRA (SRX392929, SRX393053- SRX393057) under BioProject PRJNA60277.

####  ‘CUDH2150’ x Nasik Red’ F<sub>2</sub> Progeny Pools

Propagation, phenotyping, sampling and genotyping of progeny from two
large F<sub>2</sub> families from the wide onion cross ‘CUDH2150’ x ‘Nasik Red’,
segregating for bolting (precocious flowering) were described previously
[18]. Bulbs from the two families which had been genotyped and found
homozygous at the ‘ACP267’ marker linked to the *AcBlt1* QTL were sorted
into four classes representing all combinations of the bolting phenotype
and ACP267 marker genotype (bolted and ACP267 ‘AA’; bolted and
ACP267’BB’;non bolted and ACP267’AA’;non bolted and ACP267’BB’ ; Table 1
). Tissue segments of similar size were dissected from each bulb and
stored in RNAlater (Life Technologies) at -20 °C. For RNA extraction,
the individual bulb samples from each class were divided randomly into
three separate biological replicates and ground to a powder in liquid
Nitrogen. In Teflon tubes, 1 g tissue was added to 5ml Concert Plant RNA
reagent (Invitrogen) and allowed to extract at room temperature with
shaking \> 1 hr. Samples were then spun at 12000 x g for 10 minutes at
room temperature. Supernatant (250 µL) was transferred into each of two
microcentrifuge tubes (100 mg total prep) and 500 µL Concert Plant RNA
reagent (Invitrogen) was added along with 150 µL of 5M NaCl and
inverted. Then 450 µL chloroform was added and the solution mixed
thoroughly by inversion. These were then spun at 12000 x g for 10
minutes at 4 °C and the top aqueous phase was transferred to a new
microcentrifuge tube with approximately 800 µL recovered.

An equal volume of isopropanol was added, the samples were mixed and
allowed to stand at room temperature for 10 min. This was followed by
centrifugation at 12000 x g for 10 min at 4 °C. The supernatant was very
carefully removed and 1 mL 75% ethanol added to the pellet. This was
then spun at max speed (12000 x g) for 5 minute at 4 °C. The pellet was
washed again with 500 µL of 75% ethanol and the centrifugation step
repeated. The supernatant was carefully removed and discarded. The
pellet was then spun for a further 1 min at maximum speed and any
remaining supernatant carefully removed using a pipette. The pellet was
dried for 30 min at room temperature and resuspended in 100 µL DEPC
water and combined into 1 tube. The RNA was then purified using the
RNeasy MiniElute Cleanup Kit (Qiagen) and manufacturer’s instructions
with a final eluate volume of 12 µL.

Twelve sequencing libraries were prepared by NZ Genomics Ltd using
TruSeq™ RNA Sample Preparation v2 and sequenced in two lanes of HiSeq
2000.

### Bioinformatics Analysis of RNA-Seq

#### Reference Assembly of ‘CUDH2150’

A hybrid reference assembly was produced using Trinity [27] by combining
Illumina reads with the Newbler assembly previously generated from 454
sequencing [14]. The contigs generated were then clustered using CD-Hit
[28, 29], screened for organellar and ribosomal contaminants and
submitted to Genbank as transcriptome shotgun assembly GBGJ00000000.1.

#### RNA PoolSeq Analysis

Raw reads were trimmed for quality and adapters using
fastq-mcf(https://code.google.com/p/ea-utils/wiki/FastqMcf). Reads were
mapped to the combined assembly using Bowtie2 version 2.2.1 [30] using
‘sensitive’ parameter settings. Variant detection was performed using
**samtools mpileup** [31]; version r982:295) with default settings.
Comparisons of variant frequencies among pools were performed using perl
scripts from Popoolation2 (V 1.201; [32] as follows:

1.  Variant counts were compiled from mpileup files using
    mpileup2sync.jar with minimum quality of 20.

2.  Regions around indels were removed from the sync file using
    identify-indel-regions.pl and filter-sync-by-gtf.pl

3.  Sync files were normalised by resampling to a uniform coverage of 50
    reads per site per sample using subsample-synchronized.pl with
    parameters --target-coverage 50 --max-coverage 1000 --method with
    replace

4.  Orthogonal, replicated comparisons among pool types were performed
    by Cochran-Mantel-Haenszel (CMH) tests using cmh-test.pl.

Compilation and filtering of CMH tests was performed using R.
Reproducible workflow for this is provided in Supplementary File 1
(BSA\_RNASEQ\_analysis.Rmd; Viewable online at
<http://rpubs.com/Cfljam/BSA_RNASEQ_Suppl_1>). A false discovery rate of
0.01% was used to correct for the multiple testing using the
Bioconductor qvalue package [33].

Vcf files for marker design were generated by adding read group
information to bam alignment files with Picard tools
AddOrReplaceReadGroups.jar and then calling variants with Freebayes [34]
using default settings. The final vcf file is available at
<http://dx.doi.org/10.5281/zenodo.11455>

### HRM Assay Design

#### Software Tool Development

Python scripts previously used to design PCR primers [14] were updated
to enable direct use of Primer3 V2.0 rather than the BioPython interface
to ePrimer3. Capability for predicting amplicon melting temperature of
reference and variant alleles was added to the key script
**design\_primers.py** through access to a uMelt [20] web service
provided at University of Utah
(https://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi). Python
scripts, wrappers for the Galaxy workflow environment and documentation
for command-line usage are available at GitHub
(<https://github.com/cfljam/galaxy-pcr-markers>). Tools are also
available for download from the Galaxy Toolshed
(<https://toolshed.g2.bx.psu.edu/repos/john-mccallum/pcr_markers>).

#### Assay Design for SNP Validation

Workflows for PCR assays design were developed and documented in iPython
notebooks (<http://ipython.org/notebook.html>), using the following
components:

1.  Scripts from <https://github.com/cfljam/galaxy-pcr-markers>

2.  Unix shell tools

3.  Python pandas (<http://pandas.pydata.org/>) for tabular data
    manpulation and filtering.

4.  Bedtools <http://bedtools.readthedocs.org/> and vcftools
    <http://vcftools.sourceforge.net/> for variant file manipulation

5.  ipcress from the exonerate package
    (<https://www.ebi.ac.uk/~guy/exonerate/>) for *in silico* PCR

A representative workflow is documented in [an iPython notebook
](https://github.com/cfljam/Onion_BSA_RNASEQ_Workflow/blob/master/HRM_design_analysis/02.primer_design/HRM%20Design%20from%20BSA%20RNASEQ%20Variants.ipynb)(<https://github.com/cfljam/Onion_BSA_RNASEQ_Workflow/blob/master/HRM_design_analysis/02.primer_design/HRM%20Design%20from%20BSA%20RNASEQ%20Variants.ipynb>)
in the Github repository which can also be viewed through the nbviewer
rendering website
(<http://nbviewer.ipython.org/github/cfljam/Onion_BSA_RNASEQ_Workflow/blob/master/HRM_design_analysis/02.primer_design/HRM%20Design%20from%20BSA%20RNASEQ%20Variants.ipynb>)
. We also provide a simple Vagrant (<https://www.vagrantup.com/>)
virtual machine definition for Oracle VirtualBox at
<https://github.com/cfljam/GFBToolbox>. This may be used to create a
local iPython notebook server configured to run these examples and
perform design of assays with user data.

A vcf file
(<https://github.com/cfljam/Onion_BSA_RNASEQ_Workflow/blob/master/HRM_design_analysis/02.primer_design/sig_contigs.recode.vcf>)
containing only the contigs with SNPs exhibiting frequency differences
among pools was compiled using bedtools and vcftools, and then converted
to gff3 format for primer design. Design of up to 5 primer sets
generating amplicons of 60-120 bp was specified for each SNP and the set
with highest predicted Tm difference was selected.

#### HRM analysis and mapping

Primer sets with a predicted amplicon product Tm difference greater than
0.6 between homozygotes were assessed for polymorphism using the parents
and F~1~ as described in Baldwin et al. [14] with the PCR cycling
modification described in Baldwin et al [18] (Supplementary Table 1).
Those assays that were polymorphic and where the different homozygote
classes could be resolved were then tested on 10 bin mapping lines from
the ‘2262’ population and the approximate genetic map position for the
marker located as described in Baldwin et al. [14].

Availability of supporting data
-------------------------------

The data sets supporting the results of this article are included within
the article and its additional files, or available in public
repositories at github.com, Genbank, and zenodo.org.

List of Abbreviations
---------------------

**BSA** bulked segregant analysis

**BSR-seq** bulked segregant RNA PoolSeq

**BSTA** bulked segregant transcriptome analysis

**CMH test** Cochran–Mantel–Haenszel test

**HRM** High resolution melting

**MAS** marker-assisted selection

**NGS** Next-generation sequencing

**QTL** mapping quantitative trait loci

**SSR** simple sequence repeat

Competing interests
-------------------

The authors declare no competing interests

Authors' contributions
----------------------

SB and JM conceived the study, conducted biological and *in silico*
analyses and prepared manuscript. BW,ZD and SL developed software
components. ST conducted bioinformatics analyses. RM, RL, MPJ and KW
conducted plant propagation, sampling, RNA library analysis.

Acknowledgments
---------------

This research was funded by the New Zealand Ministry for Business,
Innovation and Employment projects ‘Genes for Hybrid Seed’ and ‘Virtual
Institute for Statistical Genetics’. We gratefully acknowledge the
support of DNAture Ltd for reagents and Allium Solutions Ltd for field
trialling.

References
-----------

<span id="_ENREF_1" class="anchor"></span>1. Futschik A, Schlötterer C:
**The Next Generation of Molecular Markers From Massively Parallel
Sequencing of Pooled DNA Samples.** *Genetics* 2010, **186:**207-218.

<span id="_ENREF_2" class="anchor"></span>2. Boitard S, Schlötterer C,
Nolte V, Pandey RV, Futschik A: **Detecting Selective Sweeps from Pooled
Next-Generation Sequencing Samples.** *Molecular Biology and Evolution*
2012.

<span id="_ENREF_3" class="anchor"></span>3. Liu S, Yeh C-T, Tang HM,
Nettleton D, Schnable PS: **Gene Mapping via Bulked Segregant RNA-Seq
(BSR-Seq).** *PLoS ONE* 2012, **7:**e36406.

<span id="_ENREF_4" class="anchor"></span>4. Trick M, Adamski N, Mugford
S, Jiang C-C, Febrer M, Uauy C: **Combining SNP discovery from
next-generation sequencing data with bulked segregant analysis (BSA) to
fine-map genes in polyploid wheat.** *BMC Plant Biology* 2012,
**12:**14.

<span id="_ENREF_5" class="anchor"></span>5. Livaja M, Wang Y,
Wieckhorst S, Haseneyer G, Seidel M, Hahn V, Knapp S, Taudien S, Schon
C-C, Bauer E: **BSTA: a targeted approach combines bulked segregant
analysis with next- generation sequencing and de novo transcriptome
assembly for SNP discovery in sunflower.** *BMC Genomics* 2013,
**14:**628.

<span id="_ENREF_6" class="anchor"></span>6. Wang R, Sun L, Bao L, Zhang
J, Jiang Y, Yao J, Song L, Feng J, Liu S, Liu Z: **Bulk segregant
RNA-seq reveals expression and positional candidate genes and
allele-specific expression for disease resistance against enteric
septicemia of catfish.** *BMC Genomics* 2013, **14:**1-18.

<span id="_ENREF_7" class="anchor"></span>7. Li M, Zhou L, Palais RA,
Wittwer CT: **Genotyping Accuracy of High-Resolution DNA Melting
Instruments.** *Clinical Chemistry* 2014.

<span id="_ENREF_8" class="anchor"></span>8. Liew M, Pryor R, Palais R,
Meadows C, Erali M, Lyon E, Wittwer C: **Genotyping of single-nucleotide
polymorphisms by high-resolution melting of small amplicons. .** *Clin
Chem* 2004, **50:**1156-1164.

<span id="_ENREF_9" class="anchor"></span>9. Erali M, Wittwer CT: **High
resolution melting analysis for gene scanning.** *Methods (San Diego,
Calif)* 2010, **50:**250-261.

<span id="_ENREF_10" class="anchor"></span>10. Montgomery JL, Sanford
LN, Wittwer CT: **High-resolution DNA melting analysis in clinical
research and diagnostics.** *Expert Rev Mol Diagn* 2010, **10:**219–240.

<span id="_ENREF_11" class="anchor"></span>11. Brewster JL: *Onions and
other vegetable alliums.* CABI; 2008.

<span id="_ENREF_12" class="anchor"></span>12. Martin W, McCallum J,
Shigyo M, Jakse J, Kuhl J, Yamane N, Pither-Joyce M, Gokce A, Sink K,
Town C, Havey M: **Genetic mapping of expressed sequences in onion and
in silico comparisons with rice show scant colinearity.** *Molecular
Genetics and Genomics* 2005, **274:**197-204.

<span id="_ENREF_13" class="anchor"></span>13. Duangjit J, Bohanec B,
Chan AP, Town CD, Havey MJ: **Transcriptome sequencing to produce
SNP-based genetic maps of onion.** *Theoretical and Applied Genetics*
2013, **126:**2093-2101.

<span id="_ENREF_14" class="anchor"></span>14. Baldwin S, Revanna R,
Thomson S, Pither-Joyce M, Wright K, Crowhurst R, Fiers M, Chen L,
Macknight R, McCallum J: **A Toolkit for bulk PCR-based marker design
from next-generation sequence data: application for development of a
framework linkage map in bulb onion (Allium cepa L.).** *BMC Genomics*
2012, **13:**637.

<span id="_ENREF_15" class="anchor"></span>15. McCallum J, Thomson S,
Pither-Joyce M, Kenel F, Clarke A, Havey MJ: **Genetic Diversity
Analysis and Single-nucleotide Polymorphism Marker Development in
Cultivated Bulb Onion Based on Expressed Sequence Tag-Simple Sequence
Repeat Markers.** *Journal of the American Society for Horticultural
Science* 2008, **133:**810-818.

<span id="_ENREF_16" class="anchor"></span>16. Baldwin S, Pither-Joyce
M, Wright K, Chen L, McCallum J: **Development of robust genomic simple
sequence repeat markers for estimation of genetic diversity within and
among bulb onion (*Allium cepa* L.) populations.** *Mol Breed*
2012**:**1-11.

<span id="_ENREF_17" class="anchor"></span>17. Lee R, Baldwin S, Kenel
F, McCallum J, Macknight R: **FLOWERING LOCUS T genes control onion bulb
formation and flowering.** *Nat Commun* 2013, **4**.

<span id="_ENREF_18" class="anchor"></span>18. Baldwin S, Revanna R,
Pither-Joyce M, Shaw M, Wright K, Thomson S, Moya L, Lee R, Macknight R,
McCallum J: **Genetic analyses of bolting in bulb onion (Allium cepa
L.).** *Theoretical and Applied Genetics* 2013**:**1-13.

<span id="_ENREF_19" class="anchor"></span>19. Taylor A, Massiah AJ,
Thomas B: **Conservation of Arabidopsis thaliana Photoperiodic Flowering
Time Genes in Onion (*Allium cepa* L.).** *Plant and Cell Physiology*
2010, **51:**1638-1647.

<span id="_ENREF_20" class="anchor"></span>20. Dwight Z, Palais R,
Wittwer CT: **uMELT: prediction of high-resolution melting curves and
dynamic melting profiles of PCR products in a rich web application.**
*Bioinformatics* 2011, **27:**1019-1020.

<span id="_ENREF_21" class="anchor"></span>21. Konczal M, Koteja P,
Stuglik M, Radwan J, Babik W: **Accuracy of allele frequency estimation
using pooled RNA-Seq.** *Mol Ecol Resour* 2014, **14:**381-392.

<span id="_ENREF_22" class="anchor"></span>22. Seeb JE, Pascal CE, Grau
ED, Seeb LW, Templin WD, Harkins T, Roberts SB: **Transcriptome
sequencing and high-resolution melt analysis advance single nucleotide
polymorphism discovery in duplicated salmonids.** *Molecular Ecology
Resources* 2011, **11:**335-348.

<span id="_ENREF_23" class="anchor"></span>23. Salem M, Vallejo RL,
Leeds TD, Palti Y, Liu S, Sabbagh A, Rexroad CE, III, Yao J: **RNA-Seq
Identifies SNP Markers for Growth Traits in Rainbow Trout.** *PLoS ONE*
2012, **7:**e36264.

<span id="_ENREF_24" class="anchor"></span>24. Balasubramanian S,
Sureshkumar S, Lempe J, Weigel D: **Potent Induction of *Arabidopsis
thaliana* Flowering by Elevated Growth Temperature.** *PLoS Genet* 2006,
**2:**e106.

<span id="_ENREF_25" class="anchor"></span>25. Alan AR, Mutschler MA,
Brants A, Cobb E, Earle ED: **Production of gynogenic plants from
hybrids of Allium cepa L. and A. roylei Stearn.** *Plant Science* 2003,
**165:**1201-1211.

<span id="_ENREF_26" class="anchor"></span>26. Hyde PT, Earle ED,
Mutschler MA: **Doubled Haploid Onion (Allium cepa L.) Lines and Their
Impact on Hybrid Performance.** *Hortscience* 2012, **47:**1690-1695.

<span id="_ENREF_27" class="anchor"></span>27. Grabherr MG, Haas BJ,
Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,
Raychowdhury R, Zeng Q, et al: **Full-length transcriptome assembly from
RNA-Seq data without a reference genome.** *Nat Biotechnol* 2011,
**29:**644-652.

<span id="_ENREF_28" class="anchor"></span>28. Li W, Godzik A: **Cd-hit:
a fast program for clustering and comparing large sets of protein or
nucleotide sequences.** *Bioinformatics* 2006, **22:**1658-1659.

<span id="_ENREF_29" class="anchor"></span>29. Fu L, Niu B, Zhu Z, Wu S,
Li W: **CD-HIT: accelerated for clustering the next generation
sequencing data. .** *Bioinformatics* 2012, **28:**3150-3152.

<span id="_ENREF_30" class="anchor"></span>30. Langmead B, Salzberg SL:
**Fast gapped-read alignment with Bowtie 2.** *Nat Meth* 2012,
**9:**357-359.

<span id="_ENREF_31" class="anchor"></span>31. Li H, Handsaker B,
Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R,
Subgroup GPDP: **The Sequence Alignment/Map format and SAMtools.**
*Bioinformatics* 2009, **25:**2078-2079.

<span id="_ENREF_32" class="anchor"></span>32. Kofler R,
Orozco-terWengel P, De Maio N, Pandey RV, Nolte V, Futschik A, Kosiol C,
Schlötterer C: **PoPoolation: A Toolbox for Population Genetic Analysis
of Next Generation Sequencing Data from Pooled Individuals.** *PLoS ONE*
2011, **6:**e15925.

<span id="_ENREF_33" class="anchor"></span>33. Storey JD, Tibshirani R:
**Statistical significance for genomewide studies.** *Proceedings of the
National Academy of Sciences* 2003, **100:**9440-9445.

<span id="_ENREF_34" class="anchor"></span>34. Garrison E, Marth G:
**Haplotype-based variant detection from short-read sequencing.**
Preprint at: [arXiv:1207.3907v2
[q-bio.GN]](http://arxiv.org/abs/1207.3907v2) 

Table 1 Description of the bulks based on bolting phenotype and ACP267
genotype and number of bulbs used for each pool for RNA-seq analysis.
Bulb samples for each of the four genotype/phenotype combinations were
divided among three biological replicates.

  Bulk   Pool\_ID    Phenotype   ACP267 genotype   Number of bulb samples used for RNA extraction   M bp    Total Reads   % Q \> 30
  ------ ----------- ----------- ----------------- ------------------------------------------------ ------- ------------- -----------
  1      BoltA1      bolt        A                 113                                              4,382   43,385,360    88.6
  2      BoltA2      bolt        A                                                                  4,221   41,790,262    88.5
  3      BoltA3      bolt        A                                                                  5,046   49,957,544    88.8
  4      BoltB1      bolt        B                 36                                               5,356   53,035,938    88.7
  5      BoltB2      bolt        B                                                                  5,077   50,269,452    88.8
  6      BoltB3      bolt        B                                                                  4,681   46,346,486    89.0
  7      NonBoltA1   nonbolt     A                 25                                               3,929   38,892,038    88.9
  8      NonBoltA2   nonbolt     A                                                                  5,512   54,574,822    88.7
  9      NonBoltA3   nonbolt     A                                                                  4,352   43,082,942    88.7
  10     NonBoltB1   nonbolt     B                 78                                               4,743   46,958,184    89.3
  11     NonBoltB2   nonbolt     B                                                                  5,693   56,364,642    89.0
  12     NonBoltB3   nonbolt     B                                                                  4,493   44,489,736    88.8

Table 2 Summary of the HRM primer design and amplification results.

  **Description**                                 **Results**    ** **
  ----------------------------------------------- -------------- -------------------------------------------------
  Total number of contigs with markers designed   112             
  Duplicate SNP assays\*                          5               
  Total number of HRM markers assessed            138             
  Total amplified using standard conditions       138            100% of total assays
  Single locus and simple melt profiles           115            83% of total assays
  Polymorphic                                     95/115 (83%)   69% of total assays
  Polymorphic and homozygotes resolved            77/95 (81%)    67% of the single locus and 56% of total assays
  Bin-mapped                                                     bin map positions for 49
  resolved markers per contig                     69/112 (62%)    

\* same SNP in contig interval but different assay.

Table 3 Minor allele frequencies of SNPs assigned to chromosome 1

  **primer set**   **Genbank accession**   Minor allele frequency
  ---------------- ----------------------- ------------------------ ------- ---------- ----------
                                           BoltA
  ACPA40           GBGJ01051189            0.00
  ACPA59           GBGJ01059209            0.67
