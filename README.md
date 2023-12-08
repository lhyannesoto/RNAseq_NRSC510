# RNAseq_NRSC510
For NRSC510A course.
Lhyanne :) 

This project was made for the NRSC510A course. The main objective of this project is to learn basic R language code to be able to analyze data sets from experiments. I do not have any prior experience in R before starting this course. My R language tutor is Kaitlin.

The data I will use my R scripts to analyze are from experiments investigating the specific role and mechanisms of a transcription factor called Pu.1 in microglia. Microglia are the resident immune cells in the brain, which are important for constantly surveying the brain parenchyma for invading pathogens, maintaining brain homeostasis, and many other important functions. This data will come from a novel mouse model that has an Alzheimer's Disease associated mutation knocked in that is orthologous to a risk allele that is found to lead to early onset in humans. In an attempt to understand how this happens, I will be studying differential gene expression (RNA-seq) and eventually chromatin accessibility (ATAC-seq).

Our mouse model leads to an overexpression of Pu.1, which is a lineage determining factor in microglia. Interestingly, Pu.1 is found in all macrophages, but in the brain it is only expressed in microglia. Pu.1 binds to DNA in nucleosomes and triggers chromatin remodeling to make DNA more accessible to other transcription factors and regulatory gene expression proteins. A single nucleotide polymorphism (SNP) is knocked into a microglia-specific enhancer in the mice that leads to Pu.1 being up-regulated. It is very rare for mice and humans to have orthologous chromatin regions that also interact but it was shown by Shemer at colleagues (Shemer et al. Immunity 2020 PMID: 33049219) that this allele is orthologous and also has very similar interactions 3D interactions within the genome. It is also important that the interactions between Pu.1 and the enhancers are functional, which was shown by the Gjoneska lab (Gjoneska et al. Nature 2015 PMID: 25693568) using luciferase functional assays in vitro.

Alongside analysis of RNAseq and ATACseq, I will also be investigating changes in microglial functions to assess differences in phenotypes.

Overall research question: How does this risk allele impact gene expression and chromatin accessibility, leading to altered microglia function? By understanding how overexpression of Pu.1 leads to early onset of Alzheimer's disease and neurodegeneration, we can determine potential pharmacological and therapeutic intervention.

Experiments: Samples are collected from mouse brains and the microglia cells are isolated and sorting with flow cytometry. Next, RNA is extracted from the sorted cells and undergo RNA sequencing.

The subjects used will include experimental groups of different developmental time points such as postnatal (P) 7, P14, P30, P80 and 6 months. We aim to determine at what age the increased expression becomes significant and how this may effect microglial function.

Dataset: The data is represented in the form of a counts matrix. This is raw data as a result of bulk RNA-seq. The rows indicate gene, which is identified by the ensemble ID and the gene name. The columns represent each sample and the level of expression for each gene of that sample.
