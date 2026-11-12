# Part 3 – RNA-seq analysis of the mouse cerebellum (aging, dietary restriction, young plasma)

Jupyter notebooks (R kernel) for Part 3 of the Systems Genomics tutorial, sections 3-4 to 3-7 plus follow-ups.
The input data (RSEM tables of 27 samples) are not included; the notebooks read them from
`RNAseq_aging/for_teammates` (set `base` in the setup cell).

| Notebook | Content |
|---|---|
| `03_part3_step1b_followup.ipynb` | Sample comparison per cohort, effect of the outliers SRR21355389 (blood contamination) and SRR23914572 (AL 572), PCA on the most variable genes, how subtle the group differences are |
| `03_part3_step2.ipynb` | 3-4 to 3-7: ANCOVA and DESeq2 differential expression, grouping of DEGs, DAVID lists, GSEA (fgsea, MSigDB Hallmark + C8), cell type marker scores, isoform usage. SRR21355389 excluded |
| `03_part3_step2_clean.ipynb` | Same as step 2, with both bad samples excluded (SRR21355389 and SRR23914572) |
| `03_part3_step3_extra.ipynb` | Extra analyses: sex check, cohort-corrected overview (plots only), brain-aging genes, aging signature vs interventions, gene length, variability with age, mitochondrial/ribosomal share |

Comparisons are always made within one cohort: aging 21 mo vs 3 mo (covariate sex), aDR vs AL, YMP vs PBS.

Required R packages: tidyverse, DESeq2, ggrepel, pbapply, gplots, viridis, fgsea, msigdbr, data.table, patchwork, IRkernel.
