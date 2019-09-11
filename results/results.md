# Variables

**Importance** score is the percentage of reads identified as a specefic isomiRs divided but the total number of reads that mapped to the miRNA gene.


# Synthetic samples

## Number of synthetic miRNAs with isomiRs

All methods showed isomiRs, although the importance of each of isomiR type are different among methods (NEB) and reproducible among studies.

[figure 1](02_importance/all_samples_barplot_num_mirna.pdf)

If we filter isomiRs with importance < 5%, only a few percentage of the synthetic miRNAs showed isomiRs and there was a method bias toward different types of isomiRs. Curiosly, only one method showed no isomiRs at all (kim et al), but is confounded with the study, being impossible to know whether is reproducible among studies.

[figure 2](02_importance/all_samples_barplot_num_mirna_5pct.pdf)


## Comparing with real samples

The number of miRNAs with isomiRS after the 5% cutoff that showed improvements with the synthetic data, is 3 or 4 times more than the value for the synthetic analysis (only for some types of isomiRs). This could mean that in real samples there is isomiRs with higher abundance and don't follow the pattern showed in the synthetic data, where we assume isomiRs are artifacts. 

Samples were ordered by protocol, and real samples are in bold compared to the synthetic samples.

[figure 3](04_synth_vs_real/barplot_total100_pct5.pdf)
