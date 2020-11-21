
The file all_univariable_dens_biorxiv.txt.gz includes all summary statistics for association with T1D accompanying [Robertson, Inshaw, et al (2020)](https://www.biorxiv.org/content/10.1101/2020.06.19.158071v1)

## Explanation of columns in summary statistics download

* __alleleB__ is the effect allele, meaning the "Effect_EUR" and "Effect_phaseIandII" columns are the betas (log odds ratios) for each additional copy of the B allele
* __AF_EUR__ is the frequency of alleleB in unrelated European cases (N=13,458) and controls (N=20,143)
* __Effect_EUR__ is the log odds ratio for each additional copy of alleleB among unrelated European cases (N=13,458) and controls (N=20,143)
* __SE_EUR__ is the standard error for effect described in Effect_EUR
* __EffectphaseIandII__ is the log odds ratio for each additional copy of alleleB in meta-analysis of the full cohort (16,159 cases, 25,386 controls, and 6,143 trios)
* __SEphaseIandII__ is the standard error for effect described in Effect_phaseIandII
* __pphaseIandII__ is the p-value for the effect described in EffectphaseIandII


## FAQs

- __Which is the effect allele?__
alleleB is the effect allele, meaning the "Effect_EUR" and "Effect_phaseIandII" columns are the betas (log odds ratios) for each additional copy of the B allele.

- __Which column is used to report sample size?__
The sample sizes are not in this table. Since they do not vary by variant, we include them only once in Supplementary Table 5 in the manuscript. Additionally, for the "phaseIandII" meta-analysis results, this effect is based on a meta-analysis of case-control and family-based association data. We report the number of cases, controls, and trios used in each analysis, separately. If you need a single number as the total sample size for the meta-analysis, you can count each trio as one case/control pair (a TDT analysis of N trios is equivalent in terms of statistical power to having N cases and N controls).

- __Which effect size estimates do you suggest using your analyses?__
The "phaseIandII" summary statistics (EffectphaseIandII, SEphaseIandII, and pphaseIandII) incorporate the largest number of samples (all samples in Supplementary Table 5, except the EAS case-control samples). The total sample sizes for these columns are 16,159 cases, 25,386 controls, and 6,143 trios. If your analysis requires linkage disequilibirium patterns from a single major ancestral group, we also provide European-specific summary statistics (Effect_EUR and SE_EUR). 
