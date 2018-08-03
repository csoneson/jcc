# Junction coverage compatibility scores

This repository contains an R package aimed at facilitating the calculation of gene-wise junction coverage compatibility (JCC) scores, as described in 

- Soneson C, Love MI, Patro R, Hussain S, Malhotra D and Robinson MD: A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs. [bioRxiv doi:10.1101/378539](https://www.biorxiv.org/content/early/2018/07/28/378539).

The JCC scores are aimed at detecting genes for which estimated transcript abundances are unreliable, either because of problems in the transcript abundance estimation or because of missing or wrongly annotated reference transcripts. It does so by comparing the number of reads predicted to align across each junction, inferred from the transcript abundances and a fragment bias model, to the observed number of junction-spanning reads, obtained via alignment of the reads to the genome. A high JCC score for a gene indicates that the estimated abundances for the corresponding transcripts are unreliable and should be treated with care in downstream analyses. 