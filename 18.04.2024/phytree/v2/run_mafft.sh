#!/usr/bin/bash
mafft --thread 70 --maxiterate 1000 --auto classes_seqs_subset_pr2_asvs.fa > classes_seqs_subset_pr2_asvs.afa
mafft --thread 70 --maxiterate 1000 --localpair classes_seqs_subset_pr2_asvs.fa > localpair_classes_seqs_subset_pr2_asvs.afa
mafft --thread 70 --maxiterate 1000 --genafpair classes_seqs_subset_pr2_asvs.fa > genafpair_classes_seqs_subset_pr2_asvs.afa
mafft --thread 70 --maxiterate 1000 --globalpair classes_seqs_subset_pr2_asvs.fa > globalpair_classes_seqs_subset_pr2_asvs.afa
