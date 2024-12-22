#!/usr/bin/bash
mafft --thread 20 --maxiterate 100 --auto only_euk_classes_seqs_subset_pr2_asvs.fa > only_euk_classes_seqs_subset_pr2_asvs.afa
mafft --thread 20 --maxiterate 100 --localpair only_euk_classes_seqs_subset_pr2_asvs.fa > localpair_only_euk_classes_seqs_subset_pr2_asvs.afa
mafft --thread 20 --maxiterate 100 --genafpair only_euk_classes_seqs_subset_pr2_asvs.fa > genafpair_only_euk_classes_seqs_subset_pr2_asvs.afa
mafft --thread 20 --maxiterate 100 --globalpair only_euk_classes_seqs_subset_pr2_asvs.fa > globalpair_only_euk_classes_seqs_subset_pr2_asvs.afa
