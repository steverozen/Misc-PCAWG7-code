# Misc-PCAWG7-code
Miscellaneous code from "The Repertoire of Mutational Signatures in Human Cancer, https://www.nature.com/articles/s41586-020-1943-3

GPL-3 license

There are two code files:

* matrix_plot/matrix_plot_functions_v0.1.R, for plotting figures such as Fig 3 in https://www.nature.com/articles/s41586-020-1943-3. 
The folder also contains example input and output.

* gaddy-gram.R, for plotting figures such as the "snake plots", a.k.a. "hamburger plots" a.k.a "Gaddy grams" at 
e.g https://cancer.sanger.ac.uk/cosmic/signatures/SBS/SBS1.tt.

Some tumours were excluded from the plots:

* For all tumors and mutation types we required > 0.9 reconstruction accuracy.

* For SBS (single base substitutions) we required >= 100 total mutations.

* For DBSs (doublet base substitutions) and IDs (small insertions and deletions) we required >= 25 total mutations.

* In addition, for the snake plots for SBSs, we required >= 10 mutations per signature.

There was also some merging of cancer types.
