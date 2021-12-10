# About GFF File Parser

Generic Feature Format (GFF) is a biological sequence file format for representing features and annotations on sequences. This format depicts any experimentally varied information in the sequence, so providing simple access to specific features would help increase the accessibility of large amounts of biological sequence data. 

This file parser provides programmatic evaluation of GFF files within the R language. This program reads in the data file, converts it into a data frame format and reconciles coding sequences with their given genes, extracting the resulting matched data and storing relevant information in a final data frame for optimal accessibility. The program checks the validity of the gene matches, and accounts for common abnormalities in the inputted data.

This work was done under Nick Cooley at the Department of Biomedical Informatics at the University of Pittsburgh

# How to Use

Update ```gffFile``` variable to desired input file

Access ```parsedData``` data frame to view sequencing features
