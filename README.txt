This is just a basic Python3 text parser for Uniprot, also Swiss-protein or Tremble protein database entries.  It is set up to parse from the N-terminal starting residues, always an M followed by single letter amino acid code.  If you go longer than the first 10 amino acids, you have to use a space, as it is very basic and follows standard block 10, 60 residue per line format.  It is written very basic and modular, and currenttly set to also look for any GO, Pfam and Interpro entries, and I have it set 3 identifiers.  If you wish to expand or delete this aspect, simple comment out, or duplicate these for more identifiers as necessary, or just provide the same files for all 3x3 set using all identifiers if you wish everything.  Example input files are provided, and Uniprot .dat files can be downloaded for all entries, subsets of entries, per organisms, limited to Tremble or Swiss-protin as one wishes.  These change rapidly, especially the identifiers, such as GO or Pfam, due to constant annotation changes in the databases.  Keep in mind Tremble maintains many isoforms, unannotated and partial gene sequences, so most will lack identifiers for function or functional grouping, and only 1/3 of the human genome is currently annotated.  These files can all be downloaded from the following, and whole genome fasta can be extremely large in the gigabyte ranges;

For Uniprote .dat files organised by Kindom, or s few subgroups;

https://ftp.uniprot.org/pub/databases/uniprot/
There are several other sites with the same style data format as Uniprot

For the identifiers, these can be downloaded from the following in the .tsv file format with the columns matching properly as per the included examples;

https://www.ebi.ac.uk/interpro/entry/InterPro/#table
https://supfam.org/SUPERFAMILY/function.html
https://www.pantherdb.org/panther/prowler.jsp?
https://pfam-docs.readthedocs.io/en/latest/api.html

The script is very basic, and if you you use it or modify it, please cite;

https://doi.org/10.1101/2025.03.27.645863 

which will be replaced or alink added when published to the journal.  I need any brownie points I can get, the system in Brazil and most of the world is now running on a video game score points by references or equivalent for any job markets. 

There's also a general instruction on use at the top of any scripts, run as command line$ python3 script ...
