# mitotools

Share some scripts I usually use

- get_seq

  get_seq can be used to quickly extract annotated sequences from genbank files of mitochondrial genomes such as faa, cds, pep, trn, rrn.

  Run `get_seq --help` to show the program's usage guide.
  ```
  Usage:./get_seq -g <genbank_file> -a
  Required options:
     -g, --genbank  Intput genbank file
  Optional options:
     -pre, --prefix  Prefix of the output
     -a, --all    Flag to output all annotations
     -f, --faa    Flag to output fasta
     -p, --pep    Flag to output pep
     -c, --cds    Flag to output cds
     -t, --trn    Flag to output trn
     -r, --rrn    Flag to output rrn
     -o, --output The output path
     -h, --help      Display this help message
  
  ```
