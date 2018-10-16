Results in s3://czbiohub-mosquito/blast/*.diamond_blastx_sensitive.m8.gz

This is a Tab-separated file. Columns are described in the manual entry for `--outfmt/-f`:

```
6 BLAST tabular format (default). This format can be customized, the 6 may be followed by
a space-separated list of the following keywords, each specifying a field of the output.
     qseqid Query Seq - id
     qlen Query sequence length
     sseqid Subject Seq - id
     sallseqid All subject Seq - id(s), separated by a ’;’
     slen Subject sequence length
     qstart Start of alignment in query
     qend End of alignment in query
     sstart Start of alignment in subject
     send End of alignment in subject
     qseq Aligned part of query sequence
     sseq Aligned part of subject sequence
     full sseq Full subject sequence
     evalue Expect value
     bitscore Bit score
     score Raw score
     length Alignment length
     pident Percentage of identical matches
     nident Number of identical matches
     mismatch Number of mismatches
     positive Number of positive - scoring matches
     gapopen Number of gap openings
     gaps Total number of gaps
     ppos Percentage of positive - scoring matches
     qframe Query frame
     btop Blast traceback operations(BTOP)
     staxids Unique Subject Taxonomy ID(s), separated by a ’;’ (in numerical order). This
     field requires setting the --taxonmap parameter for makedb.
     stitle Subject Title
     salltitles All Subject Title(s), separated by a ’<>’
     qcovhsp Query Coverage Per HSP
     qtitle Query title

By default, there are 12 preconfigured fields: qseqid sseqid pident length mismatch
gapopen qstart qend sstart send evalue bitscore.
```
