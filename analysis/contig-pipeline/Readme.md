## Pipeline outline

0. Start with raw fastqs, trim adaptors
0. Map to mosquito ref
0. Subtract, map to mito ref
0. Subtract, quality filter
0. Assembly

### Alignment tools

Pick one of minimap2, STAR, or bowtie, maybe do some testing

### Quality filtering

0. PRICEseq_filter (try different settings)
0. Try another tool?

### Assembly tools

0. PRICE (contig_miner.py)
0. Spades

## Locations of references

- s3://czbiohub_mosquito/references/mosquito_genomes2

TODO: mito reference

## Samples to test

- West Nile Virus samples
  - CMS_088a
  - CMS_088b
  - CMS_089a
  - CMS_089b
- CMS_040
- Bunya ???
- Sparrow/Chicken
