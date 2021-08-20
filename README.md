# Vcf2Hgvs
This program uses the [pysam](https://pysam.readthedocs.io/en/latest/api.html) and [Biocommons HGVS](https://github.com/biocommons/hgvs) library to generate HGVS entries for all variants in a VCF file. 

The initial inspiration for this was the variantFormatter program. However, since I wanted to use actual VCF files and wanted to use latest release of the Biocommons HGVS library, I ended up rolling my own.

## Installation

This script assumes that you have both pysam and Biocommons HGVS installed (including all the accoutrements: [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta/) and [SeqRepo](https://github.com/biocommons/biocommons.seqrepo/))

## Usage
```
python3 Vcf2Hgvs.py VCF/COSMIC_chr5.vcf.gz JSON/COSMIC_chr5.json
```

The arguments are as follows:

 1. input VCF path
 2. output JSON path


### Performance

While powerful, biocommons HGVS is not super-fast. On my machine, I averaged 15 variants/s. Splitting up the VCF file into smaller pieces and using gnu parallel or a cluster is advised if you're looking to do this for several million variants.

### JSON Output

Each VCF position will be evaluated for each alternate allele and every transcript that overlaps that alternate allele. Here's the output from one such variant:

```json
{
  "variants": [
    {
      "vid": "chr17-43342141-G-C",
      "hgvsg": "NC_000017.10:g.43342141G>C",
      "transcripts": [
        {
          "hgvsc": "NR_110326.1:n.164-2339G>C",
          "hgvsp": "non-coding"
        },
        {
          "hgvsc": "NM_003954.3:c.2706C>G",
          "hgvsp": "NP_003945.2:p.(Val902=)"
        },
        {
          "hgvsc": "NM_003954.5:c.2706C>G",
          "hgvsp": "NP_003945.2:p.(Val902=)"
        },
        {
          "hgvsc": "NR_024435.2:n.265-1476G>C",
          "hgvsp": "non-coding"
        },
        {
          "hgvsc": "NM_003954.4:c.2706C>G",
          "hgvsp": "NP_003945.2:p.(Val902=)"
        },
        {
          "hgvsc": "NR_024434.2:n.80-2339G>C",
          "hgvsp": "non-coding"
        },
        {
          "hgvsc": "NR_110324.1:n.265-2339G>C",
          "hgvsp": "non-coding"
        },
        {
          "hgvsc": "NR_110325.1:n.260-2339G>C",
          "hgvsp": "non-coding"
        }
      ]
    }
  ]
}
```

