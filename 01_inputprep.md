# Preparation of the genome assemblies input files

1. First is to determine the homologous chromosomes using `scripts/mashmap.sh`
2. Determine contig N50 using `scripts/quast.sh`
3. Next is to separate each chromosomes `scripts/separate_chr_contig.sh`
4. Then remove the gaps `scripts/run_remove_gaps.sh`
5. Then concatenate the contigs per chromosome in all assemblies `scripts/concat_chr.sh`

6. The input files for PGGB would be per chromosome containing sequences from each assemblies with headers like `RVCU#0#4#piece0` from PanSN-spec (https://github.com/pangenome/PanSN-spec)

## Scaffolding the Pakistan river assemblies

To scaffold the Pakistan river assemblies, I use `scripts/ragtag.sh`. Then run `scripts/dotplot.sh` to visualise the scaffolding using dotplot.
