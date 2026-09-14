# PAN_BUF_1

This is a repository of codes for the paper: Reference-free water buffalo pangenome enables structural variant discovery and method benchmarking.

1.  [Input file preparation](01_inputprep.md)
2.  [Graph construction with PGGB](02_pggb.md)
3.  [Identifying variants of the pangenome](03_identify_variants.md)
4.  [Repeat Analysis](04_repeat.md)
5.  [Benchmarking of SVs from pangenome graph](05_SVIM.md)
6.  [Variant calling of short reads with the pangenome graph](06_variant_calling.md)
7.  [Prediction of high impact variants](07_vep.md)
8.  [Graph-based vs linear SNP calling](08_concordance.md)
9.  [Computation of graph complexity](09_complexity.md)


## HPC Hardware specifications

| Node Type       | Nodes | Node Specification                                                                                       | Total Specification                                | Performance | GoLive    |
|-----------------|-------|----------------------------------------------------------------------------------------------------------|----------------------------------------------------|-------------|-----------|
| CPU             | 105   | 2x Intel(R) Xeon(R) Platinum 8360Y CPU, 36 cores @ 2.4GHz,   255000 MiB memory(around 256 GB per node)   | 7560 CPUs, 26 TB of memory                         | N/A         | 2022-2023 |
| CPU High Memory | 8     | 2x Intel(R) Xeon(R) Platinum 8360Y CPU, 36 cores @ 2.4GHz,   1990000 MiB memory (around 2TB per node)    | 576 CPUs, 16 TB of memory                          | N/A         | 2022-2023 |
| GPU             | 50    | 2x Intel(R) Xeon(R) Platinum 8360Y CPU, 36 cores @ 2.4GHz,   515000 MiB memory, 4x Nvidia A100-SXM4-40GB | 3600 CPU, 200 GPUs, 25 TB of memory                | N/A         | 2022-2023 |
| SYSTEM TOTAL    | 163   |                                                                                                          | 11736 cores, 200 GPU   accelerators, 242 TB memory |             |           |