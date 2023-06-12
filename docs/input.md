## Input data

The input data must be an absolute path to a folder containing subfolders or files with prefix name `barcode`

It's usually the `fastq_pass` folder from an ONT run.

```
.
├── barcode04
│   ├── FAW57950_pass_barcode04_80f7ac76_c6f011b6_0.fastq.gz
│   └── FAW57950_pass_barcode04_80f7ac76_c6f011b6_1.fastq.gz
├── barcode05
│   └── FAW57950_pass_barcode05_80f7ac76_c6f011b6_0.fastq.gz
├── barcode06
│   └── FAW57950_pass_barcode06_80f7ac76_c6f011b6_0.fastq.gz
├── barcode11
│   └── FAW57950_pass_barcode11_80f7ac76_c6f011b6_0.fastq.gz
└── unclassified
    ├── FAW57950_pass_unclassified_80f7ac76_c6f011b6_0.fastq.gz
    ├── FAW57950_pass_unclassified_80f7ac76_c6f011b6_1.fastq.gz
    ├── FAW57950_pass_unclassified_80f7ac76_c6f011b6_2.fastq.gz
    ├── FAW57950_pass_unclassified_80f7ac76_c6f011b6_3.fastq.gz
    └── FAW57950_pass_unclassified_80f7ac76_c6f011b6_4.fastq.gz
```

or a folder of compressed barcode fastq files

```
.
├── barcode01.fastq.gz
├── barcode02.fastq.gz
├── barcode03.fastq.gz
├── barcode04.fastq.gz
├── barcode05.fastq.gz
├── barcode06.fastq.gz
├── barcode07.fastq.gz
├── barcode09.fastq.gz
├── barcode11.fastq.gz
└── barcode12.fastq.gz
```
