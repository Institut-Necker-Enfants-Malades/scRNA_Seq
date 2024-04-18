# BD-Rhapsody
BD Rhapsody™ Single-Cell Analysis


## Upload fastq files via Sevenbridges CLI

+ First you nedd to install the Seven Bridges Command Line Interface using the folloqing command:
```bash -c 'curl https://igor.sbgenomics.com/downloads/sb/install.sh -sSf | sudo -H sh'```

sb upload start Data/GCTACGCT_S1_R1_L001.fastq.gz Data/GCTACGCT_S1_R2_L001.fastq.gz --destination barber/barbara-pancreatic-cell-ide --tag uploadfastq
