# BD-Rhapsody
BD Rhapsody™ Single-Cell Analysis


# `Seven Bridges Command Line Interface`

### Install Seven Bridges Command Line Interface
First you nedd to install the Seven Bridges Command Line Interface using the folloqing command:

```bash -c 'curl https://igor.sbgenomics.com/downloads/sb/install.sh -sSf | sudo -H sh'```
### Configure credentials 
After installing Seven Bridges Command Line Interface, you must enter your credentials to authenticate with the Platform. You will the API andpoint and the authentification token. These information can be found in the developer panel.
Launch the command line and enter the following:

``` sb configure ```

###  Upload fastq files 

sb upload start Data/GCTACGCT_S1_R1_L001.fastq.gz Data/GCTACGCT_S1_R2_L001.fastq.gz --destination barber/barbara-pancreatic-cell-ide --tag uploadfastq
