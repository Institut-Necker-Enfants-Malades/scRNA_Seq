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
### Some commands
Get the information about the user name with the following command :

``` sb whoami```

List the cuurent projects with the following command :

```sb projects list```

To see all the available commands use the following command :

``` sb -h```
###  Upload fastq files 

To upload a large data files, the command line is recommanded by SenvenBrigdes. The command is following :

```
sb [global-parameters] upload <upload-subcommand> [command-parameters]
# The upload-subcommand :
sb upload start	# Start the upload job
sb upload status	# Check the upload status
sb upload resume	# Resume the upload job
```
Start uploading files with the following command :
``` 
sb [global-parameters] upload start [<input> | --manifest-file <manifest-file>] --destination <path> [command-parameters]
# Example :
sb upload start Data/GCTACGCT_S1_R1_L001.fastq.gz Data/GCTACGCT_S1_R2_L001.fastq.gz --destination barber/barbara-pancreatic-cell-ide --tag uploadfastq
```

