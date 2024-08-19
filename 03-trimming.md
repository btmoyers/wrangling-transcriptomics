---
title: Trimming and Filtering
teaching: 30
exercises: 25
---

::::::::::::::::::::::::::::::::::::::: objectives

- Clean FASTQ reads using cutadapt.
- Select and set multiple options for command-line bioinformatic tools.
- Write `for` loops with two variables.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How can I get rid of sequence data that does not meet my quality standards?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Cleaning reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality graphs showing the distribution of read quality at each base across all reads in a sample and extracted information about which samples have warnings or failures for which quality checks. 

It is very common to have some quality metrics fail or have some moderately concerning values, and this may or may not be a problem for your downstream application. For our RNA-Seq workflow, we can filter out reads that have remnants from library preparation and sequencing or remove some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) to filter poor quality reads and trim poor quality bases from our samples.

## Accessing the cutadapt module

Remember that we can look for modules (programs like this that might be installed but that not everyone needs all the time) with `module avail`:

```bash
$ module avail cutadapt
```

which on chimera should give you this output:

```output
-------------------------------- /share/apps/modulefiles/modules --------------------------------
   py-cutadapt-2.10-gcc-10.2.0-2x2ytr5

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching any of the
"keys".
```

As it says, if you think something **should** be there and it isn't showing up, you can try `module spider`.

In our case, we will need to also load these modules to help cutadapt work:
```
module load py-dnaio-0.4.2-gcc-10.2.0-gaqzhv4
module load py-xopen-1.1.0-gcc-10.2.0-5kpnvqq
module load py-cutadapt-2.10-gcc-10.2.0-2x2ytr5
```

## cutadapt options

cutadapt has a variety of options to trim your reads. If we run the following command, we can see some of our options.

```bash
$ cutadapt --help
```

Which will give you the following output:

```output
cutadapt version 2.10

Copyright (C) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard characters are supported. All reads from input.fastq will be written to output.fastq with the adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter sequences can be given (use further -a options), but only the best-matching adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for standard input/output. Without the -o option, output is sent to standard output.
```

Plus a lot more options.

In paired end mode, cutadapt expects the two input files (R1 and R) after the names of the output files given after the options `-o` and `-p`. These files are described below. 

| option         | meaning                                                                                                      | 
| -------------- | ------------------------------------------------------------------------------------------------------------ |
| \<outputFile1> | After -o, output file that contains trimmed or filtered reads from the first input file (typically 'R1')                        | 
| \<outputFile2> | After -p, output file that contains trimmed or filtered reads from the first input file (typically 'R2')                        | 
| \<inputFile1>   | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.                                          | 
| \<inputFile2>   | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.                         | 

Cutadapt can remove adapter sequences that are in your sequence data, trim or remove low-quality bases or reads, and do a few other useful things. We will use only a few of these options and trimming steps in our analysis. It is important to understand the steps you are using to clean your data. For more information about the cutadapt arguments and options, see [the cutadapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html).

However, a complete command for cutadapt will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

```bash
$ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25 -o example_R1_trim.fastq -p example_R2_trim.fastq exampleR1.fastq exampleR2.fastq
```

In this example, we have told cutadapt:

| code           | meaning                                                                                                      | 
| -------------- | ------------------------------------------------------------------------------------------------------------ |
| `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`               | identify and remove bases that match the Illumina adapter from each R1 read       | 
| `-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`               | identify and remove bases that match the Illumina adapter from each R2 read       | 
| `--trim-n`               | trim any N bases from the beginning or end of each read                                                  |
| `-m 25`               | discard any reads that are shorter than 25 bases after trimming                                             | 
| `-o example_R1_trim.fastq`               | the output file for trimmed and filtered reads from the R1 input file                   | 
| `-p example_R2_trim.fastq`               | the output file for trimmed and filtered reads from the R2 input file                   | 
| `exampleR1.fastq`               | the input R1 fastq file          | 
| `exampleR2.fastq`               | the input R2 fastq file            | 

:::::::::::::::::::::::::::::::::::::::::  callout

## Multi-line commands

Some of the commands we ran in this lesson are long! When typing a long command into your terminal or nano, you can use the `\` character to separate code chunks onto separate lines. This can make your code more readable.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Running cutadapt

Now we will run cutadapt on our data. To begin, navigate to your `untrimmed_fastq` data directory:

```bash
$ cd ~/1_project/data/untrimmed_fastq
```
We are going to run cutadapt on one of our paired-end samples (V1). We will identify and remove any leftover adapter sequences from the reads. We will also remove N bases from the ends of the reads and filter out any reads that are shorter than 25 bases in our trimmed sequences (like in our example above). 

```bash
$ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 35 -o V1_S1_L001_R1_001_ds_trim.fastq -p V1_S1_L001_R2_001_ds_trim.fastq V1_S1_L001_R1_001_downsampled.fastq V1_S1_L001_R2_001_downsampled.fastq
```

```output
This is cutadapt 2.10 with Python 3.8.6
Command line parameters: -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 35 -o V1_S1_L001_R1_001_ds_trim.fastq -p V1_S1_L001_R2_001_ds_trim.fastq V1_S1_L001_R1_001_downsampled.fastq V1_S1_L001_R2_001_downsampled.fastq
Processing reads on 1 core in paired-end mode ...
[         8=-] 00:00:20       962,229 reads  @     21.8 µs/read;   2.76 M reads/minute
Finished in 21.23 s (22 us/read; 2.72 M reads/minute).

=== Summary ===

Total read pairs processed:            962,229
  Read 1 with adapter:                  18,572 (1.9%)
  Read 2 with adapter:                  21,863 (2.3%)
Pairs written (passing filters):       961,614 (99.9%)

Total basepairs processed:    98,147,358 bp
  Read 1:    49,073,679 bp
  Read 2:    49,073,679 bp
Total written (filtered):     97,939,516 bp (99.8%)
  Read 1:    48,977,906 bp
  Read 2:    48,961,610 bp
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Exercise

Use the output from your cutadapt command to answer the
following questions.

1) What percent of read pairs passed our filters?
2) What percent of basepairs passed our filters?

:::::::::::::::  solution

## Solution

1) 99\.9%
2) 99\.8%

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

We can confirm that we have our output files:

```bash
$ ls V1*
```

```output
V1_S1_L001_R1_001_downsampled.fastq  V1_S1_L001_R2_001_downsampled.fastq
V1_S1_L001_R1_001_ds_trim.fastq      V1_S1_L001_R2_001_ds_trim.fastq
```

The output files are also FASTQ files. It might be smaller than our input file, if we have removed reads. We can confirm this:

```bash
$ ls -thor V1*
```

```output
-rw-rw-r-- 1 brook.moyers 146M Aug  7 00:55 V1_S1_L001_R1_001_downsampled.fastq
-rw-rw-r-- 1 brook.moyers 146M Aug  7 00:55 V1_S1_L001_R2_001_downsampled.fastq
-rw-rw-r-- 1 brook.moyers 146M Aug 16 17:22 V1_S1_L001_R1_001_ds_trim.fastq
-rw-rw-r-- 1 brook.moyers 146M Aug 16 17:22 V1_S1_L001_R2_001_ds_trim.fastq
```

Hmmm, it doesn't look that different (they are all 146 MB). Maybe they are not that different because we kept most reads! We could check the number of lines, maybe using a for loop.

```bash
$ for filename in V1*; do lines=$(wc -l ${filename}); echo $lines; done
```

```
3848916 V1_S1_L001_R1_001_downsampled.fastq
3846456 V1_S1_L001_R1_001_ds_trim.fastq
3848916 V1_S1_L001_R2_001_downsampled.fastq
3846456 V1_S1_L001_R2_001_ds_trim.fastq
```
Yes, it looks like the trimmed files are a couple thousand lines shorter.

We have just successfully run cutadapt on one pair of our FASTQ files! However, there is some bad news. cutadapt can only operate on one sample at a time and we have more than one sample. The good news is that we can use a script to iterate through our sample files quickly!

Here is the text of a script you can use to do it:

```bash
#!/bin/bash

#SBATCH --job-name=trim # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=2
#SBATCH --account=itcga # our account
#SBATCH --reservation=ITCGA_AUG2024 # this gives us special access during the workshop
#SBATCH --time=1:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --error=%x-%A.err   # a filename to save error messages into
#SBATCH --output=%x-%A.out  # a filename to save any printed output into

# Usage: sbatch cutadapt.sh path/to/input_dir path/to/output_dir
# Works for paired end reads where both end in the same *_001_downsampled.fastq

# Module load
module load py-dnaio-0.4.2-gcc-10.2.0-gaqzhv4
module load py-xopen-1.1.0-gcc-10.2.0-5kpnvqq
module load py-cutadapt-2.10-gcc-10.2.0-2x2ytr5

# Define variables
input_dir=$1 # takes this from the command line, first item after the script
output_dir=$2 # takes this from the command line, second item

for R1_fastq in ${input_dir}/*_R1_001_downsampled.fastq
 do
 # Pull basename
 name=$(basename ${R1_fastq} _R1_001_downsampled.fastq)

 # Run cutadapt
 cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25 \
 -o ${output_dir}/${name}_R1_001_ds_trim.fastq \
 -p ${output_dir}/${name}_R2_001_ds_trim.fastq \
 ${input_dir}/${R1_fastq} \
 ${input_dir}/${name}_R2_001_downsampled.fastq

 echo cutadapt is finished with $name

done
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Exercise
Take a second to look, really look, through the script and note what each line does. There are a few new elements here that you might not have encountered before: what could you do to understand them?

1. Write down in plain English what this script does, line by line.
2. How could you modify this script if your FASTQ files ended in `_1.fastq` and `_2.fastq` instead?

:::::::::::::::  solution

## Solution

1. Your answer may vary, but here is one description: This is a bash shell script that can run for up to an hour on the itcga partition and reservation with 1 node, 2 CPUs, and 4 GB of memory. First, it loads three modules necessary to run cutadapt, then stores the paths that come after the name of the script on the command line as input_dir and output_dir. Then for each file in the input_dir that ends with `_R1_001_downsampled.fastq`, it runs cutadapt with it and its paired R2 to remove adapters, trim Ns from the ends of reads, and filter out any reads shorter than 25 bases after trimming. It stores the output with the file extension `001_ds_trim.fastq` in the output_dir.

2. You will need to modify the for loop call, the basename call, and the cutadapt call as follows:

```bash
for R1_fastq in ${input_dir}/*_1.fastq
```

```bash
name=$(basename ${R1_fastq} _1.fastq)
```

```bash
 cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25 \
 -o ${output_dir}/${name}_R1_001_ds_trim.fastq \
 -p ${output_dir}/${name}_R2_001_ds_trim.fastq \
 ${input_dir}/${R1_fastq} \
 ${input_dir}/${name}_2.fastq
 ```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

It is common that you might start your work with a shared script or code from a collaborator or someone else. This kind of exercise can save you a lot of headaches when trying to adapt it to your own workflow. 

Okay, now that you really understand this script, let's run it:

```bash
$ mkdir ~/1_project/data/trimmed_fastq
$ sbatch cutadapt.sh ~/1_project/data/untrimmed_fastq ~/1_project/data/trimmed_fastq
```

If you check the trimmed_fastq directory, you can see that your trimmed fastq files have been stored there!

You might also notice that your log files are kind of annoyingly piling up in various places. Think about possible ways to manage where they are stored!

We have now completed the trimming and filtering steps of our quality control process! Before we move on, let's move our trimmed FASTQ files to a new subdirectory within our `data/` directory.

:::::::::::::::::::::::::::::::::::::::  challenge

## Bonus exercise (advanced)

Now that our samples have gone through quality control, they should perform better on the quality tests run by FastQC. Go ahead and re-run FastQC on your trimmed FASTQ files and visualize the HTML files to see whether your per base sequence quality is higher after
trimming. Hint: you might need to modify your `fastqc.sh` script.

:::::::::::::::  solution

## Solution

There are a few ways to do this, but this one will work! Modify your `fastqc.sh` script like this:

```bash
# module load
module load fastqc-0.11.9-gcc-10.2.0-osi6pqc

# Define variables
input_dir=$1 

# run fastqc
fastqc ${input_dir}/*.fastq
```

Then you can run it by specifying the input directory on the command line.

```bash
$ sbatch fastqc.sh ~/1_project/data/trimmed_fastq/
```

After trimming and filtering, our overall quality is much higher, we have a distribution of sequence lengths, and more samples pass adapter content. However, you may want to explore the other options that [cutadapt](https://cutadapt.readthedocs.io/en/stable/) includes to refine your quality filtering.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- The options you set for the command-line tools you use are important!
- Data cleaning is an essential step in a genomics workflow.

::::::::::::::::::::::::::::::::::::::::::::::::::


