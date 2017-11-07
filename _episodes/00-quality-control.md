---
title: "Assessing Read Quality"
teaching: 0
exercises: 0
questions:
- "Is my data good enough?"
objectives:
- "Describe how the FASTQ format encodes quality."
- "Evaluate a FastQC report."
keypoints:
- "First key point."
---

# Bioinformatics workflows

When working with high-throughput sequencing data, the raw reads you get off of the sequencer will need to pass
through a number of  different tools in order to generate your final desired output. The execution of this set of
tools in a specified order is commonly referred to as a *workflow* or a *pipeline*. 

An example of the workflow we will be using for our variant calling analysis is provided below with a brief
description of each step. 

![workflow](../img/variant_calling_workflow.png)


1. Quality control - Assessing quality using FastQC
2. Quality control - Trimming and/or filtering reads (if necessary)
3. Align reads to reference genome 
4. Perform post-alignment clean-up
5. Variant calling

These workflows in bioinformatics adopt a plug-and-play approach in that the output of one tool can be easily
used as input to another tool without any extensive configuration. Having standards for data formats is what 
makes this feasible. Standards ensure that data is stored in a way that is generally accepted and agreed upon 
within the community. The tools that are used to analyze data at different stages of the workflow are therefore 
built under the assumption that the data will be provided in a specific format.  


# Quality Control

The first step in the variant calling work flow is to take the FASTQ files received from the sequencing facility
and assess the quality of the sequence reads. 

![workflow_qc](../img/var_calling_workflow_qc.png)
## Details on the FASTQ format

Although it looks complicated (and maybe it is), its easy to understand the
[fastq](https://en.wikipedia.org/wiki/FASTQ_format) format with a little decoding. Some rules about the format
include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

We can view the first complete read in one of the files our dataset by using `head` to look at
the first four lines.

~~~
$ head -n4 SRR098281.fastq 
~~~
{: .bash}

~~~
@SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
~~~
{: .output}

All but one of the nucleotides in this read are unknown (`N`). This is a pretty bad read!

Line 4 shows the quality for each nucleotide in the read. Quality is interpreted as the 
probability of an incorrect base call (eg 1 in 10) or, equivalently, the base call 
accuracy (eg 90%). To make it possible to line up each individual nucleotide with its quality
score, the numerical score is converted into a code where each individual character 
represents the numerical quality score for an individual nucleotide. For example, in the line
above, the quality score line is: 

~~~
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
~~~
{: .output}

The `#` character and each of the `!` characters represent the encoded quality for an 
individual nucleotide. The numerical value assigned to each of these characters depends on the 
sequencing platform that generated the reads. The sequencing machine used to generate our data 
uses a version of this encoding called Illumina 1.8 encoding. Each character is assigned a
quality score between 0 and 40 as shown in the chart below.

~~~
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                  |         |         |         |         |
Quality score:    0........10........20........30........40                                
~~~
{: .output}

Each quality score represents the probability that the corresponding nucleotide call is
incorrect. This quality score is logarithmically based, so a quality score of 10 reflects a
base call accuracy of 90%, but a quality score of 20 reflects a base call accuracy of 99%. 
These probability values are the results from the base calling algorithm and dependent on how 
much signal was captured for the base incorporation. 

Looking back at our read: 

~~~
@SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
~~~
{: .output}

we can now see that the quality of each of the `N`s is 0 and the quality of the first 
nucleotide call (`C`) is also very poor (`#` = a quality score of 2). This is indeed a very
bad read. 


> ## Quality Encodings Vary
>
> Although we've used a particular quality encoding system to demonstrate interpretation of 
> read quality, different sequencing machines use different encoding systems. This means that, 
> depending on which sequencer you use to generate your data, a `#` may not be an indicator of 
> a poor quality base call. 
> 
> It's essential that you know which sequencing platform was
> used to generate your data, so that you can tell your quality control program which encoding
> to use. If you choose the wrong encoding, you run the risk of throwing away good reads or 
> (even worse) not throwing away bad reads!
{: .callout}

## Assessing Quality using FastQC
In real life, you won't be assessing the quality of your reads by visually inspecting your 
FASTQ files. Rather, you'll be using a software program to assess read quality and 
filter out poor quality reads. We'll first use a program called [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to visualize the quality of our reads. 
Later in our workflow, we'll use another program to filter out poor quality reads. 

FastQC has a number of features which can give you a  quick impression of any problems your
data may have, so you can take these issues into consideration before moving forward with your
analyses. Rather than looking at quality scores for each individual read, FastQC looks at
quality collectively across all reads within a sample. The image below shows a FastQC-generated plot that indicates
a very high quality sample:

![good_quality](../img/good_quality.png)

The x-axis displays the base position in the read, and the y-axis shows quality scores. In this
example, the sample contains reads that are 40 bp long. For each position, there is a 
box-and-whisker plot showing the distribution of quality scores for all reads at that position.
The horizontal red line indicates the median quality score and the yellow box shows the 2nd to
3rd quartile range. This means that 50% of reads have a quality score that falls within the
range of the yellow box at that position. The whiskers show the range to the 1st and 4th 
quartile.

For each position in this sample, the quality values do not drop much lower than 32. This 
is a high quality score. The plot background is also color-coded to identify good (green),
acceptable (yellow), and bad (red) quality scores.

Now let's take a look at a quality plot on the other end of the spectrum. 

![bad_quality](../img/bad_quality.png)

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the "bad" range, particularly on the tail end of the reads. The FastQC tool produces several other diagnostic plots to assess sample quality, in addition to the one plotted above. 

## Running FastQC  

We will be working with a set of sample data that is located in a hidden directory (`.dc_sampledata_lite`). First, we
will move some of these hidden files to the `data` directory your created at [the end of our
last lesson](http://www.datacarpentry.org/shell-genomics/06-organization/).  

~~~
$ mv ~/.dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/
~~~
{: .bash}

Navigate to your FASTQ dataset: 

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq/
~~~
{: .bash}

> ## Exercise
> 
> How many FASTQ files are in this dataset? How big are the files?  
>
>> ## Solution
>>  
> {: .solution}
{: .challenge}

To run the FastQC program, we need to tell our computer where the program is located 
(in `~/FastQC`).  FastQC can accept multiple file names as input, so we can use the *.fastq wildcard to run FastQC on all of the FASTQ files in this directory.

~~~
$ ~/FastQC/fastqc *.fastq
~~~
{: .bash}

This command created several new files within our `/data/untrimmed_fastq/` directory.

~~~
$ ls
~~~
{: .bash}

~~~
~~~
{: .output}

> ## HTML and Zip files
>
>
>
{: .callout}

However, we want to keep our data files and our results files separate, so we will move these
output files into a new directory within our `results/` directory.

~~~
$ mkdir ~/dc_workshop/results/fastqc_untrimmed_reads
$ mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
$ mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/
~~~
{: .bash}


## Results

Let's examine the results in detail

Navigate to the results and view the directory contents


    $ cd ~/dc_workshop/results/fastqc_untrimmed_reads/
    $ ls

The zip files need to be unpacked with the `unzip` program. 

Use `unzip` to unzip the FastQC results: 

    $ unzip *.zip

Did it work? No, because `unzip` expects to get only one zip file.  Welcome to the real world. We *could* do each file, one by one, but what if we have 500 files?  There is a smarter way. We can
save time by using a simple shell `for` loop to iterate through the list of files in `*.zip`. After you type the first line, you will get a special `>` prompt to type next next lines. You start with
`do`, then enter your commands, then end with `done` to execute the loop.


Build a `for` loop to unzip the files
    

    $ for zip in *.zip
    > do
    > unzip $zip
    > done

Note that, in the first line, we create a variable named `zip`.  After that, we call that variable with the syntax `$zip`.  `$zip` is assigned the value of each item (file) in the list `*.zip`, once for
each iteration of the loop.

This loop is basically a simple program.  When it runs, it will run `unzip` once for each file (whose name is stored in the `$zip` variable). The contents of each file will be unpacked into a separate
directory by the `unzip` program.

The `for` loop is interpreted as a multipart command.  If you press the up arrow on your keyboard to recall the command, it will be shown like so:
   
    $ for zip in *.zip; do echo File $zip; unzip $zip; done

When you check your history later, it will help your remember what you did!

## D. Document your work

To save a record, let's `cat` all fastqc `summary.txt`s into one file `full_report.txt` and move this to ``~/dc_workshop/docs``. You can use wildcards in paths as well as file names.  Do you remember how we
said `cat` is really meant for concatenating text files?

    cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt






