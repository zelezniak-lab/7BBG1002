---
title: "7BBG1002: Using BLAST and HMMER from the Command Line Tutorial"
author: Aleksej Zelezniak
output: html_document
editor_options: 
  chunk_output_type: inline
  markdown: 
    wrap: sentence
---


### Overview

```         
This tutorial is designed to get hands-on experience with two powerful sequence alignment tools: BLAST (Basic Local Alignment Search Tool) and HMMER. You will learn how to run these tools from the command line, understand their output, and see how to combine them for deeper sequence analysis.
```

#### Prerequisites

-   Basic understanding of sequence alignment concepts.
-   Familiarity with the command line and Linux commands.
-   Familiarity with CREATE KCL HPC system, e.g. account is required to be set

### Objectives

- Build custom BLAST databases
- Run BLAST from a command line
- Filter BLAST results using command line tools
- Build HMM profiles and use them to improve sequence detection, functional prediction to search metagenomes.

### 1. Setting Up Your Environment

1.  **Install BLAST and HMMER**:

    -   Please login to the CREATE note

        ``` sh
        srun -p cpu --time=03:00:00 --pty /bin/bash # use your assigned partition instead of cpu
        ```

    -   Navigate to your scratch space

        ``` sh
        cd /scratch/users/kYOURNUMBER/
        ```

    -   Using Conda to install BLAST and HMMER:

        ``` sh
        module load anaconda3/2022.10-gcc-13.2.0
        conda create -n alignments
        conda init bash 
        source ~/.bashrc
        conda activate alignments
        conda install -c bioconda blast hmmer seqtk mafft alv
        ```

2.  **Prepare Your Workspace**:

    -   Create a new directory to store the files for this tutorial:

        ``` sh
        mkdir alignments
        cd alignments
        ```

### 2. Building a BLAST Database

Before running a BLAST search, it's important to have a database against which to compare your query sequence.
In this section, we will build a local BLAST database using Swiss-Prot sequences, see [documentation](https://www.uniprot.org/uniprotkb).
The sequences were downloaded for you but you can access them [here](https://www.uniprot.org/help/downloads).

1.  **Obtain Sequences for the Database**:

    
    -   For this example, let's look a FASTA file named `uniref50.fasta`:
    
        ``` sh
        less /scratch/grp/msc_appbio/alignments/uniprot_sprot.fasta 
        ```

    -   Q: How many sequences are in the database?
        Answer: <font color="white"> 573661 </font>

2.  **Building the BLAST Database**:

    -   The following the `makeblastdb` command creates a BLAST database from the FASTA file.
        
        ``` sh
        makeblastdb -in /scratch/grp/msc_appbio/alignments/uniprot_sprot.fasta -dbtype prot -out uniprot_sprot_db
        ```

    -   **Options Explained**:

        -   `-in /scratch/grp/msc_appbio/alignments/uniprot_sprot.fasta`: The input FASTA file.
        -   `-dbtype prot`: Specifies that the database is protein sequences.
        -   `-out uniprot_sprot_db`: The name of the output database.

        
### 3. Running BLAST from the Command Line

BLAST is commonly used to identify sequences that are similar to a query sequence in a database.
In this section, we will conduct a BLAST search using a protein query sequence.

1.  **Obtain a Query Sequence**:

    -   For this example, let's use a sample FASTA file named `query.fasta`.
        You can create it using:

        ``` sh
        nano query.fasta
        ```

    -   Paste the following sequence (this is PETase from Ideonella sakaiensis, Yoshida et al 2016, Science), write to file (control + W) and exit nano:

        ```         
        >A0A0K8P8E7
        MQTTVTTMLLASVALAACAGGGSTPLPLPQQQPPQQEPPPPPVPLASRAACEALKDGNGDMVWPNAATVVEVAAW
        RDAAPATASAAALPEHCEVSGAIAKRTGIDGYPYEIKFRLRMPAEWNGRFFMEGGSGTNGSLSAATGSIGGGQIA
        SALSRNFATIATDGGHDNAVNDNPDALGTVAFGLDPQARLDMGYNSYDQVTQAGKAAVARFYGRAADKSYFIGCS
        EGGREGMMLSQRFPSHYDGIVAGAPGYQLPKAGISGAWTTQSLAPAAVGLDAQGVPLINKSFSDADLHLLSQAIL
        GTCDALDGLADGIVDNYRACQAAFDPATAANPANGQALQCVGAKTADCLSPVQVTAIKRAMAGPVNSAGTPLYNR
        WAWDAGMSGLSGTTYNQGWRSWWLGSFNSSANNAQRVSGFSARSWLVDFATPPEPMPMTQVAARMMKFDFDIDPL
        KIWATSGQFTQSSMDWHGATSTDLAAFRDRGGKMILYHGMSDAAFSALDTADYYERLGAAMPGAAGFARLFLVPG
        MNHCSGGPGTDRFDMLTPLVAWVERGEAPDQISAWSGTPGYFGVAARTRPLCPYPQIARYKGSGDINTEANFACA
        APP
        ```

2.  **Running a BLAST Search**:

    -   Use the BLASTP program to search against the custom database we just created:

        ``` sh
        blastp -query query.fasta -db uniprot_sprot_db -out blast_results.txt -evalue 1e-5 -outfmt 6
        ```

    -   **Options Explained**:

        -   `-query query.fasta`: The query sequence file.
        -   `-db uniprot_sprot_db`: The uniprot_sprot_db database to search against.
        -   `-out blast_results.txt`: Output file for the results.
        -   `-evalue 1e-5`: Set the E-value threshold for reporting matches.
        -   `-outfmt 6`: Output in tabular format (easy to parse).

3.  **Reviewing the Output**:

    -   Open the results file:

        ``` sh
        less blast_results.txt
        ```

    -   The output format 6 is a tab-delimited table with the following fields: query ID, subject ID, % identity, alignment length, mismatches, gap opens, query start, query end, subject start, subject end, E-value, and bit score.

    -   Q: how many significant E-value \< 1e-6 are out there?
        Answer: <font color="white"> awk '\$11 \< 1e-6' blast_results.txt \| wc -l </font>

    -   Extract significant hits with E-value \< 1e-6 a and identity \>40%

        ``` sh
        awk '$11 < 1e-6' blast_results.txt > significant_blast_hits.txt
        ```

### 4. Running HMMER from the Command Line

HMMER is a tool used for searching sequence databases with profile hidden Markov models (HMMs).
It is particularly useful for finding more remote homologs compared to BLAST.

1.  **Preparing a Multiple Sequence Alignment (MSA)**:

    -   To create a profile HMM, you need an MSA file in FASTA format.
        Let's start first from extracting identified sequences using very useful and lightweight `seqtk` utility for sequence manipulation.
        You can find more info about `seqtk` [here](https://github.com/lh3/seqtk).
        First extract the hits IDs.
        ```sh
        cat significant_blast_hits.txt | cut -f2 > significant_blast_hits.ids
        
        ```
        Then filter out from huge database file
        ``` sh
        seqtk subseq /scratch/grp/msc_appbio/alignments/uniprot_sprot.fasta significant_blast_hits.ids > blast_hits.fasta
        ```

    -   Now let's aligned the extracted sequences, using [MAFFT](https://mafft.cbrc.jp/alignment/software/).

        ``` sh
        mafft  blast_hits.fasta > alignment.fasta
        ```

    -   You can use to alignemnt viewer tools like [Jalview](https://www.jalview.org/) to visualise the alignments.
        We will use handy [`alv`](https://github.com/arvestad/alv) utility to view alignments in terminal.

        ``` sh
        alv -k alignment.fasta  | less -R
        ```

2.  **Building a Profile HMM**:

    -   Once you have the MSA, use the following command to build a profile HMM:

        ``` sh
        hmmbuild profile.hmm alignment.fasta
        ```

    -   **Options Explained**:

        -   `profile.hmm`: The output HMM profile file.
        -   `alignment.fasta`: The input MSA file.

3.  **Searching with HMMER**:

    -   Use `hmmsearch` to search for sequences that match your profile HMM in a target database.
        In our case as a database we will use metagenome sample collected as a part of [Tara Ocean Project](https://fondationtaraocean.org/en/home/).
        Let's first access the sample data, stored at [MGnify](https://www.ebi.ac.uk/metagenomics/analyses/MGYA00590504#download) database.
        We are looking at the study MGYS00002008, sample ERS494579, assembly ERZ841232.
        To safe time, the file is already dowloaded for you using the following command, **Please don't run**: 
        
        ```sh
        #wget https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00590504/file/ERZ841232_FASTA_predicted_cds.faa.gz
        ```
        
        Instead let's search the file directly using produced `profile.hmm`
        
    -   The command is:

        ``` sh
        hmmsearch --tblout hmm_results.txt profile.hmm /scratch/grp/msc_appbio/alignments/ERZ841232_FASTA_predicted_cds.faa
        ```

    -   **Options Explained**:

        -   `--tblout hmm_results.txt`: Save output in a tabular format.
        -   `profile.hmm`: The HMM profile that you made.
        -   `RZ841232_FASTA_predicted_cds.faa`: The database file to search.

4.  **Reviewing the Output**:

    -   Open the results file:

        ``` sh
        cat hmm_results.txt
        ```

    -   The output includes fields such as target name, accession, query name, E-value, and score.

### 5. Conclusion

In this tutorial, you learned how to use BLAST and HMMER from the command line for sequence similarity searches and how to combine the results to identify significant homologs.
BLAST is fast and efficient for detecting close homologs, while HMMER excels at finding more distant relationships due to its use of profile HMMs.

### 7. Further Reading

-   [BLAST Documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
-   [HMMER User Guide](http://eddylab.org/software/hmmer/Userguide.pdf)

### 8. Exercises

1.  Build an HMM from a different alignment and see how the results change when using it to search a larger database.
2.  Write a script to automate the combination of BLAST and HMMER results.
3.  Build a custom BLAST database with more sequences and test how the database size affects the search results.

Feel free to reach out with any questions or if you need further guidance on the exercises!

