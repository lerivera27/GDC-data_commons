# Task: Finding tumor mutation rates with the use of NCI Genomic Data Commons (GDC) Datasets. 
The goal is to process masked somatic mutation files in MAF format and analyze mutation rates in **TP53** and **STAT5A** genes.

## Background 

Somatic mutations are changes in DNA that occur after fertilization and may contribute to tumor development. Genes commonly mutated across multiple tumors (like **TP53**) often play crucial roles in cancer biology. TCGA provides data from 33 cancer subtypes! 


## Step 1: Install GDC Client 
### IMPORTANT MESSAGE: Please be sure to download the last version of the gdc-client and copy it into your PATH by using 'cp gdc-client ~/bin'. 

## Step 2: Search and download **Masked Somatic Mutation** MAF files for a specific cancer subtype and construct a repository for downloading the manifest file. 

##### link: https://portal.gdc.cancer.gov/

```
Project Name  = TCGA-X (X=cancer subtype of your choosing; Select under Projects tab)   #this script can be manipulated and tailored for any cancer subtype!
Then click the "Open cohort in repository" button.

Select these file types: 
Data Format  = maf (Select under Files tab)
Access = open
Data Type = Masked Somatic Mutation

Click the " Manifest" download button to retrieve the file

```
## Step 3: Download the MAF files into Linux command line: 
```  gdc-client download -m gdc-manifest.txt  ``` 

## Step 4: Move all '.maf.gz' files up a directory & uncompress them: 
```
for d in */; do cp "$d"*.maf.gz ./; done    #this will move all downloaded '.maf.gz' files to the current directory 
gunzip *.maf.gz    #uncompress all the '.maf' files

```
# Bash Script
```
#!/bin/bash

# Set the output report file path
report_file="mutation_report.txt"   #where results will be located in the output file


echo "Mutation Report for FIRSTNAME LASTNAME" > "$report_file"     #code to be used as example and needs to be replaced with your information 
echo "===================================" >> "$report_file"  
echo "" >> "$report_file"


compute_stats() {
    local values=("$@")     #function to compute mean and standard deviation
    local sum=0
    local count=${#values[@]}

    if [ "$count" -eq 0 ]; then
        echo "No data available."     #if no data is accurately downloaded this message will be echoed
        return
    fi

    for num in "${values[@]}"; do    #compute for the sum 
        sum=$((sum + num))
    done

    mean=$(echo "$sum / $count" | bc -l)     #compute for the mean

    sum_sq_diff=0
    for num in "${values[@]}"; do
        diff=$(echo "$num - $mean" | bc -l)
        sq_diff=$(echo "$diff^2" | bc -l)
        sum_sq_diff=$(echo "$sum_sq_diff + $sq_diff" | bc -l)    #compute for the standard deviation 
    done

    std_dev=$(echo "sqrt($sum_sq_diff / $count)" | bc -l)

    echo "Mean: $mean" >> "$report_file"
    echo "Standard Deviation: $std_dev" >> "$report_file"    #confirms results of all variables will be located in the output file 
    echo "" >> "$report_file"
}

tp53_counts=()     #initializes targets for storing specificed mutation counts (tp53 & stat5a) 
stat5a_counts=()

total_tp53_mutations=0
total_stat5a_mutations=0     #the variables used to track totals for summary 
total_files_processed=0

for file in *.maf; do    #processes each MAF file 
    if [ -f "$file" ]; then
        total_files_processed=$((total_files_processed + 1))
        echo "Processing file: $file"
    
        tp53_count=$(grep -c "TP53" "$file") 
        stat5a_count=$(grep -c "STAT5A" "$file")     #count mutation occurrences for tp53 & stat5a 

        total_tp53_mutations=$((total_tp53_mutations + tp53_count))       
        total_stat5a_mutations=$((total_stat5a_mutations + stat5a_count))     #used to add counts to total 

        tp53_counts+=("$tp53_count")
        stat5a_counts+=("$stat5a_count")

        # Append individual file results to the report file 
        echo "File: $file" >> "$report_file"
        echo "TP53 Mutations: $tp53_count" >> "$report_file"
        echo "STAT5A Mutations: $stat5a_count" >> "$report_file"
        echo "--------------------------------" >> "$report_file"     #lists individuals file results in the output file 
    fi
done 

echo "===== TP53 Mutation Statistics =====" >> "$report_file"   #initializes and confirms results for TP53 
compute_stats "${tp53_counts[@]}"

# Compute and append stats for STAT5A
echo "===== STAT5A Mutation Statistics =====" >> "$report_file"   #initializes and confirms results for STAT5A
compute_stats "${stat5a_counts[@]}"

echo "===== Summary of Findings =====" >> "$report_file"
echo "Total number of files processed: $total_files_processed" >> "$report_file"
echo "Total TP53 mutations: $total_tp53_mutations" >> "$report_file"
echo "Total STAT5A mutations: $total_stat5a_mutations" >> "$report_file"     #adds a fun summary of all findings at the end 

mv "$report_file" /leiarar/scratch/tcga-brca-maf    #used to ensure the report_file will be in the current working directory 
echo "Report saved to /leiarar/scratch/tcga-brca-maf/mutation_report.txt"    #the path can be tailored to your allocated path
``` 
## The Summary of findings will look like this: 

```
===== TP53 Mutation Statistics =====
Mean: .39593114241001564945
Standard Deviation: .55500203861271185816

===== STAT5A Mutation Statistics =====
Mean: .01408450704225352112
Standard Deviation: .13044548346901224885

===== Summary of Findings =====
Total number of files processed: 639
Total TP53 mutations: 253
Total STAT5A mutations: 9

```

## License - this project is licensed under the MIT License.

### Contact information: 
Leiara Rivera 
* leiarar@clemson.edu 
* lrivera@email.sc.edu 







