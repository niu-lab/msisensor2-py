MSIsensor2 (Python version)
===========
MSIsensor2 is a novel algorithm based machine learning, featuring a large upgrade in the microsatellite instability (MSI) detection for tumor only sequencing data, including Cell-Free DNA (cfDNA), Formalin-Fixed Paraffin-Embedded(FFPE) and other sample types. The original MSIsensor is specially designed for tumor/normal paired sequencing data.


Usage
-----

        Version 0.1
        Usage:  python msisensor2.py [input_sample_dis_file_name] [input_sites_list] [input_models_dir] >> [output_file]

Parameters:
-------

        input_sample_dis_file_name <input>   the name of sample that need to be predicted
        input_sites_list           <input>   list of training sites corresponding to all training site models
        input_models_dir           <input>   directory for training models
        output_file                <output>  output file of model prediction results


Install:
-------
You may already have these prerequisite packages. 
````
1. Python3.5 or above;
2. Python package:
   xgboost, pysam, multiprocessing, numpy, collections, datetime, scipy.stats, getopt, sys, etc.
````


Example
-------
1. Generate distribution files for predicted samples (that is input_sample_dis_file_name) using MSIsensor(python version)  for example:
````
chr1 30867 TCTCC 12[CT] CATTT
N: 0 0 0 0 0 1 0 0 0 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
T: 0 0 0 0 0 0 0 2 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
````
````
python msisensor.py -s test/example.microsate.sites -t test/example.tumor.bam -n test/example.normal.bam -p test/output_p_value.txt -o test/output_dis_file.txt
````

Notes: usage of msisensor.py
````
Usage:

        python msisensor.py -s sites_list -t input.tumor.bam -n input.normal.bam -p output.p.value -o output.dis

        Options:
            -s, --site-list     <input>   the MSI sites file
            -t, --tumor-bam     <input>   the tumor bam input file
            -n, --nornal-bam    <input>   the normal bam input file
            -p, --p-value-file  <output>  the output of p value
            -o, --output-dis    <output>  the output of distribution for MSI sites

            -h, --help  help

````

2. Generate sites_list.txt after deleting X, Y, and Un chromosomes:
````
cat model_dir/sites_list.txt |grep -v "chrX"|grep -v "chrY"|grep -v "chrUn"|grep -v "random" > ./sites_list.txt
echo -e "input_sample_dis_file_name \c" >> output.txt
````

3. Predict for single sample:
````
python msisensor2.py input_sample_dis_file_name ./sites_list.txt ./model_dir/models/ >> output.txt
````

4. Batch submit multiple sample jobs:
````
chmod +x submit_test.sh
sh submit_test.sh
````

Output
-------

output.txt: The result of model prediction for samples: input_sample_dis_file_name  MSIscore  uns_num  sta_num
````
TCGA-3L-AA1B_dis-MSS 0.0356259277585 72 1949
TCGA-4N-A93T_dis-MSS 0.0380085653105 71 1797
TCGA-5M-AAT4_dis-MSS 0.0495565988524 95 1822
TCGA-5M-AAT5_dis-MSI-L 0.041188738269 79 1839
TCGA-5M-AAT6_dis-MSI-H 0.683049147442 1362 632
TCGA-5M-AATA_dis-MSS 0.0334613274822 61 1762
````

Test sample
-------
We provided one small dataset (tumor only bam file) to test the msi scoring step:
````
source /your_dir/.bashrc; sh test/submit_test.sh TCGA-CRC
````
The test result TCGA-CRC.txt for TCGA-CRC is in ./test.


Contact
-------
If you have any questions, please contact one or more of the following folks: Beifang Niu <bniu@sccas.cn> or Ruilin Li <lirl@sccas.cn>
