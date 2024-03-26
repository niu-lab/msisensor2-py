import pysam
import numpy as np
from datetime import datetime
import multiprocessing
from collections import Counter
from scipy.stats import chi2_contingency
import getopt
import sys 

class Bam2Dis:
    def __init__(self,list_file,normal_bam_input,tumor_bam_input):
        self.list_file = list_file
        self.normal_bam_input = normal_bam_input
        self.tumor_bam_input = tumor_bam_input

    def read_file_in_batches(self,batch_size=100):
        # Batch processing code, batch_size is the number of bits processed each time, and how many lines are read

        with open(self.list_file,'r') as file:
            while True:
                batch = []
                for _ in range(batch_size):
                    line = file.readline()
                    if not line:
                        break
                    batch.append(line)

                if not batch:
                    break

                yield  batch


    def process_one_batch(self,one_batch):
        """
        This function processes the site information of each batch. If a single or multiple base repeat condition is not met at a certain site, it will be skipped; If the sequencing depth is less than 20, skip it;
        :param one_batch: Batch_size sites in an MS site file
        :return: Frequency matrix
        """

        normal_bam_input = self.normal_bam_input
        normal_bam_input_class = pysam.AlignmentFile(normal_bam_input, 'r')

        tumor_bam_input = self.tumor_bam_input
        tumor_bam_input_class = pysam.AlignmentFile(tumor_bam_input,'r')

        valid_chromosomes = normal_bam_input_class.references and tumor_bam_input_class.references

        tmp_output_list = []
        tmp_p_value_list = []

        for one_loc in one_batch:

            output_loc = one_loc

            one_loc = one_loc.strip().split('\t')
            chromosome = one_loc[0]

            repeat_unit_length = int(one_loc[2])
            repeat_unit_times = int(one_loc[4])

            repeat_unit = one_loc[-3]
            left_flank_bases = one_loc[-2]
            right_flank_bases = one_loc[-1]

            loc_begin = int(one_loc[1])
            # Calculate the starting and ending sites, and retrieve all reads covering this region
            loc_end = loc_begin+repeat_unit_times*repeat_unit_length

            # Filter position points based on the parameters in MSIsensor
            if (repeat_unit_length==1 and ((repeat_unit_times<10) or (repeat_unit_times>50))):
                continue

            if (repeat_unit_length>1 and ((repeat_unit_times<5) or (repeat_unit_times>40))):
                continue

            # Check if a chromosome is in the bam file, if not, skip
            if chromosome not in valid_chromosomes:
                continue

            # Used for debugging
            # if loc_begin != 11182261:
            #     continue

            # Skip directly for values less than 20
            if normal_bam_input_class.count(chromosome,loc_begin,loc_end) < 20:
                continue

            if tumor_bam_input_class.count(chromosome,loc_begin,loc_end) < 20:
                continue



            # print(normal_bam_input_class.count(chromosome,loc_begin,loc_end))

            # break
            normal_bam_iterator = normal_bam_input_class.fetch(chromosome, loc_begin, loc_end)
            tumor_bam_iterator = tumor_bam_input_class.fetch(chromosome,loc_begin,loc_end)

            tumor_frequency = [0 for i in range(100)]

            normal_frequency = [0 for i in range(100)]

            tumor_frequency = np.array(self.iterate_reads(normal_bam_iterator,tumor_frequency,repeat_unit,left_flank_bases,right_flank_bases))

            normal_frequency = np.array(self.iterate_reads(tumor_bam_iterator,normal_frequency,repeat_unit,left_flank_bases,right_flank_bases))

            if np.sum(tumor_frequency)<20 or np.sum(normal_frequency)<20:
                continue

            # print(tumor_frequency)
            # print(normal_frequency)

            tumor_index = np.nonzero(tumor_frequency)
            normal_index = np.nonzero(normal_frequency)
            total_index = np.union1d(tumor_index[0], normal_index[0])

            # print(tumor_frequency[total_index])
            # print(normal_frequency[total_index])

            df = [tumor_frequency[total_index], normal_frequency[total_index]]

            # print(df)

            stat, p, dof, expected = chi2_contingency(df, correction=False)

            # print(p)


            tmp_output_list.append(output_loc.strip('\n'))
            tmp_output_list.append('N: '+' '.join(map(str,normal_frequency)))
            tmp_output_list.append('T: '+ ' '.join(map(str,tumor_frequency)))

            tmp_p_value_list.append(p)
            counter = 0

        return tmp_output_list,tmp_p_value_list

            # print(tumor_frequency)
            # break
    def count_unit_repeats(self, read_string, repeat_unit, fbases, ebases):
        """
        The function calculates the number of repeated units in a read. If it can perfectly match the left and right wing sequences, it returns the number of repetitions in the MS region between the two sequences
        :param read_string: Enter the AGTC string for read
        :param repeat_unit: Repeating unit content
        :param fbases: left flank bases
        :param ebases: right flank bases
        :return: The size of repeating units within the read
        """
        start_pos = 0  # Initialize the starting position

        while True:
            start_pos = read_string.find(fbases, start_pos)

            if start_pos == -1:
                break

            count = 0
            # It should be noted that this does not skip the first repeating unit, but rather the left wing sequence
            tstart = start_pos + len(fbases)

            # Find duplicate units at the current position. If found, return to the current position and continue the loop. If not found, return negative one and jump out of the loop
            while tstart == read_string.find(repeat_unit, tstart):

                count = count + 1
                tstart += len(repeat_unit)

                # Once the right flank sequence is equal to the current position, it jumps out of the loop and captures the MSI information here
                tstart0 = read_string.find(ebases, tstart)
                if tstart == tstart0:
                    if ((len(repeat_unit)==1) and (count>=5)) or (len(repeat_unit)>1 and (count>=3)):
                        return count
                    else:
                        # print(count)
                        return 0

            start_pos += 1

        return 0

    def iterate_reads (self,bam_data,frequency,repeat_unit, left_flank_bases, right_flank_bases):
        """
        Iterating the reads of a certain site to return the frequency information of a certain site
        :param bam_data:
        :param frequency:
        :return:
        """

        reads_counter = 0

        for i in bam_data:
            if i.is_unmapped:
                continue

            read_frequency = self.count_unit_repeats(i.seq, repeat_unit, left_flank_bases, right_flank_bases)

            read_frequency = read_frequency - 1

            if read_frequency == -1:
                continue

            frequency[read_frequency] = frequency[read_frequency] + 1
            reads_counter = reads_counter + 1

        # print(reads_counter)

        return frequency

#def highlight_substring(input_string, substring_1,substring_2):
#    # Debugging function, input_string is the input string, substring_1 and substring_2 are the left flank sequence and the right flank sequence, respectively
#    if (substring_1 not in input_string) and (substring_2 not in input_string):
#        #print(f"Substring '{substring_1}' not found in the input string.")
#        return input_string
#
#    highlighted_string = input_string.replace(substring_1, f'\033[1;31m{substring_1}\033[0m')
#    highlighted_string = highlighted_string.replace(substring_2, f'\033[1;34m{substring_2}\033[0m')
#    return highlighted_string

def check(argv):
    ## HELP
    def Usage():
        print("\nMSIsensor python version 1.0")
        print("\nUsage:")
        print('''
        python msisensor.py -s sites_list -t input.tumor.bam -n input.normal.bam -p output.p.value -o output.dis
        
        Option:
            -s, --site-list     <file>  the MSI sites file
            -t, --tumor-bam     <file>  the tumor bam input file
            -n, --nornal-bam    <file>  the normal bam input file
            -p, --p-value-file  <file>  the output of p value
            -o, --output-dis    <file>  the output of distribution for MSI sites
            
            -h, --help  help
        ''')
    try:
        opts, args = getopt.getopt(argv, "s:t:n:p:o:h:", ["site-list", "tumor-bam", "nornal-bam", "p-value-file", "output-dis", "help"])
    except getopt.GetoptError:
        #print("\n[Error] Unknown parameter, please check.")
        Usage()
        sys.exit(-3)
    ## parameters parsing
    site_list = None
    tumor_bam = None
    normal_bam = None
    p_value_file = None
    output_dis = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
            sys.exit(-4)
        elif opt in ("-s", "--site-list"):
            site_list = arg
        elif opt in ("-t", "--tumor-bam"):
            tumor_bam = arg
        elif opt in ("-n", "--normal-bam"):
            normal_bam = arg
        elif opt in ("-p", "--p-value-file"):
            p_value_file = arg
        elif opt in ("-o", "--output-dis"):
            output_dis = arg
        elif opt in ("-h", "--help"):
            Usage()
            sys.exit(-3)

    bam2dis = Bam2Dis(site_list, tumor_bam, normal_bam)
    pool = multiprocessing.Pool(processes=40)

    i = 0
    start_time = datetime.now()
    print('Start time: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    total_test_locs = 33793502
    batch_size = 100000
    break_iterator = total_test_locs/batch_size

    lock = multiprocessing.Lock()

    this_start_time = datetime.now()

    multiprocess_results = []
    result_value_list = []
    p_value_dict = {}

    for batch in bam2dis.read_file_in_batches(batch_size=batch_size):
        i=i+1

        # if i<5000:
        #     continue
        # bam2dis.process_one_batch(batch,output_file)
        result = pool.apply_async(bam2dis.process_one_batch,args=(batch,))
        multiprocess_results.append(result)

        # if i == 20:
        #     break


    pool.close()
    pool.join()


    for result in multiprocess_results:

        result_value,p_value = result.get()
        if len(result_value) == 0:
            continue
        for one_line in result_value:
            result_value_list.append(one_line)
        #print(result_value) # Optional output of p-value result
        #print(p_value) # Optional output of p-value result

        for line_count in range(0,len(result_value),3):
            p_value_dict[result_value[line_count]] = p_value[line_count//3]
        # print(result_value)
        # break

    #print(len(result_value_list)/3) # Optional output of p-value result
    # for i in result_value_list:
    #     print(i)

    with open(output_dis,'w') as file:
        for dis_line in result_value_list:
            #file.write(str(dis_line)+'\n') # This is the standard format for MSI software, but to adapt to the input format of machine learning prediction models, the following modifications are made
            if dis_line.startswith('T') or dis_line.startswith('N'):
                file.write(str(dis_line)+'\n')
            else:
                columns = dis_line.split()
                file.write('ch' + columns[0] + '\t' + columns[1] + '\t' + columns[8] + '\t' + columns[4] + '[' + columns[7] + ']' +'\t' +columns[9] + '\n' ) # chromosome location left_flank_bases repeat_times[repeat_unit_bases] right_flank_bases

    with open(p_value_file,'w') as file:
        for key,value in p_value_dict.items():
            file.write(str(key)+": ")
            file.write(str(value)+"\n")
    print('End time: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    end_time = datetime.now()
    print('Total running time: ' + str(end_time - start_time))

def main(argv):
    check(argv)

if __name__ == '__main__' or (__name__ == 'main'):
    if len(sys.argv) > 1:
        main(sys.argv[1:])
        processes = []
    else:
        print('''
Please enter parameters and execute the following command to view help:
    python3 msisensor.py -h
        ''') 
