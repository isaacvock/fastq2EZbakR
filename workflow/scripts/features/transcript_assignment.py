'''
    File name: count_reads_only.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    Orcid: 
    Date created: 10/3/2023
    Date last modified: 10/3/2023
    Version: 1.0.0
    License: GPLv3
    Python Version: 
    Packages: pysam 
    Description: Script for creating table of transcripts read mapped to
'''

import os
import pysam
import csv
import argparse
import datetime
import sys

# # Proper logging
# with open(snakemake.log[0], "w") as f:
#     sys.stderr = sys.stdout = f

# Get input bam file
input = snakemake.input.bam
output = snakemake.output.table

#  Set .csv file for writing (simulating _counts.rds file)
header = ['qname', 'bamfile_transcripts']
myfile = open(output, 'w', newline='')
wr = csv.writer(myfile)
wr.writerow(header)

# Set .bam file for reading
samfile = pysam.AlignmentFile(input, 'rb')

print('Started analysis')

### I need to concatenate the names of all transcripts
### that a read could have mapped to and then write
### that to the file once I reach a new read name
readcount = 1
for r in samfile:

    # Name of first read
    if readcount == 1:
        
        # Save name of read to check if still on this particular read
        queryname = r.query_name

        # hack to make sure that this is only done for first read in bam file
        readcount = 2

        # Initialize list of read info
        r_info = 2*['']

        # Initialize set of transcripts read maps to
        transcripts = {str(r.reference_name)}

        # Save the read name
        r_info[0] = r.query_name
        

    # If still on the particular read
    elif queryname == r.query_name:

        # Check to see if RSEM bam file is provided
        # and thus if probability of assignment
        # can be used
        if r.has_tag("ZW"):

            if r.get_tag("ZW") != 0:

                # Add transcript to set of transcripts
                transcripts.add(str(r.reference_name))


        else:

            # Add transcript to set of transcripts
            transcripts.add(str(r.reference_name))

    # If on to the next read
    else:
        
        # Concatenate set of transcripts into single string
            # Sorted to ensure reproducibility of final string
        r_info[1] = '+'.join(sorted(transcripts))

        # Write to csv
        wr.writerow(r_info)

        # Save new read name
        queryname = r.query_name

        # Save read name to what will be another row of output table
        r_info[0] = r.query_name            
        
        # Initialize set of transcripts
        transcripts = {str(r.reference_name)}



if readcount > 1:
    r_info[1] = '+'.join(sorted(transcripts))
    wr.writerow(r_info)

samfile.close()
print('Ended analysis')

##### Close files ######
myfile.close()