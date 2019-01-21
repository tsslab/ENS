#!/usr/bin/env python
import csv
import sys

promoter_file = open(sys.argv[1], 'rb')
rna_file = open(sys.argv[2], 'rb')
promoter_reader = csv.reader(promoter_file)
rna_reader = csv.reader(rna_file)

promoter_header = next(promoter_reader, None)
rna_header = next(rna_reader, None) 

rna_list = []
promoter_list = []

# Filter down your files to only contain the rows you care about.
for promoter_row in promoter_reader:
    if float(promoter_row[12]) < 0.05 and promoter_row[29] != 'NA': # remove all padj less than 0.05
        promoter_list.append(promoter_row)

#NB: Change number based on column number in file 
for rna_row in rna_reader:
    if rna_row[7] != 'NA' and float(rna_row[7]) < 0.05 and rna_row[8] != 'NA': # filter out all Padj value below 0.05
        rna_list.append(rna_row)


# Find common genes between the two files.
result = [] # empty result list that we will add things too.
for promoter_row in promoter_list:
    for rna_row in rna_list:
        if promoter_row[29] == rna_row[8]:
            promoter_row.extend(rna_row)
            result.append(promoter_row)

# Some string manipulation to make the file name different.
existing_file_elements = sys.argv[1].split('.')
del existing_file_elements[-1]
new_file_name = '.'.join(existing_file_elements) + '_results.csv'

print "saving results to:"
print new_file_name

new_file = open(new_file_name, 'wb')
writer = csv.writer(new_file, delimiter=',')
promoter_header.extend(rna_header) #add rna_header to promoter_header
writer.writerow(promoter_header)
for row in result:
    writer.writerow(row)

promoter_file.close()
rna_file.close()
new_file.close()
