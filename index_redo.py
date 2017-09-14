import argparse
import itertools
from collections import defaultdict

#use arg parse to get args
def get_arguments():
    parser = argparse.ArgumentParser(description="bring in sequence files")
    parser.add_argument("-f2", help="set R2 file", required=True, type=str)
    parser.add_argument("-f3", help="set R3 file", required=True, type=str)
    parser.add_argument("-fi", help="set index file", required=True, type=str)
    return parser.parse_args()

args = get_arguments()

#set variables for input files
f2 = args.f2
f3 = args.f3
fi = args.fi

def rComplement(DNA):
    r = ""

    for x in DNA:
        if x == 'A':
            r = 'T' + r
        elif x == 'C':
            r = 'G' + r
        elif x == 'G':
            r = 'C' + r
        elif x == 'T':
            r = 'A' + r
        elif x == 'N':
            r = 'N'

    return r

#read in indexes to list
indexlist = []

with open(fi) as fileindex:
    for line in fileindex:
        line = line.strip()
        parts = line.split()
        indexlist.append(parts[1])


#create all combinations of indexes as keys in a dictionary
indexlist.sort()
comboindex = itertools.permutations(indexlist, 2)

#add combos to dict
indexdict = defaultdict(int)
for c in comboindex:
    indexdict[c[0] + "_" + c[1]] = 0

#add repeated indexes to dict
for i in indexlist:
    indexdict[i + "_" + i] = 0




#parse through index files and count index combos
with open(f2) as indexfile1, open(f3) as indexfile2:
    linecount = 0
    incorrect = 0
    correct = 0
    swapped = 0
    not_swapped = 0


    #intialize index variables
    i1 = ""
    i2 = ""
    i3 = ""


    #qual scores are local variables below

    for x, y in zip(indexfile1, indexfile2):


        #get indexes and seqs
        if linecount % 4 == 1:
            #indexes
            i1 = str(x.strip())
            i2 = str(y.strip())

            #convert to reverse complement
            i2 = rComplement(i2)
            #index combo
            i3 = i1 + "_" + i2

            #check to see if indexes are in dict
            if i3 in indexdict:

                #increment counter in dict
                indexdict[i3] += 1
                correct += 1

                #check if swapped
                if i1 != i2:
                    swapped = swapped + 1
                else:
                    not_swapped = not_swapped + 1
            else:
                #not in dict
                incorrect += 1

        #increment line counter
        linecount += 1

    #output results
    abspath = "/home/bnorris8/BI621/index_hopping/"
    #get total number of reads
    total_reads = incorrect + correct

    #open two output files
    with open(abspath + "ih_output_table.txt", "w") as output1, open(abspath + "swap_dist.txt", "w") as output2:
        #output not swapped index pair results first
        for k in indexdict.keys():
            if k[:8] == k[9:]:
                output1.write(k + "\t" + str(indexdict[k]) + "\t" + str((indexdict[k] / total_reads) * 100) + "\n")
        #output swapped index pairs results to 2 files one is the output table, 2 is file for plotting distribution
        for k in indexdict.keys():
            if k[:8] != k[9:]:
                output1.write(k + "\t" + str(indexdict[k]) + "\t" + str((indexdict[k] / total_reads) * 100) + "\n")
                output2.write(k + "\t" + str(indexdict[k]) + "\n")

    #output index swapping statistics to its own file
    with open(abspath + "ih_stats.txt", "w") as output:
        output.write("Total number of reads: " + str(total_reads) + "\n")
        output.write("Number of undetermined reads: " + str(incorrect) + "\n")
        output.write("Number of index swapped reads: " + str(swapped) + "\n")
        output.write("Number of correct index reads: " + str(not_swapped) + "\n")
        output.write("Percentage of reads swapped: " + str((swapped / total_reads) * 100) + "\n")
        output.write("Percentage of correct index reads: " + str((not_swapped / total_reads) * 100) + "\n")

















