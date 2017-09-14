#use arg parse to get args
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description="bring in sequence files")
    parser.add_argument("-f1", help="set R1 file", required=True, type=str)
    parser.add_argument("-f2", help="set R2 file", required=True, type=str)
    parser.add_argument("-f3", help="set R3 file", required=True, type=str)
    parser.add_argument("-f4", help="set R4 file", required=True, type=str)
    return parser.parse_args()

args = get_arguments()

f1 = args.f1
f2 = args.f2
f3 = args.f3
f4 = args.f4

def convert_to_phred(char):
    '''converts quality data to phred score'''
    return ord(char) - 33

def bp_qd(file, index, n):
    if index == True:
        bp_len = 8
    else:
        bp_len = 101

    #create empty score list
    mean_scores = []  # create empty list
    for x in range(bp_len):
        mean_scores.append(0.0)


    #add up quality scores
    with open(file) as fh:
        NR = 0  # line counter
        for line in fh:
            NR += 1
            line = line.strip()
            if NR % 4 == 0:  # picks out only quality score lines
                i = 0  # tracks number of chars converted
                for c in line:
                    phred = ord(c) - 33 # convert chars to phred score
                    mean_scores[i] += phred  # add phred score to respective element of mean scores
                    i += 1

    reads = NR / 4
    j = 0
    abspath = "/home/bnorris8/BI621/index_hopping/"

    with open(abspath + "R" + str(n) + "_qd.txt", "w") as output:
        for s in mean_scores:
            output.write(str(j) + "\t" + str(mean_scores[j] / reads) + "\n")
            j += 1

def read_freq(file, index, n):
    score_dict = {}

    if index == True:
        bp_len = 8
    else:
        bp_len = 101

    # add up quality scores
    with open(file) as fh:
        NR = 0  # line counter
        for line in fh:
            NR += 1
            line = line.strip()
            if NR % 4 == 0:  # picks out only quality score lines
                sum = 0
                for c in line:
                    sum += ord(c) - 33
                avg = int(sum / bp_len)
                if avg in score_dict:
                    score_dict[avg] += 1
                else:
                    score_dict[avg] = 1
    abspath = "/home/bnorris8/BI621/index_hopping/"
    with open(abspath + "R" + str(n) + "_rf.txt", "w") as output:
        for k in score_dict.keys():
            output.write(str(k) + "\t" + str(score_dict[k]) + "\n")


bp_qd(f1, False,1)
bp_qd(f2, True,2)
bp_qd(f3, True,3)
bp_qd(f4, False,4)

read_freq(f1, False,1)
read_freq(f2, True,2)
read_freq(f3, True,3)
read_freq(f4, False,4)