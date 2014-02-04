import sys
from collections import Counter

PBfastq = open(sys.argv[1])

def just_reads(fastq):
    global listfilelines
    listfilelines = fastq.readlines()
    seq_counter = 0
    line_counter = 0
    seqs = []
    while line_counter < len(listfilelines):
        seqs.append(listfilelines[line_counter+1])
        line_counter += 4
    return seqs

def reverse_comp(seq):
    seq = seq.replace("\n", "")
    comp_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}
    comp = ""
    for x in seq:
        comp += comp_dict[x]
    reverse_comp = comp[::-1] + "\n"       
    return reverse_comp

def freqdistcsv(sizelist):
    distribution_dictionary = (Counter(sizelist))
    csv_format = ""
    for size in distribution_dictionary:
        csv_format += str(size) + "," + str(distribution_dictionary[size])+ "\n"
    return csv_format

seq_list = just_reads(PBfastq)
PBfastq.close()


right_way_list = []
cannot_be_righted = []
for seq in seq_list:
    if "AATGATACGGCGACCACCGAGATCTACAC" in seq[:29] and "ATCTCGTATGCCGTCTTCTGCTTG" in seq[-25:]:
        right_way_list.append(seq)
    elif "CAAGCAGAAGACGGCATACGAGAT" in seq[:24]:
        if "AATGATACGGCGACCACCGAGATCTACAC" in reverse_comp(seq)[:29] and "ATCTCGTATGCCGTCTTCTGCTTG" in reverse_comp(seq)[-25:]:
            right_way_list.append(reverse_comp(seq))
        else:
            cannot_be_righted.append(seq)
    else:
        cannot_be_righted.append(seq)


barcodelist = ["ATTACTCGTATAGCCT", "ATTACTCGATAGAGGC", "ATTACTCGCCTATCCT", "ATTACTCGGGCTCTGA",
               "TCCGGAGATATAGCCT", "TCCGGAGAATAGAGGC", "TCCGGAGACCTATCCT", "TCCGGAGAGGCTCTGA",
               "CGCTCATTTATAGCCT", "CGCTCATTATAGAGGC", "CGCTCATTCCTATCCT", "CGCTCATTGGCTCTGA",
               "GAGATTCCTATAGCCT", "GAGATTCCATAGAGGC", "GAGATTCCCCTATCCT", "GAGATTCCGGCTCTGA"]

demult = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

for y in right_way_list:
    demult_list_pos = 0
    for x in barcodelist:
        if x[8:16] in y[29:37] and x[:8] in y[-33:-25]:
            demult[demult_list_pos].append(y)
        demult_list_pos += 1
        


 
print "Total reads", len(seq_list)       
print "Number of reads which are in P5 to P7 orientation or have been succesfully fliped to that direction", len(right_way_list)
print "Cannot be righted", len(cannot_be_righted), "\n"

print "           P501       P502       P503       P504"
print "P701", str(len(demult[0])).rjust(10," "),str(len(demult[1])).rjust(10," "),str(len(demult[2])).rjust(10," "),str(len(demult[3])).rjust(10," ")
print "P702", str(len(demult[4])).rjust(10," "),str(len(demult[5])).rjust(10," "),str(len(demult[6])).rjust(10," "),str(len(demult[7])).rjust(10," ")
print "P703", str(len(demult[8])).rjust(10," "),str(len(demult[9])).rjust(10," "),str(len(demult[10])).rjust(10," "),str(len(demult[11])).rjust(10," ")
print "P704", str(len(demult[12])).rjust(10," "),str(len(demult[13])).rjust(10," "),str(len(demult[14])).rjust(10," "),str(len(demult[15])).rjust(10," "),"\n"

counter = 0
distribution_csv_format = ""
poly_a_distribution_csv_format = ""
for z in demult:
    if len(z) < 100:
        try:
            if sys.argv[2] == "-lc":
                print "\nBarcode1 Barcode2 Random N Seq start -3"
                for seqs in z:
                    print seqs[29:37], seqs[-33:-25], seqs[37:45], seqs[75:130]
            pass
        except:
            pass

    else:
        print barcodelist[counter]
        print "Reads\t\t\t\t\t\t\t-\t", len(demult[counter])
        counter += 1
        Ad25 = 0
        Td15 = 0
        library_size = 0
        start_of_tail_list = []
        start_of_tail_count = 0
        straight_into_tail = 0
        GGAGATTCCCGAATAGGTTmotif = 0
        GAAGCGTCCTCAGCGACGGACmotif = 0
        list_library_size = []
        polya_list_library_size = []
        for x in z:
            list_library_size.append(len(x))
            library_size += len(x)
            if "CTTCCGATCTAAAAAAAAAA" in x:
                straight_into_tail += 1
            if "AAAAAAAAAAAAAAAAAAAAAAAAA" in x:
                start_of_tail_list.append(x.find("AAAAAAAAAAAAAAAAAAAAAAAAA"))
                start_of_tail_count += x.find("AAAAAAAAAAAAAAAAAAAAAAAAA")
                Ad25 += 1
                polya_list_library_size.append(len(x))
            if "TTCCGATCTTTTTTTTTTTTTT" in x:
                Td15 += 1
            if "GGAGATTCCCGAATAGGTT" in x:
                GGAGATTCCCGAATAGGTTmotif += 1
            if "GAAGCGTCCTCAGCGACGGAC" in x:
                GAAGCGTCCTCAGCGACGGACmotif += 1
                
        distribution_csv_format += barcodelist[counter] + "\n" + freqdistcsv(list_library_size)
        poly_a_distribution_csv_format += barcodelist[counter] + "\n" + freqdistcsv(polya_list_library_size)

        print "Average poly A start position in an Illumina read\t-\t%.1f" % ((start_of_tail_count / float(Ad25)) - 78)
        print "Percentage of reads with Ad(25) present\t\t\t-\t%.1f" % ((Ad25 / float(len(z))) * 100)
        print "Percentage that read straight into a-tail\t\t-\t%.1f" % ((straight_into_tail / float(len(z))) * 100)
        print "Average library size\t\t\t\t\t-\t%.1f" % (library_size / float(len(z)))
        print "Percentage with invert\t\t\t\t\t-\t%.1f" % (Td15 / (float(len(z))) * 100)
        print "Percentage GGAGATTCCCGAATAGGTT motif\t\t\t-\t%.1f" % (GGAGATTCCCGAATAGGTTmotif / (float(len(z))) * 100)
        print "Percentage GAAGCGTCCTCAGCGACGGAC motif\t\t\t-\t%.1f\n" % (GAAGCGTCCTCAGCGACGGACmotif / (float(len(z))) * 100)

try:
    if "-d" in sys.argv:
        print "P5 to P7 orientation/flipped dist"
        print distribution_csv_format
        print "Poly A dist"
        print poly_a_distribution_csv_format
except:
    pass
