from sys import argv
import os
import regex
import json
import collections

fq1, fq2, output_dir, output_base = argv[1:]

output1 = output_dir + output_base + ".err1.umi.txt"
output2 = output_dir + output_base + ".err1.unuse.txt"
output3 = output_dir + output_base + ".err1.umi_index.txt"


out1=open(output1, "w")
out2=open(output2, "w")
out3=open(output3, "w")


f_re = "(?e)(?P<p1>TCGCTCTTATTGACCACACC){e<=1}(?P<umi_f>[ATCG]{20}|[ATCG]{25})(?P<p2>GCTTCGGCAGCACATATACT){e<=1}"
r_re = "(?e)(?P<p3>AGTATATGTGCTGCCGAAGC){e<=1}(?P<umi_r>[ATCG]{20}|[ATCG]{25})(?P<p4>GGTGTGGTCAATAAGAGCGA){e<=1}"


with open(fq1, "r") as f1, open(fq2, "r") as f2:
    for j, element in enumerate(zip(f1, f2)):
        if j % 4 == 1:
            forward_seq = regex.search(f_re, element[0].strip(), flags=0)
            reverse_seq = regex.search(r_re, element[1].strip(), flags=0)
            if forward_seq is not None and reverse_seq is not None:
                forward_umi = forward_seq.groupdict()["umi_f"]
                reverse_umi = reverse_seq.groupdict()["umi_r"]
                reverse_umi = reverse_umi[::-1].replace('A', 't').replace('T', 'a').replace('G', "c").replace("C",'g').upper()
                if forward_umi == reverse_umi:
                    if len(forward_umi) == 25:
                        out3.write(forward_umi + "\n")
                    forward_umi = forward_umi[:20]
                    out1.write(forward_umi + "\n")
                    
            else:
                out2.write("%s %s %s\n" % (j, element[0], element[1]))
out1.close()
out2.close()
out3.close()

count_dict = {}

with open(output1, "r") as f1:
    str1 = f1.read().split("\n")
    for k in str1:
        if k in count_dict:
            count_dict[k] += 1
        else:
            count_dict[k] = 1

output = output_dir + output_base + ".txt"

with open(output, "w") as out3:
    for i in count_dict:
        out3.write('%s\t%s\n' % (i, count_dict[i]))

out3.close()

