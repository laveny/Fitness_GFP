from sys import argv
import regex as re

input_dir, name_prefix, outdir = argv[1:]
mud_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
Reverse_complemrnt = lambda x: "".join([mud_dict[i] for i in x][::-1])

geno_pattern=re.compile("AACACACATAAACAAACAAA(?P<geno>[ATCG]{804})ATAAGCGAATTTCTTATGAT")
umi_pattern=re.compile("TCGCTCTTATTGACCACACC(?P<barcode>[ATCG]{20})(?P<index>[ATCG]{5})GCTTCGGCAGCACATATACT")
outfile = open(outdir + name_prefix+".call.txt", "w")
with open(input_dir + name_prefix + ".n.ccs.sam", "r") as inputr:  # nagetive strand ccs result
    for nega in inputr:
        if nega[0] == "m":
            nega_line_list = nega.strip().split()
            nega_seq = nega_line_list[9]  # the result sequence after ccs
            zmw_nega = nega_line_list[0].split("/")[1]  # the number of ZMW
            nega_reverse = Reverse_complemrnt(nega_seq)
            geno2 = re.search(geno_pattern, nega_reverse, flags=0)
            umi2 = re.search(umi_pattern, nega_reverse, flags=0)
            if geno2 is not None and umi2 is not None:
                geno3 = geno2.groupdict()["geno"]
                umi3 = umi2.groupdict()["barcode"]
                index3 = umi2.groupdict()["index"]
                outfile.write(">" + name_prefix + "/ccs/" + str(zmw_nega) + "/n" + "\n" + geno3 + " " + umi3 + " " + index3 + "\n")
with open(input_dir + name_prefix + ".p.ccs.sam", "r") as inputf:  # positive strand ccs result
    for posi in inputf:
        if posi[0] == "m":
            posi_line_list = posi.strip().split()
            posi_seq = posi_line_list[9]
            zmw_posi = posi_line_list[0].split("/")[1]
            geno = re.search(geno_pattern, posi_seq, flags=0)
            umi = re.search(umi_pattern, posi_seq, flags=0)
            if geno is not None and umi is not None:
                geno1 = geno.groupdict()["geno"]
                umi1 = umi.groupdict()["barcode"]
                index1 = umi.groupdict()["index"]
                outfile.write(">" + name_prefix + "/ccs/" + str(zmw_posi) + "/p" + "\n" + geno1 + " " + umi1 + " " + index1 + "\n")

# os.system("rm *%s.*_ccs.bam.pbi"%(raw_bam[:-4]))
# os.system("rm *%s.*_ccs.sam"%(raw_bam[:-4]))
outfile.close()