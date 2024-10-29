import pandas as pd
from collections import defaultdict
from datetime import datetime, date, time, timedelta
import sys

def read_bar_mut(File):

    bar_mut_c_D = {}

    with open(File, 'r') as f:
        for line in f:
            value1, geno = line.strip().split("\t")
            count, bar = value1.strip().split(" ")
            if bar in bar_mut_c_D:
                bar_mut_c_D[bar].append([geno, count])
            else:
                bar_mut_c_D[bar] =[[geno, count]]
    return bar_mut_c_D

def read_cluster(File):

    cluster_D = {}

    with open(File, 'r') as f:
        for line in f:
            bar1, bar2, num = line.strip().split(",")
            if float(num) == 1:
                if bar1 in cluster_D:
                    cluster_D[bar1].append(bar2)
                else:
                    cluster_D[bar1] =[bar2]
                if bar2 in cluster_D:
                    cluster_D[bar2].append(bar1)
                else:
                    cluster_D[bar2] =[bar1]
    return cluster_D

## find families with >1 clusteration and drop this barcode
def find_drop_bar(cluster_D):
    bar_to_drop = []
    n_Family = []
    for key in cluster_D:
        Include = [key] + cluster_D[key]
        Exclude = []
        for i in Include:
            Exclude.append(i)
            Exclude += cluster_D[i]
        if set(Include) != set(Exclude):
            check = 'Fail'
            for x in set(Exclude):
                Include_c = [x] + cluster_D[x]
                Exclude_c = []
                for ix in Include_c:
                    Exclude_c.append(i)
                    Exclude_c += cluster_D[ix]
                if set(Include_c) == set(Exclude_c):
                    check = "Succ"
            if check == 'Fail':
                bar_to_drop += list(set(Exclude))
            else:
                if sorted(list(set(Exclude_c))) not in n_Family:
                    n_Family.append(sorted(list(set(Exclude_c))))
        else:
            if sorted(list(set(Exclude))) not in n_Family:
                n_Family.append(sorted(list(set(Exclude))))
    
    bar_to_drop = list(set(bar_to_drop))
    return bar_to_drop, n_Family


def find_major(dictA):
    dictA = {key: int(value) for key, value in dictA.items()}
    max_value = max(dictA.values())
    max_key = [i for i, value in dictA.items() if value == max_value]
    total_sum = sum(dictA.values())
    ## ratio
    Ratio = int(max_value) / total_sum
    return max_key, max_value, Ratio


def Filter_majority(bar_mut_dict, cluster_dict, bar_drop_list, family_list, cutoff=0.5):
    final_D = {}
    Drop_D = {}

    for key in bar_mut_dict:
        status_family = ""
        status_bar = ""
        status_geno = ""
        status_center = ""
        
        sums = defaultdict(int)

        if key not in cluster_dict:
            status_family = "uniq"
            ## the barcode is uniq
            # 1. convert it into a dict
            sums = dict(bar_mut_dict[key])
            max_key, max_value, Ratio = find_major(sums)
        
            if Ratio > 0.5:
                ## check if ins|del in major genotype
                if "ins" in max_key[0] or "del" in max_key[0] or "-" in max_key[0]:
                    status_geno = f"HasIndel"
                    Drop_D[key] = [status_family, status_geno, status_bar, status_center]
                else :
                    final_D[key] = [max_key[0], round(Ratio,2)]
            else:
                status_geno = f"Unmajor_geno({Ratio})"
                Drop_D[key] = [status_family, status_geno, status_bar, status_center]
        elif key not in bar_drop_list:
            status_family = "clustered"
            ## the barcode has similar and step = 1
            all_bar = [sub_F for sub_F in family_list if key in sub_F][0]
            ## 可能family里面的barcode已经用ccs <2 筛选掉了
            all_bar = [i for i in all_bar if i in bar_mut_dict.keys()]
            if len(all_bar) == 1:
                status_family = "AfilterUniq"
                ## the barcode is uniq
                # 1. convert it into a dict
                sums = dict(bar_mut_dict[key])
                max_key, max_value, Ratio = find_major(sums)
        
                if Ratio > 0.5:
                    ## check if ins|del in major genotype
                    if "ins" in max_key[0] or "del" in max_key[0] or "-" in max_key[0]:
                        status_geno = f"HasIndel"
                        Drop_D[key] = [status_family, status_geno, status_bar, status_center]
                    else :
                        final_D[key] = [max_key[0], round(Ratio,2)]
                else:
                    status_geno = f"Unmajor_geno({Ratio})"
                    Drop_D[key] = [status_family, status_geno, status_bar, status_center]
            else:
                all_geno = []
                for i in all_bar:
                    all_geno += bar_mut_dict[i]
                    # find major genotype
                # 1. convert it into a dict
                for geno, value in all_geno:
                    if geno not in sums:
                        sums[geno] = int(value)
                    else:
                        sums[geno] += int(value)

                max_key, max_value, Ratio = find_major(sums)

                    ## find major barcode
                bar_c_D = {}
                for i in all_bar:
                    bar_c_D[i] = sum(int(value) for _, value in bar_mut_dict[i])
                max_bar, max_bar_value, Ratio_bar = find_major(bar_c_D)


                    ## find center bar
                bar_center_D = {}
                for i in all_bar:
                    bar_center_D[i] = len(cluster_dict[i])
                max_link, max_link_value, Ratio_link = find_major(bar_center_D)

                if Ratio > 0.5 and Ratio_bar > 0.5 and max_bar[0] in max_link:
                    if max_bar[0] == key:
                        if "ins" in max_key[0] or "del" in max_key[0] or "-" in max_key[0]:
                            status_geno = f"HasIndel"
                            Drop_D[key] = [status_family, status_geno, status_bar, status_center]
                        else :
                            final_D[key] = [max_key[0], round(Ratio,2)]
                    else:
                        Drop_D[key] = ["PASS_unCenter", status_geno, status_bar, status_center]
                    ## which means only one major barcode, and it is one of the center
                else:
                    if Ratio <= 0.5:
                        status_geno = f"Unmajor_geno({Ratio})"
                    if status_family == "clustered":
                        if Ratio_bar <= 0.5:
                            status_bar = f"Unmajor_bar({Ratio_bar})"
                        if Ratio_bar > 0.5 and max_bar[0] not in max_link:
                            status_center = f"Uncentered"
                    Drop_D[key] = [status_family, status_geno, status_bar, status_center]
        else:
            Drop_D[key] = ["familyDrop", status_geno, status_bar, status_center]
    return final_D, Drop_D
def main():
    # Check if the correct number of arguments are provided
    
    # Retrieve the arguments
    Path =  sys.argv[1]
    Name = sys.argv[2]

    File = Path + Name + ".bar_mut_c.txt"
    Log = Path + Name + ".cluster.log"
    SlideSort_file = Path + Name + ".slidesort.ss.txt"
    OUTPUT = Path + Name + ".PASS.txt"
    FAIL = Path + Name + ".Fail.txt"
    
    log_file = open(Log, 'w')

        # Print the arguments
    log_file.write(f'Running barcode_geno file: {File}\n')
    log_file.write(f'Running slidesort file: {SlideSort_file}\n')
    log_file.write(f'Saving to: {OUTPUT} {FAIL}\n')
    log_file.write(f"Start Running : {datetime.now()}\n")
    log_file.flush()
    
    # 1. 读取两个表格 ： barcode_genotype 计数表 和 slidesort输出
    bar_mut_c_D = read_bar_mut(File)
    cluster_D = read_cluster(SlideSort_file)

    # 2. RAW 计数
    n_raw_bar = len(set(bar_mut_c_D.keys()))
    n_raw_geno = len(set([i[0] for key in bar_mut_c_D for i in bar_mut_c_D[key]]))

    log_file.write("# -------- 1. Before filteration and clusteration -------- #\n")
    log_file.write(f"  Number of Barcodes is : {n_raw_bar}\n")
    log_file.write(f"  Number of geno is : {n_raw_geno}\n")
    log_file.flush()

    # 3. filter ccs > 1
    keys_to_remove = []
    for key in bar_mut_c_D:
        count_L = [int(i[1]) for i in bar_mut_c_D[key]]
        SUM = sum(count_L)
        if SUM < 2:
            keys_to_remove.append(key)
    
    for key in keys_to_remove:
        del bar_mut_c_D[key]

    n_fill_bar = len(set(bar_mut_c_D.keys()))
    n_fill_geno = len(set([i[0] for key in bar_mut_c_D for i in bar_mut_c_D[key]]))

    log_file.write("# -------- 2. Filteration: CCS >= 2 -------- #\n")
    log_file.write(f"  Number of Barcodes PASS FIlTERATION is : {n_fill_bar} ( {round(n_fill_bar/n_raw_bar * 100, 2)}% )\n")
    log_file.write(f"  Number of geno PASS FILTERATION is : {n_fill_geno} ( {round(n_fill_geno/n_raw_geno * 100, 2)}% )\n")
    log_file.flush()

    # 4. clustering barcodes and drop those has a distance > 1
    bar_to_drop, n_Family = find_drop_bar(cluster_D)

    bar_uni = [i for i in bar_mut_c_D.keys() if i not in cluster_D]

    log_file.write("# -------- 3. Clustering 1bp mismatch barcodes -------- #\n")
    log_file.write(f"  Number of Barcodes to DROP is : {len(bar_to_drop)} ( {round(len(bar_to_drop)/n_fill_bar * 100, 2)}% )\n")
    log_file.write(f"  Number of Barcodes that has no similar bar is : {len(bar_uni)} ( {round(len(bar_uni)/n_fill_bar * 100, 2)}% )\n")
    log_file.write(f"  Number of Barcodes that has similar bars is : {(len(cluster_D.keys()) - len(bar_to_drop))} ( {round((len(cluster_D.keys()) - len(bar_to_drop))/n_fill_bar * 100, 2)}% )\n")
    log_file.write(f"  Number of =1 Families is : {len(n_Family)}\n")
    log_file.flush()

    # 5. 正式找Filteration
    final_dict, Drop_dict = Filter_majority(bar_mut_c_D, cluster_D, bar_to_drop, n_Family, cutoff=0.5)

    # 6. check
    if len(final_dict.keys()) + len(Drop_dict.keys()) != len(bar_mut_c_D.keys()):
        log_file.write(f"  Check Fail : SUM - {len(final_dict.keys()) + len(Drop_dict.keys())} ALL - {len(bar_mut_c_D.keys())}\n")
        log_file.flush()
    else:
        log_file.write(f"  Check FPASS\n")
        log_file.flush()

    # 7. 写出结果
    with open(OUTPUT, 'w') as f1:
        for key, values in final_dict.items():
            line = f"{key}\t{values[0]}\t{values[1]}\n"
            f1.write(line)

    with open(FAIL, 'w') as f1:
        for key, values in Drop_dict.items():
            ll = "\t".join([str(i) for i in values])
            line = f"{key}\t{ll}\n"
            f1.write(line)
    
    log_file.write(f" Success")
    log_file.write(f"End Running : {datetime.now()}")

    log_file.close()
    
if __name__ == '__main__':
    main()
