import sys
import os
os.environ['PATH']+= os.pathsep + '/mnt/data5/disk/guoziyan/anaconda3/envs/enzyme_kinetic_Python/bin/'
import psa
from Bio.Seq import Seq
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def get_first_n_items(d, n):
    subdict = {}

    for i, (header,item) in enumerate(d.items()):
        if i >= n:
            break
        subdict[header] = item
    return subdict
    
def read_call_file(file):
    seq_dict = {}

    with open(file, 'r') as f:
        seq = []
        header = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                ## Save the previous seq if it exists
                if header:
                    seq_dict[header] = seq
                header = line[1:]
            else:
                seq = line.split(" ")
        if header:
            seq_dict[header] = seq
    return seq_dict

def print_first_n_sequences(fasta_dict, n=10, index = True):
    for i, (header, sequence) in enumerate(fasta_dict.items()):
        if i >= n:
            break
        print(f"Header: {header}")
        print(f"Sequence: {sequence[0][:50]}...")  # Print first 50 characters of the sequence for brevity
        if index:
            print(f"Index: {sequence[2]}")
        print(f"Barcode: {sequence[1]} ({len(sequence[1])})")
        if len(sequence) > 2:
            print(f"nucl: {sequence[len(sequence)-2]}")
            print(f"AA: {sequence[len(sequence)-1]}")
            
def read_ref(file):
    with open(file, 'r') as f:
        Header = ''
        Seq = ''
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                Seq = line
                return [Header, Seq]
            Header = line[1:]

def find_mut_AA(qpro, spro):
    change = []
    if len(qpro) != len(spro):
        return 'ins|del'
    else:
        for i in range(len(spro)):
            if qpro[i] != spro[i]:
                change.append(f"{spro[i]}{i+1}{qpro[i]}") 

        return " ".join(change) if change else "WT"

def find_first_asterisk(protein_seq):
    # Find the index of the first occurrence of '*'
    index = protein_seq.find('*')
    return index + 1 if index != -1 else "no stop codon"

def process_sequence(header, seqL, ref):
    aln = psa.needle(moltype = 'nucl', qseq = seqL[0], sseq = ref)
    nucL = []
    for i,(qseq, sseq) in enumerate(aln):
        if qseq != sseq:
            if qseq == "-":
                nucL.append(f"{sseq}{i+1}del")
            elif sseq == "-":
                nucL.append(f"ins{i+1}{qseq}")
            else:
                nucL.append(f"{sseq}{i+1}{qseq}")
    
    nuc = " ".join(nucL) if nucL else 'WT'

    ## AA change
    qpro = Seq(seqL[0]).translate()
    spro = Seq(ref).translate()
    prot = find_mut_AA(qpro, spro)

    ## find stop codon
    sstp = find_first_asterisk(qpro)

    return header, nuc, prot, sstp

def needle_mut(fasta_dict, ref, num_threads = 10):
    results = {}
    headers = list(fasta_dict.keys())

    with ThreadPoolExecutor(max_workers= num_threads) as executor:
        with tqdm(total = len(headers), desc = "Processing sequences") as progress_bar:
           
            future_to_header = {executor.submit(process_sequence, header, fasta_dict[header], ref[1]): header for header in headers}

            for future in as_completed(future_to_header):
                header = future_to_header[future]
                try:
                    header, nuc, prot, sstp = future.result()
                    if header in fasta_dict:
                        fasta_dict[header].append(nuc)
                        fasta_dict[header].append(prot)
                        fasta_dict[header].append(sstp)
                except Exception as exc:
                    print(f'{header} generated an exception: {exc}')
                finally:
                    # Update the progress bar
                    progress_bar.update(1)
    return fasta_dict
    
def write_out(output_path, results_dict):
    with open(output_path, 'w') as file:
        for key, values in results_dict.items():
            file.write(f"{key}\t" + '\t'.join(map(str, values)) + "\n")


def main():
    # Check if the correct number of arguments are provided
    
    # Retrieve the arguments
    file =  sys.argv[1]
    ref = sys.argv[2]
    output_path = sys.argv[3]
    
    # Print the arguments
    print(f'Running call file: {file}')
    print(f'Using ref: {ref}')
    print(f'Saving to: {output_path}')

    seq_dict = read_call_file(file)
    ref = read_ref(ref)
    print(f"Header: {ref[0]} \nSeq: {ref[1]} ({len(ref[1])})")

    results = needle_mut(seq_dict, ref)
    write_out(output_path, results)
    
if __name__ == '__main__':
    main()