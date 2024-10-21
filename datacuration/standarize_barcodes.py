import re 
import os 
from multiprocessing import Pool

def process_fastq(file_path): 
    W1 = "GAGTGATTGCTTGTGACGCCTT"
    W1_REV = "AAGGCGTCACAAGCAATCACTC"
    BC2_LENGTH = 8
    UMI_LENGTH = 6
    
    base_name = os.path.basename(file_path).split('_')[0]  # Extracts SRRXXXXXXX
    output_file_path = f"{base_name}_1_bc.fastq"

    with open(file_path, "r") as f, open(output_file_path, 'w') as outf:  # the a opens it in append mode
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            match = re.search(W1, seq)
            if match != None: 
                st_w1= match.start()
                end_w1 = match.end()
                new_seq = seq[0:st_w1] + seq[end_w1:end_w1 + BC2_LENGTH + UMI_LENGTH] # remove w1 from the sequence and keep bc1, bc2, and umi
                qual = qual[0:st_w1] + qual[end_w1:end_w1 + BC2_LENGTH + UMI_LENGTH]
                # Ensure barcode1 is trimmed or padded to 11 bases
                if len(new_seq) < 25:
                    new_seq = new_seq.rjust(25, 'N')  # Pad with N
                    qual = qual.rjust(25, '!')  # Pad quality scores with lowest quality
                    header = f"{header} trimmed"
                    plus = f"{plus} trimmed" 
                else: 
                    new_seq = seq

            else: 
                new_seq = seq
            outf.write(f"{header}\n{new_seq}\n{plus}\n{qual}\n")
            #print(header)
            #print(str(len(new_seq)) + "  " + new_seq)
            #print(qual)

    print("finished trimming barcodes for " + file_path)

        
def main():
    fastq_files = ['singlecell_pancreas/data/SRR3879604_1.fastq', 'singlecell_pancreas/data/SRR3879605_1.fastq', 'singlecell_pancreas/data/SRR3879606_1.fastq']  # read file from read1 that captures the barcodes and UMI
    for ffile in fastq_files:
        process_fastq (ffile)

if __name__ == "__main__":
    main()