
import os
import sys 
import math
import subprocess
import re
from pysam import AlignmentFile, AlignedSegment

# TE_size = "/zata/zippy/boxu/annotation/hg38/ALSE.size"

user_path = os.environ.get("user_path", os.environ.get("USER_PATH", "/path/to/user"))

data_file = sys.argv[1]
TE_size = sys.argv[2]
out_file = sys.argv[3]
data_type = sys.argv[4]
caller = sys.argv[5]
depth = sys.argv[6]
out_path = sys.argv[7]
genome = sys.argv[8]
cal_div_script = sys.argv[9] if len(sys.argv) > 9 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "03_cal_div3_giab.sh"
)

te_size_dict = {}


for line in open(TE_size):
    line = line.strip().split("\t")
    te_size_dict[line[0]] = int(line[1])



te_seq_dict = {}
def read_fasta_to_dict(filename):
    """
    Read a FASTA file into a dictionary.
    Returns: {sequence_id: sequence_string}
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):  # Sequence identifier line.
                # Save the previous sequence.
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Extract the sequence ID.
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save the last sequence.
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

te_seq_dict = read_fasta_to_dict('/zata/zippy/boxu/annotation/from_temp3/GRCh38.transposon.fa')
# for seq_id, sequence in seq_dict.items():
#     print(f"ID: {seq_id}")
#     print(f"Sequence: {sequence[:50]}...")
    


def cal_divergency(read):
    """
        Function:
            Calculate divergence from supported alignment CIGAR modes.
        Parameter: 
            reads
            reference
        Returns:
            divergency
    """

    cigar_count = AlignedSegment.get_cigar_stats(read)[0]

    if cigar_count[7] >= cigar_count[0]:
        # minimap2  --eqx 
        n_PM = cigar_count[7]
        n_MM = cigar_count[8]
        n_ins = cigar_count[1]
        n_del = cigar_count[2]
        n_M = n_PM + n_MM
    else:
        n_M = cigar_count[0]
        NM = cigar_count[-1]
        n_ins = cigar_count[1]
        n_del = cigar_count[2]
        n_MM = cigar_count[10]

    return n_MM, n_ins, n_del

g = open(out_file, "a")

i = 0
for line in open(data_file):
    line_info = line.strip().split("\t")
    # print(line_info)
    gold_ins = line_info[3:9]
    caller_ins = line_info[9:15]

    gold_ins_start = int(gold_ins[1])
    gold_ins_end = int(gold_ins[2])

    gold_ins_te_info = gold_ins[3].split("~")

    gold_ins_te = gold_ins_te_info[1]

    gold_te_start = 0
    gold_te_end = te_size_dict[gold_ins_te]


    # if  gold_ins_te == "ALU":
    #     gold_te_end = 281
    # if gold_ins_te == "LINE1":
    #     gold_te_end = 6019
    #     if gold_ins_te_info[2].split(":")[0] == "TF" or gold_ins_te_info[2].split(":")[0] == "TR":
    #         trunc_te_ins = gold_ins_te_info[2].split("..")
    #         trunc_te_ins_st = int(trunc_te_ins[0].split("[")[1])
    #         trunc_te_ins_ed = int(trunc_te_ins[1].split("]")[0])

    #         if trunc_te_ins_st < 10:
    #             gold_te_start = trunc_te_ins_ed
    #             gold_te_end = 6019
    #         if 6019 - trunc_te_ins_ed < 10:
    #             gold_te_start = 0
    #             gold_te_end = trunc_te_ins_ed - 1
    #         # if trunc_te_ins_st >= 10 and 6019 - trunc_te_ins_ed >= 10:
    #         #     gold_te_start = [0,trunc_te_ins_ed - 1]
    #         #     gold_te_end = [gold_te_start - 1, 6019]
   


    caller_ins_start = int(caller_ins[1])
    caller_ins_end = int(caller_ins[2])
    


    caller_te_info = caller_ins[3].split(":")
    
    caller_te_start = caller_te_info[1]
    if caller_te_start == "NA":
        caller_te_start = 0
    else:
        caller_te_start = int(caller_te_start)

    caller_te_end = caller_te_info[2]
    if caller_te_end == "NA":
        caller_te_end = 0
    else:
        caller_te_end = int(caller_te_end)

    if caller == "xTea":
        if gold_ins_te == "SVA":
            caller_te_start = min(0, caller_te_start)

        caller_te_end = min(te_size_dict[gold_ins_te], caller_te_end)

    # gold_mid = (gold_ins_start + gold_ins_end)/2
    # caller_mid = (caller_ins_start + caller_ins_end)/ 2
    # ins_dis = abs(gold_mid - caller_mid)

    ins_dis = min([abs(gold_ins_start - caller_ins_start), abs(gold_ins_end - caller_ins_end)])


    # dis_start = abs(caller_te_start - gold_te_start )
    # dis_end = abs(caller_te_end - gold_te_end )
    # te_dis = (dis_start + dis_end ) / 2

    if caller_te_start == 0 and caller_te_end == 0:
        te_dis = "-"
    else:
        te_dis = min([ abs(caller_te_start - gold_te_start), abs(caller_te_end - gold_te_end)])


    # sequence similarity

    insert_seq = caller_ins[5]
    ref_seq = gold_ins[5]

    ### remove polyA
    if re.search("A{20,}",ref_seq):
        poly_tail = re.search("A{20,}",ref_seq).span()
        ref_seq = ref_seq[0:poly_tail[0]]
    
    if re.search("T{20,}",ref_seq):
        poly_tail = re.search("T{20,}",ref_seq).span()
        ref_seq = ref_seq[poly_tail[1]:]



    div_info = ["-", "-", "-", "-", "-"]
    seq_identity = "-"
    if caller_ins[5] != "-" and ref_seq not in ["-",""] :
        subprocess.run(
            ["bash", cal_div_script, "-s", insert_seq, "-r", ref_seq, "-o", out_path, "-g", genome],
            stderr=subprocess.DEVNULL,
            check=True,
        )
        for read in AlignmentFile(out_path + "/" +  genome + ".test3.sam", 'rb'):
            if not read.is_supplementary and not read.is_unmapped:
                align_info = cal_divergency(read)
                # print(ref_seq)
                seq_identity = min((read.query_alignment_length - align_info[1])  / len(ref_seq), 1)
                div_info = [ str(align_info[0]), str(align_info[1]), str(align_info[2]), str(read.query_alignment_length), str(read.get_tag("de")) ]

    # g.write("\t".join([data_type, caller, depth, str(ins_dis), str(te_dis), str("_".join(div_info)) ]) + "\t" + line.strip() + "\n")

    g.write("\t".join(line_info[:3]) + "\t" + depth + "\t" + "_".join(gold_ins) + "\t" + "_".join([caller, str(ins_dis), str(te_dis), str(seq_identity), str("_".join(div_info)) ]) + "\n")


g.close()
