import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import matplotlib.pyplot as plt

# 读取BCR文件
contig_annotations = pd.read_csv('filtered_contig_annotations.csv')
airr_rearrangement = pd.read_csv('airr_rearrangement.tsv', sep='\t')

# 提取BCR序列
sequences = {}
for record in SeqIO.parse('filtered_contig.fasta', 'fasta'):
    sequences[record.id] = str(record.seq)

# 参考序列，可以从IMGT下载
reference_sequences = {}
for record in SeqIO.parse('IGHV.fasta', 'fasta'):
    header = record.description.split('|')
    gene_id = header[1].split('*')[0]  # 使用 IGHV 基因名作为键
    reference_sequences[gene_id] = str(record.seq).upper()
#print(reference_sequences)

# 初始化PairwiseAligner
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")  # 使用合适的替代矩阵
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def clean_sequence(sequence):
    return ''.join([base if base in 'ACGT' else 'N' for base in sequence])

# 比对和识别突变
def calculate_shm(query_seq, reference_seq):
    query_seq = clean_sequence(query_seq)
    reference_seq=clean_sequence(reference_seq)
    alignments = aligner.align(Seq(query_seq), Seq(reference_seq))
    alignment = alignments[0]
    query_aligned = alignment.aligned[0]
    reference_aligned = alignment.aligned[1]
    mutations = sum(query_seq[start:end] != reference_seq[start:end]
                    for start, end in zip(reference_aligned[0], reference_aligned[1]))
    mutation_frequency = (mutations / len(reference_seq)) * 100
    return mutations, mutation_frequency

# 计算每个序列的SHM
results = []
for index, row in contig_annotations.iterrows():
    contig_id = row['contig_id']
    v_gene = row['v_gene']
    query_seq = sequences.get(contig_id, '')
    if v_gene in reference_sequences:
        reference_seq = reference_sequences[v_gene]
        mutations, mutation_frequency = calculate_shm(query_seq, reference_seq)
        results.append({
            'contig_id': contig_id,
            'v_gene': v_gene,
            'mutations': mutations,
            'mutation_frequency': mutation_frequency
        })

# 生成报告
shm_report = pd.DataFrame(results)
print(shm_report)

# 突变频率分布图
plt.figure(figsize=(10, 6))
plt.hist(shm_report['mutation_frequency'], bins=30, edgecolor='k', alpha=0.7)
plt.title('Distribution of Mutation Frequency')
plt.xlabel('Mutation Frequency (%)')
plt.ylabel('Number of Sequences')
plt.grid(True)
plt.show()

# 突变数目分布图
plt.figure(figsize=(10, 6))
plt.hist(shm_report['mutations'], bins=30, edgecolor='k', alpha=0.7)
plt.title('Distribution of Number of Mutations')
plt.xlabel('Number of Mutations')
plt.ylabel('Number of Sequences')
plt.grid(True)
plt.show()


# 突变频率与突变数目的散点图
plt.figure(figsize=(10, 6))
plt.scatter(shm_report['mutations'], shm_report['mutation_frequency'], alpha=0.6)
plt.title('Mutations vs Mutation Frequency')
plt.xlabel('Number of Mutations')
plt.ylabel('Mutation Frequency (%)')
plt.grid(True)
plt.show()