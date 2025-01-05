from Bio import SeqIO
from collections import Counter
import pandas as pd


def analyze_dna(file_path):
    # بارگذاری دنباله از فایل GenBank
    record = SeqIO.read(file_path, "genbank")
    dna_sequence = str(record.seq).upper()
    
    if len(dna_sequence) < 2000:
        raise ValueError("طول دنباله DNA باید حداقل 2000 نوکلئوتید باشد.")
    
    # شمارش نوکلئوتیدهای ساده
    simple_counts = Counter(dna_sequence)
    simple_percentages = {nuc: count / len(dna_sequence) * 100 for nuc, count in simple_counts.items()}
    
    # شمارش دی‌نوکلئوتیدها
    di_nucleotides = [dna_sequence[i:i+2] for i in range(len(dna_sequence) - 1)]
    di_counts = Counter(di_nucleotides)
    di_percentages = {dinuc: count / len(di_nucleotides) * 100 for dinuc, count in di_counts.items()}
    
    # شمارش تری‌نوکلئوتیدها
    tri_nucleotides = [dna_sequence[i:i+3] for i in range(len(dna_sequence) - 2)]
    tri_counts = Counter(tri_nucleotides)
    tri_percentages = {trinuc: count / len(tri_nucleotides) * 100 for trinuc, count in tri_counts.items()}
    
    # ساخت جداول
    simple_table = pd.DataFrame({
        "Nucleotide": list(simple_counts.keys()),
        "Count": list(simple_counts.values()),
        "Percentage": list(simple_percentages.values())
    })
    
    di_table = pd.DataFrame({
        "Di-nucleotide": list(di_counts.keys()),
        "Count": list(di_counts.values()),
        "Percentage": list(di_percentages.values())
    })
    
    tri_table = pd.DataFrame({
        "Tri-nucleotide": list(tri_counts.keys()),
        "Count": list(tri_counts.values()),
        "Percentage": list(tri_percentages.values())
    })
    
    # محاسبه درصد GC
    gc_content = (simple_counts.get('G', 0) + simple_counts.get('C', 0)) / len(dna_sequence) * 100
    
    return simple_table, di_table, tri_table, gc_content

# استفاده از تابع
file_path = "sequence.gb"
simple_table, di_table, tri_table, gc_content = analyze_dna(file_path)

# نمایش نتایج
print("جدول نوکلئوتیدهای ساده:")
print(simple_table)

print("\nجدول دی‌نوکلئوتیدها:")
print(di_table)

print("\nجدول تری‌نوکلئوتیدها:")
print(tri_table)

print(f"\nدرصد GC: {gc_content:.2f}%")
