def is_valid_dna(sequence: str) -> bool:
    """Check if a string represents a valid DNA sequence (ACGT only)."""
    return all(base.upper() in {'A', 'C', 'G', 'T'} for base in sequence)

def reverse_complement(dna: str) -> str:
    """Generate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna.upper()))

def gc_content(sequence: str) -> float:
    """Calculate GC content percentage of a DNA sequence."""
    gc = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc / len(sequence)) * 100 if sequence else 0.0

def translate_dna(dna: str, frame: int = 0) -> str:
    """Translate DNA sequence to protein (standard genetic code)."""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    protein = []
    for i in range(frame, len(dna)-2, 3):
        codon = dna[i:i+3].upper()
        protein.append(codon_table.get(codon, 'X'))
    return ''.join(protein)

def find_motifs(dna: str, motif: str) -> list[int]:
    """Find all starting positions of a motif in DNA sequence."""
    return [i for i in range(len(dna)-len(motif)+1) 
            if dna[i:i+len(motif)].upper() == motif.upper()]

def fasta_to_dict(fasta_str: str) -> dict:
    """Parse FASTA string into dictionary {header: sequence}."""
    entries = fasta_str.split('>')[1:]
    return {entry.split('\n')[0]: ''.join(entry.split('\n')[1:]) 
            for entry in entries}

def count_kmers(sequence: str, k: int) -> dict:
    """Count all k-mers in a DNA sequence."""
    return {sequence[i:i+k]: sequence.count(sequence[i:i+k]) 
            for i in range(len(sequence)-k+1)}

def calculate_entropy(sequence: str) -> float:
    """Calculate Shannon entropy of a DNA sequence."""
    from math import log2
    counts = {'A':0, 'C':0, 'G':0, 'T':0}
    for base in sequence.upper():
        counts[base] += 1
    total = len(sequence)
    return -sum((v/total)*log2(v/total) for v in counts.values() if v)

def vcf_to_dataframe(vcf_lines: list[str]):
    """Convert VCF file lines to pandas DataFrame (requires pandas)."""
    import pandas as pd
    headers = [line for line in vcf_lines if line.startswith('#CHROM')]
    if not headers:
        return pd.DataFrame()
    columns = headers[0].strip('#').split()
    data = [line.split() for line in vcf_lines if not line.startswith('#')]
    return pd.DataFrame(data, columns=columns)

def filter_low_quality_variants(vcf_df, min_qual: float = 20.0):
    """Filter VCF DataFrame by quality score (requires pandas)."""
    import pandas as pd
    vcf_df['QUAL'] = pd.to_numeric(vcf_df['QUAL'], errors='coerce')
    return vcf_df[vcf_df['QUAL'] >= min_qual].copy()

def codon_usage(dna: str) -> dict:
    """Calculate codon usage frequency in a DNA sequence."""
    from collections import defaultdict
    counts = defaultdict(int)
    for i in range(0, len(dna)-2, 3):
        codon = dna[i:i+3].upper()
        counts[codon] += 1
    total = sum(counts.values())
    return {k: v/total for k, v in counts.items()}

def calculate_tm(sequence: str) -> float:
    """Calculate melting temperature (Wallace rule) for primers."""
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return 2*(a_count + t_count) + 4*(g_count + c_count)
