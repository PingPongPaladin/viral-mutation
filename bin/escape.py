import pandas as pd

from utils import Seq, SeqIO

def load_russ2020(escape_criteria="escape"):
    assert escape_criteria in ("escape", "antibiotic-resistance", "combination-resistance"), \
        "Invalid escape_criteria"

    base_path = "/afs/csail.mit.edu/u/a/andytso/meng/viral-mutation/data/beta_lactamase/escape_russ2020/"
    
    # load wt sequence
    wt_fpath = base_path + "ecoli_beta_lactamase_wt.fasta"
    records = list(SeqIO.parse(wt_fpath, "fasta"))
    assert len(records) == 1, "Expecting single wt sequence"
    wt_seq = records[0].seq
    
    # load mutations
    df = pd.read_excel(
        base_path + "beta_lactamase_inhibitory_concentrations.xlsx",
        index_col="Row",
        skiprows=[1])
    df = df.rename_axis(index="Mutant")
    df = df.applymap(lambda x: float(x[1:]) if isinstance(x, str) else x)
    
    # for alignment between russ et al. 2020 indices and that of UniProt
    #   https://www.nature.com/articles/s41467-020-15666-2
    #   https://www.uniprot.org/uniprot/P00811#sequences
    offset = -15
    
    # get the single-residue mutated sequence
    def mutate(seq, mutation, offset=offset):
        """
        mutation applied to seq where
        mutation is of the form <aa_original><index><aa_mutated>
        """
        i = int(mutation[1:-1]) - offset
        mutate_from = mutation[0]
        mutate_to = mutation[-1]
        assert mutate_from == seq[i]
        
        mut_seq = seq.tomutable()
        mut_seq[i] = mutate_to
        return mut_seq.toseq()
    
    # build seqs_escape where each mutation has a list of entries with fields
    # mutation, drug, drug IC50, combination IC50, significant
    seqs_escape = {}
    for mutation, row in df.iterrows():
        if mutation == "WT":
            seq = wt_seq
            index = None
        else:
            seq = mutate(wt_seq, mutation)
            index = int(mutation[1:-1]) - offset
        seqs_escape[seq] = []
        for drug in ("PIP", "ATM", "FEP"):
            y_label = "{}_AVI".format(drug)
            wt_x = df.loc["WT", "{}".format(drug)]
            wt_y = df.loc["WT", "{}_AVI".format(drug)]
            if escape_criteria == "escape":
                is_significant = row[drug] > wt_x and row[y_label] > wt_x
            elif escape_criteria == "antibiotic-resistance":
                is_significant = row[drug] > wt_x
            elif escape_criteria == "combination-resistance":
                is_significant = row[y_label] > wt_y
            else:
                raise ValueError
            seqs_escape[seq].append({
                "mutation": mutation,
                "pos": index,
                "drug": drug,
                "drug-ic50": row[drug],
                "combination-ic50": row[y_label],
                "significant": is_significant
            })
            
    return wt_seq, seqs_escape

def load_rhee2004(drug_type="PI"):
    assert drug_type in ("PI", "NRTI", "NNRTI"), "Invalid drug_type"
    
    # load wt sequence
    base_path = "/afs/csail.mit.edu/u/a/andytso/meng/viral-mutation/data/hiv/escape_rhee2004/"
    if drug_type == "PI":
        wt_fpath =  base_path + "hiv1-pr-wt.fasta"
    else:
        wt_fpath = base_path + "hiv1-rt-wt.fasta"
    records = list(SeqIO.parse(wt_fpath, "fasta"))
    assert len(records) == 1, "Expecting single wt sequence"
    wt_seq = records[0].seq
    
    # load mutations
    df = (pd.read_csv(base_path + "drug_resistance/{}.csv".format(drug_type.lower()))
          .set_index("Mutation Patterns"))
    
    # filter out multi-mutations and ignore None
    df = df[~df.index.str.contains(",")]
    df = df[df.index != "None"]
    
    # Stanford mutation indices found at
    #   https://hivdb.stanford.edu/pages/genotype-phenotype.html
    # are 1-based not 0-based
    offset = 1
    
    # get the single-residue mutated sequence
    def mutate(seq, mutation):
        i = int(mutation[:-1])
        mutate_to = mutation[-1]
        mut_seq = seq.tomutable()
        mut_seq[i - offset] = mutate_to
        return mut_seq.toseq()

    # build seqs_escape where each mutation has a list of entries with fields
    # mutation, drug, fold_change, significant
    seqs_escape = {}
    for mutation, row in df.iterrows():
        if mutation == "None":
            seq = wt_seq
        else:
            seq = mutate(wt_seq, mutation)
        seqs_escape[seq] = []
        for drug in row.index[1:-2]:
            seqs_escape[seq].append({
                "pos": int(mutation[:-1]) - offset,
                "drug": drug,
                "resistance fold change": row[drug],
                "significant": drug in row["Drugs Escaped"]
            })
            
    return wt_seq, seqs_escape


def load_doud2018(survival_cutoff=0.05):
    pos_map = {}
    with open('data/influenza/escape_doud2018/pos_map.csv') as f:
        f.readline() # Consume header.
        for line in f:
           fields = line.rstrip().split(',')
           pos_map[fields[1]] = int(fields[0]) - 1

    fname = 'data/influenza/escape_doud2018/WSN1933_H1_HA.fa'
    seqs = []
    for record in SeqIO.parse(fname, 'fasta'):
        seq = record.seq
        seqs.append(seq)

    seqs_escape = {}
    antibodies = [
        'C179', 'FI6v3', 'H17L10', 'H17L19', 'H17L7', 'S139',
    ]
    for antibody in antibodies:
        fname = ('data/influenza/escape_doud2018/' +
                 'medianfracsurvivefiles/' +
                 'antibody_{}_median.csv'.format(antibody))
        with open(fname) as f:
            f.readline() # Consume header.
            for line in f:
                fields = line.rstrip().split(',')
                frac_survived = float(fields[3])
                pos = pos_map[fields[0]]
                if seq[pos] != fields[1]:
                    print((seq[pos], fields[1], pos))
                assert(seq[pos] == fields[1])
                escaped = seq[:pos] + fields[2] + seq[pos + 1:]
                assert(len(seq) == len(escaped))
                if escaped not in seqs_escape:
                    seqs_escape[escaped] = []
                seqs_escape[escaped].append({
                    'frac_survived': frac_survived,
                    'antibody': antibody,
                    'significant': frac_survived > survival_cutoff,
                })

    return seq, seqs_escape

def load_lee2019():
    fname = 'data/influenza/escape_lee2019/Perth2009_H3_HA.fa'
    for record in SeqIO.parse(fname, 'fasta'):
        seq = record.seq
        break

    seqs_escape = {}
    fname = 'data/influenza/escape_lee2019/avg_sel_tidy.csv'
    with open(fname) as f:
        f.readline() # Consume header.
        for line in f:
            fields = line.rstrip().split(',')
            significant = fields[14] == 'True'
            pos = int(fields[13])
            assert(seq[pos] == fields[5])
            escaped = seq[:pos] + fields[6] + seq[pos + 1:]
            assert(len(seq) == len(escaped))
            if escaped not in seqs_escape:
                seqs_escape[escaped] = []

            if '-age-' in fields[0]:
                species = 'human'
            elif 'ferret-' in fields[0]:
                species = 'ferret'
            else:
                species = 'antibody'

            seqs_escape[escaped].append({
                'abs_diff_selection': float(fields[11]),
                'antibody': fields[1],
                'species': species,
                'significant': significant,
            })

    return seq, seqs_escape

def load_dingens2019(survival_cutoff=0.11):
    pos_map = {}
    with open('data/hiv/escape_dingens2019/BG505_to_HXB2.csv') as f:
        f.readline() # Consume header.
        for line in f:
           fields = line.rstrip().split(',')
           pos_map[fields[1]] = int(fields[0]) - 1

    fname = 'data/hiv/escape_dingens2019/Env_protalign_manualeditAD.fasta'
    for record in SeqIO.parse(fname, 'fasta'):
        if record.description == 'BG505':
            seq = record.seq
            break

    seqs_escape = {}
    antibodies = [
        '101074', '10E8', '3BNC117-101074-pool', '3BNC117', 'PG9',
        'PGT121', 'PGT145', 'PGT151', 'VRC01', 'VRC34',
    ]
    for antibody in antibodies:
        fname = ('data/hiv/escape_dingens2019/FileS4/'
                 'fracsurviveaboveavg/{}.csv'.format(antibody))
        with open(fname) as f:
            f.readline() # Consume header.
            for line in f:
                fields = line.rstrip().split(',')
                frac_survived = float(fields[3])
                pos = pos_map[fields[0]]
                assert(seq[pos] == fields[1])
                escaped = seq[:pos] + fields[2] + seq[pos + 1:]
                assert(len(seq) == len(escaped))
                if escaped not in seqs_escape:
                    seqs_escape[escaped] = []
                seqs_escape[escaped].append({
                    'pos': pos,
                    'frac_survived': frac_survived,
                    'antibody': antibody,
                    'significant': frac_survived > survival_cutoff,
                })

    return seq, seqs_escape

def load_baum2020():
    seq = SeqIO.read('data/cov/cov2_spike_wt.fasta', 'fasta').seq

    muts = [
        'K417E', 'K444Q', 'V445A', 'N450D', 'Y453F', 'L455F',
        'E484K', 'G485D', 'F486V', 'F490L', 'F490S', 'Q493K',
        'H655Y', 'R682Q', 'R685S', 'V687G', 'G769E', 'Q779K',
        'V1128A',
    ]

    seqs_escape = {}
    for mut in muts:
        aa_orig = mut[0]
        aa_mut = mut[-1]
        pos = int(mut[1:-1]) - 1
        assert(seq[pos] == aa_orig)
        escaped = seq[:pos] + aa_mut + seq[pos + 1:]
        assert(len(seq) == len(escaped))
        if escaped not in seqs_escape:
            seqs_escape[escaped] = []
        seqs_escape[escaped].append({
            'mutation': mut,
            'significant': True,
        })

    return seq, seqs_escape

def load_greaney2020(survival_cutoff=0.3,
                     binding_cutoff=-0.4, expr_cutoff=-0.4):
    seq = SeqIO.read('data/cov/cov2_spike_wt.fasta', 'fasta').seq

    sig_sites = set()
    with open('data/cov/greaney2020cov2/significant_escape_sites.csv') as f:
        f.readline()
        for line in f:
            fields = line.rstrip().split(',')
            sig_sites.add(int(fields[1]) - 1)

    binding = {}
    with open('data/cov/starr2020cov2/single_mut_effects.csv') as f:
        f.readline()
        for line in f:
            fields = line.rstrip().split(',')
            pos = float(fields[1]) - 1
            aa_orig = fields[2].strip('"')
            aa_mut = fields[3].strip('"')
            if aa_mut == '*':
                continue
            if fields[8] == 'NA':
                score = float('-inf')
            else:
                score = float(fields[8])
            if fields[11] == 'NA':
                expr = float('-inf')
            else:
                expr = float(fields[11])
            binding[(pos, aa_orig, aa_mut)] = score, expr

    seqs_escape = {}
    with open('data/cov/greaney2020cov2/escape_fracs.csv') as f:
        f.readline() # Consume header.
        for line in f:
            fields = line.rstrip().split(',')
            antibody = fields[2]
            escape_frac = float(fields[10])
            aa_orig = fields[5]
            aa_mut = fields[6]
            pos = int(fields[4]) - 1
            assert(seq[pos] == aa_orig)
            escaped = seq[:pos] + aa_mut + seq[pos + 1:]
            assert(len(seq) == len(escaped))
            if escaped not in seqs_escape:
                seqs_escape[escaped] = []
            significant = (
                escape_frac > survival_cutoff and
                pos in sig_sites and
                binding[(pos, aa_orig, aa_mut)][0] > binding_cutoff and
                binding[(pos, aa_orig, aa_mut)][1] > expr_cutoff
            )
            seqs_escape[escaped].append({
                'pos': pos,
                'frac_survived': escape_frac,
                'antibody': antibody,
                'significant': significant,
            })

    return seq, seqs_escape

if __name__ == '__main__':
    load_doud2018()
    load_lee2019()
    load_dingens2019()
    load_baum2020()
    load_greaney2020()
