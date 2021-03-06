from mutation import *

np.random.seed(1)
random.seed(1)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Beta-Lactamase sequence analysis')
    parser.add_argument('model_name', type=str,
                        help='Type of language model (e.g., hmm, lstm)')
    parser.add_argument('--namespace', type=str, default='beta_lactamase',
                        help='Model namespace')
    parser.add_argument('--dim', type=int, default=512,
                        help='Embedding dimension')
    parser.add_argument('--batch-size', type=int, default=1000,
                        help='Training minibatch size')
    parser.add_argument('--n-epochs', type=int, default=4,
                        help='Number of training epochs')
    parser.add_argument('--seed', type=int, default=1,
                        help='Random seed')
    parser.add_argument('--checkpoint', type=str, default=None,
                        help='Model checkpoint')
    parser.add_argument('--train', action='store_true',
                        help='Train model')
    parser.add_argument('--train-split', action='store_true',
                        help='Train model on portion of data')
    parser.add_argument('--test', action='store_true',
                        help='Test model')
    parser.add_argument('--embed', action='store_true',
                        help='Analyze embeddings')
    parser.add_argument('--semantics', action='store_true',
                        help='Analyze mutational semantic change')
    parser.add_argument('--combfit', action='store_true',
                        help='Analyze combinatorial fitness')
    args = parser.parse_args()
    return args

def load_meta(meta_fnames):
    """
    Parse meta information from fasta headers which have the following structure:

    db|UniqueIdentifier|EntryName ProteinName OS=OrganismName \
    OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
    """
    lineage_level_keys = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus-species']
    additional_info_keys = [
        "database", "accession_id",
        "entry_name", "protein",
        "OS", "OX", "GN", "PE", "SV"
    ]
    
    def parse_lineage(lineage, lineage_levels=lineage_level_keys):
        output = dict()
        for level, classification in zip(lineage_levels, lineage.split("; ")):
            output[level] = classification
        for level in lineage_levels[len(output):]:
            output[level] = "Unknown"
        assert set(output.keys()) == set(lineage_level_keys)
        return output
    
    def parse_additional_info(additional_info, df):
        """
        Parse additional info section of fasta header which follows
        the following format:
        
        EntryName ProteinName OS=OrganismName OX=OrganismIdentifier \
        [GN=GeneName ] PE=ProteinExistence SV=SequenceVersion
        """
        split_i = additional_info.find("=") - 3
        unkeyed_info = additional_info[:split_i]
        keyed_info = additional_info[split_i:]
        
        j = unkeyed_info.find(" ")
        entry_name = unkeyed_info[:j]
        protein = unkeyed_info[j + 1:]
        output = {
            "entry_name": entry_name,
            "protein": protein
        }
        
        curr = ""
        key = ""
        for i, c in enumerate(keyed_info):
            if c == "=":
                if key != "":
                    value = curr[:-3]
                    output[key] = value
                key = curr[-2:]
                curr = ""
            else:
                curr += c
        output[key] = curr

        # gene name may not be provided
        if "GN" not in output:
            output["GN"] = "Unknown"
        
        lineage = df.loc[int(output["OX"]), "Lineage"]
        output = {**output, **parse_lineage(lineage)}
        return output
    
    # load taxonomy id lookup table
    fname = "/afs/csail.mit.edu/u/a/andytso/meng/data/bacterial-taxonomy-lookup-table.txt"
    df = pd.read_csv(fname, sep='\t')
    df = df.set_index("Taxon")

    metas = {}
    keys_to_keep = lineage_level_keys + additional_info_keys
    for fname in meta_fnames:
        with open(fname) as f:
            for line in f:
                if not line.startswith('>'):
                    continue
                accession = line[1:].rstrip()
                database, accession_id, additional_info = accession.split("|")
                metas[accession] = {
                    "database": database,
                    "accession_id": accession_id,
                    **parse_additional_info(additional_info, df)
                }
                metas[accession] = {key: metas[accession][key] for key in keys_to_keep}
    return metas

def process(args, fnames, meta_fnames):
    metas = load_meta(meta_fnames)

    seqs = {}
    for fname in fnames:
        for record in SeqIO.parse(fname, 'fasta'):
            accession = record.description
            meta = metas[accession]
            meta['seqlen'] = len(str(record.seq))
            if 'X' in record.seq:
                continue
            if record.seq not in seqs:
                seqs[record.seq] = []
            seqs[record.seq].append(meta)
    return seqs

def split_seqs(seqs, split_method='random'):
    train_seqs, test_seqs = {}, {}

    old_cutoff = 1900
    new_cutoff = 2008

    tprint('Splitting seqs...')
    for seq in seqs:
        # Pick validation set based on date.
        seq_dates = [
            meta['year'] for meta in seqs[seq]
            if meta['year'] is not None
        ]
        if len(seq_dates) == 0:
            test_seqs[seq] = seqs[seq]
            continue
        if len(seq_dates) > 0:
            oldest_date = sorted(seq_dates)[0]
            if oldest_date < old_cutoff or oldest_date >= new_cutoff:
                test_seqs[seq] = seqs[seq]
                continue
        train_seqs[seq] = seqs[seq]
    tprint('{} train seqs, {} test seqs.'
           .format(len(train_seqs), len(test_seqs)))

    return train_seqs, test_seqs

def setup(args):
    uniprot_fname = "/afs/csail.mit.edu/u/a/andytso/meng/data/uniprot_filtered_beta_lactamase.fasta"
    fnames = [uniprot_fname]
    meta_fnames = [uniprot_fname]

    seqs = process(args, fnames, meta_fnames)

    seq_len = max([ len(seq) for seq in seqs ]) + 2
    vocab_size = len(AAs) + 2

    model = get_model(args, seq_len, vocab_size,
                      inference_batch_size=1000)

    return model, seqs

def interpret_clusters(adata):
    clusters = sorted(set(adata.obs['louvain']))
    for cluster in clusters:
        tprint('Cluster {}'.format(cluster))
        adata_cluster = adata[adata.obs['louvain'] == cluster]
        for var in [ 'phylum', 'protein' ]:
            tprint('\t{}:'.format(var))
            counts = Counter(adata_cluster.obs[var])
            for val, count in counts.most_common():
                tprint('\t\t{}: {}'.format(val, count))
        tprint('')

    for attribute in [ 'phylum', 'protein' ]:
        cluster2attribute = {}
        for i in range(len(adata)):
            cluster = adata.obs['louvain'][i]
            if cluster not in cluster2attribute:
                cluster2attribute[cluster] = []
            cluster2attribute[cluster].append(adata.obs['subtype'][i])
        largest_pct_attribute = []
        for cluster in cluster2attribute:
            count = Counter(cluster2attribute[cluster]).most_common(1)[0][1]
            pct_attribute = float(count) / len(cluster2attribute[cluster])
            largest_pct_attribute.append(pct_attribute)
            tprint('\tCluster {}, largest {} % = {}'
                   .format(cluster, attribute, pct_attribute))
        tprint('Purity, Louvain and {}: {}'
               .format(attribute, np.mean(largest_pct_attribute)))

def plot_umap(adata):
    sc.tl.umap(adata, min_dist=1.)
    sc.pl.umap(adata, color='louvain', save='_beta_lactamase_louvain.png')
    
    # plot umap of top 10 most common phylum
    most_common_phylum = [
        'Proteobacteria',
        'Firmicutes',
        'Actinobacteria',
        'Bacteroidetes',
        'environmental samples',
        'Cyanobacteria',
        'Chloroflexi',
        'Acidobacteria',
        'Unknown',
        'Spirochaetes'
    ]
    mask = np.array([False] * len(adata.obs["phylum"]))
    for phylum in most_common_phylum:
        mask |= (adata.obs["phylum"] == phylum)
    filtered_phylum_adata = adata[mask, :]
    sc.pl.umap(filtered_phylum_adata, color='phylum', save='_beta_lactamase_phylum.png')

    # plot umap of top 10 most common protein
    most_common_protein = [
        'Beta-lactamase',
        'Ribonuclease J',
        'Hydroxyacylglutathione hydrolase',
        'Beta-lactamase (Fragment)',
        'TPR_REGION domain-containing protein',
        'MBL fold metallo-hydrolase',
        'Uncharacterized protein',
        'B2 metallo-beta-lactamase',
        'Tetratricopeptide repeat protein',
        'Hydroxyacylglutathione hydrolase (Fragment)'
    ]
    mask = np.array([False] * len(adata.obs["protein"]))
    for protein in most_common_protein:
        mask |= (adata.obs["protein"] == protein)
    filtered_protein_adata = adata[mask, :]
    sc.pl.umap(filtered_protein_adata, color='protein', save='_beta_lactamase_protein.png')


def populate_embedding(args, model, seqs, vocabulary,
                       use_cache=False, namespace=None):
    if namespace is None:
        namespace = args.namespace

    if use_cache:
        mkdir_p('target/{}/embedding'.format(namespace))
        embed_prefix = ('target/{}/embedding/{}_{}'
                        .format(namespace, args.model_name, args.dim))

    sorted_seqs = np.array([ str(s) for s in sorted(seqs.keys()) ])
    batch_size = 3000
    n_batches = math.ceil(len(sorted_seqs) / float(batch_size))
    for batchi in range(n_batches):
        # Identify the batch.
        start = batchi * batch_size
        end = (batchi + 1) * batch_size
        sorted_seqs_batch = sorted_seqs[start:end]
        seqs_batch = { seq: seqs[seq] for seq in sorted_seqs_batch }

        # Load from cache if available.
        if use_cache:
            embed_fname = embed_prefix + '.{}.npy'.format(batchi)
            if os.path.exists(embed_fname):
                X_embed = np.load(embed_fname, allow_pickle=True)
                if X_embed.shape[0] == len(sorted_seqs_batch):
                    for seq_idx, seq in enumerate(sorted_seqs_batch):
                        for meta in seqs[seq]:
                            meta['embedding'] = X_embed[seq_idx]
                    continue

        # Embed the sequences.
        seqs_batch = embed_seqs(args, model, seqs_batch, vocabulary,
                                use_cache=False)
        if use_cache:
            X_embed = []
        for seq in sorted_seqs_batch:
            for meta in seqs[seq]:
                meta['embedding'] = seqs_batch[seq][0]['embedding'].mean(0)
            if use_cache:
                X_embed.append(seqs[seq][0]['embedding'].ravel())
        del seqs_batch

        if use_cache:
            np.save(embed_fname, np.array(X_embed))

    return seqs

def analyze_embedding(args, model, seqs, vocabulary):
    seqs = populate_embedding(args, model, seqs, vocabulary,
                              use_cache=True)

    X, obs = [], {}
    obs['n_seq'] = []
    obs['seq'] = []
    for seq in seqs:
        meta = seqs[seq][0]
        X.append(meta['embedding'])
        for key in meta:
            if key == 'embedding':
                continue
            if key not in obs:
                obs[key] = []
            obs[key].append(Counter([
                meta[key] for meta in seqs[seq]
            ]).most_common(1)[0][0])
        obs['n_seq'].append(len(seqs[seq]))
        obs['seq'].append(str(seq))
    X = np.array(X)

    adata = AnnData(X)
    for key in obs:
        adata.obs[key] = obs[key]

    sc.pp.neighbors(adata, n_neighbors=200, use_rep='X')
    sc.tl.louvain(adata, resolution=1.)

    sc.set_figure_params(dpi_save=500)
    plot_umap(adata)

    interpret_clusters(adata)

if __name__ == '__main__':
    args = parse_args()

    AAs = [
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
        'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
        'Y', 'V', 'X', 'Z', 'J', 'U', 'B',
    ]
    vocabulary = { aa: idx + 1 for idx, aa in enumerate(sorted(AAs)) }

    model, seqs = setup(args)

    if args.checkpoint is not None:
        model.model_.load_weights(args.checkpoint)
        tprint('Model summary:')
        tprint(model.model_.summary())

    if args.train:
        batch_train(args, model, seqs, vocabulary, batch_size=5000)

    if args.train_split or args.test:
        train_test(args, model, seqs, vocabulary, split_seqs)

    if args.embed:
        if args.checkpoint is None and not args.train:
            raise ValueError('Model must be trained or loaded '
                             'from checkpoint.')
        no_embed = { 'hmm' }
        if args.model_name in no_embed:
            raise ValueError('Embeddings not available for models: {}'
                             .format(', '.join(no_embed)))
        analyze_embedding(args, model, seqs, vocabulary)

    if args.semantics:
        if args.checkpoint is None and not args.train:
            raise ValueError('Model must be trained or loaded '
                             'from checkpoint.')

        from escape import load_russ2020
        tprint('Russ et al. 2020...')
        seq_to_mutate, escape_seqs = load_russ2020(escape_criteria="combination-resistance")
        min_pos = 0
        max_pos = len(seq_to_mutate) - 1
        analyze_semantics(
            args, model, vocabulary, seq_to_mutate, escape_seqs,
            min_pos=min_pos, max_pos=max_pos,
            prob_cutoff=0., beta=1., plot_acquisition=True,
            plot_namespace="beta_lactamase"
        )

    if args.combfit:
        from combinatorial_fitness import load_haddox2018
        tprint('Haddox et al. 2018...')
        wt_seqs, seqs_fitness = load_haddox2018()
        strains = sorted(wt_seqs.keys())
        for strain in strains:
            analyze_comb_fitness(args, model, vocabulary,
                                 strain, wt_seqs[strain], seqs_fitness,
                                 prob_cutoff=0., beta=1.)
