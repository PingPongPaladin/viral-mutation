from mutation import *

np.random.seed(1)
random.seed(1)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='HIV sequence analysis')
    parser.add_argument('model_name', type=str,
                        help='Type of language model (e.g., hmm, lstm)')
    parser.add_argument('--namespace', type=str, default='hiv_pr',
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
    metas = {}
    for fname in meta_fnames:
        with open(fname) as f:
            for line in f:
                if not line.startswith('>'):
                    continue
                accession = line[1:].rstrip()
                fields = line.rstrip().split('.')
                subtype, country, year, strain = fields[:4]
                drug_treatment = fields[5] if "stanford" in fname else "LANL"

                if year == '-':
                    year = None
                else:
                    year = int(year)

                subtype = subtype.split('_')[-1]
                subtype = subtype.lstrip('>0123')

                keep_subtypes = {
                    'A', 'A1', 'A1A2', 'A1C', 'A1D', 'A2', 'A3', 'A6',
                    'AE', 'AG', 'B', 'C', 'BC', 'D',
                    'F', 'F1', 'F2', 'G', 'H', 'J',
                    'K', 'L', 'N', 'O', 'P', 'U',
                }
                if subtype not in keep_subtypes:
                    subtype = 'Other'

                metas[accession] = {
                    'subtype': subtype,
                    'country': country,
                    'year': year,
                    'strain': strain,
                    'drug_treatment': drug_treatment
                }
    return metas

def process(args, fnames, meta_fnames):
    metas = load_meta(meta_fnames)

    seqs = {}
    for fname in fnames:
        for record in SeqIO.parse(fname, 'fasta'):
            accession = record.description
            meta = metas[accession]
            meta['seqlen'] = len(str(record.seq))
            if args.namespace == 'hiva' and \
               (not meta['subtype'].startswith('A')):
                continue
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
    lanl_fname = "/afs/csail.mit.edu/u/a/andytso/meng/data/HIV-1_pr.fasta"
    stanford_fname = "/afs/csail.mit.edu/u/a/andytso/meng/data/stanford_pr.fasta"
    fnames = [lanl_fname, stanford_fname]
    meta_fnames = [lanl_fname, stanford_fname]

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
        for var in [ 'year', 'country', 'subtype' ]:
            tprint('\t{}:'.format(var))
            counts = Counter(adata_cluster.obs[var])
            for val, count in counts.most_common():
                tprint('\t\t{}: {}'.format(val, count))
        tprint('')

    cluster2subtype = {}
    for i in range(len(adata)):
        cluster = adata.obs['louvain'][i]
        if cluster not in cluster2subtype:
            cluster2subtype[cluster] = []
        cluster2subtype[cluster].append(adata.obs['subtype'][i])
    largest_pct_subtype = []
    for cluster in cluster2subtype:
        count = Counter(cluster2subtype[cluster]).most_common(1)[0][1]
        pct_subtype = float(count) / len(cluster2subtype[cluster])
        largest_pct_subtype.append(pct_subtype)
        tprint('\tCluster {}, largest subtype % = {}'
               .format(cluster, pct_subtype))
    tprint('Purity, Louvain and subtype: {}'
           .format(np.mean(largest_pct_subtype)))

def plot_umap(adata):
    sc.tl.umap(adata, min_dist=1.)
    single_drug_adata = adata[(adata.obs["drug_treatment"] != "None") &
                              (adata.obs["drug_treatment"] != "LANL"), :]
    sc.pl.umap(single_drug_adata, color='drug_treatment', save='_hiv_drug_treatment.png')
    sc.pl.umap(adata, color='louvain', save='_hiv_louvain.png')
    sc.pl.umap(adata, color='subtype', save='_hiv_subtype.png')


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

        from escape import load_rhee2004
        tprint('Rhee et al. 2004...')
        # protease
        pol_wt, pi_escape_seqs = load_rhee2004("PI")
        analyze_semantics(
            args, model, vocabulary, pol_wt, pi_escape_seqs,
            min_pos=54, max_pos=145,
            prob_cutoff=0., beta=1., plot_acquisition=True,
            plot_namespace="hiv_pol_pr",
        )

        # reverse transcriptase
        pol_wt, nrti_escape_seqs = load_rhee2004("NRTI")
        pol_wt, nnrti_escape_seqs = load_rhee2004("NNRTI")
        rt_escape_seqs = {**nrti_escape_seqs, **nnrti_escape_seqs}
        assert len(nrti_escape_seqs) + len(nnrti_escape_seqs) == \
            len(rt_escape_seqs)
        analyze_semantics(
            args, model, vocabulary, pol_wt, rt_escape_seqs,
            min_pos=209, max_pos=380,
            prob_cutoff=0., beta=1., plot_acquisition=True,
            plot_namespace="hiv_pol_rt",
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
