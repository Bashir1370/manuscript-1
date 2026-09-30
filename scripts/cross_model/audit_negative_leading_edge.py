"""Independent reconstruction of original ranks and negative leading-edge CSVs.

Run from repository root with Python 3 (standard library only). This does not run
GSEA/GSVA, generate plots, or claim R execution. Original per-study FDR is retained.
"""
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STUDIES = ('OIPN_GSE160543', 'NC_GSE246156', 'CCI_GSE212311')
BASES = ('results/GSE160543_Oxaliplatin_vs_Vehicle',
         'results/GSE246156_NC_L5_day7', 'results/GSE212311_CCI_L4L6_day11')
SUBS = ('Pathway_analysis', 'pathway_analysis', 'pathway_analysis_source_aware')
INPUTS = []

def read(path):
    INPUTS.append(path)
    with (ROOT / path).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def number(value):
    if value in ('', 'NA', None):
        return None
    result = float(value)
    assert math.isfinite(result), value
    return result

def evidence(row, symbol, feature, source):
    fdr = number(row['padj'])
    assert fdr is None or 0 <= fdr <= 1
    return dict(symbol=symbol, feature_id=feature, statistic=number(row['stat']),
                log2FC=number(row['log2FoldChange']), gene_FDR=fdr, mapping_source=source)

def write(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({k: ('TRUE' if v else 'FALSE') if isinstance(v, bool)
                        else '' if v is None else v for k, v in row.items()})

def load_studies():
    tables, gseas = {}, {}
    for study, base, sub in zip(STUDIES, BASES, SUBS):
        gsea = read(f'{base}/{sub}/GSEA_Hallmark_results.csv')
        assert len(gsea) == 50 and len({x['pathway'] for x in gsea}) == 50
        gseas[study] = {x['pathway']: dict(NES=number(x['NES']), FDR=number(x['padj']),
                           edges=set(x['leadingEdge'].split(';'))) for x in gsea}
        if study == STUDIES[0]:
            raw = read(f'{base}/DE_all_genes.csv')
            raw = [x for x in raw if x['symbol'] and number(x['stat']) is not None]
            raw.sort(key=lambda x: -number(x['stat']))  # stable source tie order
            candidates = [evidence(x, x['symbol'], x['gene_id'], 'original_DE_symbol') for x in raw]
        elif study == STUDIES[1]:
            raw = read(f'{base}/DESeq2/GSE246156_Compression_vs_Sham_DE_all_genes.csv')
            mapping = read(f'{base}/{sub}/Ensembl_annotation_all_mappings.csv')
            pairs = defaultdict(set)
            for x in mapping:
                if x['SYMBOL'] and x['SYMBOL'] != 'NA':
                    pairs[x['ENSEMBL']].add(x['SYMBOL'])
            candidates = []
            for x in raw:
                ens = re.sub(r'\.[0-9]+$', '', x['gene'])
                if len(pairs[ens]) == 1 and number(x['stat']) is not None:
                    candidates.append(evidence(x, next(iter(pairs[ens])), x['gene'],
                                               'original_unambiguous_Ensembl_symbol'))
            candidates.sort(key=lambda x: (-abs(x['statistic']), re.sub(r'\.[0-9]+$', '', x['feature_id'])))
        else:
            raw = read(f'{base}/{sub}/feature_to_symbol_audit.csv')
            candidates = [evidence(x, x['symbol'], x['gene_id'], x['source']) for x in raw]
            assert len({x['symbol'] for x in candidates}) == len(candidates)
        table = {}
        for x in candidates:
            table.setdefault(x['symbol'], x)
        rank = read(f'{base}/{sub}/ranked_statistics.csv')
        assert len(rank) == len(table) and {x['symbol'] for x in rank} == set(table)
        for x in rank:
            assert math.isclose(number(x['statistic']), table[x['symbol']]['statistic'], rel_tol=1e-8, abs_tol=1e-8)
        for g in gseas[study].values():
            assert g['edges'] <= table.keys() and '' not in g['edges']
        tables[study] = table
    return tables, gseas

def reconstruct(tables, gseas, selected):
    long, summary, totals = [], [], []
    for p in selected:
        union = sorted(set.union(*(gseas[s][p]['edges'] for s in STUDIES)))
        current = []
        for symbol in union:
            rows, member = [], []
            for study in STUDIES:
                gene = tables[study].get(symbol)
                inside = symbol in gseas[study][p]['edges']
                member.append(inside)
                row = dict(pathway=p, symbol=symbol, study=study, rank_available=gene is not None,
                           leading_edge=inside if gene else None,
                           status='unavailable_to_GSEA' if gene is None else
                                  'leading_edge' if inside else 'ranked_not_leading_edge',
                           feature_id=gene['feature_id'] if gene else None,
                           mapping_source=gene['mapping_source'] if gene else None,
                           rank_statistic=gene['statistic'] if gene else None,
                           log2FC=gene['log2FC'] if gene else None,
                           gene_FDR=gene['gene_FDR'] if gene else None,
                           pathway_NES=gseas[study][p]['NES'], pathway_FDR=gseas[study][p]['FDR'])
                long.append(row); rows.append(row)
            ranked = all(x['rank_available'] for x in rows)
            lfc = [x['log2FC'] for x in rows]; fdr = [x['gene_FDR'] for x in rows]
            result = dict(pathway=p, symbol=symbol, n_ranked=sum(x['rank_available'] for x in rows),
                          n_leading_edge=sum(member),
                          n_LE_in_significant_pathways=sum(v and gseas[s][p]['FDR'] < .05 for s, v in zip(STUDIES, member)),
                          n_positive_log2FC=sum(x is not None and x > 0 for x in lfc),
                          n_negative_log2FC=sum(x is not None and x < 0 for x in lfc),
                          n_gene_FDR_lt_0_05=sum(x is not None and x < .05 for x in fdr),
                          all_three_ranked=ranked, shared_ge3=all(member), shared_all3=all(member),
                          positive_all3=ranked and all(x > 0 for x in lfc),
                          negative_all3=ranked and all(x < 0 for x in lfc))
            for s, x in zip(STUDIES, rows):
                result.update({s+'_LE': x['leading_edge'], s+'_log2FC': x['log2FC'], s+'_gene_FDR': x['gene_FDR']})
            summary.append(result); current.append(result)
        totals.append(dict(pathway=p, union_LE_genes=len(union),
                           shared_ge3=sum(x['shared_ge3'] for x in current),
                           shared_ge3_all3_ranked=sum(x['shared_ge3'] and x['all_three_ranked'] for x in current),
                           shared_all3=sum(x['shared_all3'] for x in current),
                           shared_ge3_positive_all3=sum(x['shared_ge3'] and x['positive_all3'] for x in current),
                           shared_ge3_negative_all3=sum(x['shared_ge3'] and x['negative_all3'] for x in current)))
    return long, summary, totals

def priorities(shared):
    genes = []
    for symbol in sorted({x['symbol'] for x in shared}):
        rows = [x for x in shared if x['symbol'] == symbol]
        cols = [s+'_log2FC' for s in STUDIES] + [s+'_gene_FDR' for s in STUDIES]
        for col in cols:
            assert len({x[col] for x in rows}) == 1
        g = {'symbol': symbol, **{col: rows[0][col] for col in cols}}
        sig = [g[s+'_gene_FDR'] is not None and g[s+'_gene_FDR'] < .05 for s in STUDIES]
        g.update(n_shared_pathways=len(rows), shared_pathways=';'.join(sorted(x['pathway'] for x in rows)),
                 negative_all3=all(g[s+'_log2FC'] < 0 for s in STUDIES),
                 n_gene_FDR_lt_0_05=sum(sig), OIPN_significant=sig[0], NC_significant=sig[1], CCI_significant=sig[2])
        g.update(selected_OIPN_plus_physical=g['negative_all3'] and sig[0] and (sig[1] or sig[2]),
                 strict_significant_all3=g['negative_all3'] and all(sig),
                 physical_only_significant=g['negative_all3'] and not sig[0] and sig[1] and sig[2],
                 support_pattern='OIPN_NC_CCI' if all(sig) else 'OIPN_NC' if sig[0] and sig[1] else
                                 'OIPN_CCI' if sig[0] and sig[2] else 'NC_CCI' if sig[1] and sig[2] else 'fewer_than_two_significant')
        genes.append(g)
    return genes

def main():
    selected = sorted(x['pathway'] for x in read('results/manuscript_hallmark_three_dataset/shared_negative.csv'))
    assert selected == ['HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION']
    tables, gseas = load_studies()
    assert all(gseas[s][p]['NES'] < 0 and gseas[s][p]['FDR'] < .05 for p in selected for s in STUDIES)
    long, summary, totals = reconstruct(tables, gseas, selected)
    shared = [x for x in summary if x['shared_all3']]
    genes = priorities(shared)
    assert genes, 'No shared genes: handle empty priority schema before exporting.'
    priority = [x for x in genes if x['selected_OIPN_plus_physical']]
    strict = [x for x in genes if x['strict_significant_all3']]
    two = [x for x in priority if x['n_gene_FDR_lt_0_05'] == 2]
    physical = [x for x in genes if x['physical_only_significant']]
    # Reconstruct archived positive findings as a regression control on rank/joins.
    positive = sorted(p for p in gseas[STUDIES[0]] if all(gseas[s][p]['NES'] > 0 and gseas[s][p]['FDR'] < .05 for s in STUDIES))
    _, pos, _ = reconstruct(tables, gseas, positive)
    pos_shared = [x for x in pos if x['shared_all3']]
    assert len(positive) == 10 and len(pos_shared) == 109 and len({x['symbol'] for x in pos_shared}) == 78
    out = ROOT / 'results/shared_negative_Hallmark_leading_edge_three_dataset'
    for name, rows in [('gene_evidence_long.csv', long), ('gene_membership_summary.csv', summary),
                       ('shared_ge3.csv', shared), ('shared_all3.csv', shared), ('pathway_overlap_summary.csv', totals)]:
        write(out / name, rows)
    write(out / 'pathway_evidence.csv', [dict(study=s, pathway=p, NES=gseas[s][p]['NES'],
          pathway_FDR=gseas[s][p]['FDR'], leading_edge_genes=len(gseas[s][p]['edges']), ranked_symbols=len(tables[s]))
          for s in STUDIES for p in selected])
    write(out / 'shared_gene_pathway_counts.csv', [dict(symbol=x['symbol'], n_selected_pathways=x['n_shared_pathways']) for x in genes])
    for name, rows in [('all_shared_genes_with_priority_flags.csv', genes), ('priority_OIPN_plus_physical.csv', priority),
                       ('strict_significant_all3.csv', strict), ('priority_significant_exactly2.csv', two),
                       ('excluded_physical_only_significant.csv', physical)]:
        write(out / 'gene_prioritization' / name, rows, list(genes[0]))
    memberships = sorted([x for x in shared if x['symbol'] in {g['symbol'] for g in priority}], key=lambda x:(x['symbol'],x['pathway']))
    write(out / 'gene_prioritization/priority_pathway_memberships.csv', memberships, list(summary[0]))
    for name, rows in [('priority_genes_STRING.txt', priority), ('strict_genes_STRING.txt', strict)]:
        (out / 'gene_prioritization' / name).write_text(''.join(x['symbol']+'\n' for x in rows))
    counts = dict(shared_unique_genes=len(genes), negative_all3=sum(x['negative_all3'] for x in genes),
                  FDR_ge2_any_pair_negative_all3=sum(x['negative_all3'] and x['n_gene_FDR_lt_0_05'] >= 2 for x in genes),
                  priority_OIPN_plus_physical=len(priority), strict_significant_all3=len(strict),
                  priority_significant_exactly2=len(two), excluded_physical_only=len(physical), priority_pathway_memberships=len(memberships))
    write(out / 'gene_prioritization/selection_summary.csv', [dict(metric=k, count=v) for k,v in counts.items()])
    write(out / 'input_checksums.csv', [dict(input=p, md5=hashlib.md5((ROOT/p).read_bytes()).hexdigest()) for p in INPUTS])
    write(out / 'gene_prioritization/input_checksums.csv', [dict(input=str((out/'shared_all3.csv').relative_to(ROOT)),
          md5=hashlib.md5((out/'shared_all3.csv').read_bytes()).hexdigest())])
    write(out / 'plot_settings.csv', [dict(color_scale='linear_log2FC', lower_limit=-6, upper_limit=6,
          shared_across_studies=True, shared_across_pathways=True, display_clipping_only=True)])
    checks = dict(source='Python reconstruction of archived inputs; no R run or plot generation',
                  pathway_counts=totals, selection_counts=counts, priority_genes=[x['symbol'] for x in priority],
                  strict_genes=[x['symbol'] for x in strict], original_rank_matches=True,
                  positive_regression=dict(pathways=10, memberships=109, genes=78),
                  union_memberships=len(summary), evidence_rows=len(long), shared_memberships=len(shared))
    (out/'validation.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(checks, indent=2))

if __name__ == '__main__':
    main()
