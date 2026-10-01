#!/usr/bin/env python3
"""Extract an explicitly requested gene panel from archived DE, without refitting."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/targeted_pain_gene_panel_three_dataset'
GENES = 'ADORA1 ANO1 ASIC3 BACE1 CARTPT COMT CXCL12 DISC1 EPHB1 FYN HTR2A ITGA2 KCNA1 LXN MMP24 NR2F6 NTRK1 NTSR1 PHF24 PRDM12 SCN11A SCN1A SCN9A SCRN3 TAC1 TAC4 TMEM120A TNF TRPA1 TRPV1'.split()
PATHS = {
    'OIPN': 'results/GSE160543_Oxaliplatin_vs_Vehicle/DE_all_genes.csv',
    'NC': 'results/GSE246156_NC_L5_day7/DESeq2/GSE246156_Compression_vs_Sham_DE_all_genes.csv',
    'NC_mapping': 'results/GSE246156_NC_L5_day7/pathway_analysis/Ensembl_annotation_all_mappings.csv',
    'CCI_mapping': 'results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/feature_to_symbol_audit.csv',
    'CCI': 'results/GSE212311_CCI_L4L6_day11/DESeq2/GSE212311_CCI_vs_Sham_DE_all_genes.csv',
}

def read(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None

def write(rows, name):
    with (OUT / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def run():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {key: ROOT / value for key, value in PATHS.items()}
    # A verified exact subset supports the hosted audit when the full CCI file
    # is not downloaded locally. Full checkouts always use the original file.
    if not paths['CCI'].exists():
        paths['CCI'] = OUT / 'CCI_selected_original_DE.csv'
    tables = {key: read(path) for key, path in paths.items()}
    mapping = {}
    for r in tables['NC_mapping']:
        if r['SYMBOL'] not in ('', 'NA'):
            mapping.setdefault(r['ENSEMBL'], set()).add(r['SYMBOL'])
    nc = {}
    for r in tables['NC']:
        symbols = mapping.get(r['gene'].split('.')[0], set())
        if len(symbols) == 1:
            nc.setdefault(next(iter(symbols)).upper(), []).append(r)
    oipn = {}
    for r in tables['OIPN']:
        oipn.setdefault(r['symbol'].upper(), []).append(r)
    cci_de = {r['gene_id']: r for r in tables['CCI']}
    assert len(cci_de) == len(tables['CCI']), 'Duplicate CCI feature IDs'
    cci = {}
    for r in tables['CCI_mapping']:
        if r['symbol'].upper() not in GENES:
            continue
        assert r['gene_id'] in cci_de, 'Mapped CCI feature absent from DE'
        d = cci_de[r['gene_id']]
        for field in ('log2FoldChange', 'padj', 'stat'):
            a, b = number(r[field]), number(d[field])
            assert (a is None and b is None) or (a is not None and b is not None and math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)), 'CCI mapping/DE mismatch'
        cci.setdefault(r['symbol'].upper(), []).append({**d, 'symbol': r['symbol'], 'mapping_source': r['source']})
    long, wide = [], []
    for gene in GENES:
        row = {'requested_gene': gene}
        for study, data in (('OIPN', oipn), ('NC', nc), ('CCI', cci)):
            matches = data.get(gene, [])
            if len(matches) > 1:
                raise ValueError(f'Ambiguous tested feature: {gene}, {study}; do not choose by P value')
            d = matches[0] if matches else None
            lfc = number(d['log2FoldChange']) if d else None
            fdr = number(d['padj']) if d else None
            status = 'not_in_archived_symbol_mapping' if d is None else 'FDR_unavailable' if fdr is None else 'significant' if fdr < .05 else 'not_significant'
            result = {'requested_gene': gene, 'study': study, 'feature_id': (d.get('gene_id', d.get('gene', '')) if d else ''), 'mapping_source': (d.get('mapping_source', 'original_DE_symbol' if study == 'OIPN' else 'archived_unambiguous_Ensembl') if d else ''), 'baseMean': number(d['baseMean']) if d else None, 'log2FC': lfc, 'lfcSE': number(d['lfcSE']) if d else None, 'stat': number(d['stat']) if d else None, 'pvalue': number(d['pvalue']) if d else None, 'FDR': fdr, 'direction': ('up' if lfc > 0 else 'down' if lfc < 0 else 'zero') if lfc is not None else '', 'status': status}
            long.append(result)
            for field in ('log2FC', 'FDR', 'status'):
                row[f'{study}_{field}'] = result[field]
        wide.append(row)
    assert len(long) == 90 and len(wide) == 30
    write(long, 'gene_panel_long.csv')
    write(wide, 'gene_panel_comparison.csv')
    counts = {s: {status: sum(r['study'] == s and r['status'] == status for r in long) for status in ('significant', 'not_significant', 'FDR_unavailable', 'not_in_archived_symbol_mapping')} for s in ('OIPN', 'NC', 'CCI')}
    checksums = []
    for key, path in paths.items():
        b = path.read_bytes()
        checksums.append({'input': str(path.relative_to(ROOT)), 'MD5': hashlib.md5(b).hexdigest(), 'git_blob_SHA1': hashlib.sha1(b'blob ' + str(len(b)).encode() + b'\0' + b).hexdigest()})
    write(checksums, 'input_checksums.csv')
    (OUT / 'validation.json').write_text(json.dumps({'genes': 30, 'study_gene_rows': 90, 'counts': counts, 'threshold': 'original genome-wide gene FDR < 0.05; no panel-specific recalculation', 'CCI_mapping_matches_original_DE': True, 'CCI_full_original_input_used': paths['CCI'] == ROOT / PATHS['CCI']}, indent=2) + '\n')
    print(json.dumps(counts, indent=2))
    return wide

if __name__ == '__main__':
    run()
