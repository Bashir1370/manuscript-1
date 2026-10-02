"""Audit archived STRING associations and redraw the author's network layout.

Run from repository root: python scripts/gobp/09_shared_STRING_network.py
Requires matplotlib for figure export; all tabular calculations use stdlib.
No web requests, inference of missing interactions, or expression refitting.
"""
import csv
import hashlib
import json
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/String'


def read(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write(name, rows):
    with (OUT / name).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def components(adj):
    unseen = set(adj)
    result = []
    while unseen:
        pending = [min(unseen)]
        found = set()
        while pending:
            node = pending.pop()
            if node in found:
                continue
            found.add(node)
            pending.extend(adj[node] - found)
        unseen -= found
        result.append(sorted(found))
    return sorted(result, key=lambda group: (-len(group), group))


def main():
    edges_path = OUT / 'string_interactions_short_37_nodes.tsv'
    nodes_path = OUT / 'string_interactions_short_37_nodes.tsv default node.csv'
    shared_path = ROOT / 'results/GO_BP_three_dataset/representative_11/displayed_gene_evidence.csv'
    opposite_path = ROOT / 'results/GO_BP_three_dataset/divergent_representative_9/priority_final_genes.csv'
    sensitivity_path = ROOT / 'results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/gene_sensitivity_summary.csv'
    layout_path = OUT / 'network_layout_37.csv'
    shared = read(shared_path)
    genes = {row['symbol'] for row in shared}
    assert len(shared) == len(genes) == 37
    edges = read(edges_path, '\t')
    pairs = {tuple(sorted((row['#node1'], row['node2']))) for row in edges}
    assert len(edges) == len(pairs) == 55
    for row in edges:
        assert row['#node1'] in genes and row['node2'] in genes
        assert row['#node1'] != row['node2']
        assert all(row[key].startswith('10116.') for key in ('node1_string_id', 'node2_string_id'))
        assert 0.4 <= float(row['combined_score']) <= 1.0
    node_output, threshold_output = [], []
    adjacency = {}
    for threshold in (0.4, 0.7, 0.9):
        adj = {gene: set() for gene in genes}
        kept = [row for row in edges if float(row['combined_score']) >= threshold]
        for row in kept:
            adj[row['#node1']].add(row['node2'])
            adj[row['node2']].add(row['#node1'])
        groups = components(adj)
        threshold_output.append(dict(minimum_score=threshold, n_query_genes=37,
            n_edges=len(kept), n_connected_genes=sum(bool(value) for value in adj.values()),
            n_isolates=sum(not value for value in adj.values()),
            component_sizes=';'.join(str(len(group)) for group in groups)))
        for component_id, group in enumerate(groups, 1):
            for gene in group:
                node_output.append(dict(symbol=gene, minimum_score=threshold,
                    degree=len(adj[gene]), component_id=component_id,
                    component_size=len(group), isolated=not bool(adj[gene])))
        adjacency[threshold] = adj
    archived_nodes = read(nodes_path)
    assert len(archived_nodes) == 30
    assert {row['name'] for row in archived_nodes} == {g for g in genes if adjacency[0.4][g]}
    assert all(int(row['Degree']) == len(adjacency[0.4][row['name']]) for row in archived_nodes)
    write('network_node_audit_37.csv', node_output)
    write('network_threshold_sensitivity.csv', threshold_output)
    evidence = {row['symbol']: row for row in shared + read(opposite_path)}
    sensitivity = {(row['id'], row['study']): row for row in read(sensitivity_path)}
    panel = [('Cdkn1a', 'shared_positive_injury_response'),
             ('Cdk1', 'shared_positive_injury_response'),
             ('Csf1', 'shared_positive_immune_response'),
             ('Hmgcs1', 'shared_negative_sterol_biosynthesis'),
             ('Cav1', 'opposite_WNT_annotation'),
             ('Col4a2', 'opposite_collagen_annotation')]
    studies = ('OIPN_GSE160543', 'NC_GSE246156', 'CCI_GSE212311')
    panel_rows = []
    for gene, role in panel:
        row = evidence[gene]
        result = dict(symbol=gene, panel_role=role,
            OIPN_RNAseq_fold_change=2 ** float(row['OIPN_GSE160543_log2FC']))
        for study in studies:
            result[study + '_log2FC'] = row[study + '_log2FC']
            result[study + '_gene_FDR'] = row[study + '_gene_FDR']
            result[study + '_descriptive_direction_stable'] = sensitivity[(gene, study)]['all_omissions_preserve_direction']
        result['shared_network_degree_0_4'] = len(adjacency[0.4][gene]) if gene in genes else 'not_in_shared_network'
        panel_rows.append(result)
    write('qpcr_selected_panel_6.csv', panel_rows)
    positions = {row['symbol']: (float(row['x']), float(row['y'])) for row in read(layout_path)}
    assert len(positions) == 37 and set(positions) == genes
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from matplotlib.lines import Line2D
    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                                'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(10, 8))
    for row in edges:
        x1, y1 = positions[row['#node1']]
        x2, y2 = positions[row['node2']]
        score = float(row['combined_score'])
        ax.plot([x1, x2], [y1, y2], color='#858585',
                linewidth=0.6 + 1.8 * (score - 0.4) / 0.6, zorder=1)
    sterol = {'Dhcr7', 'Hmgcs1', 'Insig1', 'Sqle'}
    colors = {'main': '#A1D99B', 'sterol': '#47C7E8', 'isolated': '#FDAE6B'}
    for gene, (x, y) in positions.items():
        category = 'isolated' if not adjacency[0.4][gene] else 'sterol' if gene in sterol else 'main'
        ax.add_patch(Ellipse((x, y), 110, 50, facecolor=colors[category], edgecolor='none', zorder=2))
        ax.text(x, y, gene, fontsize=10, ha='center', va='center', zorder=3)
    ax.set(xlim=(30, 1190), ylim=(855, 15))
    ax.set_aspect('equal')
    ax.axis('off')
    handles = [Line2D([], [], marker='o', linestyle='none', markerfacecolor=colors[k],
               markeredgecolor='none', markersize=9, label=label) for k, label in
               [('main', 'Main connected component'), ('sterol', 'Sterol-associated component'),
                ('isolated', 'No association at score >= 0.4')]]
    legend = fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.04, 0.015),
                        fontsize=9, frameon=False)
    score_handles = [Line2D([], [], color='#858585', linewidth=0.6 + 1.8 * (score - 0.4) / 0.6,
                    label=str(score)) for score in (0.4, 0.7, 1.0)]
    fig.legend(handles=score_handles, title='STRING confidence score', loc='lower right',
               bbox_to_anchor=(0.97, 0.015), fontsize=9, title_fontsize=9, frameon=False)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.15)
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(OUT / ('Fig_shared_37_STRING_network.' + extension), dpi=600, facecolor='white')
    plt.close(fig)
    provenance = dict(source_commit='fa278673ab02ec4eb023326352d802e80ae14783',
        species='Rattus norvegicus', taxon_id=10116, primary_score=0.4,
        additional_interactors=0, network_type='full functional association network',
        active_sources='all channels plus evidence transfer, as documented by author settings screenshot',
        STRING_version='not recorded in supplied export', Cytoscape_version='not recorded in supplied export',
        ppi_enrichment_p_value='not supplied; not inferred from edge list',
        layout='Manually transcribed from author Cytoscape PNG; no algorithmic clustering',
        figure_line_width_points='linear: score 0.4 -> 0.6 pt; score 1.0 -> 2.4 pt',
        python_version=platform.python_version(), matplotlib_version=matplotlib.__version__,
        input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in (edges_path, nodes_path, shared_path, opposite_path, sensitivity_path, layout_path)})
    (OUT / 'network_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'thresholds': threshold_output, 'panel': [gene for gene, role in panel],
                      'archived_degrees_match': True}, indent=2))


if __name__ == '__main__':
    main()
