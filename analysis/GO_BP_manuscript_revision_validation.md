# GO BP manuscript revision validation

## Source and scope

Reviewed branch: `three-dataset-gobp-exploration`, source commit `b08296108e276225ce168f36c4bd53ba903f86e1`. The full repository tree and relevant GO:BP scripts, Hallmark protocol, study scope, gene-set manifest, representative selections, local production tables and existing supplementary document were cross-referenced. This update prepares the GO:BP computational section for integration; it does not replace unseen experimental or full-manuscript sections.

## Scientific checks

- The uploaded production expression, diagnostic, sensitivity and concordance tables match the earlier independent calculation within 1e-10 numeric tolerance. Original gene/pathway statistics and locked selections are retained.
- Independently recomputed all 840 gene and 400 pathway sample-omission contrasts, expression log transforms and within-study z-scores. Verified summary ranges and direction-preservation counts.
- All 33 input MD5 checksums match source files with LF/Windows CRLF accounted for. The production R session records R 4.5.2 and DESeq2 1.50.2.
- Manuscript counts, significance distinctions and unstable genes/pathways are tied to the selection and stage-08 tables. The stable CCI Schwann differentiation GSVA effect is explicitly distinguished from its opposing GSEA direction.
- No new inferential P values, subset FDR adjustments, exclusions or selections were introduced.

## Script and presentation checks

- Parsed the original and revised stage-08 R syntax trees: only `gobp08_plot` changes; all analysis functions and the entry point are identical.
- All 40 related R scripts parse successfully in R 4.3.3.
- Executed the revised plotting function on uploaded production tables and inspected both figures. PNGs are complete and readable, with 300-dpi dimensions of 4500 x 4500 and 3600 x 3000; the PDF versions are single-page vector plots.
- Separate figure input checksums and rendering session information identify the presentation update. Original production checksums, session information and scientific CSVs remain unchanged.
- Verified the supplementary source document by Git blob SHA `85cdb33750707ca91896cf1e208e49f5e1238192`. Preserved all original paragraphs, images and both gene tables; appended Table S3 with eight summary rows.
- Rendered the updated supplementary document and inspected all nine pages for clipping, overlaps, table splitting and legibility. The full opposite-direction pathway list was not appended.

## Integration limits

Final main-figure numbering, journal-specific style, software-reference bibliography and integration with the full manuscript, network analysis and in vivo sections remain manuscript integration tasks. Those sections were not rewritten without their current source documents. No claim of cell-specific mechanism, OIPN specificity or independent GSEA-GSVA validation is supported by this update.
