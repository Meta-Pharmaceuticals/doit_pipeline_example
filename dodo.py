from tasks.Rtask import RtaskBuilder
from taskUtils import download_geo_sup_file
from pathlib import Path
import tarfile


gsenum = 'GSE197023'
gse_raw_data_dir = Path.cwd() / "data/raw"

if not gse_raw_data_dir.exists():
    download_geo_sup_file(gsenum, gse_raw_data_dir)

meta_data = RtaskBuilder().setScript('R/get_meta_data.R') \
    .setOutputs(gse_raw_data_dir / 'meta.csv') \
    .setExtraParam(params=[gsenum]).build()

raw_seurat = RtaskBuilder().setScript('R/create_seurat_object.R') \
    .setInputs(gse_raw_data_dir) \
    .setOutputs('data/raw_seurat.rds') \
    .setExtraParam(params=[gsenum]).build()

seurat_with_meta = RtaskBuilder().setScript('R/load_meta_data_and_combine_into_seurat.R') \
    .setInputs(meta_data, raw_seurat) \
    .setOutputs('data/seurat_with_meta.rds') \
    .setExtraParam(params=[gsenum]).build()

seurat_after_qc_res = RtaskBuilder().setScript('R/seurat_qc.R') \
    .setInputs(seurat_with_meta) \
    .setOutputs('data/seurat_after_qc.rds', 'data/QC_metrics.png', 'data/seurat_umap.png') \
    .setExtraParam(
        params={'minGene': 100, 'max_mt': 25}, paramOrder=['minGene', 'max_mt']
).build()

seurat_after_qc = seurat_after_qc_res.outputs[0]

annotated_seurat_res = RtaskBuilder().setScript('R/annotate_clusters.R') \
    .setInputs(seurat_after_qc) \
    .setOutputs('data/annotated_seurat.rds', 'data/plots/annotated_seurat_umap.png').build()
annotated_seurat = annotated_seurat_res.outputs[0]

recluster_res = RtaskBuilder().setScript('R/recluster_T_cells.R') \
    .setInputs(annotated_seurat) \
    .setOutputs('data/recluster_Keratinocytes_cells.rds', 'data/plots/recluster_Keratinocytes_cells.umap.png') \
    .setExtraParam(['Keratinocytes']).build()
recluster = recluster_res.outputs[0]

gene_of_interest_plots = RtaskBuilder().setScript('R/genes_of_interest.R') \
    .setInputs(recluster) \
    .setOutputs('data/plots') \
    .build()