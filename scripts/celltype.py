import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
import celltypist
print(f'celltypist version: {celltypist.__version__}')
from celltypist import models
import sys

file_path = sys.argv[1]
usemodel = sys.argv[2]
adata = sc.read(file_path)
adata = adata.T
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,target_sum=1e4)
models.models_path
model = models.Model.load(model = usemodel)
predictions = celltypist.annotate(adata, model = usemodel, majority_voting = True, mode = 'prob match', p_thres = 0.5)
adata_new = predictions.to_adata()
majority_voting = adata_new.obs['majority_voting']
output_file = f"out_{sys.argv[1]}.csv"
majority_voting.to_csv(output_file, header=True)
