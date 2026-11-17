#!/usr/bin/env python
# coding: utf-8


import spatialdata as sd
from spatialdata_io import xenium
import sopa
import scanpy as sc


# preparation of file path
input_path = "path/to/the/input_xenium_bundle"

# Please make your file name like:
# "0076861__Normal_lung.explorer", "0076861__tkbadeno1.explorer"
explorer_path = "path/to/the/output_xenium_bundle"

# 
sdata = sopa.io.xenium(input_path, cells_boundaries=True, cells_table=True, nucleus_labels=True, cells_labels=True, cells_as_circles=False, nucleus_boundaries=True)
sc.pp.filter_cells(sdata.tables["table"], min_counts=10)
sc.pp.filter_cells(sdata.tables["table"], min_genes=5)    
sopa.make_transcript_patches(sdata, patch_width=None, prior_shapes_key="auto", min_points_per_patch=1)
sopa.segmentation.proseg(sdata, command_line_suffix="--nthreads 16")
sopa.aggregate(sdata)
sdata.tables["table"].layers["counts"] = sdata.tables["table"].X.copy()
sopa.io.explorer.write(explorer_path, sdata, save_h5ad=True)