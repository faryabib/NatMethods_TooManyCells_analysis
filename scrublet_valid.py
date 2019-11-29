import os
import os.path
import random
import sys

import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.sparse

import scrublet as scr

random.seed(0)
np.random.seed(0)

inputs = sys.argv[1:]

# Load matrix
rawMat = scipy.io.mmread(os.path.join(inputs[0], "matrix.mtx"))
labels = np.array([inputs[0].split("_")[-1]] * rawMat.shape[1])
items = pd.read_csv(os.path.join(inputs[0], "barcodes.tsv"), sep = "\t", header
                    = None).iloc[:,0].to_numpy()
for i in range(1, len(inputs)):
    rawMat_i = scipy.io.mmread(os.path.join(inputs[i], "matrix.mtx"))
    rawMat = scipy.sparse.hstack([rawMat, rawMat_i])
    labels = np.concatenate(
        [labels, np.array([inputs[i].split("_")[-1]] * rawMat_i.shape[1])])
    items = np.concatenate([items, pd.read_csv(os.path.join(inputs[i],
                                                            "barcodes.tsv"),
                                               sep = "\t", header
                                               =
                                               None).iloc[:,0].to_numpy()])
mat2 = rawMat.T
mat = mat2.toarray()

# Cluster matrix
scrub = scr.Scrublet(mat)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                        min_cells=3,
                                                        min_gene_variability_pctl=85,
                                                        n_prin_comps=30)

# Label matrix
validBarcodes = items[np.invert(predicted_doublets)]

# Output results
np.savetxt(sys.stdout.buffer, validBarcodes, delimiter=",", fmt="%s")
