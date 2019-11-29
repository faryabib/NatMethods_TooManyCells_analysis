import csv
import os
import os.path
import random
import sys
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import altair as alt
import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.sparse
from sklearn.manifold import TSNE

import altairThemes
import backspinpy as bs


# From https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
# Ignore std output from a function.
@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


random.seed(0)
np.random.seed(0)

output = sys.argv[1]
inputs = sys.argv[2:]

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")

# Defaults from original code (https://github.com/linnarsson-lab/BackSPIN/blob/master/backspinpy/backspin)
input_path = None
outfiles_path = None
numLevels = 4  # -d 2
feature_fit = False  # -f
feature_genes = 2000
first_run_iters = 10  # -t
first_run_step = 0.1  # -s
runs_iters = 8  # -T
runs_step = 0.3  # -S
split_limit_g = 2  # -g
split_limit_c = 2  # -c
stop_const = 1.15  # -k
low_thrs = 0.2  # -r
normal_spin = False  # -b
normal_spin_axis = "both"
verbose = False  # -v

# Load matrix
matRaw = scipy.io.mmread(os.path.join(inputs[0], "matrix.mtx"))
labels = np.array([inputs[0].split("_")[-1]] * matRaw.shape[1])
items = pd.read_csv(
    os.path.join(inputs[0], "barcodes.tsv"), sep="\t",
    header=None).iloc[:, 0].to_numpy()
for i in range(1, len(inputs)):
    mat_i = scipy.io.mmread(os.path.join(inputs[i], "matrix.mtx"))
    matRaw = scipy.sparse.hstack([matRaw, mat_i])
    labels = np.concatenate(
        [labels, np.array([inputs[i].split("_")[-1]] * mat_i.shape[1])])
    items = np.concatenate([
        items,
        pd.read_csv(
            os.path.join(inputs[i], "barcodes.tsv"), sep="\t",
            header=None).iloc[:, 0].to_numpy()
    ])
mat = matRaw.toarray()
matSc = mat

# Feature selection
ixFeatures = bs.feature_selection(matSc, feature_genes)
matSc = matSc[ixFeatures, :]

# Normalization
matSc = np.log2(matSc + 1)
matSc = matSc - matSc.mean(1)[:, np.newaxis]

# Cluster matrix
with suppress_stdout_stderr():
    res = bs.backSPIN(matSc, numLevels, first_run_iters, first_run_step,
                      runs_iters, runs_step, split_limit_g, split_limit_c,
                      stop_const, low_thrs, verbose)
    clusters = res.cells_gr_level[:, -1].astype(int)

# Label matrix
almost = np.column_stack([labels, clusters, items])
colnames = np.array(["label", "cluster", "item"])
final = np.row_stack([colnames, almost])

# Output results
np.savetxt(sys.stdout.buffer, final, delimiter=",", fmt="%s")

# Dimensionality reduction
with suppress_stdout_stderr():
    coordinates = TSNE(
        n_components=2, verbose=0, perplexity=30).fit_transform(mat.T)

df = pd.DataFrame({
    "TSNE 1": coordinates[:, 0],
    "TSNE 2": coordinates[:, 1],
    "label": pd.Series(labels),
    "cluster": pd.Series(clusters)
})

# Create output directory
if not os.path.exists("./backspin_img"):
    os.makedirs("./backspin_img")

# Plot results
clusterChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="TSNE 1",
    y="TSNE 2",
    color=alt.Color("cluster:N",
                    scale=alt.Scale(scheme="tableau20"))).properties(
                        width=180, height=180)
clusterChart.save("./backspin_img/" + output + "_backspin_cluster.svg")

domain = sorted(set(labels.tolist()))
range_ = ["#999999", "#e41a1c", "#377eb8"]
labelChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="TSNE 1",
    y="TSNE 2",
    color=alt.Color("label", scale=alt.Scale(domain=domain,
                                             range=range_))).properties(
                                                 width=180, height=180)
labelChart.save("./backspin_img/" + output + "_backspin_label.svg")
