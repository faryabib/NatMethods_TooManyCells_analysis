import os
import os.path
import random
import sys
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import altair as alt
import numpy as np
import pandas as pd
import phenograph
import scipy
import scipy.io
import scipy.sparse
from phenograph.core import find_neighbors
from sklearn.manifold import TSNE

import altairThemes
import umap


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

# Load matrix
mat = scipy.io.mmread(os.path.join(inputs[0], "matrix.mtx"))
labels = np.array([inputs[0].split("_")[-1]] * mat.shape[1])
items = pd.read_csv(os.path.join(inputs[0], "barcodes.tsv"), sep = "\t", header
                    = None).iloc[:,0].to_numpy()
for i in range(1, len(inputs)):
    mat_i = scipy.io.mmread(os.path.join(inputs[i], "matrix.mtx"))
    mat = scipy.sparse.hstack([mat, mat_i])
    labels = np.concatenate(
        [labels, np.array([inputs[i].split("_")[-1]] * mat_i.shape[1])])
    items = np.concatenate([items, pd.read_csv(os.path.join(inputs[i],
                                                            "barcodes.tsv"),
                                               sep = "\t", header
                                               =
                                               None).iloc[:,0].to_numpy()])
mat2 = mat.T
mat3 = mat2.toarray()

# Cluster matrix
with suppress_stdout_stderr():
    communities, graph, Q = phenograph.cluster(mat3, n_jobs = 2)

# Label matrix
print(items.shape)
almost = np.column_stack([labels, communities, items])
colnames = np.array(["label", "cluster", "item"])
final = np.row_stack([colnames, almost])

# Output results
np.savetxt(sys.stdout.buffer, final, delimiter=",", fmt="%s")

# Dimensionality reduction
with suppress_stdout_stderr():
    coordinates = TSNE(
        n_components=2, verbose=1, perplexity=30).fit_transform(mat3)

df = pd.DataFrame({
    "TSNE 1": coordinates[:, 0],
    "TSNE 2": coordinates[:, 1],
    "label": pd.Series(labels),
    "cluster": pd.Series(communities)
})

# Create output directory
if not os.path.exists("./phenograph_img"):
    os.makedirs("./phenograph_img")

# Plot results
clusterChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="TSNE 1", y="TSNE 2", color="cluster:N").properties(
        width=180, height=180)
clusterChart.save("./phenograph_img/" + output + "_phenograph_cluster.svg")

domain = sorted(set(labels.tolist()))
range_ = ["#999999", "#e41a1c", "#377eb8"]
labelChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="TSNE 1",
    y="TSNE 2",
    color=alt.Color("label", scale=alt.Scale(domain=domain,
                                             range=range_))).properties(
                                                 width=180, height=180)
labelChart.save("./phenograph_img/" + output + "_phenograph_label.svg")
