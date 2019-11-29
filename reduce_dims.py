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
import sklearn.manifold

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
labelFile = sys.argv[2]
inputs = sys.argv[3:]

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")

# Load matrix
labelsDf = pd.read_csv(labelFile)
mat = scipy.io.mmread(os.path.join(inputs[0], "matrix.mtx"))
items = pd.read_csv(os.path.join(inputs[0], "barcodes.tsv"), sep = "\t", header
                    = None).iloc[:,0].to_numpy()
for i in range(1, len(inputs)):
    mat_i = scipy.io.mmread(os.path.join(inputs[i], "matrix.mtx"))
    mat = scipy.sparse.hstack([mat, mat_i])
    items = np.concatenate([items, pd.read_csv(os.path.join(inputs[i],
                                                            "barcodes.tsv"),
                                               sep = "\t", header
                                               =
                                               None).iloc[:,0].to_numpy()])
mat2 = mat.T
mat3 = mat2.toarray()

# Label matrix
itemsDf = pd.DataFrame({"item": pd.Series(items)})
itemsDf = pd.merge(itemsDf, labelsDf, on = "item", how = "left", sort = False)

# Dimensionality reduction
with suppress_stdout_stderr():
    coordinates = sklearn.manifold.TSNE(
        n_components=2, verbose=1, perplexity=30).fit_transform(mat3)
    umapCoordinates = umap.UMAP().fit_transform(mat3)

df = pd.DataFrame({
    "barcode": itemsDf["item"],
    "TSNE 1": coordinates[:, 0],
    "TSNE 2": coordinates[:, 1],
    "UMAP 1": umapCoordinates[:, 0],
    "UMAP 2": umapCoordinates[:, 1],
    "label": itemsDf["label"]
})

# Output results
# df.to_csv(sys.stdout.buffer, index=False)
df.to_csv(output + ".csv", index=False)

# Plot results
labelChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="TSNE 1",
    y="TSNE 2",
    color="label").properties(width=180, height=180)
labelChart.save(output + "_tsne_label.svg")

umapLabelChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="UMAP 1",
    y="UMAP 2",
    color="label").properties(width=180, height=180)
umapLabelChart.save(output + "_umap_label.svg")
