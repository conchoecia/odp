#!/usr/bin/env python3

"""
This script's function is to test plotting the coo matrices as pixel-based images.
The goal is to keep the species in the same row, and sort the columns by similarity.
The colors will be based on the value in the matrix.
"""
from PIL import Image
import numpy as np
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz, csr_matrix
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy
import sys

def plot_coo_matrix(coomatrixpath) -> Image:
    """
    Required Loadings:
      from PIL import Image
      import numpy as np

    Description:
      - Takes a coo matrix, plots it as an image.
    """
    coo = load_npz(coomatrixpath)
    #change all the 0 values of coo to -1
    coo.data[coo.data == 0] = -1
    # Convert COO matrix to dense
    dense_matrix = coo.toarray()
    # change all the coo_matrix 0 values to 999999999
    dense_matrix[dense_matrix == 0] = 999999999
    # print all the -1s to 0s
    dense_matrix[dense_matrix == -1] = 0
    # Find columns with all 999999999
    nan_columns = np.all(dense_matrix == 999999999, axis=0)
    # remove the columns with all 999999999
    dense_matrix = dense_matrix[:, ~nan_columns]
    # group the columns by similarity

    similarity_matrix = cosine_similarity(dense_matrix)
    print("Similarity matrix shape:", similarity_matrix.shape)
    # Perform hierarchical clustering
    linkage_matrix = hierarchy.linkage(similarity_matrix, method='average')
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    ## Dendrogram for visualization
    #dendrogram = hierarchy.dendrogram(linkage_matrix, labels=similarity_matrix.columns, orientation='top')
    #plt.show()
    # Get the clustered column order
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    print("the clustered order is:", clustered_order)
    ## Rearrange columns based on clustering
    #coloc_df = coloc_df.iloc[:, clustered_order]

    #print("Sorted dense matrix based on row similarity:")
    #print(sorted_dense_matrix)
    #print("Sorted shape:", sorted_dense_matrix.shape)
    #print("Original shape:", dense_matrix.shape)

    sys.exit()

    # Just get a subdf of the columns that have the unique entries.
    # These are the ALGs.
    filtdf = ALG_pres_abs_dataframe[unique_entries]
    matrix = filtdf.values
    # make a PIL image in the same coordinates as the filtdf
    image = Image.new('RGB', (len(filtdf.columns), len(filtdf)),
                      color = (0, 0, 0))
    # iterate through the matrix i,j, and
    # if the value is 0, black, any other value, white
    for i in range(len(filtdf)):
        for j in range(len(filtdf.columns)):
            if matrix[i][j] == 0:
                #image.putpixel((j,i), (0,0,0))
                pass
            else:
                color = (255,255,255)
                #if color_dict is not None:
                #    color = hex_to_rgb(color_dict[filtdf.columns[j]])
                image.putpixel((j,i), color)

    if color_dict is not None:
        # make a ten-pixel high colorbar on the top.
        colorbar = Image.new('RGB', (len(filtdf.columns), 10), color = (0, 0, 0))
        # The colorbar encodes the ALG colocalization pairs.
        # The top 8 pixels are colored, and the bottom two pixels are black
        for i in range(len(filtdf.columns)):
            color = hex_to_rgb(color_dict[filtdf.columns[i]])
            for j in range(8):
                colorbar.putpixel((i,j), color)
        # concatenate the colorbar to the top of the image
        image = Image.fromarray(np.concatenate((np.array(colorbar), np.array(image)), axis=0))

    return image

coopath = "allsamples.coo.npz"
plot_coo_matrix(coopath)