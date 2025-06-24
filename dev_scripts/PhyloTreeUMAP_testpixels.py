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

def log_scale(value, max_value):
    # Add 1 to avoid logarithm of zero
    scaled_value = np.log(value + 1) / np.log(max_value + 1)
    return scaled_value

def resize_image(original_image, resize_factor):
    """
    Resize a PIL image based on a resize factor.
    """
    # make sure that the resive factor is an even number
    if resize_factor % 2 != 0:
        raise IOError("The resize factor must be an even number.")

    rf = resize_factor
    af = int(resize_factor/2)
    # Resize the image to 1/4 of its original size
    resized_image = original_image.resize((original_image.width // rf, original_image.height // rf))
    # Average every four pixels
    resized_image = resized_image.resize((resized_image.width // af, resized_image.height // af), Image.ANTIALIAS)
    return resized_image

def plot_coo_matrix(coomatrixpath) -> Image:
    """
    Required Loadings:
      from PIL import Image
      import numpy as np

    Description:
      - Takes a coo matrix, plots it as an image.
    """
    coo = load_npz(coomatrixpath)
    # Get the first 10000 columns of the coo matrix
    coo = coo.tocsc()[:, :1000]
    coo = coo.tocoo()
    ##change all the 0 values of coo to -1
    #coo.data[coo.data == 0] = -1
    ## Convert COO matrix to dense
    #dense_matrix = coo.toarray()
    ## change all the coo_matrix 0 values to 999999999
    #dense_matrix[dense_matrix == 0] = 999999999
    ## print all the -1s to 0s
    #dense_matrix[dense_matrix == -1] = 0
    ## Find columns with all 999999999
    #nan_columns = np.all(dense_matrix == 999999999, axis=0)
    ## remove the columns with all 999999999
    #dense_matrix = dense_matrix[:, ~nan_columns]
    ## group the columns by similarity


    # SIMILARITY OF COLUMNS
    similarity_matrix = cosine_similarity(coo.T)
    print("Col similarity matrix shape:", similarity_matrix.shape)
    # Perform hierarchical clustering
    linkage_matrix = hierarchy.linkage(similarity_matrix, method='average')
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    ## Dendrogram for visualization
    #dendrogram = hierarchy.dendrogram(linkage_matrix, labels=similarity_matrix.columns, orientation='top')
    #plt.show()
    # Get the clustered column order
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    col_old_ix_to_new_ix = {sorted_ix: old_ix
                        for old_ix, sorted_ix in enumerate(clustered_order)}

    # SIMILARITY OF ROWS
    similarity_matrix = cosine_similarity(coo)
    print("Row similarity matrix shape:", similarity_matrix.shape)
    # Perform hierarchical clustering
    linkage_matrix = hierarchy.linkage(similarity_matrix, method='average')
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    ## Dendrogram for visualization
    #dendrogram = hierarchy.dendrogram(linkage_matrix, labels=similarity_matrix.columns, orientation='top')
    #plt.show()
    # Get the clustered column order
    clustered_order = hierarchy.leaves_list(linkage_matrix)
    row_old_ix_to_new_ix = {sorted_ix: old_ix
                        for old_ix, sorted_ix in enumerate(clustered_order)}
    #print("The row order is:", row_old_ix_to_new_ix)
    #counter = 0
    #for k in row_old_ix_to_new_ix:
    #    print(k, row_old_ix_to_new_ix[k])
    #    if counter == 10:
    #        break
    #    counter += 1

    # get the max value of the coo matrix that is less than 999999999
    maxval = np.max(coo.data[coo.data != 0])
    # make a PIL image in the same coordinates as the filtdf
    img_rowUn_colUn   = Image.new('RGB', (coo.shape[1], coo.shape[0]),
                            color = (255, 255, 255))
    img_rowUn_colSort = Image.new('RGB', (coo.shape[1], coo.shape[0]),
                            color = (255, 255, 255))
    img_rowSort_colUn   = Image.new('RGB', (coo.shape[1], coo.shape[0]),
                            color = (255, 255, 255))
    img_rowSort_colSort = Image.new('RGB', (coo.shape[1], coo.shape[0]),
                             color = (255, 255, 255))

    # iterate through the coo matrix
    for i, j, v in zip(coo.row, coo.col, coo.data):
        log_scaled = log_scale(v, maxval)
        cv = 255 - int(log_scaled * 255)
        #cv = 255 - int((v / maxval) * 255)
        unsort_row = i
        unsort_col = j
        new_col = col_old_ix_to_new_ix[j]
        new_row = row_old_ix_to_new_ix[i]
        # if the value is 0, black, any other value, white
        img_rowUn_colUn.putpixel(    (unsort_col, unsort_row), (255, cv, cv))
        img_rowUn_colSort.putpixel(  (   new_col, unsort_row), (255, cv, cv))
        img_rowSort_colUn.putpixel(  (unsort_col, new_row),    (255, cv, cv))
        img_rowSort_colSort.putpixel((   new_col, new_row),    (255, cv, cv))

    img_rowUn_colUn     =resize_image(img_rowUn_colUn    , 4)
    img_rowUn_colSorto  =resize_image(img_rowUn_colSorto , 4)
    img_rowSort_colUn   =resize_image(img_rowSort_colUn  , 4)
    img_rowSort_colSort =resize_image(img_rowSort_colSort, 4)

    # save the image
    img_rowUn_colUn.save(    "coo_matrix_rowUn_colUn.png")
    img_rowUn_colSort.save(  "coo_matrix_rowUn_colSort.png")
    img_rowSort_colUn.save(  "coo_matrix_rowSort_colUn.png")
    img_rowSort_colSort.save("coo_matrix_rowSort_colSort.png")


coopath = "allsamples.coo.npz"
plot_coo_matrix(coopath)