#!/usr/bin/env python
# coding: utf-8

# In[10]:


get_ipython().system('jupyter nbconvert --to python my_functions.ipynb')


# In[9]:


def create_z_map(adata, tma_num, punch_num, heatmap_printout):
    """
    Parameters:
    - adata: AnnData object
    - tma_num: The TMA sample (1=A, 2=B, 3=C)
    - punch_num: The number of the punch (1-25)
    - sample_string: str, name of the cell (must be in adata.obs_names)
    - heatmap_printout: if set to 1 it will 

    Returns:
    - 
    """
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Makes subset of the adata object
    sample_string=f"c_{tma_num}_{punch_num}_"
    # Get's Patient ID and response status from metadata. Pulls from first entry from each test
    j=1
    end_point=False
    while end_point==False:
        try:
            patient_ID=adata.obs.loc[f"{sample_string}{j}", "Patient.ID"]
            response_status=adata.obs.loc[f"{sample_string}{j}", "Response"]
            end_point=True
        except:
            j+=1
            if j>5000:
                end_point=True #prevent forever loop
    
    adata_sub = adata[adata.obs.cell_ID.str.contains(sample_string)].copy()
    
    #Compute spatial matrix 
    import squidpy as sq
    sq.gr.spatial_neighbors(adata_sub, coord_type="generic",delaunay=True)
    
    # Creates Matrix to Display
    dense_df = pd.DataFrame(
    adata_sub.obsp["spatial_connectivities"].toarray(),
    index=adata_sub.obs_names,
    columns=adata_sub.obs_names
    )

    # Ensure 'cell type' is categorical 
    adata_sub.obs['merged_annot_cluster'] = adata_sub.obs['merged_annot_cluster'].astype('category')

    # Run neighborhood enrichment on the subset
    sq.gr.nhood_enrichment(adata_sub, cluster_key="merged_annot_cluster")
    
    # Extract Z-scores and labels
    adata_sub.uns
    z = adata_sub.uns["merged_annot_cluster_nhood_enrichment"]["zscore"]
    RawNums = adata_sub.uns["merged_annot_cluster_nhood_enrichment"]["count"]
    
    # Convert to labeled DataFrame
    z_df = pd.DataFrame(z)
    RawNums_df = pd.DataFrame(RawNums)
    
    # Get the Z-score matrix
    z = adata_sub.uns["merged_annot_cluster_nhood_enrichment"]["zscore"]
    
    # Manually get cluster names from adata_sub.obs["merged_annot_cluster"]
    cluster_labels = np.unique(adata_sub.obs["merged_annot_cluster"])
    
    # Sanity check that dimensions match
    assert z.shape == (len(cluster_labels), len(cluster_labels)), "Mismatch in shape"
    
    # Make labeled matrix
    z_df = pd.DataFrame(z, index=cluster_labels, columns=cluster_labels)
    RawNums_df = pd.DataFrame(RawNums, index=cluster_labels, columns=cluster_labels)

    # PMN Count 
    pmn_count=pmn_counter(adata_sub,sample_string)["PMN Count"]

    #Creates Heatmap
    if heatmap_printout==1:     
        create_heatmap(z_df,pmn_count, tma_num, punch_num, patient_ID ,response_status)

    return {
        "Cell Spatial Matrix": dense_df,
        "Raw Interaction Matrix":RawNums_df,
        "Z-Score Interaction Matrix":z_df,
        "Patient ID": patient_ID,
        "Response Status": response_status,
        "PMN Count": pmn_count
    }

def create_heatmap(z_matrix,pmn_count, tma_num, punch_num, patient_ID ,response_status):
    """
    Parameters:
    - z_matrix: matrix of z_scores
    - pmn_count: the number of nuetrophils
    Returns:
    - 
    """
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        z_matrix,
        cmap="coolwarm",     # diverging color map: red = enriched, blue = depleted
        center=0,            # center at 0 for symmetric color scale
        annot=True,          # show numbers in each cell
        fmt=".1f",           # format the Z-scores
        linewidths=0.5,      # grid lines
        cbar_kws={"label": "Z-score"}
    )
    plt.title(f"TMA:{tma_num}-Sample:{punch_num} Patient: {patient_ID} Response Status: {response_status} - NE Z-score Heatmap - PMN Count: {pmn_count}")
    plt.xlabel("Neighbor Cell Type")
    plt.ylabel("Anchor Cell Type")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

def pmn_counter(adata, name_string):
    """
    Retrieve the number of PMN cells for a given TMA punch 
    Parameters:
    - adata: AnnData object
    - name_string: The name of the punch which without the cell number. Ex "c_1_1" 
    Returns:
    - int with the number of PMNs in that punch
    """
    pmn_count=0
    end=False
    cell_type=""
    i=0
     # This variable sets the breakpoint requires >100 failed requests in a row. Note that there are rondom cells throughout the data sets
    j=0

    while end==False:
        i+=1
        test_string=f"{name_string}{i}"
        try: 
            cell_type=adata.obs.loc[test_string, "merged_annot_cluster"]
            j=0
            if "Neutrophil" in cell_type:
                pmn_count+=1
        except:
           j+=1
        if j>100: 
            end=True
    return {
        "PMN Count": pmn_count
    }


# In[ ]:




