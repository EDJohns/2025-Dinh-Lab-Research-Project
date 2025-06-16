#!/usr/bin/env python
# coding: utf-8

# In[30]:


get_ipython().system('jupyter nbconvert --to python my_functions.ipynb')


# In[29]:


def create_z_map(adata, tma_num, punch_num, heatmap_printout, radius=-1.0):
    """
    Parameters:
    - adata: AnnData object
    - tma_num: The TMA sample (1=A, 2=B, 3=C)
    - punch_num: The number of the punch (1-25)
    - sample_string: str, name of the cell (must be in adata.obs_names)
    - heatmap_printout: if set to 1 it will 
    - Radius: OPT argument to set radius
    Returns:
    - 
    """
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import squidpy as sq

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
    if radius<1:
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic",delaunay=True)
    else:
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic",radius=radius)
    
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
    - 'PMN Count' int with the number of PMNs in that punch
    - 'PMN Names' list of all the cell names that are PMNs
    """
    pmn_count=0
    end=False
    cell_type=""
    pmn_names=[]
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
                pmn_names.append(f"{name_string}{i}")
        except:
           j+=1
        if j>100: 
            end=True
    return {
        "PMN Count": pmn_count,
        "PMN Names": pmn_names
    }
    
def cluster_response(matrix, run_col, cluster_col, response_col):
    # Group by cluster and response status
    import pandas as pd
    grouped = matrix.groupby([cluster_col, response_col]).size().unstack(fill_value=0)

    # Rename columns for clarity
    grouped.columns = ['Num_NonResponder', 'Non_Responder'] if 0 in grouped.columns else ['Non_Responder', 'Num_NonResponder']

    # Ensure both columns exist (in case one type is missing)
    if 'Num_Responder' not in grouped.columns:
        grouped['Num_Responder'] = 0
    if 'Num_NonResponder' not in grouped.columns:
        grouped['Num_NonResponder'] = 0

    # Calculate percent responders
    grouped['Percent_Responder'] = (grouped['Non_Responder'] / (grouped['Non_Responder'] + grouped['Num_NonResponder'])) * 100

    # Reset index to turn cluster into a column
    cluster_output = grouped.reset_index()

    return {
        "Output Matrix": cluster_output
    }
    
def type_neighboors(adata, tma_num, punch_num, spatial_method="n_neighboors", radius=135):
    """
    Analyze the spatial neighbors of PMN cells within a specified TMA punch.

    Parameters:
    - adata: AnnData object containing spatial transcriptomics data
    - tma_num: Integer indicating the TMA identifier
    - punch_num: Integer indicating the punch (core) within the TMA
    - spatial_method: Programs how the spatial matrix is calculated (radius, delaunay, n_neighboors) OPT
    - radius: Specifies radius for radial method OPT

    Returns:
    - 'num_neigh': pandas DataFrame indexed by PMN cell ID, with the number of neighbors for each PMN
    - 'type_neigh': pandas DataFrame with neighbor cell type breakdowns per PMN
    """ 
    import pandas as pd
    import importlib
    import squidpy as sq
    import my_functions
    pmn_results = []
    pmn_results2 = []
    pmn_count=0 
    # Makes subset of the adata object
    adata_sub = adata[adata.obs.cell_ID.str.contains(f"c_{tma_num}_{punch_num}_")].copy()
    
    # Compute spatial matrix 
    if spatial_method == "radius":
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic", radius=135)
    elif spatial_method == "delaunay":
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic", delaunay=True)
    else:
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic")
        
    
    # Display the results nicely
    dense_df = pd.DataFrame(
        adata_sub.obsp["spatial_connectivities"].toarray(),
        index=adata_sub.obs_names,
        columns=adata_sub.obs_names
    )

    dense_summed = dense_df.sum(axis=1) #Sum each row

    #Counts of number of close cell for each PMN
    pmn_cell_counts=pmn_counter(adata,f"c_{tma_num}_{punch_num}_")['PMN Names']

    #Creates matrix of PMN cells
    for z in pmn_cell_counts:
        # print(f"Cell ID: {z} Number Neighboors: {dense_summed.loc[z]} ")
        new_row=[z, tma_num, punch_num, dense_summed.loc[z]]
        pmn_results.append(new_row)

    # Creates matrix of PMN cells
    for z in pmn_cell_counts:
        # Get neighbors of PMN `z` from the spatial connectivities matrix
        neighbors = dense_df.loc[z]
        
        # Keep only neighbors with a non-zero connection
        neighbor_ids = neighbors[neighbors > 0].index.tolist()
        
        # Get their cell types from adata_sub.obs
        neighbor_types = adata_sub.obs.loc[neighbor_ids, 'cell_type']  # replace 'cell_type' with your actual column
        
        # Count frequency of each neighbor cell type
        neighbor_type_counts = neighbor_types.value_counts().to_dict()
        
        # Store the data
        new_row = [z, tma_num, punch_num, dense_summed.loc[z], neighbor_type_counts]
        pmn_results2.append(new_row)

    dense_df = pd.DataFrame(
    pmn_results,
    columns=["Cell Id", "TMA", "Sample", "Number of Cells"]
    )
    # set_index() can set columns 
    dense_df.set_index("Cell Id", inplace=True)

    pmn_results2 = pd.DataFrame(pmn_results2, columns=[
        "Cell Id", "TMA_num", "punch_num", "Number of Neighbors", "Neighbor Cell Types"
    ])
    
    return{
        "num_neigh": dense_df,
        "type_neigh": pmn_results2
    }

# def average_cell_counts(adata,tma, punch):
#     """
#     Calculate the total and percentage counts of neighbor cell types for a given TMA punch.

#     Parameters:
#     - adata: matrix of cosmx protein data
#     - tma: Identifier for the tissue microarray (TMA)
#     - punch: Identifier for the punch (core) within the TMA

#     Returns:
#     - 'Total Cell Counts': pandas DataFrame with raw counts of each neighbor cell type,
#       including a total count column and indexed by sample ID.
#     - 'Cell Type Percentages': pandas DataFrame with percentages of each cell type 
#       relative to the total, indexed by sample ID.
#     """
#     import pandas as pd
#     from collections import Counter
#     import my_functions

#     total_counts = Counter()
#     res = my_functions.type_neighboors(adata, tma, punch) 
#     cell_matrix = res["type_neigh"]

#     for index, row in cell_matrix.iterrows():
#         total_counts.update(row['Neighbor Cell Types'])

#     # Create the total counts matrix
#     total_counts = dict(total_counts)
#     total_cell_counts = pd.DataFrame([total_counts])
#     total_cell_counts['Total Cell Count'] = total_cell_counts.sum(axis=1)
#     total_cell_counts['Sample ID'] = f"c_{tma}_{punch}"
#     total_cell_counts.set_index('Sample ID', inplace=True)

#     # Create the percentages matrix (excluding the 'Total Cell Count')
#     cell_percentages = total_cell_counts.drop(columns=['Total Cell Count']).div(
#         total_cell_counts['Total Cell Count'], axis=0
#     ) * 100
#     cell_percentages['Sample ID'] = f"c_{tma}_{punch}"
#     cell_percentages.set_index('Sample ID', inplace=True)

#     return {
#         "Total Cell Counts": total_cell_counts,
#         "Cell Type Percentages": cell_percentages
#     }

def average_cell_counts(adata, tma, punch, spatial_method, radius=135):
    """
    Calculate total and percentage counts of neighbor cell types for a given TMA punch,
    using a specified spatial method for neighbor calculation.

    Parameters:
    - adata: AnnData object containing spatial transcriptomics data
    - tma: Identifier for the tissue microarray (TMA), e.g., integer or string
    - punch: Identifier for the punch (core) within the TMA, e.g., integer or string
    - spatial_method: String indicating method to compute spatial neighbors; expected values:
        "radius" or "delaunay" or others handled by default method in type_neighboors
    - radius: Numeric value specifying the radius for neighbor detection when using "radius" method (default=135)

    Returns:
    - dict with keys:
        "Total Cell Counts": pandas DataFrame, one row with raw counts of neighbor cell types,
            includes 'Total Cell Count' column and indexed by sample ID.
        "Cell Type Percentages": pandas DataFrame, one row with percentages of each cell type relative
            to the total, indexed by sample ID.
    """
    import pandas as pd
    from collections import Counter
    import my_functions

    total_counts = Counter()
    res = my_functions.type_neighboors(adata, tma, punch, spatial_method) 
    cell_matrix = res["type_neigh"]

    for index, row in cell_matrix.iterrows():
        total_counts.update(row['Neighbor Cell Types'])

    # Create the total counts matrix
    total_counts = dict(total_counts)
    total_cell_counts = pd.DataFrame([total_counts])
    total_cell_counts['Total Cell Count'] = total_cell_counts.sum(axis=1)
    total_cell_counts['Sample ID'] = f"c_{tma}_{punch}"
    total_cell_counts.set_index('Sample ID', inplace=True)

    # Create the percentages matrix (excluding the 'Total Cell Count')
    cell_percentages = total_cell_counts.drop(columns=['Total Cell Count']).div(
        total_cell_counts['Total Cell Count'], axis=0
    ) * 100
    cell_percentages['Sample ID'] = f"c_{tma}_{punch}"
    cell_percentages.set_index('Sample ID', inplace=True)

    return {
        "Total Cell Counts": total_cell_counts,
        "Cell Type Percentages": cell_percentages
    }
    
def plot_cell_type_pie(sample_id, all_percentages, width=800, height=800, donut=False):
    """
    Display an interactive pie (or donut) chart of cell type percentages for a given sample.

    Parameters:
    - sample_id (str): Sample ID (e.g., "c_1_5")
    - all_percentages (pd.DataFrame): DataFrame containing percentage data (output from your analysis)
    - width (int): Width of the chart in pixels (default: 800)
    - height (int): Height of the chart in pixels (default: 800)
    - donut (bool): Whether to render as a donut chart (default: False)

    Returns:
    - Plotly interactive pie chart
    """

    import plotly.express as px
    
    if sample_id not in all_percentages.index:
        raise ValueError(f"Sample ID '{sample_id}' not found in the dataframe.")

    # Extract row and convert to DataFrame
    percentages = all_percentages.loc[sample_id]
    df = percentages.reset_index()
    df.columns = ['Cell Type', 'Percentage']

    # Create pie chart
    fig = px.pie(
        df,
        values='Percentage',
        names='Cell Type',
        title=f"Cell Type Percentages for {sample_id}",
        hover_data=['Percentage'],
        hole=0.3 if donut else 0
    )

    # Update layout for size and readability
    fig.update_traces(textinfo='percent+label')
    fig.update_layout(
        width=width,
        height=height,
        legend_title="Cell Types",
        title_font_size=24
    )

    fig.show()

def plot_cell_type_pie(sample_id, all_percentages, width=800, height=800, donut=False):
    """
    Display an interactive pie (or donut) chart of cell type percentages for a given sample.
    
    Parameters:
    - sample_id (str): Sample ID (e.g., "c_1_5")
    - all_percentages (pd.DataFrame): DataFrame with percentage values
    - width (int): Width of the chart (default: 800)
    - height (int): Height of the chart (default: 800)
    - donut (bool): If True, renders a donut chart (default: False)

    Returns:
    - Displays an interactive Plotly pie chart.
    """
    
    import plotly.express as px
    
    # Predefined consistent colors for all possible cell types
    color_map = {
        "CD4+ T cell": "#1f77b4",
        "Endothelial": "#ff7f0e",
        "B cell": "#2ca02c",
        "Unknown": "#d62728",
        "Other immune": "#9467bd",
        "Macrophage/monocyte": "#8c564b",
        "Fibroblast": "#e377c2",
        "CD8+ T cell": "#7f7f7f",
        "Vascular smooth muscle": "#bcbd22",
        "neutrophil": "#17becf",
        "Adipocyte": "#f7b6d2",
        "NK cell": "#c49c94",
        "Epithelial": "#aec7e8",
        "Dendritic cell": "#98df8a"
    }

    if sample_id not in all_percentages.index:
        raise ValueError(f"Sample ID '{sample_id}' not found in the dataframe.")

    # Extract and format the data
    percentages = all_percentages.loc[sample_id]
    df = percentages.reset_index()
    df.columns = ['Cell Type', 'Percentage']

    # Filter out zero values (optional)
    df = df[df['Percentage'] > 0]

    # Generate the figure
    fig = px.pie(
        df,
        values='Percentage',
        names='Cell Type',
        title=f"Cell Type Percentages for {sample_id}",
        hole=0.3 if donut else 0,
        color='Cell Type',
        color_discrete_map=color_map
    )

    fig.update_traces(textinfo='percent+label')
    fig.update_layout(
        width=width,
        height=height,
        legend_title="Cell Types",
        title_font_size=24
    )

    fig.show()


# In[ ]:




