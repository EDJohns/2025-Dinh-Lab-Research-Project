#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().system('jupyter nbconvert --to python my_functions.ipynb')


# In[2]:


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
    # Ensure name_string ends with "_"
    if not name_string.endswith("_"):
        name_string += "_"
    pmn_count=0
    end=False
    cell_type=""
    pmn_names=[]
    overall_cell_count=0
    i=0
     # This variable sets the breakpoint requires >100 failed requests in a row. Note that there are rondom cells throughout the data sets
    j=0

    while end==False:
        i+=1
        test_string=f"{name_string}{i}"
        try: 
            cell_type=adata.obs.loc[test_string, "merged_annot_cluster"]
            j=0
            overall_cell_count+=1
            if "Neutrophil" in cell_type:
                pmn_count+=1
                pmn_names.append(f"{name_string}{i}")
        except:
           j+=1
        if j>100: 
            end=True

    # Calculate the overall percent of PMN cells
    if overall_cell_count>0:
        decimal_pmn=pmn_count/overall_cell_count
    else:
        decimal_pmn=0

    
    return {
        "PMN Count": pmn_count,
        "PMN Names": pmn_names,
        "PMN decimal": decimal_pmn
    }
    
def cluster_response(matrix, run_col, cluster_col, response_col):
    """
    Analyze response distribution and tumor cell information by cluster, including average immune cell values.
    Parameters:
    ----------
    matrix : pandas.DataFrame
        The input data containing at least the columns:
        - cluster_col
        - response_col (assumed binary: 0 = Non-Responder, 1 = Responder)
        - Tumor_cells
        - CD4+T_cells
        - CD8+T_cells
        - Tregs
    run_col : str
        Not used in this function but kept for compatibility with existing interfaces.
    cluster_col : str
        Name of the column indicating cluster assignments.
    response_col : str
        Name of the column indicating response status (assumed binary: 0 = Non-Responder, 1 = Responder).
    Returns:
    -------
    dict
            "Output Matrix": pandas.DataFrame
                Contains cluster-level counts of Responders, Non-Responders, and Percent_Responder.
            "cluster_output_with_tumor_info": pandas.DataFrame
                Same as Output Matrix with additional columns:
                - 'Avg Tumor Cell Delta'
                - 'Avg CD4+T_cells'
                - 'Avg CD8+T_cells'
                - 'Avg Tregs'
    """
    import pandas as pd

    # Group by cluster and response status
    grouped = matrix.groupby([cluster_col, response_col]).size().unstack(fill_value=0)

    # Rename columns for clarity
    grouped.columns = ['Non_Responder', 'Responder']

    # Ensure both columns exist
    if 'Responder' not in grouped.columns:
        grouped['Responder'] = 0
    if 'Non_Responder' not in grouped.columns:
        grouped['Non_Responder'] = 0

    # Calculate percent responders
    grouped['Percent_Responder'] = (grouped['Responder'] / (grouped['Non_Responder'] + grouped['Responder'])) * 100

    # Reset index to turn cluster into a column
    cluster_output = grouped.reset_index()

    # Calculate average Tumor_cells and immune cell values per cluster
    avg_metrics = matrix.groupby(cluster_col)[['Tumor_cells', 'CD4+T_cells', 'CD8+T_cells', 'Tregs']].mean().reset_index()
    avg_metrics.rename(columns={
        'Tumor_cells': 'Avg Tumor Cell Delta',
        'CD4+T_cells': 'Avg CD4+T_cells',
        'CD8+T_cells': 'Avg CD8+T_cells',
        'Tregs': 'Avg Tregs'
    }, inplace=True)

    # Merge averages with cluster output
    cluster_output_with_tumor_info = pd.merge(cluster_output, avg_metrics, on=cluster_col, how='left')

    return {
        "Output Matrix": cluster_output,
        "cluster_output_with_tumor_info": cluster_output_with_tumor_info
    }
    
def type_neighboors(adata, tma_num, punch_num, spatial_method="n_neighboors", radius=15):
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
    pmn_results = []
    pmn_results2 = []
    pmn_count=0 
    # Makes subset of the adata object
    adata_sub = adata[adata.obs.cell_ID.str.contains(f"c_{tma_num}_{punch_num}_")].copy()

    # Changes Radius to use pixels instead of microns
    radius=radius/0.168
    
    # Compute spatial matrix 
    if spatial_method == "radius" or spatial_method == "radial":
        sq.gr.spatial_neighbors(adata_sub, coord_type="generic", radius=radius)
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
        neighbor_types = adata_sub.obs.loc[neighbor_ids, 'merged_annot_cluster']  # replace 'cell_type' with your actual column
        
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
def average_cell_counts(adata, tma, punch, spatial_method, radius=15):
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

    total_counts = Counter()
    res = type_neighboors(adata, tma, punch, spatial_method, radius) 
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

def plot_cell_type_pie(sample_id, all_percentages, title, width=800, height=800, donut=False, legend=True):
    """
    Display an interactive pie (or donut) chart of cell type percentages for a given sample.
    
    Parameters:
    - sample_id (str): Sample ID (e.g., "c_1_5")
    - all_percentages (pd.DataFrame): DataFrame with percentage values
    - width (int): Width of the chart (default: 800)
    - height (int): Height of the chart (default: 800)
    - donut (bool): If True, renders a donut chart (default: False)
    - legend (bool): True by deafault, False removes legend

    Returns:
    - Displays an interactive Plotly pie chart.
    """
    
    import plotly.express as px
    
    # Predefined consistent colors for all possible cell types
    color_map = {
        "Plasma_cells": "#1f77b4",           # blue
        "Fibroblasts/SMCs": "#ff7f0e",       # orange
        "Endothelial_cells": "#2ca02c",      # green
        "Monocytes": "#d62728",              # red
        "Tumor_cells": "#9467bd",            # purple
        "Neutrophils": "#8c564b",            # brown
        "CD8+T_cells": "#e377c2",            # pink
        "CD4+T_cells": "#7f7f7f",            # gray
        "Tregs": "#bcbd22",                  # yellow-green
        "Macrophages": "#17becf",            # cyan
        "DCs": "#f7b6d2",                    # light pink
        "B_cells": "#c49c94",                # beige
        "NK_cells": "#aec7e8"                # light blue
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
        title=title,
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
    if legend==False:
        fig.update_layout(
        width=width,
        height=height,
        showlegend=False,  # disables legend display
        title_font_size=24
        )
        


    fig.show()
    
def get_sample_info(adata,sample_name):
    """
    Retrieve the response status and patient ID for a given sample prefix by
    scanning through possible cell ID entries in the AnnData object (`adata`).

    Parameters:
    - adata (AnnDataOject): The data object 
    - sample_name (str): The prefix of the sample’s cell IDs (e.g., "c_1_5_").
                         The function appends increasing integers (starting at 1)
                         to this prefix while searching the AnnData index.

    Returns:
    - dict:
        "Response"    : The clinical response status for the first matching cell ID.
                        Returns "NA" if no valid ID is found after 5,000 attempts.
        "Patient ID"  : The patient identifier for the matching cell.
                        Returns "NA" if no valid ID is found.
    """
    # Ensure sample_id ends with "_"
    if not sample_name.endswith("_"):
        sample_name += "_"
    
    j=1
    end_point=False
    patient_ID="NA" # Deafault Value
    response_status="NA"
    while end_point==False:
        try:
            patient_ID=adata.obs.loc[f"{sample_name}{j}", "Patient.ID"]
            response_status=adata.obs.loc[f"{sample_name}{j}", "Response"]
            end_point=True
        except:
            j+=1
            if j>5000:
                end_point=True #prevent forever loop
    return{
        "Response":response_status,
        "Patient ID":patient_ID
    }
    
def get_sample_cell_percentages(adata, sample_id):
    """
    Calculate total cell counts and percentage breakdown of cell types for a given sample.

    Parameters:
    ----------
    adata : AnnData
        Annotated data matrix containing single-cell observations.
    sample_id : str
        Sample identifier prefix (e.g., "c_1_1" or "c_1_1_"). Underscore will be added if missing.

    Returns:
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - First DataFrame: one row with total counts of each cell type.
        - Second DataFrame: one row with the percentage of each cell type (summing to ~100).
    """
     # Ensure sample_id ends with "_"
    if not sample_id.endswith("_"):
        sample_id += "_"
    
    # Filter the AnnData object to only the cells from this sample
    sample_cells = adata[adata.obs.index.str.startswith(sample_id)]
    
    # Count occurrences of each cell type
    counts = sample_cells.obs["merged_annot_cluster"].value_counts()
    
    # Create a row-style DataFrame (1 row, multiple columns)
    count_row = counts.T.to_frame().T
    count_row.index = [sample_id]
    
    # Compute percentages
    total = counts.sum()
    percent_row = (counts / total * 100).T.to_frame().T
    percent_row.index = [sample_id]
    
    return {
        "Total Counts": count_row,
        "Percentage Counts": percent_row
    }
    
def get_xy_coords(adata,cell_ID):
    """
    Retrieve the XY pixel coordinates for a given cell.
    Parameters:
    ----------
    adata : AnnData
    cell_ID : str
        Unique identifier for the cell.
    Returns:
    -------
    dict
        Dictionary containing:
        - "X": X-coordinate in pixels (int)
        - "Y": Y-coordinate in pixels (int)
    """
    x_coord=adata.obs.loc[cell_ID, "x_FOV_px"]
    y_coord=adata.obs.loc[cell_ID, "y_FOV_px"]
    return{
        "X":x_coord,
        "Y":y_coord
    }
    
def find_distance(adata,cell_ID_1,cell_ID_2,unit_conversion=0.168):
    """
    Calculate the Euclidean distance between two cells based on their XY coordinates.
    Parameters:
    ----------
    cell_ID_1 : str
        Unique identifier for the first cell.
    cell_ID_2 : str
        Unique identifier for the second cell.
    unit_conversion : float, optional
        Conversion factor from pixels to microns. Default is 0.168 µm/px.

    Returns:
    -------
    dict
        Dictionary containing:
        - "Distance px": Distance between cells in pixels (float)
        - "Distance um": Distance between cells in microns (float)
    """
    from math import sqrt
    # Get Cell Locations
    cell_1_xy=get_xy_coords(adata,cell_ID_1)
    cell_2_xy=get_xy_coords(adata,cell_ID_2)

    # Distance Function
    Distance_px = sqrt((cell_2_xy["X"] - cell_1_xy["X"])**2 + (cell_2_xy["Y"] - cell_1_xy["Y"])**2)

    # Convert pixel distance to micron
    Distance_um = Distance_px*unit_conversion

    return{
            "Distance px":Distance_px,
            "Distance um": Distance_um
        }


def nearest_cells_of_particular_type(adata,sample_ID,cell_distance_type):
    """
    Calculates distance matrix between PMNs and target cells of a particular type,
    along with summary statistics.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing spatial information for cells.
    sample_ID : str
        Identifier for the specific sample to analyze.
    cell_distance_type : str
        Cell type to use as the target for distance calculations (e.g., "CD4+T_cells").

    Returns:
    --------
    dict
        A dictionary containing:
        - 'Distance Matrix' (pd.DataFrame): Matrix of distances between PMNs and target cells,
          including additional columns for 'Minimum Distance' and 'Average Distance' per PMN.
        - 'Average Distance' (float): Overall average distance across the entire matrix.
        - 'Average Mininum Distance' (float): Average of each PMN's minimum distance to target cells.
    """
    import pandas as pd

    # Creates lists of cell types for PMN and the target comparison
    pmn_list=pmn_counter(adata,sample_ID)['PMN Names']
    target_cell_list=matching_cell_list(adata,sample_ID,cell_distance_type)['Cell Names']
    
    # Initialize results storage: one column for this cell type
    distance_matrix = pd.DataFrame(index=pmn_list, columns=target_cell_list)
    
    for pmn in pmn_list:
        for target_cell in target_cell_list:
            distance_matrix.loc[pmn,target_cell]=find_distance(adata,pmn, target_cell)["Distance um"]
    
    # Ensure numeric types
    distance_matrix = distance_matrix.apply(pd.to_numeric, errors='coerce')
    
    # Calculate overall average distance across entire matrix
    overall_avg_distance = distance_matrix.mean().mean()

    # Calculate average mininum distance
    average_min_distance = distance_matrix.min(axis=1).mean()
            
    # Add columns first (unrounded)
    distance_matrix['Minimum Distance'] = distance_matrix.min(axis=1).round(2)
    distance_matrix['Average Distance'] = distance_matrix.mean(axis=1).round(2)
    
    return{
        "Distance Matrix":distance_matrix,
        "Average Distance":overall_avg_distance,
        "Average Mininum Distance": average_min_distance
    }

def get_interaction_status(adata,sample_id, distance_threshold, cell_type_1="CD8+T_cells",cell_type_2="Tumor_cells"):
    """
    Classifies PMN cell interactions with two specified cell types (e.g., CD8+ T cells and Tumor cells)
    based on spatial proximity from a distance matrix.
    Parameters
    ----------
    adata : AnnData
        The annotated data object containing spatial transcriptomics or single-cell data.
    sample_id : str
        Identifier for the sample to analyze within `adata`.
    distance_threshold : float
        Distance in microns to define interaction (e.g., 20 means a cell is interacting if it's ≤ 20 microns away).
    cell_type_1 : str, optional
        First target cell type to assess interaction with (default is "CD8+T_cells").
    cell_type_2 : str, optional
        Second target cell type to assess interaction with (default is "Tumor_cells").
    Returns
    -------
    dict
        A dictionary with one key:
        - "Interaction Matrix": pandas.DataFrame with:
            * Boolean columns: 
              - '{cell_type_1} Interacting'
              - '{cell_type_2} Interacting'
            * Categorical column:
              - 'Interaction_Status' (values: 'Both', '{cell_type_1} only', '{cell_type_2} only', 'Neither')
    """
    import pandas as pd
    
    # Load distance matrices
    results_CD8 = nearest_cells_of_particular_type(adata,sample_id, cell_type_1)["Distance Matrix"]
    results_Tumor = nearest_cells_of_particular_type(adata,sample_id, cell_type_2)["Distance Matrix"]
    
    # Set interaction threshold (in microns)
    threshold = distance_threshold
    
    # Initialize result matrix
    output = pd.DataFrame(index=results_CD8.index)
    output[f'{cell_type_1} Interacting'] = results_CD8['Minimum Distance'] <= threshold
    output[f'{cell_type_2} Interacting'] = results_Tumor['Minimum Distance'] <= threshold
    
    # Compute interaction status
    def classify(row):
        if row[f'{cell_type_1} Interacting'] and row[f'{cell_type_2} Interacting']:
            return 'Both'
        elif row[f'{cell_type_1} Interacting']:
            return f'{cell_type_1} only'
        elif row[f'{cell_type_2} Interacting']:
            return f'{cell_type_2} only'
        else:
            return 'Neither'
    
    output['Interaction_Status'] = output.apply(classify, axis=1)

    return{
        "Interaction Matrix":output
    }

def plot_interaction_status_pie(interaction_matrix, title = "PMN Interaction Status Distribution"):
    """
    Plots a pie chart of PMN cell interaction status counts.
    Parameters
    ----------
    interaction_matrix : pandas.DataFrame
        The output DataFrame from `get_interaction_status()["Interaction Matrix"]` that contains 
        the 'Interaction_Status' column.
    title : str, optional
        Title for the pie chart (default is "PMN Interaction Status Distribution").
    Returns
    -------
    None
        Displays the pie chart using matplotlib.
    """
    import matplotlib.pyplot as plt
    
    # Count the frequency of each interaction status
    status_counts = interaction_matrix['Interaction_Status'].value_counts()

    # Define colors (optional customization)
    colors = {
        'Both': '#1f77b4',
        'CD8+T_cells only': '#2ca02c',
        'Tumor_cells only': '#d62728',
        'Neither': '#7f7f7f'
    }

    # Match colors to present statuses
    status_colors = [colors.get(status, '#cccccc') for status in status_counts.index]

    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(
        status_counts,
        labels=status_counts.index,
        autopct='%1.1f%%',
        startangle=140,
        colors=status_colors
    )
    plt.title(title)
    plt.axis('equal')  # Equal aspect ratio for a perfect circle
    plt.show()

def get_interaction_status_v2(adata, sample_id, distance_threshold, target_cell_type, cell_type_1="Neutrophil", cell_type_2="NA"):
    """
    Determine whether a particular cell group is interacting with one or two other cell groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial transcriptomics or single-cell data.
    sample_id : str
        Identifier for the sample to analyze within `adata`.
    distance_threshold : float
        Distance in microns to define interaction.
    target_cell_type : str
        Cell type to check interactions for.
    cell_type_1 : str, optional
        First target cell type to assess interaction with (default: "Neutrophil").
    cell_type_2 : str, optional
        Second target cell type to assess interaction with (default: "NA").

    Returns
    -------
    dict
        Contains:
        - "Interaction Matrix": pandas.DataFrame with interaction booleans and status labels.
    """
    import pandas as pd

    # Get distance matrix for cell_type_1
    results_1 = nearest_cells_of_particular_type_v2(adata, sample_id, target_cell_type, cell_type_1)["Distance Matrix"]
    threshold = distance_threshold

    output = pd.DataFrame(index=results_1.index)
    output[f'{cell_type_1} Interacting'] = results_1['Minimum Distance'] <= threshold

    if cell_type_2 != "NA":
        # Get distance matrix for cell_type_2
        results_2 = nearest_cells_of_particular_type_v2(adata, sample_id, target_cell_type, cell_type_2)["Distance Matrix"]
        output[f'{cell_type_2} Interacting'] = results_2['Minimum Distance'] <= threshold

        # Classify into four categories
        def classify(row):
            if row[f'{cell_type_1} Interacting'] and row[f'{cell_type_2} Interacting']:
                return 'Both'
            elif row[f'{cell_type_1} Interacting']:
                return f'{cell_type_1} only'
            elif row[f'{cell_type_2} Interacting']:
                return f'{cell_type_2} only'
            else:
                return 'Neither'
    else:
        # Classify into two categories
        def classify(row):
            if row[f'{cell_type_1} Interacting']:
                return 'Interacting'
            else:
                return 'Not Interacting'

    output['Interaction_Status'] = output.apply(classify, axis=1)

    return {
        "Interaction Matrix": output
    }
def nearest_cells_of_particular_type_v2(adata,sample_ID,orgin_cell_type,cell_distance_type):
    """
    Calculates distance matrix between PMNs and target cells of a particular type,
    along with summary statistics.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing spatial information for cells.
    sample_ID : str
        Identifier for the specific sample to analyze.
    orgin_cell_type : str
        The string name of the orgin cell specified. This will mostly commonly be neutrophils
    cell_distance_type : str
        Cell type to use as the target for distance calculations (e.g., "CD4+T_cells").

    Returns:
    --------
    dict
        A dictionary containing:
        - 'Distance Matrix' (pd.DataFrame): Matrix of distances between orgin and target cells,
          including additional columns for 'Minimum Distance' and 'Average Distance' per orging cell.
        - 'Average Distance' (float): Overall average distance across the entire matrix.
        - 'Average Mininum Distance' (float): Average of each orgin cell's minimum distance to target cells.
    """
    import pandas as pd

    # Creates lists of cell types for PMN and the target comparison
    orgin_cell_type_list=matching_cell_list(adata,sample_ID,orgin_cell_type)['Cell Names']
    target_cell_list=matching_cell_list(adata,sample_ID,cell_distance_type)['Cell Names']
    
    # Initialize results storage: one column for this cell type
    distance_matrix = pd.DataFrame(index=orgin_cell_type_list, columns=target_cell_list)
    
    for pmn in orgin_cell_type_list:
        for target_cell in target_cell_list:
            distance_matrix.loc[pmn,target_cell]=find_distance(adata,pmn, target_cell)["Distance um"]
    
    # Ensure numeric types
    distance_matrix = distance_matrix.apply(pd.to_numeric, errors='coerce')
    
    # Calculate overall average distance across entire matrix
    overall_avg_distance = distance_matrix.mean().mean()

    # Calculate average mininum distance
    average_min_distance = distance_matrix.min(axis=1).mean()
            
    # Add columns first (unrounded)
    distance_matrix['Minimum Distance'] = distance_matrix.min(axis=1).round(2)
    distance_matrix['Average Distance'] = distance_matrix.mean(axis=1).round(2)
    
    return{
        "Distance Matrix":distance_matrix,
        "Average Distance":overall_avg_distance,
        "Average Mininum Distance": average_min_distance
    }

def plot_interaction_status_pie_v2(interaction_matrix, title = "PMN Interaction Status Distribution"):
    """
    Plots a pie chart of PMN cell interaction status counts.
    Parameters
    ----------
    interaction_matrix : pandas.DataFrame
        The output DataFrame from `get_interaction_status()["Interaction Matrix"]` that contains 
        the 'Interaction_Status' column.
    title : str, optional
        Title for the pie chart (default is "PMN Interaction Status Distribution").
    Returns
    -------
    None
        Displays the pie chart using matplotlib.
    """
    import matplotlib.pyplot as plt
    
    # Count the frequency of each interaction status
    status_counts = interaction_matrix['Interaction_Status'].value_counts()

    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(
        status_counts,
        labels=status_counts.index,
        autopct='%1.1f%%',
        startangle=140
    )
    plt.title(title)
    plt.axis('equal')  # Equal aspect ratio for a perfect circle
    plt.show()

def plot_grouped_protein_heatmap(adata, interaction_matrix, protein_layer=None,
                                  title="Protein Expression by Interaction Group",
                                  group_order=None, vmax=8):
    """
    Plots a heatmap of average protein expression across interaction groups,
    including cell counts in row labels.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing protein expression data (in `.X` or a `.layers[]` slot).

    interaction_matrix : pandas.DataFrame
        Output of get_interaction_status()["Interaction Matrix"], must have 'Interaction_Status'.

    protein_layer : str or None, optional
        If provided, uses adata.layers[protein_layer] for protein expression.
        If None, uses adata.X.

    title : str
        Title of the heatmap plot.

    group_order : list of str, optional
        Custom ordering of interaction groups. If None, groups will be sorted alphabetically.

    vmax : float, optional
        Max value for color scale (default: 8).

    Returns
    -------
    None
        Displays the heatmap.
    """

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Align AnnData with interaction matrix
    adata_subset = adata[interaction_matrix.index].copy()
    adata_subset.obs['Interaction_Status'] = interaction_matrix['Interaction_Status']

    # Get expression matrix
    expr = adata_subset.layers[protein_layer] if protein_layer else adata_subset.X
    expr_df = pd.DataFrame(
        expr.toarray() if hasattr(expr, 'toarray') else expr,
        index=adata_subset.obs_names,
        columns=adata_subset.var_names
    )
    expr_df['Interaction_Status'] = adata_subset.obs['Interaction_Status']

    # Count cells per group
    group_counts = expr_df['Interaction_Status'].value_counts()

    # Group and average expression
    grouped_expr = expr_df.groupby('Interaction_Status').mean()

    # Rename rows to include counts
    new_index = [
        f"{group} (n={group_counts[group]})"
        for group in grouped_expr.index
    ]
    grouped_expr.index = new_index

    # Apply custom order (with counts)
    if group_order:
        # Map raw group names to their labeled form
        order_with_counts = [
            f"{g} (n={group_counts[g]})" for g in group_order if g in group_counts
        ]
        grouped_expr = grouped_expr.reindex(order_with_counts)
    else:
        grouped_expr = grouped_expr.sort_index()

    # Plot heatmap
    plt.figure(figsize=(max(10, len(grouped_expr.columns) * 0.5), len(grouped_expr) * 0.6 + 2))
    sns.heatmap(grouped_expr, cmap="viridis", annot=True, vmin=0, vmax=vmax)
    plt.title(title)
    plt.ylabel('Interaction Group')
    plt.xlabel('Protein Marker')
    plt.tight_layout()
    plt.show()
    
def matching_cell_list(adata,name_string,cell_type_string):
    
    """
    Retrieve the number of tartget cells for a given TMA punch 
    Parameters:
    - adata: AnnData object
    - name_string: The name of the punch which without the cell number. Ex "c_1_1" 
    - cell_type_string: string with the cell type which you want counts for.
    Returns:
    - 'Cell Count' int with the number of target cells in that punch
    - 'Cell Names' list of all the cell names that are PMNs
    - 'Overall Cell Count' The overall number of cells in the sample
    """
    # Ensure name_string ends with "_"
    if not name_string.endswith("_"):
        name_string += "_"
    cell_count=0
    end=False
    cell_type=""
    cell_names=[]
    i=0
    overall_cell_count=0
    
     # This variable sets the breakpoint requires >100 failed requests in a row. Note that there are rondom cells throughout the data sets
    j=0

    while end==False:
        i+=1
        test_string=f"{name_string}{i}"
        try: 
            cell_type=adata.obs.loc[test_string, "merged_annot_cluster"]
            j=0
            overall_cell_count+=1
            if cell_type_string in cell_type:
                cell_count+=1
                cell_names.append(f"{name_string}{i}")
        except:
           j+=1
        if j>100: 
            end=True

    # Calculate the overall percent of targeted cells
    if overall_cell_count>0:
        decimal_count=cell_count/overall_cell_count
    else:
        decimal_count=0
    
    return {
        "Cell Count": cell_count,
        "Cell Names": cell_names,
        "Overall Cell Count":overall_cell_count,
        "Decimal Count": decimal_count
    }

