from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Allen Atlas connectivity and structure tree data pulls
mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
structure_tree = mcc.get_structure_tree()
root = structure_tree.get_structures_by_name(['root'])[0]
root_children = structure_tree.child_ids([root["id"]])[0]
name_map = structure_tree.get_name_map()
acro_map = structure_tree.get_id_acronym_map()
acro_map = {v: k for k, v in acro_map.items()} # flip the native mapping
hemisphere_name_map = {1: "left", 2: "right", 3: ""}
hemisphere_acro_map = {1: "-L", 2: "-R", 3: ""}


# get dataframe of experiments from wildtype (non cre) mice
all_experiments = mcc.get_experiments(dataframe=True, cre=False)


# get all structure unionizes from wildtype experiments
unz = pd.DataFrame()
for row in all_experiments.iterrows():
    structure_unionizes = mcc.get_structure_unionizes([row[1].id], is_injection=None, structure_ids=[root['id']], hemisphere_ids=[1, 2], include_descendants=True).reset_index(drop=True)
    structure_unionizes["inj_id"] = row[1].structure_id
    structure_unionizes["inj_hemi"] = (structure_unionizes["hemisphere_id"][structure_unionizes["is_injection"]]).unique()[0]
    unz = pd.concat([unz, structure_unionizes], axis=0)

# adding/changing name columns to match LifeCanvas
unz["name"] = unz["structure_id"].map(name_map)
unz["hemisphere_name"] = unz["hemisphere_id"].map(hemisphere_name_map)
unz["hemisphere_acro"] = unz["hemisphere_id"].map(hemisphere_acro_map)
unz["acronym"] = unz["structure_id"].map(acro_map)
unz["name_LR"] = unz["hemisphere_name"] + " " + unz["name"]
unz["acronym_LR"] = unz["acronym"] + unz["hemisphere_acro"]

unz["inj_name"] = unz["inj_id"].map(name_map)
unz["inj_name_LR"] = unz["inj_hemi"].map(hemisphere_name_map) + " " + unz["inj_name"]
unz["inj_acronym_LR"] = unz["inj_id"].map(acro_map) + unz["hemisphere_acro"]

unz_inj_only = unz[unz["name_LR"].isin(unz["inj_name_LR"])]

pivot_pd = unz.groupby(["inj_acronym_LR","acronym_LR"]).agg({"projection_density": "mean"}).reset_index().pivot(index='inj_acronym_LR', columns='acronym_LR')['projection_density']

pivot_pd_inj = unz_inj_only.groupby(["inj_acronym_LR","acronym_LR"]).agg({"projection_density": "mean"}).reset_index().pivot(index='inj_acronym_LR', columns='acronym_LR')['projection_density']

pivot_en_inj = unz_inj_only.groupby(["inj_acronym_LR","acronym_LR"]).agg({"projection_energy": "mean"}).reset_index().pivot(index='inj_acronym_LR', columns='acronym_LR')['projection_energy']



X = pivot_en_inj.T
distorsions = []
for k in range(2, 20):
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(X)
    distorsions.append(kmeans.inertia_)

fig = plt.figure(figsize=(15, 5))
plt.plot(range(2, 20), distorsions)
plt.grid(True)
plt.title('Elbow curve')

# n_clust = 10
# kmeans = KMeans(n_clusters=n_clust)
# kmeans.fit(X)

# labels = pd.DataFrame({"name_LR": pivot_pd_inj.index, "label": kmeans.labels_})

# labels_sorted = labels.sort_values("label")

# row_sort = pivot_pd_inj.iloc[labels_sorted.index]
# col_sort = row_sort.reindex(columns=labels_sorted["name_LR"])
# plt.imshow(col_sort)

# plt.show()
\

# https://wil.yegelwel.com/cluster-correlation-matrix/

import scipy.cluster.hierarchy as sch
import seaborn as sns
plt.gcf().subplots_adjust(bottom=0.5, left=0.5)

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]

mask_yfp = np.triu(yfp_region_corr_cluster)
mask_chr2 = np.triu(chr2_region_corr_cluster)

# output pattern correlation by injection site
# try penalizing this by structure centroid distance (are there any regions separated by a large distance that nonetheless project to similar regions)
pivot_pd_output = pivot_pd.T.fillna(0)
pivot_pd_output_corr = pivot_pd_output.corr()
pivot_pd_output_clust = cluster_corr(pivot_pd_output_corr)
output_mask = np.triu(pivot_pd_output_clust)
plt.figure(figsize = (15,15))
sns.heatmap(pivot_pd_output_clust, cmap = "icefire", vmax=1, vmin=-1, mask=output_mask, cbar_kws={'label': 'correlation'}, xticklabels=1, yticklabels=1)
plt.title("Output correlation")
plt.show()


# iput pattern correlation by injection site
pivot_pd_input = pivot_pd.fillna(0)
pivot_pd_input_corr = pivot_pd_input.corr()
pivot_pd_input_clust = cluster_corr(pivot_pd_input_corr)
input_mask = np.triu(pivot_pd_input_clust)
plt.figure(figsize = (15,15))
sns.heatmap(pivot_pd_input_clust, cmap = "icefire", vmax=1, vmin=-1, cbar_kws={'label': 'correlation'}, xticklabels=1, yticklabels=1)
plt.ylim([500,600])
plt.xlim([500,600])
plt.title("Input correlation")
plt.show()


pivot_npv = unz.groupby(["inj_name_LR","name_LR"]).agg({"normalized_projection_volume": "mean"}).reset_index().pivot(index='inj_name_LR', columns='name_LR')['normalized_projection_volume']


structure_unionizes = mcc.get_structure_unionizes([pl_exp_ids[0]], is_injection=None, structure_ids=[root['id']], include_descendants=True).reset_index(drop=True)


pm = mcc.get_projection_matrix(experiment_ids = exp_ids, 
                               projection_structure_ids = root_children,
                               hemisphere_ids= [1,2,3], # right hemisphere, ipsilateral
                               parameter = 'projection_density')




structure_tree = mcc.get_structure_tree()

structures = structure_tree.get_structures_by_name(["Prelimbic area"])
pd.DataFrame(structures)


pl_experiment = mcc.get_experiments(cre=False, injection_structure_ids=[pl['id']])


pl_exp_ids = [i['id'] for i in pl_experiment]
pm = mcc.get_projection_matrix(experiment_ids = pl_exp_ids, 
                               projection_structure_ids = root_children,
                               hemisphere_ids= [2], # right hemisphere, ipsilateral
                               parameter = 'projection_density')
projd, pd_info = mcc.get_projection_density(pl_exp_ids[0])

ind, ind_info = mcc.get_injection_density(pl_exp_ids[0])

