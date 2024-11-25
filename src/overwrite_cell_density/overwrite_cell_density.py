
#!/usr/bin/python
# -*- coding: latin-1 -*-


# =============================================================================================
# Libraries importation
import nrrd
import numpy as np
import os
from voxcell import RegionMap
import sys
import pandas
import json
import shutil
from pathlib import Path
import glob

from src.literature import voxel_volume


def absolute_nrrd_file_paths(directory):
    """
    Return the absolute path of all the files in a directory
    :param str directory: The target directory
    """

    for dir_path, _, filenames in os.walk(directory):
        for f in filenames:
            if ".nrrd" in f:
                yield os.path.abspath(os.path.join(dir_path, f))


def read_and_sum_files(cell_type_files):
    """
    Return the sum of all the nrrd files in a list.
    :param list of str cell_type_files: The list of files to be averaged
    """
    dens_matrix, hd_density_matrix = nrrd.read(cell_type_files[0])
    dens_matrix = 0 * dens_matrix

    for filename in cell_type_files:
        local_matrix, hd_density_matrix = nrrd.read(filename)
        dens_matrix = local_matrix + dens_matrix
    return dens_matrix, hd_density_matrix


def overwrite_cell_density_constrain_me(out_dict="", cell_t="", exp_list=[], reg_id=0, density_value=0., distrib="homogeneous"):
    """
    Perform the cell density scaling.
    :param dict out_dict: Where to put the files
    :param str cell_t: which cell type to scale
    :param str exp_list: The list of files to scale
    :param int reg_id: where to scale the density
    :param float density_value: the new density
    :param str distrib: can be either homogeneous or heterogeneous
    """

    def homogenous_processing(cell_type_files, density_val, mask, coord):
        """
        Scale with uniform distribution.
        :param list of str cell_type_files: the files to be scaled
        :param float density_val: the new density
        :param int array mask: region_mask in matrix for region
        :param int array coord: indices corresponding to region_mask
        """
        density_matr, header = read_and_sum_files(cell_type_files)

        total_mean = np.mean(density_matr[mask])

        for filename in cell_type_files:
            local_matrix, hd_local_matrix = nrrd.read(filename)
            if total_mean == 0:
                local_matrix[coord] = density_val / len(cell_type_files)
            else:
                current_mean = np.mean(density_matr[mask])
                local_matrix[coord] = local_matrix[coord] * current_mean / total_mean * density_val
            if np.sum(np.isnan(local_matrix.flat)):
                print("Nan detected homogenous")
                raise NaNError("Nan detected homogenous")

#            output_density_matrix_path = os.path.join(output_folder, os.path.basename(filename))
            nrrd.write(filename, local_matrix, header=hd_local_matrix)



    def heterogeneous_processing(cell_type_files, density_val, mask, coord):
        """
        Scale with heterogeneous distribution.
        :param list of str cell_type_files: the files to be scaled
        :param float density_val: the new density
        :param int array mask: region_mask in matrix for region
        :param int array coord: indices corresponding to region_mask

        """
        density_matr, header = read_and_sum_files(cell_type_files)

        current_mean = np.mean(density_matr[mask])
        ratio = density_val / current_mean

        for filename in cell_type_files:
            cell_matrix, hd_cell_matrix = nrrd.read(filename)
            if current_mean==0:
                cell_matrix[coord] = density_val / len(cell_type_files) #go for homogenous if no choice
            else:
                cell_matrix[coord] = cell_matrix[coord] * ratio
            if np.sum(np.isnan(cell_matrix.flat)):
                print("Nan detected heterogenous")
                raise NaNError("Nan detected homogenous")

            # save it
            nrrd.write(filename, cell_matrix, header=hd_cell_matrix)



    expand_list2 = []
    for cell in exp_list:
        full_path_n = os.path.join(out_dict[cell], (cell + ".nrrd"))
        expand_list2.append(full_path_n)

    # Input hard-coded paths: to be adapted for new versions of for any ME-type or m-type
    #annotation_path = "/gpfs/bbp.cscs.ch/home/piluso/cell_atlas/03_warped_annotation_fix_last/blue_brain_atlas_pipeline/leaves_only/annotation_ccfv2_l23split_barrelsplit.nrrd"
    #annotation_path = "/gpfs/bbp.cscs.ch/data/project/proj162/Model_Data/Brain_atlas/Mouse/resolution_25_um/version_1.1.0/Annotation_volume/annotation_ccfv3_l23split_barrelsplit_validated.nrrd"
    #annotation_path = "/gpfs/bbp.cscs.ch/home/dakeller/annotation_ccfv3_l23split_barrelsplit_validated.nrrd"    
    # Initialization of the hierarchy
    #hierarchy_input_path = "/gpfs/bbp.cscs.ch/home/piluso/cell_atlas/03_warped_annotation_fix_last/blue_brain_atlas_pipeline/leaves_only/hierarchy_ccfv2_l23split_barrelsplit.json"
    #hierarchy_input_path = "/gpfs/bbp.cscs.ch/data/project/proj162/Model_Data/Brain_atlas/Mouse/resolution_25_um/version_1.1.0/Parcellation_ontology/mba_hierarchy.json"
    #hierarchy_input_path = "/gpfs/bbp.cscs.ch/home/dakeller/mba_hierarchy.json"
    region_map = RegionMap.load_json(hierarchy_input_path)
    print("    Done: All files read")

    # Reading the input files
    print("\n1. Reading the input files...")
    print("    Reading the annotation file...")
    annotation, hd_annot = nrrd.read(annotation_path)

    # Creating the region_mask
    print("\n2. Creating the region_mask...")
    labels = region_map.find(reg_id, "id", with_descendants=True)
    region_name = region_map.get(reg_id, "name")
    print("    Region selected:", region_name)
    print("    Labels concerned:", labels)
    region_mask = np.isin(annotation, list(labels))
    coordinates = np.where(region_mask)
    nb_values_in_mask = np.sum(region_mask)
    region_volume = round(nb_values_in_mask * voxel_volume, 3)
    print("    Support: " + str(nb_values_in_mask) + " voxels")
    print("    Region_volume =", str(region_volume), "mm^3")
    print("    Done: Mask created")

    # Setting the artificial densities
    print("\n3. Setting the density...")
    print("    Setting density of the region to", density_value)
    print("    Cell type:", cell_t)
    print("    Cell expansion: ",exp_list)
    if distrib == "homogeneous":
        print("    Homogeneous distribution")
        homogenous_processing(expand_list2, density_value, region_mask, coordinates)
    elif distrib == "heterogeneous":
        print("    Heterogeneous distribution")
        heterogeneous_processing(expand_list2, density_value, region_mask, coordinates)

    print("    Done: Density overwritten")

    # Writing the output volumes
    print("\n4. Writing the output volumes...")
    print("    Writing cell density volume...")

    return 

class NaNError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def expand_list(key, cell_gr):
    """
    Returns a full list of cells corresponding to a group name
    :param str key : the name to be expanded
    :param dict of str cell_gr: the groups to be expanded
    """
    exp_list = []
    list1 = cell_gr[key]
    for val in list1:
        if val in cell_gr.keys():
            if val != key:
                exp_list = exp_list + expand_list(val, cell_gr)
            else:
                exp_list.append(val)
        else:
            exp_list.append(val)
    return exp_list


def get_directory_mapping(input_folder, in_dict, out_dict, output_folder):
    """
    Find where to put files.
    :param str input_folder: the input folder
    :param dict in_dict: the input file dictionary
    :param dict out_dict: the output file dictionary
    :param str output_folder: the output folder
    """

    files = glob.iglob(os.path.join(os.path.normpath(input_folder), "*.nrrd"))
    for file in files:
        if os.path.isfile(file):
            shutil.copy2(file, output_folder)

    only_files = [f for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, f))]
    for file1 in only_files:
        if ".nrrd" in file1:
            in_dict[file1.replace(".nrrd", "")] = input_folder
            out_dict[file1.replace(".nrrd", "")] = output_folder

    return in_dict, out_dict


"""
    Setting cell densities to a given region (does not handle glia cell subtypes) in a constrained way: retrieving the cell density proportion rules
    @method overwrite_cell_density_constrain
    @param {String} input_folder0 first folder to look for files
    @param {String} input_folder1 second folder to look for files
    @param {String} input_folder2 third folder to look for files
    @param {String} cell_groups_json The names of the nrrd files corresponding to each group of cells.
    @param {String} csv_file The master overwrite file.
    @param {String} output_folder0 The output folder where to write all volumes with the overwritten density corresponding to those in input_folder0.
    @param {String} output_folder1 The output folder where to write all volumes with the overwritten density corresponding to those in input_folder1.
    @param {String} output_folder2 The output folder where to write all volumes with the overwritten density corresponding to those in input_folder2. 
    @param {String} annotation_path where the nrrd annotation is
    @param {String} hierarchy_input_path where the json hierarchy file is
"""

# Executions

# Command to system


input_folder0 = os.path.abspath(sys.argv[1])
input_folder1 = os.path.abspath(sys.argv[2])
input_folder2 = os.path.abspath(sys.argv[3])
cell_groups_json = os.path.abspath(sys.argv[4])  # 'cell_groups.json'
overwrite_csv = os.path.abspath(sys.argv[5])  # 'overwrite_csv.csv'
output_folder0 = os.path.abspath(sys.argv[6])
output_folder1 = os.path.abspath(sys.argv[7])
output_folder2 = os.path.abspath(sys.argv[8])

annotation_path = os.path.abspath(sys.argv[9])

hierarchy_input_path = os.path.abspath(sys.argv[10])

input_dict = {}
output_dict = {}

for i, in_out in enumerate(zip([input_folder0, input_folder1, input_folder2], [output_folder0, output_folder1, output_folder2]), 1):
    # make output folders if they do not exist
    try:
        print("making ", in_out[1])
        os.makedirs(in_out[1], exist_ok=True)
    except:
        print("ERROR: dir creation failed")
        exit()

    # remove existing nrrds in output folder
    for path in Path(in_out[1]).glob("*.nrrd"):
        path.unlink()

    if i >= 2:
        if not sys.argv[i]:
            continue
    input_dict, output_dict = get_directory_mapping(in_out[0], input_dict, output_dict, in_out[1])


print("Starting overwrite_cell_density...")
# read the csv file
df = pandas.read_csv(overwrite_csv)
df = df.reset_index()  # make sure indexes pair with number of rows


# read the json file
with open(cell_groups_json) as json_data:
    cell_groups = json.load(json_data)
    json_data.close()


# parse the csv file and perform the scaling
results = {}
for index, row in df.iterrows():
    print(row['cell_type'], row['region'], row['density'], row['distribution'])

    cell_type = row['cell_type']
    region_id = int(row['region'])
    density = float(row['density'])
    distribution = row['distribution']

    expanded_list = []
    if cell_type in cell_groups.keys():
        expanded_list = (expand_list(cell_type, cell_groups))
    else:
        pathname = os.path.join(output_dict[cell_type], (cell_type+".nrrd"))
        if os.path.isfile(pathname):
            expanded_list = [cell_type]

    if len(expanded_list):
        overwrite_cell_density_constrain_me(output_dict, cell_type, expanded_list, region_id, density, distribution)
    else:
        print("No cell found")


print("\n5. Writing groups")

new_folder=os.path.join(output_folder2,"groups")
print(new_folder)
os.makedirs(new_folder, exist_ok=True)
# run the constraints on other cell groups that result from the modifications
for cell_type in cell_groups.keys():
    print("    Writing "+cell_type)
    expanded_list = []
    if cell_type in cell_groups.keys():
        expanded_list = (expand_list(cell_type, cell_groups))

    expanded_list2 = []
    for file_name in expanded_list:
        full_path_name = os.path.join(output_dict[file_name], (file_name + ".nrrd"))
        expanded_list2.append(full_path_name)
    pathname = os.path.join(new_folder, (cell_type+".nrrd"))
    density_matrix, header1 = read_and_sum_files(expanded_list2)
    nrrd.write(pathname, density_matrix, header=header1)


print("\nDensity successfully overwritten")
# =============================================================================================
