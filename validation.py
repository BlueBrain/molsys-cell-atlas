#!/usr/bin/python
# -*- coding: latin-1 -*-



# =============================================================================================
# Import librairies

import numpy as np
from voxcell import VoxelData
from voxcell.nexus.voxelbrain import RegionMap
from argparse import ArgumentParser
import os, sys
import warnings
import inspect
# =============================================================================================



# =============================================================================================
# Parse arguments

parser = ArgumentParser()
parser.add_argument("annotation")
parser.add_argument("hierarchy")
parser.add_argument("density_folder")
parser.add_argument("cell_density")
parser.add_argument("inhibitory_density_folder")
parser.add_argument("excitatory_ME_types_folder")
parser.add_argument("inhibitory_ME_types_folder")
parser.add_argument("excitatory_ME_types_transplant_folder")
parser.add_argument("inhibitory_ME_types_transplant_folder")
# parser.add_argument("output_log_file_path")
# parser.add_argument("literature_values")

args = parser.parse_args()
brain_regions = args.annotation
hierarchy_json = args.hierarchy
density_folder = args.density_folder
cell_density = args.cell_density
inhibitory_density_folder = args.inhibitory_density_folder
excitatroy_ME_types_folder = args.excitatory_ME_types_folder
inhibitory_ME_types_folder = args.inhibitory_ME_types_folder
excitatory_ME_types_transplant_folder = args.excitatory_ME_types_transplant_folder
inhibitory_ME_types_transplant_folder = args.inhibitory_ME_types_transplant_folder
# output_log_file_path = args.output_log_file_path
# literature_values = args.literature_values
# =============================================================================================



# =============================================================================================
# Load files and filters initialization

# Reading files
print("\nReading the input files...")
annotation = VoxelData.load_nrrd(brain_regions).raw
cell = VoxelData.load_nrrd(cell_density).raw
neuron = VoxelData.load_nrrd(os.path.join(density_folder, "neuron_density.nrrd")).raw
glia = VoxelData.load_nrrd(os.path.join(density_folder, "glia_density.nrrd")).raw
astrocyte = VoxelData.load_nrrd(os.path.join(density_folder, "astrocyte_density.nrrd")).raw
microglia = VoxelData.load_nrrd(os.path.join(density_folder, "microglia_density.nrrd")).raw
oligodendrocyte = VoxelData.load_nrrd(os.path.join(density_folder, "oligodendrocyte_density.nrrd")).raw
gad = VoxelData.load_nrrd(os.path.join(inhibitory_density_folder, "gad67+_density.nrrd")).raw
excitatory_path_list = [os.path.join(excitatroy_ME_types_folder,f) for f in os.listdir(excitatroy_ME_types_folder) if f.endswith(".nrrd")]
excitatory_transplant_path_list = [os.path.join(excitatory_ME_types_transplant_folder,f) for f in os.listdir(excitatory_ME_types_transplant_folder) if f.endswith(".nrrd")]
inhibitory_path_list = [os.path.join(inhibitory_ME_types_folder,f) for f in os.listdir(inhibitory_ME_types_folder) if f.endswith(".nrrd")]
inhibitory_transplant_path_list = [os.path.join(inhibitory_ME_types_transplant_folder,f) for f in os.listdir(inhibitory_ME_types_transplant_folder) if f.endswith(".nrrd")]
pv = VoxelData.load_nrrd(os.path.join(inhibitory_density_folder, "pv+_density.nrrd")).raw
sst = VoxelData.load_nrrd(os.path.join(inhibitory_density_folder, "sst+_density.nrrd")).raw
vip = VoxelData.load_nrrd(os.path.join(inhibitory_density_folder, "vip+_density.nrrd")).raw
generic_excitatory = VoxelData.load_nrrd(os.path.join(excitatroy_ME_types_folder, "Generic_Excitatory_Neuron_MType_Generic_Excitatory_Neuron_EType.nrrd")).raw
generic_inhibitory = VoxelData.load_nrrd(os.path.join(excitatroy_ME_types_folder, "Generic_Inhibitory_Neuron_MType_Generic_Inhibitory_Neuron_EType.nrrd")).raw

print("Done")
# =============================================================================================



# =============================================================================================
# Region filters

print("\nRegion filters setting")
region_map = RegionMap.load_json(hierarchy_json)

# filter_isocortex = list(rmap.find("Isocortex", "acronym", with_descendants=True))
# filter_isocortex = np.isin(annotation, filter_isocortex)
# filter_layers = list(rmap.find("@.*[1-6][ab]?$", "acronym", with_descendants=True))
# filter_layers = np.isin(annotation, filter_layers)
# filter_layers = np.logical_and(filter_isocortex, filter_layers)
cerebellum = region_map.find(
        "Cerebellum", attr="name", with_descendants=True
    ) | region_map.find("arbor vitae", attr="name", with_descendants=True)
isocortex = (
    region_map.find("Isocortex", attr="name", with_descendants=True)
    | region_map.find("Entorhinal area", attr="name", with_descendants=True)
    | region_map.find("Piriform area", attr="name", with_descendants=True)
)
fiber_tracts_ids = (
    region_map.find("fiber tracts", attr="name", with_descendants=True)
    | region_map.find("grooves", attr="name", with_descendants=True)
    | region_map.find("ventricular systems", attr="name", with_descendants=True)
)
hippocampus = (
    region_map.find("Hippocampal region", attr="name", with_descendants=True)
)
hippocampal_formation = (
    region_map.find("Hippocampal formation", attr="name", with_descendants=True)
)
thalamus = (
    region_map.find("Thalamus", attr="name", with_descendants=True)
)
striatum = (
    region_map.find("Striatum", attr="name", with_descendants=True)
)
VPL = (
    region_map.find("Ventral posterolateral nucleus of the thalamus", attr="name", with_descendants=True)
)
LGd = (
    region_map.find("Dorsal part of the lateral geniculate complex", attr="name", with_descendants=True)
)
VPM = (
    region_map.find("Ventral posteromedial nucleus of the thalamus", attr="name", with_descendants=True)
)
MOB = (
    region_map.find("Main olfactory bulb", attr="name", with_descendants=True)
)
# rest_ids = region_map.find("root", attr="name", with_descendants=True)
# rest_ids -= cerebellum_group_ids | isocortex_group_ids

if cerebellum == [] or \
isocortex == [] or \
fiber_tracts_ids == [] or \
hippocampus == [] or \
striatum == []:
    raise ValueError("ERROR: some region filters return empty sets")

print("Done")
# =============================================================================================



# =============================================================================================
# Literature values setting

# # Literature values and initialization of paramters (to be moved into a configuration file)
print("\nLiterature values setting...")
voxel_volume = (25 * 1.0e-3)**3
nb_voxels_ccfv2_2011 = 31984720 # 499.8 mm^3
nb_voxels_ccfv2_2015 = 31992356 # = 499.9 mm^3
nb_voxels_ccfv2_2015_fiber_tracts = 3525860 # = 55.1 mm^3
nb_voxels_ccfv2_2015_whithout_fiber_tracts = nb_voxels_ccfv2_2015 - nb_voxels_ccfv2_2015_fiber_tracts # 444.8 mm^3
nb_voxels_ccfv3_2015 = 32387385 # = 506.1 mm^3
nb_voxels_ccfv3_2022 = 32261391 # = 504.1 mm^3
nb_voxels_ccfv3_2022_fiber_tracts = 2980469 # = 46.6 mm^3
nb_voxels_ccfv2_2022_whithout_fiber_tracts = nb_voxels_ccfv3_2022 - nb_voxels_ccfv3_2022_fiber_tracts # 457.5 mm^3
nb_voxels_ccfv3_2022_augmented = None # Not yet available
default_neuron_proportion = 0.16 # Default neuron density tolerance inheritated from total neuron tolerance when no data is available
default_glia_proportion = 0.20 # Default glia density tolerance inheritated from total glia tolerance when no data is available
default_cell_proportion = 0.18 # Default cell density tolerance inheritated from total cell tolerance when no data is available

# Literature values for the whole brain (to be moved into a configuration file)
wh_mouse_brain_vol_litt_m = 32570240 # = 508.91 mm^3 in Badea et al., 2007
wh_mouse_brain_vol_tolerance_m = 1498880 # = 23.42 mm^3 (5%) in Badea et al., 2007wh_mouse_brain_vol_litt_m = 32570240 # = 508.91 mm^3 in Badea et al., 2007
wh_mouse_brain_vol_litt = 508.91 # = 508.91 mm^3 in Badea et al., 2007
wh_mouse_brain_vol_tolerance = 23.42 # = 23.42 mm^3 (5%) in Badea et al., 2007
neuron_dens_fiber_tracts_litt = 0 # Rodarie et al., 2022
neuron_dens_fiber_tracts_tolerance = 0 # Rodarie et al., 2022
neuron_dens_litt = 71760000/wh_mouse_brain_vol_litt # = 67,870,000 + 3,890,000 Table 1 in Herculano-Houzel et al., 2011
neuron_dens_tolerance = 11660000/wh_mouse_brain_vol_litt # = 10,410,000 + 1,250,000 (16%) Table 1 in Herculano-Houzel et al., 2011
glia_dens_litt = 39320000/wh_mouse_brain_vol_litt # = 33,860,000 + 5,460,000 Table 1 in Herculano-Houzel et al., 2011
glia_dens_tolerance = 7810000/wh_mouse_brain_vol_litt # = 6,660,000 + 1,150,000 (20%) Table 1 in Herculano-Houzel et al., 2011
cell_dens_litt = 111080000/wh_mouse_brain_vol_litt # = 71,760,000 + 39,320,000 Table 1 in Herculano-Houzel et al., 2011 after summing neuron + glia
cell_dens_tolerance = 19470000/wh_mouse_brain_vol_litt # = 11,660,000 + 7,810,000 (18%) Table 1 in Herculano-Houzel et al., 2011 after summing neuron + glia
astrocyte_dens_litt = None # Not available
astrocyte_dens_tolerance = None # Not available
microglia_dens_litt = None # Not available
microglia_dens_tolerance = None # Not available
oligodendrocyte_dens_litt = None # Not available
oligodendrocyte_dens_tolerance = None # Not available
inhibitory_neuron_dens_litt = 14550000/wh_mouse_brain_vol_litt # Table 3 in Rodarie et al., 2022
inhibitory_neuron_dens_tolerance = round(inhibitory_neuron_dens_litt * default_neuron_proportion) # default Not available /!\
excitatory_neurons_dens_litt = None # Not available
excitatory_neurons_dens_tolerance = None # Not available
pv_dens_litt = 2631372/wh_mouse_brain_vol_litt # = 5916 * 445 in measurements.csv (~17.6% * Inhibitory in Table 3) in Rodarie et al., 2022
pv_dens_tolerance = 237072/wh_mouse_brain_vol_litt # = 533 (9%) in measurements.csv in Rodarie et al., 2022
sst_dens_litt = 2253658/wh_mouse_brain_vol_litt# = 5067 * 445 (~15.8% * Inhibitory in Table 3) in Rodarie et al., 2022
sst_dens_tolerance = 234498/wh_mouse_brain_vol_litt # = 527 (10%) in measurements.csv in Rodarie et al., 2022
vip_dens_litt = 434928/wh_mouse_brain_vol_litt # = 978 * 445 (~3.1% * Inhibitory in Table 3) in Rodarie et al., 2022
vip_dens_tolerance = 30616/wh_mouse_brain_vol_litt # 69 (7%) in measurements.csv in Rodarie et al., 2022
rest_inhi_dens_litt = (inhibitory_neuron_dens_litt - (pv_dens_litt + sst_dens_litt + vip_dens_litt))/wh_mouse_brain_vol_litt # (~63.5% * Inhibitory in Table 3) in Rodarie et al., 2022
rest_inhi_dens_tolerance = (inhibitory_neuron_dens_tolerance - (pv_dens_tolerance + sst_dens_tolerance + vip_dens_tolerance))/wh_mouse_brain_vol_litt # Rodarie et al., 2022

# Region support in number of voxels
isocortex_nb_vox = len(np.where(np.isin(annotation, list(isocortex)) != 0)[0])
hippocampus_nb_vox = len(np.where(np.isin(annotation, list(hippocampus)) != 0)[0])
striatum_nb_vox = len(np.where(np.isin(annotation, list(striatum)) != 0)[0])
thalamus_nb_vox = len(np.where(np.isin(annotation, list(thalamus)) != 0)[0])
VPL_nb_vox = len(np.where(np.isin(annotation, list(VPL)) != 0)[0])
LGd_nb_vox = len(np.where(np.isin(annotation, list(LGd)) != 0)[0])
VPM_nb_vox = len(np.where(np.isin(annotation, list(VPM)) != 0)[0])
MOB_nb_vox = len(np.where(np.isin(annotation, list(MOB)) != 0)[0])
hippocampal_formation_nb_vox = len(np.where(np.isin(annotation, list(hippocampal_formation)) != 0)[0])

# Literature values for sub-regions (to be moved into a configuration file)
isocortex_neuron_dens_litt = 2*5048837/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
isocortex_neuron_dens_tolerance = 2*412123/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
isocortex_glia_dens_litt = 2*6640234/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
isocortex_glia_dens_tolerance = 2*244643/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
isocortex_cell_dens_litt = (isocortex_neuron_dens_litt + isocortex_glia_dens_litt) # Table 1 in Herculano_Houzel et al., 2013
isocortex_cell_dens_tolerance = (isocortex_neuron_dens_tolerance + isocortex_glia_dens_tolerance) # Table 1 in Herculano_Houzel et al., 2013
isocortex_oligo_dens_litt = 12500 # Table 1 in Erö et al., 2018
isocortex_oligo_dens_tolerance = round(isocortex_oligo_dens_litt * default_glia_proportion) # default Not available /!\
isocortex_astro_dens_litt = 15696 # Table 1 in Erö et al., 2018
isocortex_astro_dens_tolerance = round(isocortex_astro_dens_litt * default_glia_proportion) # default Not available /!\
isocortex_microglia_dens_litt = 6500 # Table 1 in Erö et al., 2018
isocortex_microglia_dens_tolerance = round(isocortex_microglia_dens_litt * default_glia_proportion) # default Not available /!\
cerebellum_volume = 59.65 # mm^3 from Table 1 in Zhang et al., 2011
cerebellum_volume_torlerance = 3.65 # mm^3 from Table 1 in Zhang et al., 2011
cerebellum_neuron_dens_litt = 42220000/cerebellum_volume # Table 1 in Herculano_Houzel et al., 2011
cerebellum_neuron_dens_tolerance = 9280000/cerebellum_volume # Table 1 in Herculano_Houzel et al., 2011
cerebellum_glia_dens_litt = 6950000/cerebellum_volume # Table 1 in Herculano_Houzel et al., 2011
cerebellum_glia_dens_tolerance = 1500000/cerebellum_volume # Table 1 in Herculano_Houzel et al., 2011
cerebellum_cell_dens_litt = (cerebellum_neuron_dens_litt + cerebellum_glia_dens_litt) # Table 1 in Herculano_Houzel et al., 2011
cerebellum_cell_dens_tolerance = (cerebellum_neuron_dens_tolerance + cerebellum_glia_dens_tolerance) # Table 1 in Herculano_Houzel et al., 2011
cerebellum_oligo_dens_litt = 13750 # = Table 1 in Erö et al., 2018
cerebellum_oligo_dens_tolerance = 1768 # Table 1 in Erö et al., 2018
cerebellum_astro_dens_litt = 1512 # Table 1 in Erö et al., 2018
cerebellum_astro_dens_tolerance = round(cerebellum_astro_dens_litt * default_glia_proportion) # default Not available /!\
cerebellum_microglia_dens_litt = 8624 # Table 1 in Erö et al., 2018
cerebellum_microglia_dens_tolerance = 659 # Table 1 in Erö et al., 2018
hippocampus_neuron_dens_litt = 20848 # from Table 3 in Keller et al., 2018
hippocampus_neuron_dens_tolerance = hippocampus_neuron_dens_litt * default_neuron_proportion # default Not available /!\
hippocampus_oligo_dens_litt = 9425 # Table 1 in Erö et al., 2018
hippocampus_oligo_dens_tolerance = round(hippocampus_oligo_dens_litt * default_glia_proportion) # default Not available /!\
hippocampus_astro_dens_litt = 16737 # Table 1 in Erö et al., 2018
hippocampus_astro_dens_tolerance = 10496 # Table 1 in Erö et al., 2018
hippocampus_microglia_dens_litt = 3248 # Table 1 in Erö et al., 2018
hippocampus_microglia_dens_tolerance = 1563 # Table 1 in Erö et al., 2018
# striatum_neuron_dens_litt = 802679/(striatum_nb_vox * voxel_volume)  # Andsberg et al., 2001
# striatum_neuron_dens_tolerance = 31665/(striatum_nb_vox * voxel_volume)# Andsberg et al., 2001
striatum_neuron_dens_litt = 120110 # from Table 5 in Keller et al., 2018
striatum_neuron_dens_tolerance = 31800 # from Table 5 in Keller et al., 2018
striatum_oligo_dens_litt = 9950 # Table 1 in Erö et al., 2018
striatum_oligo_dens_tolerance = 4036 # Table 1 in Erö et al., 2018
striatum_astro_dens_litt = 9867 # Table 1 in Erö et al., 2018
striatum_astro_dens_tolerance = 5547 # Table 1 in Erö et al., 2018
striatum_microglia_dens_litt = 12101 # Table 1 in Erö et al., 2018
striatum_microglia_dens_tolerance = 1930 # Table 1 in Erö et al., 2018
LGd_neuron_dens_litt = 141000 # from Table 5 in Keller et al., 2018
LGd_neuron_dens_tolerance = 24000 # from Table 5 in Keller et al., 2018
thalamus_glia_dens_litt = 80966 # from Table 5 in Keller et al., 2018
thalamus_glia_dens_tolerance = thalamus_glia_dens_litt * default_glia_proportion # default Not available /!\
thalamus_cell_dens_litt = 122450 # from Table 5 in Keller et al., 2018
thalamus_cell_dens_tolerance = thalamus_cell_dens_litt * default_cell_proportion # default Not available /!\
VPM_neuron_dens_litt = 83100 # from Table 5 in Keller et al., 2018
VPM_neuron_dens_tolerance = 7900 # from Table 5 in Keller et al., 2018
VPL_neuron_dens_litt = 57466.92 # LNMC data from 2019 Thalamic Release Report pages 21-23 https://docs.google.com/document/d/1maQ8VIwaFeyOQpfUy6PoLcamp48NW9NUy6BPN_oLRQk/edit#heading=h.otpvbup5qa45
# VPL_neuron_dens_tolerance = 5201.4 # LNMC data from 2019 Thalamic Release Report pages 21-23 https://docs.google.com/document/d/1maQ8VIwaFeyOQpfUy6PoLcamp48NW9NUy6BPN_oLRQk/edit#heading=h.otpvbup5qa45
VPL_neuron_dens_tolerance = VPL_neuron_dens_litt * default_neuron_proportion # default Not available /!\
VPL_pv_dens_litt = 1238.528635	# from Kim et al., 2017
VPL_pv_dens_tolerance = 575.6900057 # from Kim et al., 2017
VPL_sst_dens_litt = 1609.664213 # from Kim et al., 2017
VPL_sst_dens_tolerance = 1186.589608 # from Kim et al., 2017
VPL_vip_dens_litt = 0.746468649 # from Kim et al., 2017
VPL_vip_dens_tolerance = 1.296087826 # from Kim et al., 2017
MOB_cell_dens_litt = 383148 # from Tables in Parrish-Aungst et al., 2007
MOB_cell_dens_tolerance = 27027 # from Tables in Parrish-Aungst et al., 2007
MOB_neuron_dens_litt = 246422 # from Tables in Parrish-Aungst et al., 2007
MOB_neuron_dens_tolerance = 17488 # from Tables in Parrish-Aungst et al., 2007
MOB_glia_dens_litt = 136725 # from Tables in Parrish-Aungst et al., 2007
MOB_glia_dens_tolerance = 9539 # from Tables in Parrish-Aungst et al., 2007



print("Done")
# =============================================================================================



# =============================================================================================
# Funtions

print("\nLaunching the assertions and writing ouptut result in the log file...")

# with open(output_log_file_path, "a") as log_file:

# Catching all prints to write them into a log file
# sys.stdout = log_file

# Functions
def print_range_bar(value, min_value, max_value, bar_length=40):
    """
    Printing the range bar of the given density compared to literature plus some statisics elements (std, z-score)
    @method print_range_bar
    @param {Float} value The input dentity value to assert
    @param {Float} min_value The min density value set by literature
    @param {Float} max_value The max density value set by literature
    @param {Integer} bar_length The length of the bar to plot
    @return {None}
    """
    progress = (value - min_value) / (max_value - min_value)
    progress = max(0, min(1, progress))
    arrow = ' ' * int(round(bar_length * progress)) + '*'
    spaces = ' ' * (bar_length - len(arrow))
    # Format the values as integers and display without decimal places
    value = int(round(value))
    min_value = int(round(min_value))
    max_value = int(round(max_value))
    mean_val = (min_value + max_value)/2
    std = mean_val - min_value
    z_score = round((value - mean_val)/std, 2)
    std_percentage = round((mean_val - min_value)/mean_val*100, 1)
    range_to_print = f'Range: [{arrow}{spaces}] {value}  | z-score = {z_score} |  -{std_percentage}% [{min_value}, {max_value}] +{std_percentage}%'
    if value < min_value:
        range_to_print = range_to_print.replace("[*", "*[")
        print(range_to_print)
        # try:
        #     print("diff = " + str(mean_val - value) + "  |   -" + str(round(mean_val - value)/mean_val*100,1) + "%")
        # except:
        #     print("EXCEPTION: low limit")
        #     pass
    elif value > max_value:
        # try:
        #     print("diff = " + str(mean_val - value) + "  |   +" + str(round(mean_val - value)/mean_val*100,1) + "%")
        # except:
        #     print("EXCEPTION: high limit")
        #     pass
        range_to_print = range_to_print.replace("*]", "]*")
        print(range_to_print)
    else:
        print(range_to_print)
    return


class DensityError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def z_score_assertion(value = 0, min_value = 0, max_value = 0, assertion_message = ""):
    """
    Asserting the z-score for a given assertion is in the right range from literature:
        - |z| <= 1: VALIDATED
        - 1 < |z| <= 2: WARNING
        - |z| > 2: ERROR
    @method z_score_assertion
    @param {Float} value The input dentity value to assert
    @param {Float} min_value The min density value set by literature
    @param {Float} max_value The max density value set by literature
    @param {String} assertion_message The assertion message to print if a Warning or an Error is raised
    @return {None}
    """
    mean_val = (min_value + max_value)/2
    std = mean_val - min_value
    z_score = round((value - mean_val)/std, 2)
    if abs(z_score) <= 1:
        print("Validated")
    elif (abs(z_score) > 1) and (abs(z_score) <= 2):
        warnings.warn(assertion_message, UserWarning)
        print("WARNING:", assertion_message)
    elif abs(z_score) > 2:
        print("ERROR:", assertion_message)
        raise DensityError(assertion_message)
    else:
        raise ValueError("Uknown value")
    return


def z_score_assertion_sub_regions(value = 0, min_value = 0, max_value = 0, assertion_message = ""):
    """
    Asserting the z-score for a given assertion is in the right range from literature:
        - |z| <= 1: VALIDATED
        - |z| > 1: WARNING
    @method z_score_assertion_sub_regions
    @param {Float} value The input dentity value to assert
    @param {Float} min_value The min density value set by literature
    @param {Float} max_value The max density value set by literature
    @param {String} assertion_message The assertion message to print if a Warning is raised
    @return {None}
    """
    mean_val = (min_value + max_value)/2
    std = mean_val - min_value
    z_score = round((value - mean_val)/std, 2)
    if abs(z_score) <= 1:
        print("Validated")
    elif abs(z_score) > 1:
        warnings.warn(assertion_message, UserWarning)
        print("WARNING:", assertion_message)
    else:
        raise ValueError("Uknown value")
    return


def z_score_assertion_after_transplant(value = 0, min_value = 0, max_value = 0, assertion_message = ""):
    """
    Asserting the z-score for a given assertion after transplant is in the right range from literature:
        - |z| <= 1: VALIDATED
        - |z| > 1: ERROR
    @method z_score_assertion
    @param {Float} value The input dentity value to assert
    @param {Float} min_value The min density value set by literature
    @param {Float} max_value The max density value set by literature
    @param {String} assertion_message The assertion message to print if a Warning or an Error is raised
    @return {None}
    """
    mean_val = (min_value + max_value)/2
    std = mean_val - min_value
    z_score = round((value - mean_val)/std, 2)
    if abs(z_score) <= 1:
        print("Validated")
    elif (abs(z_score) > 1):
        print("ERROR:", assertion_message)
        raise DensityError(assertion_message)
    else:
        raise ValueError("Uknown value")
    return

# =============================================================================================



# =============================================================================================
# Assertions


# 1/ TECHNICAL VALIDATION

# To be filled if need be


# 2/ SCIENTIFIC VALIDATION

# 2.1/ Assertion on whole brain volumetry and densities
print("Assertion on whole brain volumetry and densities...")

# Assertion on the total volumetry of the annotation
whole_brain_annotation_dens = len(np.where(annotation != 0)[0])
annotation_dens_diff = abs(wh_mouse_brain_vol_litt_m - whole_brain_annotation_dens)
print("\nAssertion on total annotation volumetry (mm^3)")
print_range_bar(whole_brain_annotation_dens, wh_mouse_brain_vol_litt_m - wh_mouse_brain_vol_tolerance_m, wh_mouse_brain_vol_litt_m + wh_mouse_brain_vol_tolerance_m)
assertion_message = "total annotation volumetry out of literature range"
z_score_assertion(whole_brain_annotation_dens, wh_mouse_brain_vol_litt_m - wh_mouse_brain_vol_tolerance_m, wh_mouse_brain_vol_litt_m + wh_mouse_brain_vol_tolerance_m, assertion_message)

# Assertion on total cell densities
cell_dens = np.sum(cell) / whole_brain_annotation_dens# * voxel_volume
cell_dens_diff = abs(cell_dens_litt - cell_dens)
print("\nAssertion on total cell densities (/mm^3)")
print_range_bar(cell_dens, cell_dens_litt - cell_dens_tolerance, cell_dens_litt + cell_dens_tolerance)
assertion_message = "total cell densities out of literature range"
z_score_assertion(cell_dens, cell_dens_litt - cell_dens_tolerance, cell_dens_litt + cell_dens_tolerance, assertion_message)

# Assertion on total neuron densities
neuron_dens = np.sum(neuron) / whole_brain_annotation_dens # * voxel_volume
neuron_dens_diff = abs(neuron_dens_litt - neuron_dens)
print("\nAssertion on total neuron densities (/mm^3)")
print_range_bar(neuron_dens, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance)
assertion_message = "total neuron densities out of literature range"
z_score_assertion(neuron_dens, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance, assertion_message)

# Assertion on total glia densities
glia_dens = np.sum(glia) / whole_brain_annotation_dens # * voxel_volume
glia_dens_diff = abs(glia_dens_litt - glia_dens)
print("\nAssertion on total glia densities (/mm^3)")
print_range_bar(glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance)
assertion_message = "total glia densities out of literature range"
z_score_assertion(glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance, assertion_message)

# Assertion on total sum of inhibitory + excitatory neuron densities
exci_inhib_sum = VoxelData.load_nrrd(excitatory_path_list[0]).raw
print("\nadding initialization 1/" + str(len(excitatory_path_list)))
for i in range (1, len(excitatory_path_list)):
    print("+ adding file " + str(i+1) + "/" + str(len(excitatory_path_list)))
    exci_inhib = VoxelData.load_nrrd(excitatory_path_list[i]).raw
    exci_inhib_sum += exci_inhib
exci_inhib_sum_sum = np.sum(exci_inhib_sum) / whole_brain_annotation_dens # * voxel_volume
diff_exci_inhib_sum = abs(neuron_dens_litt - exci_inhib_sum_sum)
print("\nAssertion on sum of inhibitory + excitatory neuron densities (/mm^3)")
print_range_bar(exci_inhib_sum_sum, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance)
# print("/!\ Data not available")
assertion_message = "sum of inhibitory + excitatory neuron densities out of literature range"
z_score_assertion(exci_inhib_sum_sum, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance, assertion_message)

# Assertion on total sum of inhibitory + excitatory neuron densities after transplant
exci_inhib_transplant_sum = VoxelData.load_nrrd(excitatory_transplant_path_list[0]).raw
print("\nadding initialization 1/" + str(len(excitatory_transplant_path_list)))
for i in range (1, len(excitatory_transplant_path_list)):
    print("+ adding file " + str(i+1) + "/" + str(len(excitatory_transplant_path_list)))
    exci_inhib_transplant = VoxelData.load_nrrd(excitatory_transplant_path_list[i]).raw
    exci_inhib_transplant_sum += exci_inhib_transplant
exci_inhib_transplant_sum_sum = np.sum(exci_inhib_transplant_sum) / whole_brain_annotation_dens # * voxel_volume
diff_exci_inhib_transplant_sum = abs(neuron_dens_litt - exci_inhib_transplant_sum_sum)
print("\nAssertion on sum of inhibitory + excitatory neuron densities after transplant (/mm^3)")
print_range_bar(exci_inhib_transplant_sum_sum, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance)
# print("/!\ Data not available")
assertion_message = "sum of inhibitory + excitatory neuron densities after transplant out of literature range"
z_score_assertion_after_transplant(exci_inhib_transplant_sum_sum, neuron_dens_litt - neuron_dens_tolerance, neuron_dens_litt + neuron_dens_tolerance, assertion_message)

# Assertion on total excitatory neuron density
print("\nAssertion on total excitatory neuron densities")
print("/!\ Literature data not available")

# Assertion on total sum of inhibitory neuron density subtypes pv + sst + vip + rest_inhi
print("\nAssertion on sum of inhibitory neuron density subtypes pv + sst + vip + rest_inhi")
print("/!\ Part of the data not available")

# Assertion on total inhibitory neuron densities
inhi_dens = np.sum(gad) / whole_brain_annotation_dens # * voxel_volume
diff_inhi = abs(inhibitory_neuron_dens_litt - inhi_dens)
print("\nAssertion on total inhibitory neuron densities (/mm^3)")
print("/!\ Tolerance set to default for neuron density")
# print("/!\ Default tolerance increased by a factor of 2.5")
print_range_bar(inhi_dens, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_tolerance, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance)
# assert (diff_inhi <= inhibitory_neuron_dens_tolerance * 2.5) #, f"diff_ini = " + str(diff_inhi) + "   |   tolerence = " + str(inhibitory_neuron_dens_tolerance)
assertion_message = "total inhibitory neuron densities out of literature range"
z_score_assertion(inhi_dens, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_tolerance, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance, assertion_message)

# Assertion on total sum of inhibitory ME-type neuron densities which should be inferior or equal to the total inhibitory neuron densities
inhib_sum = VoxelData.load_nrrd(inhibitory_path_list[0]).raw
print("\nadding initialization 1/" + str(len(inhibitory_path_list)))
for i in range (1, len(inhibitory_path_list)):
    print("+ adding file " + str(i+1) + "/" + str(len(inhibitory_path_list)))
    inhib = VoxelData.load_nrrd(inhibitory_path_list[i]).raw
    inhib_sum += inhib
inhib_sum_sum = np.sum(inhib_sum) / whole_brain_annotation_dens # * voxel_volume
# diff_inhib_sum = abs(inhibitory_neuron_dens_litt - inhib_sum_sum)
print("\nAssertion on total sum of inhibitory ME-type neuron densities which should be inferior or equal to the total inhibitory neuron densities")
print_range_bar(inhib_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_litt, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance)
# assert(diff_sum_neuron <= neuron_dens_tolerance)
assertion_message = "sum of inhibitory ME-type neuron densities out of literature range"
z_score_assertion(inhib_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_litt, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance, assertion_message)

# Assertion on total inhibitory sum density
tot_inhib_sum = inhib_sum + generic_inhibitory
tot_inhib_sum_sum = np.sum(tot_inhib_sum) / whole_brain_annotation_dens
diff_tot_inhib = abs(inhibitory_neuron_dens_litt - tot_inhib_sum_sum)
print("\nAssertion total inhibitory sum density (/mm^3)")
print_range_bar(tot_inhib_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_tolerance, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance)
assertion_message = "total inhibitory sum density out of literature range"
z_score_assertion(tot_inhib_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_tolerance, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance, assertion_message)

# Assertion on total sum of inhibitory ME-type neuron densities which should be inferior or equal to the total inhibitory neurons after transplant
inhib_transplant_sum = VoxelData.load_nrrd(inhibitory_transplant_path_list[0]).raw
print("\nadding initialization 1/" + str(len(inhibitory_transplant_path_list)))
for i in range (1, len(inhibitory_transplant_path_list)):
    print("+ adding file " + str(i+1) + "/" + str(len(inhibitory_transplant_path_list)))
    inhib_transplant = VoxelData.load_nrrd(inhibitory_transplant_path_list[i]).raw
    inhib_transplant_sum += inhib_transplant
inhib_transplant_sum_sum = np.sum(inhib_transplant_sum) / whole_brain_annotation_dens # * voxel_volume
# diff_inhib_transplant_sum = abs(inhibitory_neuron_dens_litt - inhib_transplant_sum_sum)
print("\nAssertion on total sum of inhibitory ME-type neuron densities after transplant which should be inferior or equal to the total inhibitory neurons")
print_range_bar(inhib_transplant_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_litt, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance)
assertion_message = "sum of inhibitory ME-type neuron densities after transplant out of literature range"
z_score_assertion_after_transplant(inhib_transplant_sum_sum, inhibitory_neuron_dens_litt - inhibitory_neuron_dens_litt, inhibitory_neuron_dens_litt + inhibitory_neuron_dens_tolerance, assertion_message)

# Assertion on total sum of glia density subtypes
sum_glia = astrocyte + microglia + oligodendrocyte
sum_glia_dens = int(np.sum(sum_glia)) / whole_brain_annotation_dens # * voxel_volume)
diff_sum_glia = abs(glia_dens_litt - sum_glia_dens)
print("\nAssertion on sum of glia density subtypes (/mm^3)")
print_range_bar(sum_glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance)
assertion_message = "total sum of glia density subtypes out of literature range"
z_score_assertion(sum_glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance, assertion_message)

# Assertion on total sum of glia density subtypes
sum_glia = astrocyte + microglia + oligodendrocyte
sum_glia_dens = int(np.sum(sum_glia)) / whole_brain_annotation_dens # * voxel_volume)
diff_sum_glia = abs(glia_dens_litt - sum_glia_dens)
print("\nAssertion on sum of glia density subtypes (/mm^3)")
print_range_bar(sum_glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance)
assertion_message = "total sum of glia density subtypes out of literature range"
z_score_assertion(sum_glia_dens, glia_dens_litt - glia_dens_tolerance, glia_dens_litt + glia_dens_tolerance, assertion_message)

# Assertion on total astrocyte densities
print("\nAssertion on total astrocyte densities")
print("/!\ Litterature data not available")

# Assertion on total microglia densities
print("\nAssertion on total microglia densities")
print("/!\ Litterature data not available")

# Assertion on total oligodendrocyte density
print("\nAssertion on total oligodendrocyte densities")
print("/!\ Data not available")

# Assertion on total sst densities
sst_dens = np.sum(sst) / whole_brain_annotation_dens # * voxel_volume
diff_sst_dens = abs(sst_dens_litt - sst_dens)
print("\nAssertion on total sst densities (/mm^3)")
print_range_bar(sst_dens, sst_dens_litt - sst_dens_tolerance, sst_dens_litt + sst_dens_tolerance)
assertion_message = "total sst densities out of literature range"
z_score_assertion(sst_dens, sst_dens_litt - sst_dens_tolerance, sst_dens_litt + sst_dens_tolerance, assertion_message)

# Assertion on total pv densities
pv_dens = np.sum(pv) / whole_brain_annotation_dens # * voxel_volume
diff_pv_dens = abs(pv_dens_litt - pv_dens)
print("\nAssertion on total pv densities (/mm^3)")
print_range_bar(pv_dens, pv_dens_litt - pv_dens_tolerance, pv_dens_litt + pv_dens_tolerance)
assertion_message = "total pv densities out of literature range"
z_score_assertion(pv_dens, pv_dens_litt - pv_dens_tolerance, pv_dens_litt + pv_dens_tolerance, assertion_message)

# Assertion on total vip density
vip_dens = np.sum(vip) / whole_brain_annotation_dens # * voxel_volume
diff_vip_dens = abs(vip_dens_litt - vip_dens)
print("\nAssertion on total vip densities (/mm^3)")
print_range_bar(vip_dens, vip_dens_litt - vip_dens_tolerance, vip_dens_litt + vip_dens_tolerance)
assertion_message = "total vip densities out of literature range"
z_score_assertion(vip_dens, vip_dens_litt - vip_dens_tolerance, vip_dens_litt + vip_dens_tolerance, assertion_message)

# Assertion on total rest_inhib density
print("\nAssertion on total rest_inhib densities")
print("/!\ Data not available")


# 2.2 Assertion on sub regions densities
print("\n\n==================================================")
print("\nAssertion on sub-region densities...")

# Assertion on fiber tracks + grooves + ventricular_system where no neuron should be found
fiber_tracts_neuron = neuron[np.isin(annotation, list(fiber_tracts_ids))]
fiber_tracts_neuron_sum = np.sum(fiber_tracts_neuron) * voxel_volume
diff_fiber_tracts_neuron_sum = abs(neuron_dens_fiber_tracts_litt - fiber_tracts_neuron_sum)
print("\nAssertion on fiber tracks + grooves + ventricular_system where no neuron should be found")
# assert(diff_fiber_tracts_neuron_sum <= neuron_dens_fiber_tracts_tolerance)
if not diff_fiber_tracts_neuron_sum <= neuron_dens_fiber_tracts_tolerance:
    warning_message = "fiber tracks + grooves + ventricular_system where no neuron should be found not consistent with literature"
    warnings.warn(warning_message, UserWarning)
    print("WARNING:", warning_message)
else:
    print("Validated")

# ---------------------------------------------------------------------------------------------
# ISOCORTEX
print("\n\n----------ISOCORTEX----------")

# Assertion on isocortex cell densities
isocortex_cell_dens = cell[np.isin(annotation, list(isocortex))]
isocortex_cell_dens_sum = np.sum(isocortex_cell_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_cell_dens = abs(isocortex_cell_dens_litt - isocortex_cell_dens_sum)
print("\nAssertion on isocortex cell densities (/mm^3)")
print_range_bar(isocortex_cell_dens_sum, isocortex_cell_dens_litt - isocortex_cell_dens_tolerance, isocortex_cell_dens_litt + isocortex_cell_dens_tolerance)
assertion_message = "isocortex cell densities out of literature range"
z_score_assertion_sub_regions(isocortex_cell_dens_sum, isocortex_cell_dens_litt - isocortex_cell_dens_tolerance, isocortex_cell_dens_litt + isocortex_cell_dens_tolerance, assertion_message)

# Assertion on isocortex neuron densities
isocortex_neuron_dens = neuron[np.isin(annotation, list(isocortex))]
isocortex_neuron_dens_sum = np.sum(isocortex_neuron_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_neuron_dens = abs(isocortex_neuron_dens_litt - isocortex_neuron_dens_sum)
print("\nAssertion on isocortex neuron densities (/mm^3)")
print("/!\ Tolerance set to default for neuron density")
isocortex_neuron_dens_default_tolerance = isocortex_neuron_dens_litt * default_neuron_proportion
print_range_bar(isocortex_neuron_dens_sum, isocortex_neuron_dens_litt - isocortex_neuron_dens_default_tolerance, isocortex_neuron_dens_litt + isocortex_neuron_dens_default_tolerance)
assertion_message = "isocortex neuron densities out of literature range"
z_score_assertion_sub_regions(isocortex_neuron_dens_sum, isocortex_neuron_dens_litt - isocortex_neuron_dens_default_tolerance, isocortex_neuron_dens_litt + isocortex_neuron_dens_default_tolerance, assertion_message)

# Assertion on isocortex glia densities
isocortex_glia_dens = glia[np.isin(annotation, list(isocortex))]
isocortex_glia_dens_sum = np.sum(isocortex_glia_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_glia_dens = abs(isocortex_glia_dens_litt - isocortex_glia_dens_sum)
print("\nAssertion on isocortex glia densities (/mm^3)")
print("/!\ Tolerance set to default for glia density")
isocortex_glia_dens_default_tolerance = isocortex_glia_dens_litt * default_glia_proportion
print_range_bar(isocortex_glia_dens_sum, isocortex_glia_dens_litt - isocortex_glia_dens_default_tolerance, isocortex_glia_dens_litt + isocortex_glia_dens_default_tolerance)
assertion_message = "isocortex glia densities out of literature range"
z_score_assertion_sub_regions(isocortex_glia_dens_sum, isocortex_glia_dens_litt - isocortex_glia_dens_default_tolerance, isocortex_glia_dens_litt + isocortex_glia_dens_default_tolerance, assertion_message)

# Assertion on isocortex oligodendrocyte densities
isocortex_oligo_dens = oligodendrocyte[np.isin(annotation, list(isocortex))]
isocortex_oligo_dens_sum = np.sum(isocortex_oligo_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_oligo_dens = abs(isocortex_oligo_dens_litt - isocortex_oligo_dens_sum)
print("\nAssertion on isocortex oligodendrocyte densities (/mm^3)")
print("/!\ Literature figures not consistent + Tolerance not available, set by default")
print_range_bar(isocortex_oligo_dens_sum, isocortex_oligo_dens_litt - isocortex_oligo_dens_tolerance, isocortex_oligo_dens_litt + isocortex_oligo_dens_tolerance)
assertion_message = "isocortex oligodendrocyte densities out of literature range"
z_score_assertion_sub_regions(isocortex_oligo_dens_sum, isocortex_oligo_dens_litt - isocortex_oligo_dens_tolerance, isocortex_oligo_dens_litt + isocortex_oligo_dens_tolerance, assertion_message)

# Assertion on isocortex astrocyte densities
isocortex_astro_dens = astrocyte[np.isin(annotation, list(isocortex))]
isocortex_astro_dens_sum = np.sum(isocortex_astro_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_astro_dens = abs(isocortex_astro_dens_litt - isocortex_astro_dens_sum)
print("\nAssertion on isocortex astrocyte densities (/mm^3)")
print("/!\ Literature figures not consistent + Tolerance not available, set by default")
print_range_bar(isocortex_astro_dens_sum, isocortex_astro_dens_litt - isocortex_astro_dens_tolerance, isocortex_astro_dens_litt + isocortex_astro_dens_tolerance)
assertion_message = "isocortex astrocyte densities out of literature range"
z_score_assertion_sub_regions(isocortex_astro_dens_sum, isocortex_astro_dens_litt - isocortex_astro_dens_tolerance, isocortex_astro_dens_litt + isocortex_astro_dens_tolerance, assertion_message)

# Assertion on isocortex microglia densities
isocortex_microglia_dens = microglia[np.isin(annotation, list(isocortex))]
isocortex_microglia_dens_sum = np.sum(isocortex_microglia_dens) / isocortex_nb_vox # * voxel_volume
diff_isocortex_microglia_dens = abs(isocortex_microglia_dens_litt - isocortex_microglia_dens_sum)
print("\nAssertion on isocortex microglia densities (/mm^3)")
print("/!\ Literature figures not consistent + Tolerance not available, set by default")
print_range_bar(isocortex_microglia_dens_sum, isocortex_microglia_dens_litt - isocortex_microglia_dens_tolerance, isocortex_microglia_dens_litt + isocortex_microglia_dens_tolerance)
assertion_message = "isocortex microglia densities out of literature range"
z_score_assertion_sub_regions(isocortex_microglia_dens_sum, isocortex_microglia_dens_litt - isocortex_microglia_dens_tolerance, isocortex_microglia_dens_litt + isocortex_microglia_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
# CEREBELLUM
print("\n\n----------CEREBELLUM----------")

# Assertion on cerebellum cell densities
cerebellum_nb_vox = len(np.where(np.isin(annotation, list(cerebellum)) != 0)[0])
cerebellum_cell_dens = cell[np.isin(annotation, list(cerebellum))]
cerebellum_cell_dens_sum = np.sum(cerebellum_cell_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_cell_dens = abs(cerebellum_cell_dens_litt - cerebellum_cell_dens_sum)
print("\nAssertion on cerebellum cell densities (/mm^3)")
print_range_bar(cerebellum_cell_dens_sum, cerebellum_cell_dens_litt - cerebellum_cell_dens_tolerance, cerebellum_cell_dens_litt + cerebellum_cell_dens_tolerance)
assertion_message = "cerebellum cell densities out of literature range"
z_score_assertion_sub_regions(cerebellum_cell_dens_sum, cerebellum_cell_dens_litt - cerebellum_cell_dens_tolerance, cerebellum_cell_dens_litt + cerebellum_cell_dens_tolerance, assertion_message)

# Assertion on cerebellum neuron densities
cerebellum_neuron_dens = neuron[np.isin(annotation, list(cerebellum))]
cerebellum_neuron_dens_sum = np.sum(cerebellum_neuron_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_neuron_dens = abs(cerebellum_neuron_dens_litt - cerebellum_neuron_dens_sum)
print("\nAssertion on cerebellum neuron densities (/mm^3)")
print_range_bar(cerebellum_neuron_dens_sum, cerebellum_neuron_dens_litt - cerebellum_neuron_dens_tolerance, cerebellum_neuron_dens_litt + cerebellum_neuron_dens_tolerance)
assertion_message = "cerebellum neuron densities out of literature range"
z_score_assertion_sub_regions(cerebellum_neuron_dens_sum, cerebellum_neuron_dens_litt - cerebellum_neuron_dens_tolerance, cerebellum_neuron_dens_litt + cerebellum_neuron_dens_tolerance, assertion_message)

# Assertion on cerebellum glia densities
cerebellum_glia_dens = glia[np.isin(annotation, list(cerebellum))]
cerebellum_glia_dens_sum = np.sum(cerebellum_glia_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_glia_dens = abs(cerebellum_glia_dens_litt - cerebellum_glia_dens_sum)
print("\nAssertion on cerebellum glia densities (/mm^3)")
print_range_bar(cerebellum_glia_dens_sum, cerebellum_glia_dens_litt - cerebellum_glia_dens_tolerance, cerebellum_glia_dens_litt + cerebellum_glia_dens_tolerance)
assertion_message = "cerebellum glia densities out of literature range"
z_score_assertion_sub_regions(cerebellum_glia_dens_sum, cerebellum_glia_dens_litt - cerebellum_glia_dens_tolerance, cerebellum_glia_dens_litt + cerebellum_glia_dens_tolerance, assertion_message)

# Assertion on cerebellum oligodendrocyte densities
cerebellum_oligo_dens = oligodendrocyte[np.isin(annotation, list(cerebellum))]
cerebellum_oligo_dens_sum = np.sum(cerebellum_oligo_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_oligo_dens = abs(cerebellum_oligo_dens_litt - cerebellum_oligo_dens_sum)
print("\nAssertion on cerebellum oligodendrocyte densities (/mm^3)")
print("/!\ Literature figures not consistent")
print_range_bar(cerebellum_oligo_dens_sum, cerebellum_oligo_dens_litt - cerebellum_oligo_dens_tolerance, cerebellum_oligo_dens_litt + cerebellum_oligo_dens_tolerance)
assertion_message = "cerebellum oligodendrocyte densities out of literature range"
z_score_assertion_sub_regions(cerebellum_oligo_dens_sum, cerebellum_oligo_dens_litt - cerebellum_oligo_dens_tolerance, cerebellum_oligo_dens_litt + cerebellum_oligo_dens_tolerance, assertion_message)

# Assertion on cerebellum astrocyte densities
cerebellum_astro_dens = astrocyte[np.isin(annotation, list(cerebellum))]
cerebellum_astro_dens_sum = np.sum(cerebellum_astro_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_astro_dens = abs(cerebellum_astro_dens_litt - cerebellum_astro_dens_sum)
print("\nAssertion on cerebellum astrocyte densities (/mm^3)")
print("/!\ Literature figures not consistent + Tolerance not available, set by default")
print_range_bar(cerebellum_astro_dens_sum, cerebellum_astro_dens_litt - cerebellum_astro_dens_tolerance, cerebellum_astro_dens_litt + cerebellum_astro_dens_tolerance)
assertion_message = "cerebellum astrocyte densities out of literature range"
z_score_assertion_sub_regions(cerebellum_astro_dens_sum, cerebellum_astro_dens_litt - cerebellum_astro_dens_tolerance, cerebellum_astro_dens_litt + cerebellum_astro_dens_tolerance, assertion_message)

# Assertion on cerebellum microglia densities
cerebellum_microglia_dens = microglia[np.isin(annotation, list(cerebellum))]
cerebellum_microglia_dens_sum = np.sum(cerebellum_microglia_dens) / cerebellum_nb_vox # * voxel_volume
diff_cerebellum_microglia_dens = abs(cerebellum_microglia_dens_litt - cerebellum_microglia_dens_sum)
print("\nAssertion on cerebellum microglia densities (/mm^3)")
print("/!\ Literature figures not consistent")
print_range_bar(cerebellum_microglia_dens_sum, cerebellum_microglia_dens_litt - cerebellum_microglia_dens_tolerance, cerebellum_microglia_dens_litt + cerebellum_microglia_dens_tolerance)
assertion_message = "cerebellum microglia densities out of literature range"
z_score_assertion_sub_regions(cerebellum_microglia_dens_sum, cerebellum_microglia_dens_litt - cerebellum_microglia_dens_tolerance, cerebellum_microglia_dens_litt + cerebellum_microglia_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# STRIATUM

print("\n\n----------STRIATUM----------")

# Assertion on striatum neuron densities
striatum_neuron_dens = neuron[np.isin(annotation, list(striatum))]
striatum_neuron_dens_sum = np.sum(striatum_neuron_dens) / striatum_nb_vox # * voxel_volume
diff_striatum_neuron_dens = abs(striatum_neuron_dens_litt - striatum_neuron_dens_sum)
print("\nAssertion on striatum neuron densities (/mm^3)")
print_range_bar(striatum_neuron_dens_sum, striatum_neuron_dens_litt - striatum_neuron_dens_tolerance, striatum_neuron_dens_litt + striatum_neuron_dens_tolerance)
assertion_message = "striatum neuron densities out of literature range"
z_score_assertion_sub_regions(striatum_neuron_dens_sum, striatum_neuron_dens_litt - striatum_neuron_dens_tolerance, striatum_neuron_dens_litt + striatum_neuron_dens_tolerance, assertion_message)

# Assertion on striatum oligodendrocyte densities
striatum_oligo_dens = oligodendrocyte[np.isin(annotation, list(striatum))]
striatum_oligo_dens_sum = np.sum(striatum_oligo_dens) / striatum_nb_vox # * voxel_volume
diff_striatum_oligo_dens = abs(striatum_oligo_dens_litt - striatum_oligo_dens_sum)
print("\nAssertion on striatum oligodendrocyte densities (/mm^3)")
print_range_bar(striatum_oligo_dens_sum, striatum_oligo_dens_litt - striatum_oligo_dens_tolerance, striatum_oligo_dens_litt + striatum_oligo_dens_tolerance)
assertion_message = "striatum oligodendrocyte densities out of literature range"
z_score_assertion_sub_regions(striatum_oligo_dens_sum, striatum_oligo_dens_litt - striatum_oligo_dens_tolerance, striatum_oligo_dens_litt + striatum_oligo_dens_tolerance, assertion_message)

# Assertion on striatum astrocyte densities
striatum_astro_dens = astrocyte[np.isin(annotation, list(striatum))]
striatum_astro_dens_sum = np.sum(striatum_astro_dens) / striatum_nb_vox # * voxel_volume
diff_striatum_astro_dens = abs(striatum_astro_dens_litt - striatum_astro_dens_sum)
print("\nAssertion on striatum astrocyte densities (/mm^3)")
print_range_bar(striatum_astro_dens_sum, striatum_astro_dens_litt - striatum_astro_dens_tolerance, striatum_astro_dens_litt + striatum_astro_dens_tolerance)
assertion_message = "striatum astrocyte densities out of literature range"
z_score_assertion_sub_regions(striatum_astro_dens_sum, striatum_astro_dens_litt - striatum_astro_dens_tolerance, striatum_astro_dens_litt + striatum_astro_dens_tolerance, assertion_message)

# Assertion on striatum microglia densities
striatum_microglia_dens = microglia[np.isin(annotation, list(striatum))]
striatum_microglia_dens_sum = np.sum(striatum_microglia_dens) / striatum_nb_vox # * voxel_volume
diff_striatum_microglia_dens = abs(striatum_microglia_dens_litt - striatum_microglia_dens_sum)
print("\nAssertion on striatum microglia densities (/mm^3)")
print_range_bar(striatum_microglia_dens_sum, striatum_microglia_dens_litt - striatum_microglia_dens_tolerance, striatum_microglia_dens_litt + striatum_microglia_dens_tolerance)
warning_message = "striatum microglia densities out of literature range"
z_score_assertion_sub_regions(striatum_microglia_dens_sum, striatum_microglia_dens_litt - striatum_microglia_dens_tolerance, striatum_microglia_dens_litt + striatum_microglia_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# HIPPOCAMPUS

print("\n\n----------HIPPOCAMPUS----------")

# Assertion on hippocampus neuron densities
hippocampus_neuron_dens = neuron[np.isin(annotation, list(hippocampal_formation))]
hippocampus_neuron_dens_sum = np.sum(hippocampus_neuron_dens) / hippocampal_formation_nb_vox # * voxel_volume
diff_hippocampus_neuron_dens = abs(hippocampus_neuron_dens_litt - hippocampus_neuron_dens_sum)
print("\nAssertion on hippocampus neuron densities (/mm^3)")
print("/!\ Tolerance not available, set by default")
print_range_bar(hippocampus_neuron_dens_sum, hippocampus_neuron_dens_litt - hippocampus_neuron_dens_tolerance, hippocampus_neuron_dens_litt + hippocampus_neuron_dens_tolerance)
assertion_message = "hippocampus neuron densities out of literature range"
z_score_assertion_sub_regions(hippocampus_neuron_dens_sum, hippocampus_neuron_dens_litt - hippocampus_neuron_dens_tolerance, hippocampus_neuron_dens_litt + hippocampus_neuron_dens_tolerance, assertion_message)

# Assertion on hippocampus oligodendrocyte densities
hippocampus_oligo_dens = oligodendrocyte[np.isin(annotation, list(hippocampus))]
hippocampus_oligo_dens_sum = np.sum(hippocampus_oligo_dens) / hippocampus_nb_vox # * voxel_volume
diff_hippocampus_oligo_dens = abs(hippocampus_oligo_dens_litt - hippocampus_oligo_dens_sum)
print("\nAssertion on hippocampus oligodendrocyte densities (/mm^3)")
print("/!\ Tolerance not available, set by default")
print_range_bar(hippocampus_oligo_dens_sum, hippocampus_oligo_dens_litt - hippocampus_oligo_dens_tolerance, hippocampus_oligo_dens_litt + hippocampus_oligo_dens_tolerance)
assertion_message = "hippocampus oligodendrocyte densities out of literature range"
z_score_assertion_sub_regions(hippocampus_oligo_dens_sum, hippocampus_oligo_dens_litt - hippocampus_oligo_dens_tolerance, hippocampus_oligo_dens_litt + hippocampus_oligo_dens_tolerance, assertion_message)

# Assertion on hippocampus astrocyte densities
hippocampus_astro_dens = astrocyte[np.isin(annotation, list(hippocampus))]
hippocampus_astro_dens_sum = np.sum(hippocampus_astro_dens) / hippocampus_nb_vox # * voxel_volume
diff_hippocampus_astro_dens = abs(hippocampus_astro_dens_litt - hippocampus_astro_dens_sum)
print("\nAssertion on hippocampus astrocyte densities (/mm^3)")
print_range_bar(hippocampus_astro_dens_sum, hippocampus_astro_dens_litt - hippocampus_astro_dens_tolerance, hippocampus_astro_dens_litt + hippocampus_astro_dens_tolerance)
assertion_message = "hippocampus astrocyte densities out of literature range"
z_score_assertion_sub_regions(hippocampus_astro_dens_sum, hippocampus_astro_dens_litt - hippocampus_astro_dens_tolerance, hippocampus_astro_dens_litt + hippocampus_astro_dens_tolerance, assertion_message)

# Assertion on hippocampus microglia densities
hippocampus_microglia_dens = microglia[np.isin(annotation, list(hippocampus))]
hippocampus_microglia_dens_sum = np.sum(hippocampus_microglia_dens) / hippocampus_nb_vox # * voxel_volume
diff_hippocampus_microglia_dens = abs(hippocampus_microglia_dens_litt - hippocampus_microglia_dens_sum)
print("\nAssertion on hippocampus microglia densities (/mm^3)")
print_range_bar(hippocampus_microglia_dens_sum, hippocampus_microglia_dens_litt - hippocampus_microglia_dens_tolerance, hippocampus_microglia_dens_litt + hippocampus_microglia_dens_tolerance)
assertion_message = "hippocampus microglia densities out of literature range"
z_score_assertion_sub_regions(hippocampus_microglia_dens_sum, hippocampus_microglia_dens_litt - hippocampus_microglia_dens_tolerance, hippocampus_microglia_dens_litt + hippocampus_microglia_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# THALAMUS

print("\n\n----------THALAMUS----------")

# Assertion on Talamus cell densities
thalamus_cell_dens = cell[np.isin(annotation, list(thalamus))]
thalamus_cell_dens_sum = np.sum(thalamus_cell_dens) / thalamus_nb_vox # * voxel_volume
diff_thalamus_cell_dens = abs(thalamus_cell_dens_litt - thalamus_cell_dens_sum)
print("\nAssertion on Thalamus cell densities (/mm^3)")
print("/!\ Tolerance set to default for cell density")
print_range_bar(thalamus_cell_dens_sum, thalamus_cell_dens_litt - thalamus_cell_dens_tolerance, thalamus_cell_dens_litt + thalamus_cell_dens_tolerance)
assertion_message = "Thalamus cell densities out of literature range"
z_score_assertion_sub_regions(thalamus_cell_dens_sum, thalamus_cell_dens_litt - thalamus_cell_dens_tolerance, thalamus_cell_dens_litt + thalamus_cell_dens_tolerance, assertion_message)

# Assertion on Thalamus glia densities
thalamus_glia_dens = cell[np.isin(annotation, list(thalamus))]
thalamus_glia_dens_sum = np.sum(thalamus_glia_dens) / thalamus_nb_vox # * voxel_volume
diff_thalamus_glia_dens = abs(thalamus_glia_dens_litt - thalamus_glia_dens_sum)
print("\nAssertion on Thalamus glia densities (/mm^3)")
print("/!\ Tolerance set to default for glia density")
print_range_bar(thalamus_glia_dens_sum, thalamus_glia_dens_litt - thalamus_glia_dens_tolerance, thalamus_glia_dens_litt + thalamus_glia_dens_tolerance)
assertion_message = "Thalamus cell densities out of literature range"
z_score_assertion_sub_regions(thalamus_glia_dens_sum, thalamus_glia_dens_litt - thalamus_glia_dens_tolerance, thalamus_glia_dens_litt + thalamus_glia_dens_tolerance, assertion_message)

# Assertion on LGd neuron densities
LGd_neuron_dens = neuron[np.isin(annotation, list(LGd))]
LGd_neuron_dens_sum = np.sum(LGd_neuron_dens) / LGd_nb_vox # * voxel_volume
diff_LGd_neuron_dens = abs(LGd_neuron_dens_litt - LGd_neuron_dens_sum)
print("\nAssertion on LGd neuron densities (/mm^3)")
print_range_bar(LGd_neuron_dens_sum, LGd_neuron_dens_litt - LGd_neuron_dens_tolerance, LGd_neuron_dens_litt + LGd_neuron_dens_tolerance)
assertion_message = "LGd neuron densities out of literature range"
z_score_assertion_sub_regions(LGd_neuron_dens_sum, LGd_neuron_dens_litt - LGd_neuron_dens_tolerance, LGd_neuron_dens_litt + LGd_neuron_dens_tolerance, assertion_message)

# Assertion on VPM neuron densities
VPM_neuron_dens = neuron[np.isin(annotation, list(VPM))]
VPM_neuron_dens_sum = np.sum(VPM_neuron_dens) / VPM_nb_vox # * voxel_volume
diff_VPM_neuron_dens = abs(VPM_neuron_dens_litt - VPM_neuron_dens_sum)
print("\nAssertion on VPM neuron densities (/mm^3)")
print_range_bar(VPM_neuron_dens_sum, VPL_neuron_dens_litt - VPM_neuron_dens_tolerance, VPM_neuron_dens_litt + VPM_neuron_dens_tolerance)
assertion_message = "VPM neuron densities out of literature range"
z_score_assertion_sub_regions(VPM_neuron_dens_sum, VPM_neuron_dens_litt - VPM_neuron_dens_tolerance, VPM_neuron_dens_litt + VPM_neuron_dens_tolerance, assertion_message)

# Assertion on VPL neuron densities
VPL_neuron_dens = neuron[np.isin(annotation, list(VPL))]
VPL_neuron_dens_sum = np.sum(VPL_neuron_dens) / VPL_nb_vox # * voxel_volume
diff_VPL_neuron_dens = abs(VPL_neuron_dens_litt - VPL_neuron_dens_sum)
print("\nAssertion on VPL neuron densities (/mm^3)")
print("/!\ Tolerance set to default for neuron density")
print_range_bar(VPL_neuron_dens_sum, VPL_neuron_dens_litt - VPL_neuron_dens_tolerance, VPL_neuron_dens_litt + VPL_neuron_dens_tolerance)
assertion_message = "VPL neuron densities out of literature range"
z_score_assertion_sub_regions(VPL_neuron_dens_sum, VPL_neuron_dens_litt - VPL_neuron_dens_tolerance, VPL_neuron_dens_litt + VPL_neuron_dens_tolerance, assertion_message)

# Assertion on VPL pv densities
VPL_pv_dens = pv[np.isin(annotation, list(VPL))]
VPL_pv_dens_sum = np.sum(VPL_pv_dens) / VPL_nb_vox # * voxel_volume
diff_VPL_pv_dens = abs(VPL_pv_dens_litt - VPL_pv_dens_sum)
print("\nAssertion on VPL pv densities (/mm^3)")
print_range_bar(VPL_pv_dens_sum, VPL_pv_dens_litt - VPL_pv_dens_tolerance, VPL_pv_dens_litt + VPL_pv_dens_tolerance)
assertion_message = "VPL pv densities out of literature range"
z_score_assertion_sub_regions(VPL_pv_dens_sum, VPL_pv_dens_litt - VPL_pv_dens_tolerance, VPL_pv_dens_litt + VPL_pv_dens_tolerance, assertion_message)

# Assertion on VPL sst densities
VPL_sst_dens = sst[np.isin(annotation, list(VPL))]
VPL_sst_dens_sum = np.sum(VPL_sst_dens) / VPL_nb_vox # * voxel_volume
diff_VPL_sst_dens = abs(VPL_sst_dens_litt - VPL_sst_dens_sum)
print("\nAssertion on VPL sst densities (/mm^3)")
print_range_bar(VPL_sst_dens_sum, VPL_sst_dens_litt - VPL_sst_dens_tolerance, VPL_sst_dens_litt + VPL_sst_dens_tolerance)
assertion_message = "VPL sst densities out of literature range"
z_score_assertion_sub_regions(VPL_sst_dens_sum, VPL_sst_dens_litt - VPL_sst_dens_tolerance, VPL_sst_dens_litt + VPL_sst_dens_tolerance, assertion_message)

# Assertion on VPL vip densities
VPL_vip_dens = vip[np.isin(annotation, list(VPL))]
VPL_vip_dens_sum = np.sum(VPL_vip_dens) / VPL_nb_vox # * voxel_volume
diff_VPL_vip_dens = abs(VPL_vip_dens_litt - VPL_vip_dens_sum)
print("\nAssertion on VPL vip densities (/mm^3)")
print_range_bar(VPL_vip_dens_sum, VPL_vip_dens_litt - VPL_vip_dens_tolerance, VPL_vip_dens_litt + VPL_vip_dens_tolerance)
assertion_message = "VPL vip densities out of literature range"
z_score_assertion_sub_regions(VPL_vip_dens_sum, VPL_vip_dens_litt - VPL_vip_dens_tolerance, VPL_vip_dens_litt + VPL_vip_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
# MAIN OLFACTORY BULB

print("\n\n----------MAIN OLFACTORY BULB----------")

# Assertion on MOB cell densities
MOB_cell_dens = cell[np.isin(annotation, list(MOB))]
MOB_cell_dens_sum = np.sum(MOB_cell_dens) / MOB_nb_vox # * voxel_volume
diff_MOB_cell_dens = abs(MOB_cell_dens_litt - MOB_cell_dens_sum)
print("\nAssertion on MOB cell densities (/mm^3)")
print_range_bar(MOB_cell_dens_sum, MOB_cell_dens_litt - MOB_cell_dens_tolerance, MOB_cell_dens_litt + MOB_cell_dens_tolerance)
assertion_message = "MOB_cell densities out of literature range"
z_score_assertion_sub_regions(MOB_cell_dens_sum, MOB_cell_dens_litt - MOB_cell_dens_tolerance, MOB_cell_dens_litt + MOB_cell_dens_tolerance, assertion_message)

# Assertion on MOB neuron densities
MOB_neuron_dens = neuron[np.isin(annotation, list(MOB))]
MOB_neuron_dens_sum = np.sum(MOB_neuron_dens) / MOB_nb_vox # * voxel_volume
diff_MOB_neuron_dens = abs(MOB_neuron_dens_litt - MOB_neuron_dens_sum)
print("\nAssertion on MOB neuron densities (/mm^3)")
print_range_bar(MOB_neuron_dens_sum, MOB_neuron_dens_litt - MOB_neuron_dens_tolerance, MOB_neuron_dens_litt + MOB_neuron_dens_tolerance)
assertion_message = "MOB_neuron densities out of literature range"
z_score_assertion_sub_regions(MOB_neuron_dens_sum, MOB_neuron_dens_litt - MOB_neuron_dens_tolerance, MOB_neuron_dens_litt + MOB_neuron_dens_tolerance, assertion_message)

# Assertion on MOB glia densities
MOB_glia_dens = glia[np.isin(annotation, list(MOB))]
MOB_glia_dens_sum = np.sum(MOB_glia_dens) / MOB_nb_vox # * voxel_volume
diff_MOB_glia_dens = abs(MOB_glia_dens_litt - MOB_glia_dens_sum)
print("\nAssertion on MOB glia densities (/mm^3)")
print_range_bar(MOB_glia_dens_sum, MOB_glia_dens_litt - MOB_glia_dens_tolerance, MOB_glia_dens_litt + MOB_glia_dens_tolerance)
assertion_message = "MOB_glia densities out of literature range"
z_score_assertion_sub_regions(MOB_glia_dens_sum, MOB_glia_dens_litt - MOB_glia_dens_tolerance, MOB_glia_dens_litt + MOB_glia_dens_tolerance, assertion_message)

# ---------------------------------------------------------------------------------------------



# Count the number of warnings raised ans tests passed
# print("\n" + str(successfully_passed) + " tests successfully passed")
# print("\nThe validation script raised {} warnings.".format(number_of_warnings))
print("\n\n==================================")
print("\nAll assertions successfully tested")

# # Writing all prints into a log file
# sys.stdout = sys.__stdout__
# print("Ouptut log file successfully written")
# =============================================================================================
