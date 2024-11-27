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
import re
from src.literature import *
# =============================================================================================



# =============================================================================================
# Parse arguments

def parse_args(args):
    """Parse command line parameters
    Args:
      args ([str]): command line parameters as list of strings
    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """

    parser = ArgumentParser()

    parser.add_argument(
        "--annotation",
        dest = "brain_regions",
        required=True,
        metavar="<FILE PATH>",
        help="The annotation file incuding labels")

    parser.add_argument(
        "--error_fatal",
        dest = "error_fatal",
        required=False,
        metavar="int",
        default=1,
        help="Whether or not to catastrophically fail if an error is found.")

    parser.add_argument(
        "--hierarchy",
        dest="hierarchy_json",
        required=False,
        metavar="<FILE PATH>",
        help="The hierarchy JSON file, sometimes called 1.json")

    parser.add_argument(
        "--neuron_glia_density_folder",
        dest="neuron_glia_density_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find the density files calculated for neuron, glia, microglia, oligodensdrocyte and astrocyte")

    parser.add_argument(
        "--cell_density",
        dest="cell_density",
        required=False,
        metavar="<FILE PATH>",
        help="The cell density file calculated for the whole brain")

    parser.add_argument(
        "--inhibitory_density_folder",
        dest="inhibitory_density_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find all the inhibitory density files gad67, pv, sst and vip")

    parser.add_argument(
        "--excitatory_ME_types_folder",
        dest="excitatory_ME_types_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find all the excitatory split ME-types density files")

    parser.add_argument(
        "--inhibitory_ME_types_folder",
        dest="inhibitory_ME_types_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find all the inhibitory split ME-types density files")

    parser.add_argument(
        "--excitatory_ME_types_transplant_folder",
        dest="excitatory_ME_types_transplant_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find all the excitatory split ME-types density files after transplant")

    parser.add_argument(
        "--inhibitory_ME_types_transplant_folder",
        dest="inhibitory_ME_types_transplant_folder",
        required=False,
        metavar="<DIRECTORY PATH>",
        help="The folder where to find all the inhibitory split ME-types density files after transplant")

    return parser.parse_args(args)

# =============================================================================================



# =============================================================================================

# print("\nLaunching the assertions and writing ouptut result in the log file...")

# with open(output_log_file_path, "a") as log_file:

# Catching all prints to write them into a log file
# sys.stdout = log_file

# Functions
def flatten_list(ll):
    return [i for sublist in ll for i in sublist]

def assert_densities(annotation,flattened_region_list, dens_type, gad,error_fatal, neuron=[]):
    """
    Performs assertion on a specific density
    @method assert_density
    @param {np.array} annotation
    @param {list} flattened list of region names
    @param {str} type of density, e.g. "excitatory neuron density"
    @param {np.array} gad
    @param {np.array/None} neuron array. If provided, uses neuron - gad, otherwise only uses gad
    @return {None}

    """
    for region_ids in flattened_region_list:
        assert_density(annotation, region_ids, dens_type, gad,error_fatal, neuron=neuron, region_label=None)


def assert_density(annotation, region_ids, dens_type, gad, error_fatal, neuron=[], region_label=None):
    """
    Performs assertion on a specific density
    @method assert_density
    @param {np.array} annotation
    @param {str} region names
    @param {str} type of density, e.g. "excitatory neuron density"
    @param {np.array} gad
    @param {np.array/None} neuron array. If provided, uses neuron - gad, otherwise only uses gad
    @return {None}

    """
    if len(neuron):
        dens = neuron[np.isin(annotation, list(region_ids))] - gad[np.isin(annotation, list(region_ids))]
    else:
        dens = gad[np.isin(annotation, list(region_ids))]
    region_label = region_label if region_label else region_ids
    print(f"Assertion on {dens_type} for {region_label} for error flag {error_fatal}")
    assertion_message = f"ERROR: {dens_type} is zero for {region_label} with id(s) {region_ids}"
    if dens.sum() > 0:
        print(f"Validated: {dens_type} is not zero for {region_label}")
    else:
        print(assertion_message)
        if error_fatal:
            raise DensityError(assertion_message)

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
    elif value > max_value:
        range_to_print = range_to_print.replace("*]", "]*")
        print(range_to_print)
    else:
        print(range_to_print)
    return


class DensityError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def z_score_assertion(value = 0, min_value = 0, max_value = 0, assertion_message = "",error_fatal = 1, z_error=None):
    """
    Asserting the z-score for a given assertion is in the right range from literature:
        - |z| <= 1: VALIDATED
        - 1 < |z| <= z_error: WARNING
        - |z| > z_error: ERROR
    @method z_score_assertion
    @param {Float} value The input dentity value to assert
    @param {Float} min_value The min density value set by literature
    @param {Float} max_value The max density value set by literature
    @param {String} assertion_message The assertion message to print if a Warning or an Error is raised
    @param {Int} z_error The z_score threshold to raise an error
    @return {None}
    """
    mean_val = (min_value + max_value)/2
    std = mean_val - min_value
    z_score = round((value - mean_val)/std, 2)

    z_error = z_error if z_error else 99
    if abs(z_score) <= 1:
        print("Validated")
    elif (abs(z_score) > 1) and (abs(z_score) <= z_error):
        warnings.warn(assertion_message, UserWarning)
        print("WARNING:", assertion_message)
    elif abs(z_score) > z_error:
        print("ERROR:", assertion_message)
        if error_fatal:
            raise DensityError(assertion_message)
    else:
        if error_fatal:
            raise ValueError("Unknown value")
    return


def validate_density_volume(density_name, density_volume, min_value,
                            max_value, z_error,error_fatal, tolerance_default_msg=None):
    print(f"\nAssertion on average {density_name} (/mm^3)")
    if tolerance_default_msg:
        print(f"/!\ Tolerance set to default for {tolerance_default_msg}")
    print_range_bar(density_volume, min_value, max_value)
    assertion_message = f"Average {density_name} out of literature range"
    z_score_assertion(density_volume, min_value, max_value, assertion_message,error_fatal,
                      z_error=z_error)


def validate_average_density(density_name, density_volume, voxel_number, min_value,
                             max_value, z_error,error_fatal, tolerance_default_msg=None):
    average_density = np.sum(density_volume) / voxel_number
    validate_density_volume(density_name, average_density, min_value, max_value,
                            z_error,error_fatal, tolerance_default_msg)


def validate_density_path(density_name, density_path, voxel_number, min_value,
                          max_value, z_error,error_fatal, tolerance_default_msg=None):
    # Reading the input file
    density_volume = VoxelData.load_nrrd(density_path).raw
    # Assertion on average density
    validate_average_density(density_name, density_volume, voxel_number, min_value,
                            max_value, z_error, tolerance_default_msg,error_fatal)
    return density_volume

# =============================================================================================
# Load files and perform assertions

def main():
    """
    Main entry point allowing external calls
    """

    # Parse args
    args = parse_args(sys.argv[1:])
    error_fatal=bool(int(args.error_fatal))
 
    # 2.1/ Assertion on whole brain volumetry and densities
    print("Assertion on whole brain volumetry and densities...")

    if args.brain_regions is not None:

        # Reading the input file
        annotation = VoxelData.load_nrrd(args.brain_regions).raw

        # Assertion on the total volumetry of the annotation
        whole_brain_annotation_nb_voxels = len(np.where(annotation != 0)[0])
        whole_brain_annotation_vol = whole_brain_annotation_nb_voxels * voxel_volume
        validate_density_volume(density_name="total annotation volumetry",
            density_volume=whole_brain_annotation_vol,
            min_value=wh_mouse_brain_nb_vox_lit*voxel_volume - wh_mouse_brain_nb_vox_tolerance*voxel_volume,
            max_value=wh_mouse_brain_nb_vox_lit*voxel_volume + wh_mouse_brain_nb_vox_tolerance*voxel_volume,
            z_error=2,error_fatal=error_fatal)

    if args.cell_density is not None:
        # Assertion on average cell density
        cell = validate_density_path(density_name="cell density", density_path=args.cell_density,
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=cell_dens_lit - cell_dens_tolerance,
            max_value=cell_dens_lit + cell_dens_tolerance, z_error=2,error_fatal=error_fatal)

    if args.neuron_glia_density_folder is not None:
        # Assertion on average neuron density
        neuron = validate_density_path(density_name="neuron density", density_path=os.path.join(args.neuron_glia_density_folder, "neuron_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=neuron_dens_lit - neuron_dens_tolerance,
            max_value=neuron_dens_lit + neuron_dens_tolerance, z_error=2,error_fatal=error_fatal)

        # Assertion on average glia density
        glia = validate_density_path(density_name="glia density", density_path=os.path.join(args.neuron_glia_density_folder, "glia_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=glia_dens_lit - glia_dens_tolerance,
            max_value=glia_dens_lit + glia_dens_tolerance, z_error=2,error_fatal=error_fatal)


        astrocyte = VoxelData.load_nrrd(os.path.join(args.neuron_glia_density_folder, "astrocyte_density.nrrd")).raw
        microglia = VoxelData.load_nrrd(os.path.join(args.neuron_glia_density_folder, "microglia_density.nrrd")).raw
        oligodendrocyte = VoxelData.load_nrrd(os.path.join(args.neuron_glia_density_folder, "oligodendrocyte_density.nrrd")).raw
        # Assertion on average sum of glia density subtypes
        sum_glia = astrocyte + microglia + oligodendrocyte
        validate_average_density("sum of glia subtypes densities", sum_glia,
            whole_brain_annotation_nb_voxels, glia_dens_lit - glia_dens_tolerance,
            glia_dens_lit + glia_dens_tolerance, z_error=2,error_fatal=error_fatal)

        # Assertion on average astrocyte density
        print("\nAssertion on average astrocyte density")
        print("/!\ Litterature data not available")

        # Assertion on average microglia density
        print("\nAssertion on average microglia density")
        print("/!\ Litterature data not available")

        # Assertion on average oligodendrocyte density
        print("\nAssertion on average oligodendrocyte density")
        print("/!\ Data not available")



    if args.inhibitory_density_folder is not None:

        # Reading the input files

        # Assertion on average inhibitory neuron density
        gad = validate_density_path(density_name="inhibitory neuron density", density_path=os.path.join(args.inhibitory_density_folder, "gad67+_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=inhibitory_neuron_dens_lit - inhibitory_neuron_dens_tolerance,
            max_value=inhibitory_neuron_dens_lit + inhibitory_neuron_dens_tolerance,
            z_error=None,error_fatal=error_fatal, tolerance_default_msg="neuron density")

        # Assertion on average sst density
        sst = validate_density_path(density_name="sst density", density_path=os.path.join(args.inhibitory_density_folder, "sst+_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=sst_dens_lit - sst_dens_tolerance,
            max_value=sst_dens_lit + sst_dens_tolerance, z_error=2,error_fatal=error_fatal)

        # Assertion on average pv density
        pv = validate_density_path(density_name="pv density", density_path=os.path.join(args.inhibitory_density_folder, "pv+_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=pv_dens_lit - pv_dens_tolerance,
            max_value=pv_dens_lit + pv_dens_tolerance, z_error=2,error_fatal=error_fatal)

        # Assertion on average vip density
        vip = validate_density_path(density_name="vip density", density_path=os.path.join(args.inhibitory_density_folder, "vip+_density.nrrd"),
            voxel_number=whole_brain_annotation_nb_voxels,
            min_value=vip_dens_lit - vip_dens_tolerance,
            max_value=vip_dens_lit + vip_dens_tolerance, z_error=2,error_fatal=error_fatal)

        # Assertion on average rest_inhib density
        print("\nAssertion on average rest_inhib density")
        print("/!\ Data not available")

        # Assertion on average sum of inhibitory neuron density subtypes pv + sst + vip + rest_inhi
        print("\nAssertion on sum of inhibitory neuron density subtypes pv + sst + vip + rest_inhi")
        print("/!\ Part of the data not available")


    if args.inhibitory_ME_types_folder is not None:

        # Reading the input file
        inhibitory_path_list = [os.path.join(args.inhibitory_ME_types_folder,f) for f in os.listdir(args.inhibitory_ME_types_folder) if f.endswith(".nrrd")]

        # Assertion on average sum of inhibitory ME-type neuron densities which should be inferior or equal to the average inhibitory neuron density
        inhib_sum = VoxelData.load_nrrd(inhibitory_path_list[0]).raw
        print("\nadding initialization 1/" + str(len(inhibitory_path_list)))
        for i in range (1, len(inhibitory_path_list)):
            print("+ adding file " + str(i+1) + "/" + str(len(inhibitory_path_list)))
            inhib = VoxelData.load_nrrd(inhibitory_path_list[i]).raw
            inhib_sum += inhib
        validate_average_density("sum of inhibitory ME-type neuron densities which should be inferior or equal to the average inhibitory neuron density", inhib_sum,
            whole_brain_annotation_nb_voxels,
            inhibitory_neuron_dens_lit - inhibitory_neuron_dens_lit,
            inhibitory_neuron_dens_lit + inhibitory_neuron_dens_tolerance, z_error=2,error_fatal=error_fatal)


    if args.excitatory_ME_types_folder is not None:

        # Reading the input files
        excitatory_path_list = [os.path.join(args.excitatory_ME_types_folder,f) for f in os.listdir(args.excitatory_ME_types_folder) if f.endswith(".nrrd")]
        # generic_excitatory = VoxelData.load_nrrd(os.path.join(args.excitatory_ME_types_folder, "Generic_Excitatory_Neuron_MType_Generic_Excitatory_Neuron_EType.nrrd")).raw
        generic_inhibitory = VoxelData.load_nrrd(os.path.join(args.excitatory_ME_types_folder, "Generic_Inhibitory_Neuron_MType_Generic_Inhibitory_Neuron_EType.nrrd")).raw

        # Assertion on average sum of inhibitory + excitatory neuron densities
        exci_inhib_sum = VoxelData.load_nrrd(excitatory_path_list[0]).raw
        print("\nadding initialization 1/" + str(len(excitatory_path_list)))
        for i in range (1, len(excitatory_path_list)):
            print("+ adding file " + str(i+1) + "/" + str(len(excitatory_path_list)))
            exci_inhib = VoxelData.load_nrrd(excitatory_path_list[i]).raw
            exci_inhib_sum += exci_inhib
        validate_average_density("sum of excitatory + generic inhibitory ME-type neuron densities",
            exci_inhib_sum, whole_brain_annotation_nb_voxels,
            neuron_dens_lit - neuron_dens_tolerance,
            neuron_dens_lit + neuron_dens_tolerance, z_error=2,error_fatal=error_fatal)


    if args.inhibitory_ME_types_folder is not None and args.excitatory_ME_types_folder is not None:

        # Assertion on average inhibitory sum density
        tot_inhib_sum = inhib_sum + generic_inhibitory
        validate_average_density("sum of inhibitory ME-type neuron densities",
            tot_inhib_sum, whole_brain_annotation_nb_voxels,
            inhibitory_neuron_dens_lit - inhibitory_neuron_dens_tolerance,
            inhibitory_neuron_dens_lit + inhibitory_neuron_dens_tolerance, None,error_fatal=error_fatal)

        # Assertion on average excitatory neuron density
        print("\nAssertion on average excitatory neuron density")
        print("/!\ Literature data not available")


    if args.excitatory_ME_types_transplant_folder is not None:

        # Reading the input files
        excitatory_transplant_path_list = [os.path.join(args.excitatory_ME_types_transplant_folder,f) for f in os.listdir(args.excitatory_ME_types_transplant_folder) if f.endswith(".nrrd")]

        # Assertion on average sum of inhibitory + excitatory neuron density after transplant
        exci_inhib_transplant_sum = VoxelData.load_nrrd(excitatory_transplant_path_list[0]).raw
        print("\nadding initialization 1/" + str(len(excitatory_transplant_path_list)))
        for i in range (1, len(excitatory_transplant_path_list)):
            print("+ adding file " + str(i+1) + "/" + str(len(excitatory_transplant_path_list)))
            exci_inhib_transplant = VoxelData.load_nrrd(excitatory_transplant_path_list[i]).raw
            exci_inhib_transplant_sum += exci_inhib_transplant
        validate_average_density("sum of excitatory + generic inhibitory ME-type neuron densities after transplant",
            exci_inhib_transplant_sum, whole_brain_annotation_nb_voxels,
            neuron_dens_lit - neuron_dens_tolerance,
            neuron_dens_lit + neuron_dens_tolerance, 1,error_fatal=error_fatal)


    if args.inhibitory_ME_types_transplant_folder is not None:
        # Reading the input files
        inhibitory_transplant_path_list = [os.path.join(args.inhibitory_ME_types_transplant_folder,f) for f in os.listdir(args.inhibitory_ME_types_transplant_folder) if f.endswith(".nrrd")]

        # Assertion on average sum of inhibitory ME-type neuron density which should be inferior or equal to the average inhibitory neurons after transplant
        inhib_transplant_sum = VoxelData.load_nrrd(inhibitory_transplant_path_list[0]).raw
        print("\nadding initialization 1/" + str(len(inhibitory_transplant_path_list)))
        for i in range (1, len(inhibitory_transplant_path_list)):
            print("+ adding file " + str(i+1) + "/" + str(len(inhibitory_transplant_path_list)))
            inhib_transplant = VoxelData.load_nrrd(inhibitory_transplant_path_list[i]).raw
            inhib_transplant_sum += inhib_transplant
        validate_average_density("sum of inhibitory ME-type neuron density after transplant which should be inferior or equal to the average inhibitory neurons",
            inhib_transplant_sum, whole_brain_annotation_nb_voxels,
            inhibitory_neuron_dens_lit - inhibitory_neuron_dens_lit,
            inhibitory_neuron_dens_lit + inhibitory_neuron_dens_tolerance, 1,error_fatal=error_fatal)


    if args.hierarchy_json:

        # 2.2 Assertion on sub regions densities
        print("\n\n==================================================")
        print("\nAssertion on sub-region densities...")

        # Region filters
        print("\nRegion filters setting")
        region_map = RegionMap.load_json(args.hierarchy_json)

        # Brain region filters

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
        Field_CA1 = (
            region_map.find("Field CA1", attr="name", with_descendants=True)
        )
        Field_CA2 = (
            region_map.find("Field CA2", attr="name", with_descendants=True)
        )
        Field_CA3 = (
            region_map.find("Field CA3", attr="name", with_descendants=True)
        )

        # barrels regions identification
        pattern = r"^SSp-bfd-[A-Za-z0-9]+$"
        region_ids = np.array(list(region_map.find("SSp-bfd", attr="acronym", with_descendants=True)))
        barrel_hierarchy_names = []
        for rid in region_ids:
            racronym = region_map.get(rid, attr="acronym")
            if re.match(pattern, racronym):
                barrel_hierarchy_names.append(racronym)
        children_barrel_name_list = []
        barrel_hierarchy_names = sorted(barrel_hierarchy_names)
        # for i in range(len(barrel_hierarchy_names)):
        for acronym in barrel_hierarchy_names:
            ids = list(region_map.find(acronym, "acronym", with_descendants=True))
            name_list = []
            for id_ in ids:
                name = region_map.get(id_, "name")
                name_list.append(name)
            name_list.sort()
            children_barrel_name_list.append(name_list)
        flattened_children_barrel_name_list = flatten_list(children_barrel_name_list)

        # rest_ids = region_map.find("root", attr="name", with_descendants=True)
        # rest_ids -= cerebellum_group_ids | isocortex_group_ids

        if cerebellum == [] or \
        isocortex == [] or \
        fiber_tracts_ids == [] or \
        hippocampus == [] or \
        striatum == []:
            raise ValueError("ERROR: some region filters return empty sets")
        # Volumes and literature definition for subregions
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
        isocortex_neuron_dens_lit = 2*5048837/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_neuron_dens_tolerance = 2*412123/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_glia_dens_lit = 2*6640234/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_glia_dens_tolerance = 2*244643/(isocortex_nb_vox * voxel_volume) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_cell_dens_lit = (isocortex_neuron_dens_lit + isocortex_glia_dens_lit) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_cell_dens_tolerance = (isocortex_neuron_dens_tolerance + isocortex_glia_dens_tolerance) # Table 1 in Herculano_Houzel et al., 2013
        isocortex_oligo_dens_lit = 12500 # Table 1 in Erö et al., 2018
        isocortex_oligo_dens_tolerance = round(isocortex_oligo_dens_lit * default_glia_proportion) # default Not available /!\
        isocortex_astro_dens_lit = 15696 # Table 1 in Erö et al., 2018
        isocortex_astro_dens_tolerance = round(isocortex_astro_dens_lit * default_glia_proportion) # default Not available /!\
        isocortex_microglia_dens_lit = 6500 # Table 1 in Erö et al., 2018
        isocortex_microglia_dens_tolerance = round(isocortex_microglia_dens_lit * default_glia_proportion) # default Not available /!\


        if args.neuron_glia_density_folder is not None and args.cell_density is not None:

            # Assertion on fiber tracks + grooves + ventricular_system where no neuron should be found
            fiber_tracts_neuron = neuron[np.isin(annotation, list(fiber_tracts_ids))]
            fiber_tracts_neuron_sum = np.sum(fiber_tracts_neuron) * voxel_volume
            diff_fiber_tracts_neuron_sum = abs(neuron_dens_fiber_tracts_lit - fiber_tracts_neuron_sum)
            print("\nAssertion on fiber tracks + grooves + ventricular_system where no neuron should be found")
            if not diff_fiber_tracts_neuron_sum <= neuron_dens_fiber_tracts_tolerance:
                warning_message = "fiber tracks + grooves + ventricular_system where no neuron should be found not consistent with literature"
                warnings.warn(warning_message, UserWarning)
                print("WARNING:", warning_message)
            else:
                print("Validated")

            # ---------------------------------------------------------------------------------------------
            # ISOCORTEX
            print("\n\n----------ISOCORTEX----------")

            # Assertion on isocortex cell density
            isocortex_cell_dens = cell[np.isin(annotation, list(isocortex))]
            isocortex_cell_dens_sum = np.sum(isocortex_cell_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex cell density (/mm^3)")
            print_range_bar(isocortex_cell_dens_sum, isocortex_cell_dens_lit - isocortex_cell_dens_tolerance, isocortex_cell_dens_lit + isocortex_cell_dens_tolerance)
            assertion_message = "Average isocortex cell density out of literature range"
            z_score_assertion(isocortex_cell_dens_sum, isocortex_cell_dens_lit - isocortex_cell_dens_tolerance, isocortex_cell_dens_lit + isocortex_cell_dens_tolerance, assertion_message,error_fatal)

            # Assertion on isocortex neuron density
            isocortex_neuron_dens = neuron[np.isin(annotation, list(isocortex))]
            isocortex_neuron_dens_sum = np.sum(isocortex_neuron_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex neuron density (/mm^3)")
            print("/!\ Tolerance set to default for neuron density")
            isocortex_neuron_dens_default_tolerance = isocortex_neuron_dens_lit * default_neuron_proportion
            print_range_bar(isocortex_neuron_dens_sum, isocortex_neuron_dens_lit - isocortex_neuron_dens_default_tolerance, isocortex_neuron_dens_lit + isocortex_neuron_dens_default_tolerance)
            assertion_message = "Average isocortex neuron density out of literature range"
            z_score_assertion(isocortex_neuron_dens_sum, isocortex_neuron_dens_lit - isocortex_neuron_dens_default_tolerance, isocortex_neuron_dens_lit + isocortex_neuron_dens_default_tolerance, assertion_message,error_fatal)

            # Assertion on isocortex glia density
            isocortex_glia_dens = glia[np.isin(annotation, list(isocortex))]
            isocortex_glia_dens_sum = np.sum(isocortex_glia_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex glia density (/mm^3)")
            print("/!\ Tolerance set to default for glia density")
            isocortex_glia_dens_default_tolerance = isocortex_glia_dens_lit * default_glia_proportion
            print_range_bar(isocortex_glia_dens_sum, isocortex_glia_dens_lit - isocortex_glia_dens_default_tolerance, isocortex_glia_dens_lit + isocortex_glia_dens_default_tolerance)
            assertion_message = "Average isocortex glia density out of literature range"
            z_score_assertion(isocortex_glia_dens_sum, isocortex_glia_dens_lit - isocortex_glia_dens_default_tolerance, isocortex_glia_dens_lit + isocortex_glia_dens_default_tolerance, assertion_message,error_fatal)

            # Assertion on isocortex oligodendrocyte density
            isocortex_oligo_dens = oligodendrocyte[np.isin(annotation, list(isocortex))]
            isocortex_oligo_dens_sum = np.sum(isocortex_oligo_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex oligodendrocyte density (/mm^3)")
            print("/!\ Literature figures not consistent + Tolerance not available, set by default")
            print_range_bar(isocortex_oligo_dens_sum, isocortex_oligo_dens_lit - isocortex_oligo_dens_tolerance, isocortex_oligo_dens_lit + isocortex_oligo_dens_tolerance)
            assertion_message = "Average isocortex oligodendrocyte density out of literature range"
            z_score_assertion(isocortex_oligo_dens_sum, isocortex_oligo_dens_lit - isocortex_oligo_dens_tolerance, isocortex_oligo_dens_lit + isocortex_oligo_dens_tolerance, assertion_message,error_fatal)

            # Assertion on isocortex astrocyte density
            isocortex_astro_dens = astrocyte[np.isin(annotation, list(isocortex))]
            isocortex_astro_dens_sum = np.sum(isocortex_astro_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex astrocyte density (/mm^3)")
            print("/!\ Literature figures not consistent + Tolerance not available, set by default")
            print_range_bar(isocortex_astro_dens_sum, isocortex_astro_dens_lit - isocortex_astro_dens_tolerance, isocortex_astro_dens_lit + isocortex_astro_dens_tolerance)
            assertion_message = "Average isocortex astrocyte density out of literature range"
            z_score_assertion(isocortex_astro_dens_sum, isocortex_astro_dens_lit - isocortex_astro_dens_tolerance, isocortex_astro_dens_lit + isocortex_astro_dens_tolerance, assertion_message,error_fatal)

            # Assertion on isocortex microglia density
            isocortex_microglia_dens = microglia[np.isin(annotation, list(isocortex))]
            isocortex_microglia_dens_sum = np.sum(isocortex_microglia_dens) / isocortex_nb_vox # * voxel_volume
            print("\nAssertion on isocortex microglia density (/mm^3)")
            print("/!\ Literature figures not consistent + Tolerance not available, set by default")
            print_range_bar(isocortex_microglia_dens_sum, isocortex_microglia_dens_lit - isocortex_microglia_dens_tolerance, isocortex_microglia_dens_lit + isocortex_microglia_dens_tolerance)
            assertion_message = "Average isocortex microglia density out of literature range"
            z_score_assertion(isocortex_microglia_dens_sum, isocortex_microglia_dens_lit - isocortex_microglia_dens_tolerance, isocortex_microglia_dens_lit + isocortex_microglia_dens_tolerance, assertion_message,error_fatal)

            # Assertion on barrel inhibitory neuron densities
            assert_densities(annotation, flattened_children_barrel_name_list, 'inhibitory neuron density',gad,error_fatal)
            
            # Assertion on barrel excitatory neuron densities (except layer 1)
            filtered_list = [item for item in flattened_children_barrel_name_list if "layer 1" not in item]
            assert_densities(annotation, filtered_list,
                           'excitatory neuron density', gad,error_fatal)

            # ---------------------------------------------------------------------------------------------


            # ---------------------------------------------------------------------------------------------
            # CEREBELLUM
            print("\n\n----------CEREBELLUM----------")

            # Assertion on cerebellum cell density
            cerebellum_nb_vox = len(np.where(np.isin(annotation, list(cerebellum)) != 0)[0])
            cerebellum_cell_dens = cell[np.isin(annotation, list(cerebellum))]
            cerebellum_cell_dens_sum = np.sum(cerebellum_cell_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum cell density (/mm^3)")
            print_range_bar(cerebellum_cell_dens_sum, cerebellum_cell_dens_lit - cerebellum_cell_dens_tolerance, cerebellum_cell_dens_lit + cerebellum_cell_dens_tolerance)
            assertion_message = "Average cerebellum cell density out of literature range"
            z_score_assertion(cerebellum_cell_dens_sum, cerebellum_cell_dens_lit - cerebellum_cell_dens_tolerance, cerebellum_cell_dens_lit + cerebellum_cell_dens_tolerance, assertion_message,error_fatal)

            # Assertion on cerebellum neuron density
            cerebellum_neuron_dens = neuron[np.isin(annotation, list(cerebellum))]
            cerebellum_neuron_dens_sum = np.sum(cerebellum_neuron_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum neuron density (/mm^3)")
            print_range_bar(cerebellum_neuron_dens_sum, cerebellum_neuron_dens_lit - cerebellum_neuron_dens_tolerance, cerebellum_neuron_dens_lit + cerebellum_neuron_dens_tolerance)
            assertion_message = "Average cerebellum neuron density out of literature range"
            z_score_assertion(cerebellum_neuron_dens_sum, cerebellum_neuron_dens_lit - cerebellum_neuron_dens_tolerance, cerebellum_neuron_dens_lit + cerebellum_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on cerebellum glia density
            cerebellum_glia_dens = glia[np.isin(annotation, list(cerebellum))]
            cerebellum_glia_dens_sum = np.sum(cerebellum_glia_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum glia density (/mm^3)")
            print_range_bar(cerebellum_glia_dens_sum, cerebellum_glia_dens_lit - cerebellum_glia_dens_tolerance, cerebellum_glia_dens_lit + cerebellum_glia_dens_tolerance)
            assertion_message = "Average cerebellum glia density out of literature range"
            z_score_assertion(cerebellum_glia_dens_sum, cerebellum_glia_dens_lit - cerebellum_glia_dens_tolerance, cerebellum_glia_dens_lit + cerebellum_glia_dens_tolerance, assertion_message,error_fatal)

            # Assertion on cerebellum oligodendrocyte density
            cerebellum_oligo_dens = oligodendrocyte[np.isin(annotation, list(cerebellum))]
            cerebellum_oligo_dens_sum = np.sum(cerebellum_oligo_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum oligodendrocyte density (/mm^3)")
            print("/!\ Literature figures not consistent")
            print_range_bar(cerebellum_oligo_dens_sum, cerebellum_oligo_dens_lit - cerebellum_oligo_dens_tolerance, cerebellum_oligo_dens_lit + cerebellum_oligo_dens_tolerance)
            assertion_message = "Average cerebellum oligodendrocyte density out of literature range"
            z_score_assertion(cerebellum_oligo_dens_sum, cerebellum_oligo_dens_lit - cerebellum_oligo_dens_tolerance, cerebellum_oligo_dens_lit + cerebellum_oligo_dens_tolerance, assertion_message,error_fatal)

            # Assertion on cerebellum astrocyte density
            cerebellum_astro_dens = astrocyte[np.isin(annotation, list(cerebellum))]
            cerebellum_astro_dens_sum = np.sum(cerebellum_astro_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum astrocyte density (/mm^3)")
            print("/!\ Literature figures not consistent + Tolerance not available, set by default")
            print_range_bar(cerebellum_astro_dens_sum, cerebellum_astro_dens_lit - cerebellum_astro_dens_tolerance, cerebellum_astro_dens_lit + cerebellum_astro_dens_tolerance)
            assertion_message = "Average cerebellum astrocyte density out of literature range"
            z_score_assertion(cerebellum_astro_dens_sum, cerebellum_astro_dens_lit - cerebellum_astro_dens_tolerance, cerebellum_astro_dens_lit + cerebellum_astro_dens_tolerance, assertion_message,error_fatal)

            # Assertion on cerebellum microglia density
            cerebellum_microglia_dens = microglia[np.isin(annotation, list(cerebellum))]
            cerebellum_microglia_dens_sum = np.sum(cerebellum_microglia_dens) / cerebellum_nb_vox # * voxel_volume
            print("\nAssertion on cerebellum microglia density (/mm^3)")
            print("/!\ Literature figures not consistent")
            print_range_bar(cerebellum_microglia_dens_sum, cerebellum_microglia_dens_lit - cerebellum_microglia_dens_tolerance, cerebellum_microglia_dens_lit + cerebellum_microglia_dens_tolerance)
            assertion_message = "Average cerebellum microglia density out of literature range"
            z_score_assertion(cerebellum_microglia_dens_sum, cerebellum_microglia_dens_lit - cerebellum_microglia_dens_tolerance, cerebellum_microglia_dens_lit + cerebellum_microglia_dens_tolerance, assertion_message,error_fatal)

            # ---------------------------------------------------------------------------------------------



            # ---------------------------------------------------------------------------------------------
            # STRIATUM

            print("\n\n----------STRIATUM----------")

            # Assertion on striatum neuron density
            striatum_neuron_dens = neuron[np.isin(annotation, list(striatum))]
            striatum_neuron_dens_sum = np.sum(striatum_neuron_dens) / striatum_nb_vox # * voxel_volume
            print("\nAssertion on striatum neuron density (/mm^3)")
            print_range_bar(striatum_neuron_dens_sum, striatum_neuron_dens_lit - striatum_neuron_dens_tolerance, striatum_neuron_dens_lit + striatum_neuron_dens_tolerance)
            assertion_message = "Average striatum neuron density out of literature range"
            z_score_assertion(striatum_neuron_dens_sum, striatum_neuron_dens_lit - striatum_neuron_dens_tolerance, striatum_neuron_dens_lit + striatum_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on striatum oligodendrocyte density
            striatum_oligo_dens = oligodendrocyte[np.isin(annotation, list(striatum))]
            striatum_oligo_dens_sum = np.sum(striatum_oligo_dens) / striatum_nb_vox # * voxel_volume
            print("\nAssertion on striatum oligodendrocyte density (/mm^3)")
            print_range_bar(striatum_oligo_dens_sum, striatum_oligo_dens_lit - striatum_oligo_dens_tolerance, striatum_oligo_dens_lit + striatum_oligo_dens_tolerance)
            assertion_message = "Average striatum oligodendrocyte density out of literature range"
            z_score_assertion(striatum_oligo_dens_sum, striatum_oligo_dens_lit - striatum_oligo_dens_tolerance, striatum_oligo_dens_lit + striatum_oligo_dens_tolerance, assertion_message,error_fatal)

            # Assertion on striatum astrocyte density
            striatum_astro_dens = astrocyte[np.isin(annotation, list(striatum))]
            striatum_astro_dens_sum = np.sum(striatum_astro_dens) / striatum_nb_vox # * voxel_volume
            print("\nAssertion on striatum astrocyte density (/mm^3)")
            print_range_bar(striatum_astro_dens_sum, striatum_astro_dens_lit - striatum_astro_dens_tolerance, striatum_astro_dens_lit + striatum_astro_dens_tolerance)
            assertion_message = "Average striatum astrocyte density out of literature range"
            z_score_assertion(striatum_astro_dens_sum, striatum_astro_dens_lit - striatum_astro_dens_tolerance, striatum_astro_dens_lit + striatum_astro_dens_tolerance, assertion_message,error_fatal)

            # Assertion on striatum microglia density
            striatum_microglia_dens = microglia[np.isin(annotation, list(striatum))]
            striatum_microglia_dens_sum = np.sum(striatum_microglia_dens) / striatum_nb_vox # * voxel_volume
            print("\nAssertion on striatum microglia density (/mm^3)")
            print_range_bar(striatum_microglia_dens_sum, striatum_microglia_dens_lit - striatum_microglia_dens_tolerance, striatum_microglia_dens_lit + striatum_microglia_dens_tolerance)
            assertion_message = "Average striatum microglia density out of literature range"
            z_score_assertion(striatum_microglia_dens_sum, striatum_microglia_dens_lit - striatum_microglia_dens_tolerance, striatum_microglia_dens_lit + striatum_microglia_dens_tolerance, assertion_message,error_fatal)

            # ---------------------------------------------------------------------------------------------



            # ---------------------------------------------------------------------------------------------
            # HIPPOCAMPUS

            print("\n\n----------HIPPOCAMPUS----------")

            # Assertion on hippocampus neuron density
            hippocampus_neuron_dens = neuron[np.isin(annotation, list(hippocampal_formation))]
            hippocampus_neuron_dens_sum = np.sum(hippocampus_neuron_dens) / hippocampal_formation_nb_vox # * voxel_volume
            print("\nAssertion on hippocampus neuron density (/mm^3)")
            print("/!\ Tolerance not available, set by default")
            print_range_bar(hippocampus_neuron_dens_sum, hippocampus_neuron_dens_lit - hippocampus_neuron_dens_tolerance, hippocampus_neuron_dens_lit + hippocampus_neuron_dens_tolerance)
            assertion_message = "Average hippocampus neuron density out of literature range"
            z_score_assertion(hippocampus_neuron_dens_sum, hippocampus_neuron_dens_lit - hippocampus_neuron_dens_tolerance, hippocampus_neuron_dens_lit + hippocampus_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on hippocampus oligodendrocyte density
            hippocampus_oligo_dens = oligodendrocyte[np.isin(annotation, list(hippocampus))]
            hippocampus_oligo_dens_sum = np.sum(hippocampus_oligo_dens) / hippocampus_nb_vox # * voxel_volume
            print("\nAssertion on hippocampus oligodendrocyte density (/mm^3)")
            print("/!\ Tolerance not available, set by default")
            print_range_bar(hippocampus_oligo_dens_sum, hippocampus_oligo_dens_lit - hippocampus_oligo_dens_tolerance, hippocampus_oligo_dens_lit + hippocampus_oligo_dens_tolerance)
            assertion_message = "Average hippocampus oligodendrocyte density out of literature range"
            z_score_assertion(hippocampus_oligo_dens_sum, hippocampus_oligo_dens_lit - hippocampus_oligo_dens_tolerance, hippocampus_oligo_dens_lit + hippocampus_oligo_dens_tolerance, assertion_message,error_fatal)

            # Assertion on hippocampus astrocyte density
            hippocampus_astro_dens = astrocyte[np.isin(annotation, list(hippocampus))]
            hippocampus_astro_dens_sum = np.sum(hippocampus_astro_dens) / hippocampus_nb_vox # * voxel_volume
            print("\nAssertion on hippocampus astrocyte density (/mm^3)")
            print_range_bar(hippocampus_astro_dens_sum, hippocampus_astro_dens_lit - hippocampus_astro_dens_tolerance, hippocampus_astro_dens_lit + hippocampus_astro_dens_tolerance)
            assertion_message = "Average hippocampus astrocyte density out of literature range"
            z_score_assertion(hippocampus_astro_dens_sum, hippocampus_astro_dens_lit - hippocampus_astro_dens_tolerance, hippocampus_astro_dens_lit + hippocampus_astro_dens_tolerance, assertion_message,error_fatal)

            # Assertion on hippocampus microglia density
            hippocampus_microglia_dens = microglia[np.isin(annotation, list(hippocampus))]
            hippocampus_microglia_dens_sum = np.sum(hippocampus_microglia_dens) / hippocampus_nb_vox # * voxel_volume
            print("\nAssertion on hippocampus microglia density (/mm^3)")
            print_range_bar(hippocampus_microglia_dens_sum, hippocampus_microglia_dens_lit - hippocampus_microglia_dens_tolerance, hippocampus_microglia_dens_lit + hippocampus_microglia_dens_tolerance)
            assertion_message = "Average hippocampus microglia density out of literature range"
            z_score_assertion(hippocampus_microglia_dens_sum, hippocampus_microglia_dens_lit - hippocampus_microglia_dens_tolerance, hippocampus_microglia_dens_lit + hippocampus_microglia_dens_tolerance, assertion_message,error_fatal)

            # Assertions on inhibitory and excitatory neuron for Field CA1/2/3
            field_CAn = {'Field CA1': Field_CA1,
                         'Field CA2': Field_CA2,
                         'Field CA3': Field_CA3,}
            for label, region in field_CAn.items():
                # Assertion on inhibitory neuron density for Field CAn
                assert_density(annotation, region, 'inhibitory neuron density',
                               gad,error_fatal, region_label=label)
                # Assertion on excitatory neuron density for Field CAn
                assert_density(annotation, region, 'excitatory neuron density',
                               gad,error_fatal, neuron=neuron, region_label=label)

            # ---------------------------------------------------------------------------------------------
            # THALAMUS

            print("\n\n----------THALAMUS----------")

            # Assertion on Talamus cell density
            thalamus_cell_dens = cell[np.isin(annotation, list(thalamus))]
            thalamus_cell_dens_sum = np.sum(thalamus_cell_dens) / thalamus_nb_vox # * voxel_volume
            print("\nAssertion on Thalamus cell density (/mm^3)")
            print("/!\ Tolerance set to default for cell density")
            print_range_bar(thalamus_cell_dens_sum, thalamus_cell_dens_lit - thalamus_cell_dens_tolerance, thalamus_cell_dens_lit + thalamus_cell_dens_tolerance)
            assertion_message = "Average thalamus cell density out of literature range"
            z_score_assertion(thalamus_cell_dens_sum, thalamus_cell_dens_lit - thalamus_cell_dens_tolerance, thalamus_cell_dens_lit + thalamus_cell_dens_tolerance, assertion_message,error_fatal)

            # Assertion on Thalamus glia density
            thalamus_glia_dens = glia[np.isin(annotation, list(thalamus))]
            thalamus_glia_dens_sum = np.sum(thalamus_glia_dens) / thalamus_nb_vox # * voxel_volume
            print("\nAssertion on Thalamus glia density (/mm^3)")
            print("/!\ Tolerance set to default for glia density")
            print_range_bar(thalamus_glia_dens_sum, thalamus_glia_dens_lit - thalamus_glia_dens_tolerance, thalamus_glia_dens_lit + thalamus_glia_dens_tolerance)
            assertion_message = "Average thalamus glia density out of literature range"
            z_score_assertion(thalamus_glia_dens_sum, thalamus_glia_dens_lit - thalamus_glia_dens_tolerance, thalamus_glia_dens_lit + thalamus_glia_dens_tolerance, assertion_message,error_fatal)

            # Assertion on LGd neuron density
            LGd_neuron_dens = neuron[np.isin(annotation, list(LGd))]
            LGd_neuron_dens_sum = np.sum(LGd_neuron_dens) / LGd_nb_vox # * voxel_volume
            print("\nAssertion on LGd neuron density (/mm^3)")
            print_range_bar(LGd_neuron_dens_sum, LGd_neuron_dens_lit - LGd_neuron_dens_tolerance, LGd_neuron_dens_lit + LGd_neuron_dens_tolerance)
            assertion_message = "Average LGd neuron density out of literature range"
            z_score_assertion(LGd_neuron_dens_sum, LGd_neuron_dens_lit - LGd_neuron_dens_tolerance, LGd_neuron_dens_lit + LGd_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on VPM neuron density
            VPM_neuron_dens = neuron[np.isin(annotation, list(VPM))]
            VPM_neuron_dens_sum = np.sum(VPM_neuron_dens) / VPM_nb_vox # * voxel_volume
            print("\nAssertion on VPM neuron density (/mm^3)")
            print_range_bar(VPM_neuron_dens_sum, VPL_neuron_dens_lit - VPM_neuron_dens_tolerance, VPM_neuron_dens_lit + VPM_neuron_dens_tolerance)
            assertion_message = "Average VPM neuron density out of literature range"
            z_score_assertion(VPM_neuron_dens_sum, VPM_neuron_dens_lit - VPM_neuron_dens_tolerance, VPM_neuron_dens_lit + VPM_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on VPL neuron density
            VPL_neuron_dens = neuron[np.isin(annotation, list(VPL))]
            VPL_neuron_dens_sum = np.sum(VPL_neuron_dens) / VPL_nb_vox # * voxel_volume
            print("\nAssertion on VPL neuron density (/mm^3)")
            print("/!\ Tolerance set to default for neuron density")
            print_range_bar(VPL_neuron_dens_sum, VPL_neuron_dens_lit - VPL_neuron_dens_tolerance, VPL_neuron_dens_lit + VPL_neuron_dens_tolerance)
            assertion_message = "Average VPL neuron density out of literature range"
            z_score_assertion(VPL_neuron_dens_sum, VPL_neuron_dens_lit - VPL_neuron_dens_tolerance, VPL_neuron_dens_lit + VPL_neuron_dens_tolerance, assertion_message,error_fatal)

        if args.inhibitory_density_folder is not None:

            # Assertion on VPL pv density
            VPL_pv_dens = pv[np.isin(annotation, list(VPL))]
            VPL_pv_dens_sum = np.sum(VPL_pv_dens) / VPL_nb_vox # * voxel_volume
            print("\nAssertion on VPL pv density (/mm^3)")
            print_range_bar(VPL_pv_dens_sum, VPL_pv_dens_lit - VPL_pv_dens_tolerance, VPL_pv_dens_lit + VPL_pv_dens_tolerance)
            assertion_message = "Average VPL pv density out of literature range"
            z_score_assertion(VPL_pv_dens_sum, VPL_pv_dens_lit - VPL_pv_dens_tolerance, VPL_pv_dens_lit + VPL_pv_dens_tolerance, assertion_message,error_fatal)

            # Assertion on VPL sst density
            VPL_sst_dens = sst[np.isin(annotation, list(VPL))]
            VPL_sst_dens_sum = np.sum(VPL_sst_dens) / VPL_nb_vox # * voxel_volume
            print("\nAssertion on VPL sst density (/mm^3)")
            print_range_bar(VPL_sst_dens_sum, VPL_sst_dens_lit - VPL_sst_dens_tolerance, VPL_sst_dens_lit + VPL_sst_dens_tolerance)
            assertion_message = "Average VPL sst density out of literature range"
            z_score_assertion(VPL_sst_dens_sum, VPL_sst_dens_lit - VPL_sst_dens_tolerance, VPL_sst_dens_lit + VPL_sst_dens_tolerance, assertion_message,error_fatal)

            # Assertion on VPL vip density
            VPL_vip_dens = vip[np.isin(annotation, list(VPL))]
            VPL_vip_dens_sum = np.sum(VPL_vip_dens) / VPL_nb_vox # * voxel_volume
            print("\nAssertion on VPL vip density (/mm^3)")
            print_range_bar(VPL_vip_dens_sum, VPL_vip_dens_lit - VPL_vip_dens_tolerance, VPL_vip_dens_lit + VPL_vip_dens_tolerance)
            assertion_message = "Average VPL vip density out of literature range"
            z_score_assertion(VPL_vip_dens_sum, VPL_vip_dens_lit - VPL_vip_dens_tolerance, VPL_vip_dens_lit + VPL_vip_dens_tolerance, assertion_message,error_fatal)

            # ---------------------------------------------------------------------------------------------


            # ---------------------------------------------------------------------------------------------
            # MAIN OLFACTORY BULB

            print("\n\n----------MAIN OLFACTORY BULB----------")

            # Assertion on MOB cell density
            MOB_cell_dens = cell[np.isin(annotation, list(MOB))]
            MOB_cell_dens_sum = np.sum(MOB_cell_dens) / MOB_nb_vox # * voxel_volume
            print("\nAssertion on MOB cell density (/mm^3)")
            print_range_bar(MOB_cell_dens_sum, MOB_cell_dens_lit - MOB_cell_dens_tolerance, MOB_cell_dens_lit + MOB_cell_dens_tolerance)
            assertion_message = "Average MOB_cell density out of literature range"
            z_score_assertion(MOB_cell_dens_sum, MOB_cell_dens_lit - MOB_cell_dens_tolerance, MOB_cell_dens_lit + MOB_cell_dens_tolerance, assertion_message,error_fatal)

            # Assertion on MOB neuron density
            MOB_neuron_dens = neuron[np.isin(annotation, list(MOB))]
            MOB_neuron_dens_sum = np.sum(MOB_neuron_dens) / MOB_nb_vox # * voxel_volume
            print("\nAssertion on MOB neuron density (/mm^3)")
            print_range_bar(MOB_neuron_dens_sum, MOB_neuron_dens_lit - MOB_neuron_dens_tolerance, MOB_neuron_dens_lit + MOB_neuron_dens_tolerance)
            assertion_message = "Average MOB_neuron density out of literature range"
            z_score_assertion(MOB_neuron_dens_sum, MOB_neuron_dens_lit - MOB_neuron_dens_tolerance, MOB_neuron_dens_lit + MOB_neuron_dens_tolerance, assertion_message,error_fatal)

            # Assertion on MOB glia density
            MOB_glia_dens = glia[np.isin(annotation, list(MOB))]
            MOB_glia_dens_sum = np.sum(MOB_glia_dens) / MOB_nb_vox # * voxel_volume
            print("\nAssertion on MOB glia density (/mm^3)")
            print_range_bar(MOB_glia_dens_sum, MOB_glia_dens_lit - MOB_glia_dens_tolerance, MOB_glia_dens_lit + MOB_glia_dens_tolerance)
            assertion_message = "Average MOB_glia density out of literature range"
            z_score_assertion(MOB_glia_dens_sum, MOB_glia_dens_lit - MOB_glia_dens_tolerance, MOB_glia_dens_lit + MOB_glia_dens_tolerance, assertion_message,error_fatal)


    print("\n==================================")
    print("All assertions successfully tested")

    return

# =============================================================================================
