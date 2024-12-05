The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government’s ETH Board of the Swiss Federal Institutes of Technology.

Copyright (c) 2024 Blue Brain Project/EPFL

Example of command to be run:

python validation.py \
    /home/piluso/data/00_allen_brain_atlas/ccfv2/annotation_25_ccfv2_combined.nrrd \
    /home/piluso/data/00_allen_brain_atlas/1.json \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/cell_densities_correctednissl/  \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/overall_cell_density_correctednissl.nrrd \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/inhibitory_neuron_densities_linprog_correctednissl/ \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/excitatory_split \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/mtypes_densities_probability_map \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/excitatory_split_transplant/ \
    /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/mtypes_densities_probability_map_transplant/


pip install setup.py

densities_validation \
    --annotation /home/piluso/data/00_allen_brain_atlas/ccfv2/annotation_25_ccfv2_combined.nrrd \
    --hierarchy /home/piluso/data/00_allen_brain_atlas/1.json \
    --density_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/cell_densities_correctednissl/ \
    --cell_density /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/overall_cell_density_correctednissl.nrrd \
    --inhibitory_density_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/inhibitory_neuron_densities_linprog_correctednissl/ \
    --excitatory_ME_types_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/excitatory_split \
    --inhibitory_ME_types_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/mtypes_densities_probability_map \
    --excitatory_ME_types_transplant_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/excitatory_split_transplant/ \
    --inhibitory_ME_types_transplant_folder /gpfs/bbp.cscs.ch/home/lcristel/fix_leaves_only_rootLabel/mtypes_densities_probability_map_transplant/
