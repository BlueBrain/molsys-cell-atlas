# Literature values setting

# Literature values and initialization of paramters
voxel_volume = (25 * 1.0e-3)**3 # mm^3 
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

# Literature values for the whole brain
wh_mouse_brain_vol_litt = 508.91 # = 508.91 mm^3 in Badea et al., 2007
wh_mouse_brain_vol_tolerance = 23.42 # = 23.42 mm^3 (5%) in Badea et al., 2007
wh_mouse_brain_vol_n_litt = int(wh_mouse_brain_vol_litt/voxel_volume) # = 508.91 mm^3 (count) in Badea et al., 2007
wh_mouse_brain_vol_n_tolerance = int(wh_mouse_brain_vol_tolerance/voxel_volume) # = 23.42 mm^3 (5%) (count) in Badea et al., 2007
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
