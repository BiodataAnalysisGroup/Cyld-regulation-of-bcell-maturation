#!/bin/sh
#
#

#sc CYLD

#CYLD
python src/scFEA.py --data_dir data --input_dir input/Our_data/case/ --test_file Data_CYLD.csv --moduleGene_file module_gene_complete_mouse_m168.csv --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv --sc_imputation True --output_flux_file output/mouse_flux.csv --output_balance_file output/mouse_balance.csv

#control
python src/scFEA.py --data_dir data --input_dir input/Our_data/control/ --test_file Data.csv --moduleGene_file module_gene_complete_mouse_m168.csv --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv --sc_imputation True --output_flux_file output/mouse_control_flux.csv --output_balance_file output/mouse_control_balance.csv
