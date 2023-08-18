#!/bin/bash
set -eo pipefail

cd notebooks
wget --quiet -O escapecalculator.py https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/escapecalculator.py
wget --quiet -O cleaned_Kds_RBD_ACE2_withx.tsv https://github.com/desai-lab/compensatory_epistasis_omicron/raw/main/Final_Figures/data/cleaned_Kds_RBD_ACE2_withx.tsv
