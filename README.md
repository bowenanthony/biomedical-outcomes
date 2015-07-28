# biomedical-outcomes

This repository contains raw data and Perl code used in the 2015 manuscript by Bowen and Casadevall entitled "Increasing disparities between resource inputs and outcomes, as measured by certain health deliverables, in biomedical research." Publication is currently in press at PNAS.

All PubMed data gathered in the manuscript were obtained using "pubmed_scriptV1-7.pl." The output files from this script that were used in the manuscript are in the folder "raw_pubmed_data." In this folder, you will find two sub-directories corresponding to PubMed data with a "USA" affiliation and corresponding to PubMed data with no affiliation restriction. The affiliation restriction constraint can be modified in the Perl script. Output files consist of a raw data file for each analyzed year ("YEAR_pubmed_authors.tab"), a single file showing certain data fields for all analyzed years ("ss_data.tab"), and another file showing the distribution of authors per publication for each year ("app_data.tab"). The "data_modes.tab" file is only used to compare script output from CSV and Medline sources when both data types are downloaded from PubMed. A "Perl_Temp" folder is needed to store output files while running the Perl script.

Other raw data used in the manuscript is in the "other_raw_data" folder and contains files obtained from Eurostat (1), and The World Bank (2).

For questions about this code or manuscript, please contact Anthony Bowen at anthony.bowen@med.einstein.yu.edu.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18788.svg)](http://dx.doi.org/10.5281/zenodo.18788)

## References
1. European Commission Eurostat (October 29, 2014) Life expectancy by age and sex. Available at http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_mlexpec&lang=en.
2. The World Bank. World Bank open data. Available at http://data.worldbank.org.
