# What is OSP?
The OSP (Oxidation State Probability) filter is a simple tool designed to predict the stability of compounds through the analysis of elemental oxidation states. The process begins with the selection of a ternary system for examination. The Materials Project database is then consulted to identify compounds composed of the chosen elements, which are subsequently verified against the Inorganic Crystal Structure Database to confirm their authenticity.

The OSP filtering method calculates the likelihood of each element adopting various oxidation states within these compounds by analyzing the frequency of known oxidation states across the authenticated compounds. A fundamental criterion for the OSP filter is that each compound must exhibit a net oxidation state of zero to be considered stable. This is determined by the oxidation states of all elements in the compound, factoring in their quantities. The Pauling electronegativity values guide the assignment of oxidation states, ensuring that the most electronegative element is assigned negative states, while the least electronegative is given positive states.

By focusing exclusively on compounds that satisfy the charge neutrality condition, the OSP filtering method streamlines the analysis, making it an invaluable resource for researchers in the material design field. It aids in the efficient identification of viable compounds for further exploration or synthesis, serving as a predictive tool in the material design process.

# How to use the filter
Compitable with python3.10. Install all required packages:

```bash
pip install -r requirement.txt
```

# How to use the OSP filter on a tenary system
Update the Section 1 (define variables) in element_focused_CS_analyzer.py and run the python file.

## database 
df = pd.read_hdf("all_known_materials.h5")
#You may query it through Material Project. The format of the dataframe is as it is downloaded from Material Project with the following essential column names: {'formation_energy_per_atom', 'unit_cell_formula', 'pretty_formula', 'lements', 'nelements', 'e_above_hull', 'icsd_ids'}
        
## maximum number of atoms in a compound to be generated
max_atom_in_compound = 10

## threshold to ignore low frequency charge states
cut_off_freq_percentage = 0.1 

## the tenary system that you are interested 
elements_2find_original = ['Cu','V', 'Nb']

## file generation 
After running the element_focused_CS_analyzer.py until the last Section 14, a csv file that include the OSP will be generated. 

## typical Run Time
Typically for elements that has been in the database, the python file will finish running in less than a minute. 

## updating the database
If the CP base (charge_probability_database.csv) is outdated, manually delete the entries in the csv file and run "element_focused_CS_analyzer.py". Query new compounds through material project or from a similar source with the datastructure specified above. The codes in Section 1 to 6 will try to analyze the new database file (dataframe) and update the csv file with a new charge probabilty record. 
