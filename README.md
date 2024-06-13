# What is OSP?
The OSP (Oxidation State Probability) filter is a simple tool designed to predict the stability of compounds through the analysis of elemental oxidation states. The process begins with the selection of a ternary system for examination. The Materials Project database is then consulted to identify compounds composed of the chosen elements, which are subsequently verified against the Inorganic Crystal Structure Database to confirm their authenticity.

The OSP filtering method calculates the likelihood of each element adopting various oxidation states within these compounds by analyzing the frequency of known oxidation states across the authenticated compounds. A fundamental criterion for the OSP filter is that each compound must exhibit a net oxidation state of zero to be considered stable. This is determined by the oxidation states of all elements in the compound, factoring in their quantities. The Pauling electronegativity values guide the assignment of oxidation states, ensuring that the most electronegative element is assigned negative states, while the least electronegative is given positive states.

By focusing exclusively on compounds that satisfy the charge neutrality condition, the OSP filtering method streamlines the analysis, making it an invaluable resource for researchers in the material design field. It aids in the efficient identification of viable compounds for further exploration or synthesis, serving as a predictive tool in the material design process.

# How to use the model
Install all required packages:

```bash
pip install -r requirement.txt
```
