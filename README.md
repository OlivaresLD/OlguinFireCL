# Raw data, processing scripts, and references for the 2011 Olguín Fire (Torres del Paine, Chile)

## Manuscript reference

**Post-Fire Soil Organic Matter Stabilization in Fire-Sensitive Patagonian Ecosystems**

Olivares-Martinez, Luis D. <sup>1 \*</sup>, Garcia-Franco, Noelia <sup>2,3 \*</sup>, Colocho Hurtarte, Luis C. <sup>4</sup>, García-Carmona, Minerva <sup>1</sup>,
Jaña, Jorge E. <sup>5</sup>, Arellano, Eduardo, C. <sup>5,6</sup>, Kögel-Knabner, Ingrid <sup>2</sup>, Wiesmeier, Martin <sup>2,7</sup>, Mataix-Solera, Jorge <sup>1</sup>

<sup>1</sup> Soil Science and Environmental Technologies Group (GETECMA), Department of Agrochemistry and Environment, Miguel Hernández University, Elche, Spain,
<sup>2</sup> Chair of Soil Science, TUM School of Life Sciences Weihenstephan, Technical University of Munich, Freising, Germany,
<sup>3</sup> Agroscope, Field-Crop Systems and Plant Nutrition Group, Route de Duillier 60, Nyon 1260, Switzerland
<sup>4</sup> Institute of Earth Surface Dynamics, University of Lausanne, Lausanne, Switzerland
<sup>5</sup> Facultad de Agronomía y Sistemas Naturales, Pontificia Universidad Católica de Chile, Santiago, Chile
<sup>6</sup> Center of Applied Ecology and Sustainability (CAPES), Pontificia Universidad Católica de Chile, Santiago, Chile
<sup>7</sup> Bavarian State Research Center for Agriculture, Institute for Agroecology and Organic Farming, Freising, Germany

## Supplementary data description

This repository contains the raw data, processing scripts, and references associated with the manuscript **"Post-Fire Soil Organic Matter Stabilization in Fire-Sensitive Patagonian Ecosystems"**, which investigates soil physical, chemical, and biochemical properties following the 2011 Olguín Fire in Torres del Paine National Park, Chile.

The supplementary information provides access to the datasets underlying the published results, enabling reproducibility and further analysis of post-fire soil dynamics in Patagonian ecosystems.

### Repository structure

The repository contains the following files:

- **`harmonized_DB_Chile_v3.4.xlsx`**: Contains the raw data of physical, chemical, biochemical, and fractionated soil organic matter properties, along with fire severity metrics.

- **`OM_stats.R`**: Script for data processing, statistical analysis, and plot generation of organic matter fractions for burned sites. The script includes instructions for data loading, exploratory analysis, statistical modeling, and post-hoc analysis to evaluate the effects of fire occurrence and severity.

- **`emm_abs.csv`**: Contains emmeans summaries for the models with the random effects of lowest AIC, showing the effect of fire presence/absence on each soil organic matter fraction.

- **`emm_sev.csv`**: Contains emmeans summaries for the models with the random effects of lowest AIC, showing the effect of fire severity on each soil organic matter fraction.

### Data content

The datasets include:

- **Soil physical properties**: Bulk density, particle size distribution, hydrophobicity, 
- **Soil chemical properties**: Electrical conductivity, pH, organic carbon, total nitrogen, C/N ratio, as well as free, occluded and mineral-associated organic matter fractions.
- **Soil organic matter properties by pool**: Organic carbon, total nitrogen, C/N ratio, functional groups and derived biomolecules within free, occluded and mineral-associated organic matter fractions.
- **Soil biochemical properties**: Enzymatic activities, microbial biomass carbon.
- **Fire severity metrics**: Burn severity indices derived from remote sensing data (i.e. dNBR) and derived for the Soil Burn Severity Index from USDA.

### Data availability

All data are openly available at:
https://github.com/OlivaresLD/OlguinFireCL

### Citation

If you use the data or scripts from this repository, please cite the associated manuscript:

Olivares-Martinez, L.D., Garcia-Franco, N., Colocho Hurtarte, L.C., García-Carmona, M., Jaña, J.E., Arellano, E.C., Kögel-Knabner, I., Wiesmeier, M., & Mataix-Solera, J. (2024). Post-Fire Soil Organic Matter Stabilization in Fire-Sensitive Patagonian Ecosystems. *Geoderma*.
