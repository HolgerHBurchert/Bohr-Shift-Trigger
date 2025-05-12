# Bohr-Shift-Trigger
This repository contains the entire data set, analyses code and figures for “Oxygen Dissociation Curve Inflection Point During Incremental Exercise: A Trigger for the Bohr Effect”
The data in Stringer_1994.csv were extracted from Figure 1 from the article:

Stringer W, Wasserman K, Casaburi R, Porszasz J, Maehara K, French W. 
"Lactic acidosis as a facilitator of oxyhemoglobin dissociation during exercise" 
Journal of Applied Physiology. 1994 Apr;76(4):1462–7.
https://journals.physiology.org/doi/abs/10.1152/jappl.1994.76.4.1462

## Dependencies:
The R code requires that the following packages are installed, which are available at https://cran.r-project.org/ and may be directly installed by the following command:

```{r}
install.packages(c('ggplot2', 'dplyr', 'tidyr', 'pracma', 'patchwork', 'cowplot', 'ggpubr'))
```
