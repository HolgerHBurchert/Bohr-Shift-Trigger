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
install.packages(c('ggplot2', 'dplyr', 'tidyr', 'pracma', 'patchwork', 'cowplot', 'ggpubr', 'magick'))
```

## Content and how to use:
There are two R scripts in this repository. They are executing the same analysis, 
but with different inputs as described in the corresponding manuscript:

1. Script_T_37_to_38_DPG_00290: This script executes the analysis with a linear
temperature increase from 37 to 38 degrees Celsius and [2,3-DPG] set to 0.00290 M. 

2. Script_T_37_DPG_00350: Is the same code, but temperature is fixed to 37 degrees 
Celsius and [2,3-DPG] set to 0.00350 M. 

NOTE: Upon execution of either script, Figure_3, Figure_6_pcVO2max, Figure_4,
Figure_5 will be created according to that analysis. At the end of each script, 
Figure_4 and Figure_3 are combined to Combined_Figure3_AB in script 1. and to
Combined_Figure3_CD in script 2. Script 2. will finally combine Combined_Figure3_AB, 
Combined_Figure3_CD to the full Combined_Figure_3ABCD. For this reason, Figure_5 and 
Figure_6_pcVO2max are Fig. 4 and Fig. 5 respectively in the manuscript.

