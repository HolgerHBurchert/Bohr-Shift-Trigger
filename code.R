# Developed by Holger Burchert (DPhil) - Department of Sport, Exercise and Health 
# University of Basel, CH-4052 Basel, Switzerland
# Email: holger.burchert@unibas.ch

# Many thanks to Dr. Denis Infanger for reviewing the code and pointing out more 
# efficient computations and helping to implement them

# Libraries required to run the code if not already installed run first this line
# install.packages(c('ggplot2', 'dplyr', 'tidyr', 'patchwork', 'cowplot', 'ggpubr'))
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(ggpubr)

#---------------------------------------------------------------------------------------
# Simulation of oxyhemoglobin (HbO2) and carbaminohemoglobin (HbCO2) dissociation
# curves and computation of total O2 and CO2 contents in whole blood, revised from
# the original model of Dash and Bassingthwaighte, ABME 38(4):1683-1701, 2010. The
# revision makes the model further simplified, as it bypasses the computations of
# the indices n1, n2, n3 and n4, which are complex expressions. Rather the revision
# necessiates the computations of K4p in terms of P50. Also the calculations of P50
# in terms of pH is enhanced based on a 3rd degree polynomial interpolation. 
#---------------------------------------------------------------------------------------
# Developed by: Ranjan Dash, PhD (Last modified: 2/16/2017; structured input and output)
# Department of Physiology and Biotechnology and Bioengineering Center
# Medical College of Wisconsin, Milwaukee, WI-53226
#---------------------------------------------------------------------------------------

# Input physiological variables for calculations of SHbO2CO2 and blood O2CO2 contents
SHbO2CO2 <- function(PO2 = 100, PCO2 = 40, pHrbc = 7.24, DPGrbc = 0.00465, Temp = 37, Hbrbc = Hbrbc_global, Hct = 0.45) {
  
  # This model can provide unrealistic P50 values if pHrbc > 10 (highly unphysiological).
  if (pHrbc > 10) {
    #disp('nonrealistic pHrbc > 10; so set it as 10'); disp(pHrbc);
    pHrbc <- 10
  }
  
  # Parameters those are fixed in the model (i.e. water fractions, RBCs hemoglobin
  # concentration, equilibrium constants, and Hill coefficient)
  Wpl <- 0.94                           # fractional water space of plasma; unitless
  Wrbc <- 0.65                          # fractional water space of RBCs; unitless
  Wbl <- (1 - Hct) * Wpl + Hct * Wrbc   # fractional water space of blood; unitless
  K1 <- 10^(-6.12)                      # CO2 hydration reaction equilibrium constant K
  K2 <- 21.5e-6                         # CO2 + HbNH2 equilibrium constant; unitless
  K2dp <- 1e-6                          # HbNHCOOH dissociation constant; M
  K2p <- K2 / K2dp                      # kf2p/kb2p; 1/M
  K3 <- 11.3e-6                         # CO2 + O2HbNH2 equilibrium constant; unitless
  K3dp <- 1e-6                          # O2HbNHCOOH dissociation constant; M
  K3p <- K3 / K3dp                      # kf3p/kb3p; 1/M
  K5dp <- 2.4e-8                        # HbNH3+ dissociation constant; M
  K6dp <- 1.2e-8                        # O2HbNH3+ dissociation constant; M
  Rrbc <- 0.69                          # Gibbs-Donnan ratio across the RBC membrane
  mol2ml <- 22267.4                     # Conversion factor from mol of gas to ml of gas at STP
  
  # Variables those are fixed in the model with values at standard physiological conditions
  # (i.e. PO20, PCO20, pHpl0, pHrbc0, DPGrbc0, Temp0, P500,aO20, aCO20)
  PO20 <- 100                           # standard O2 partial pressure in blood; mmHg
  PCO20 <- 40                           # standard CO2 partial pressure in blood; mmHg
  pHrbc0 <- 7.24                        # standard pH in RBCs; unitless
  pHpl0 <- pHrbc0 - log10(Rrbc)         # standard pH in plsama; unitless
  DPGrbc0 <- 4.65e-3                    # standard 2,3-DPG concentration in RBCs; M
  Temp0 <- 37                           # standard temperature in blood; degC
  P500 <- 26.8                          # standard PO2 for 50% SHbO2; mmHg
  aO20 <- 1.46e-6                       # solubility of O2 in water at 37 C; M/mmHg
  aCO20 <- 32.66e-6                     # solubility of CO2 in water at 37 C; M/mmHg
  # aO20 = 1.1129e-6;                   # solubility of O2 in water at 37 C; M/mmHg
  # aCO20 = 2.6715e-5;                  # solubility of CO2 in water at 37 C; M/mmHg
  
  # Calculation of intermediate variables in the computations of SHbO2 and SHbCO2
  pHpl <- pHrbc - log10(Rrbc)           # Rrbc = Hpl/Hrbc = 10^-(pHpl-pHrbc)
  delpHrbc <- pHrbc - pHrbc0; delpHpl <- pHpl - pHpl0
  delPCO2 <- PCO2 - PCO20; delDPGrbc <- DPGrbc - DPGrbc0; delTemp <- Temp - Temp0
  aO2 <- aO20 * (1 - 1e-2 * delTemp + 4.234e-4 * delTemp^2) # Corrected solubility of O2
  aCO2 <- aCO20 * (1 - 1.86e-2 * delTemp + 6.515e-4 * delTemp^2) # Corrected solubility of CO2
  O2 <- aO2 * PO2
  CO2 <- aCO2 * PCO2
  Hrbc <- 10^(-pHrbc)
  Hpl <- 10^(-pHpl)
  
  P501 <- P500 + 1.2 * (-21.279 * delpHrbc + 8.872 * delpHrbc^2 - 1.47 * delpHrbc^3) # all standard conditions, except pH
  P502 <- P500 + 1.7 * (4.28e-2 * delPCO2 + 3.64e-5 * delPCO2^2) # all standard conditions, except CO2
  P503 <- P500 + 1.0 * (795.633533 * delDPGrbc - 19660.8947 * delDPGrbc^2) # all standard conditions, except DPG
  P504 <- P500 + 0.98 * (1.4945 * delTemp + 4.335e-2 * delTemp^2 + 7e-4 * delTemp^3) # all standard conditions, except T
  P50 <- P500 * (P501 / P500) * (P502 / P500) * (P503 / P500) * (P504 / P500)
  C50 <- aO2 * P50
  
  # NEW IN THIS VERSION OF THE CODE: PO2 dependent variable Hill coefficient
  # nH = 2.7; Hill coefficient; unitless (redefined as a function PO2)
  alpha <- 2.8    #  Roughton et al data
  beta  <- 1.2    #  Roughton et al data
  gamma <- 29.2   #  Roughton et al data
  nH <- alpha - beta * 10^(-PO2 / gamma)    #  PO2 dependent variable nH
  
  # Compute the apparent equilibrium constant of Hb with O2 and CO2 (KHbO2 and KHbCO2); 
  # O2 and CO2 saturations of Hb (SHbO2 and SHbCO2); and O2 and CO2 contents in blood. 
  BPH1 <- 1 + K2dp / Hrbc # Binding polynomial involving K2dp and Hrbc 
  BPH2 <- 1 + K3dp / Hrbc # Binding polynomial involving K3dp and Hrbc 
  BPH3 <- 1 + Hrbc / K5dp # Binding polynomial involving K5dp and Hrbc 
  BPH4 <- 1 + Hrbc / K6dp # Binding polynomial involving K6dp and Hrbc 
  K4p <- (O2^(nH - 1) * (K2p * BPH1 * CO2 + BPH3)) / (C50^nH * (K3p * BPH2 * CO2 + BPH4))
  KHbO2 <- K4p * (K3p * BPH2 * CO2 + BPH4) / (K2p * BPH1 * CO2 + BPH3)
  KHbCO2 <- (K2p * BPH1 + K3p * K4p * BPH2 * O2) / (BPH3 + K4p * BPH4 * O2)
  SHbO2 <- KHbO2 * O2 / (1 + KHbO2 * O2)
  SHbCO2 <- KHbCO2 * CO2 / (1 + KHbCO2 * CO2)
  
  O2free <- Wbl * O2 #  M (mol O2 per L blood)
  O2bound <- 4 * Hct * Hbrbc * SHbO2 # M (mol O2 per L blood)
  O2tot <- O2free + O2bound # M (mol O2 per L blood)
  O2cont <- mol2ml * O2tot / 10 # mL O2/100 mL blood
  CO2free <- Wbl * CO2 # M (mol CO2 per L blood)
  CO2bicarb <- ((1 - Hct) * Wpl + Hct * Wrbc * Rrbc) * (K1 * CO2 / Hpl) # M (mol CO2 per L blood)
  CO2bound <- 4 * Hct * Hbrbc * SHbCO2 # M (mol CO2 per L blood)
  CO2totfb <- CO2free + CO2bound # M (mol CO2 per L blood)
  CO2contfb <- mol2ml * CO2totfb / 10 # mL CO2/100 mL blood
  CO2tot <- CO2free + CO2bicarb + CO2bound # M (mol CO2 per L blood)
  CO2cont <- mol2ml * CO2tot / 10 # mL CO2/100 mL blood
  
  # ALTERNATIVE CALCULATIONS FOR SHBO2 AND SHBCO2 BASED ON DIFFERENT HB-BOUND SPECIES
  HbNH2 <- Hbrbc / ((K2p * CO2 * BPH1 + BPH3) + K4p * O2 * (K3p * CO2 * BPH2 + BPH4))
  HbNH3p <- HbNH2 * Hrbc / K5dp
  O2HbNH2 <- K4p * O2 * HbNH2
  O2HbNH3p <- O2HbNH2 * Hrbc / K6dp
  HbNHCOOH <- K2p * CO2 * HbNH2
  HbNHCOOm <- K2dp * HbNHCOOH / Hrbc
  O2HbNHCOOH <- K3p * CO2 * O2HbNH2
  O2HbNHCOOm <- K3dp * O2HbNHCOOH / Hrbc
  SHbO2kin <- (O2HbNH2 + O2HbNH3p + O2HbNHCOOH + O2HbNHCOOm) / Hbrbc
  SHbCO2kin <- (HbNHCOOH + HbNHCOOm + O2HbNHCOOH + O2HbNHCOOm) / Hbrbc
  
  # Output physiological variables that are computed (e.g. SHbO2CO2 and blood O2CO2 contents)
  return(list(aO2 = aO2, 
              aCO2 = aCO2, 
              nH = nH, 
              P50 = P50, 
              K4p = K4p, 
              KHbO2 = KHbO2, 
              KHbCO2 = KHbCO2, 
              SHbO2 = SHbO2, 
              SHbCO2 = SHbCO2, 
              O2tot = O2tot, 
              O2cont = O2cont, 
              CO2tot = CO2tot, 
              CO2cont = CO2cont, 
              CO2totfb = CO2totfb,
              CO2contfb = CO2contfb, 
              HbNH2 = HbNH2, 
              HbNH3p = HbNH3p, 
              O2HbNH2 = O2HbNH2, 
              O2HbNH3p = O2HbNH3p, 
              HbNHCOOH = HbNHCOOH, 
              HbNHCOOm = HbNHCOOm, 
              O2HbNHCOOH = O2HbNHCOOH, 
              O2HbNHCOOm = O2HbNHCOOm, 
              SHbO2kin = SHbO2kin, 
              SHbCO2kin = SHbCO2kin,
              CO2bound = CO2bound,
              O2bound = O2bound))
  
}

#---------------------------------------------------------------------------------------
# BIOCHEMICAL REACTIONS FOR DERIVATION OF THE SHBO2 AND SHBCO2 EQUATIONS----------------
# The equations for O2 and CO2 saturations of hemoglobin (SHbO2 and SHbCO2) are  
# derived by considering the various kinetic reactions involving the binding of
# O2 and CO2 with hemoglobin in RBCs:
#
#            kf1p       K1dp
# 1. CO2+H2O <--> H2CO3 <--> HCO3- + H+;  K1=(kf1p/kb1p)*K1dp
#            kb1p		K1 = 7.43e-7 M; K1dp = 5.5e-4 M
#
#              kf2p          K2dp
# 2. CO2+HbNH2 <--> HbNHCOOH <--> HbNHCOO- + H+;  K2=(kf2p/kb2p)*K2dp
#              kb2p		K2 = 21.5e-6; K2dp = 1.0e-6 M
#
#                kf3p            K3dp
# 3. CO2+O2HbNH2 <--> O2HbNHCOOH <--> O2HbNHCOO- + H+; K3=(kf3p/kb3p)*K3dp
#                kb3p		K3 = 11.3e-6; K3dp = 1.0e-6 M
#
#              kf4p          
# 4. O2+HbNH2 <--> O2HbNH2;  K4p=func([O2];[H+];[CO2];[DPG];T)
#              kb4p
#
#           K5dp
# 5. HbNH3+ <--> HbNH2 + H+; K5dp = 2.4e-8 M
#
#             K6dp
# 6. O2HbNH3+ <--> O2HbNH2 + H+; K6dp = 1.2e-8 M
#---------------------------------------------------------------------------------------


#                   VALIDATING ABOVE VERSION of Dash's MODEL 
# ##############################################################################
# Verifying that this R implementation matches Prof. Dash's original Matlab code,
# using the default inputs where Matlab outputs SHbO2kin = 0.97554965890423.
Output <- SHbO2CO2(
  PO2 = 100,
  PCO2 = 40,
  pHrbc = 7.24,
  DPGrbc = 0.00465,
  Temp = 37,
  Hbrbc = 0.00528,
  Hct = 0.45
)

# SHbO2kin from R and Matlab will be formatted to character class so both 
# outputs have the same decimal places and class to compare if they are identical
SHbO2kin_R        <- format(Output$SHbO2kin,  digits = 14)
SHbO2kin_Matlab   <- format(0.97554965890423, digits = 14)
comparison_result <- identical(SHbO2kin_Matlab, SHbO2kin_R)
print(comparison_result) 


#                    LOADING and PREPARING STRINGER's DATA
# #############################################################################
# Loading PCO2, pH, PO2 values from Figure 4 of Stringer's 1994 article:
# Stringer W, Wasserman K, Casaburi R, Porszasz J, Maehara K, French W. 
# "Lactic acidosis as a facilitator of oxyhemoglobin dissociation during exercise" 
# Journal of Applied Physiology. 1994 Apr;76(4):1462–7.
# https://journals.physiology.org/doi/abs/10.1152/jappl.1994.76.4.1462
data <- read.csv("Stringer_1994.csv", header = TRUE, sep = ",")

# Convert Stringer's plasma pH (pHpl) to red blood cell pH (pHrbc); required by Dash's model.
# Hence, rearranging Gibbs-Donnan Equilibrium from Dash's model above: pHpl <- pHrbc - log10(Rrbc)
data$pH <- data$pH + log10(0.69)
names(data)[names(data) == "pH"] <- "pHrbc" # Now pH resembles pHrbc, renamed accordingly
names(data)[names(data) == "pco2"] <- "PCO2" # rename pco2 column to PCO2 for consistency
data$Sat_measured <- data$Sat_measured / 100 # Convert Sat_measured to fractional saturation (%saturation/100)

# Set global values for Temp, Hbrbc, Hct
Temp_global  <- 37 # Temperature of blood in Dash's model
Hbrbc_global <- 0.00528 # Hemoglobin concentration in red blood cells in M
Hct_global   <- 0.45 # Hematocrit 
max_bounded  <- 4*Hbrbc_global # The maximum of the bounded variables as PO2->100, to create right hand y-axis later in Figure 4


#                             FINDING 2,3-DPG VALUES
################################################################################
# Dash's model calculates O2 saturation from PO2, PCO2, pH, Hbrbc, Hct, Temp, and 2,3-DPG.
# Given Stringer's PO2, PCO2, pH, and measured O2 saturation (Sat_measured), with Temp=37, 
# Hbrbc=0.00528, and Hct=0.45, 2,3-DPG is the only unknown variable and can be back-calculated.
# Works as Hbrbc and Hct do not impact O2 saturation, and Temp is constant during the 10min test.

# Returns difference between calculated and observed O2 saturation for a given DPGrbc,
# used to find the DPGrbc that gives a root (zero) for this difference.
find_dpg <- function(dpg, i, data) {
  # Extract measured values for Stringer's Sample at index i
  PO2   <- data$PO2[i]
  PCO2  <- data$PCO2[i]
  pHrbc <- data$pHrbc[i]
  Sat_measured_fraction <- data$Sat_measured[i]
  # Return the difference between the modeled and measured O2 saturation
  SHbO2CO2(PO2,
           PCO2,
           pHrbc,
           dpg,
           Temp = Temp_global,
           Hbrbc = Hbrbc_global,
           Hct = Hct_global)$SHbO2kin - Sat_measured_fraction
}

# For each row in 'data', uniroot finds the DPGrbc value within the specified interval 
# that makes find_dpg return zero. uniroot adjusts DPGrbc until modeled & observed O2 saturation match.
data$DPGrbc <-
  sapply(seq_len(nrow(data)), FUN = \(i) (uniroot(
    find_dpg,
    interval = c(1e-10, 0.01),
    i = i,
    data = data,
    tol = 1e-20 
  )$root))

# Calculate saturation based on the found DPGrbc values
data$SHbO2kin <-
  unlist(
    mapply(
      SHbO2CO2,
      PO2    = data$PO2,
      PCO2   = data$PCO2,
      pHrbc  = data$pHrbc,
      DPGrbc = data$DPGrbc,
      Temp   = Temp_global,
      Hbrbc  = Hbrbc_global,
      Hct    = Hct_global
    )["SHbO2kin",]
  )


#                COMPARE MEASURED VS CALCULATED O2 SATURATION
################################################################################
# Creating Bland-Altman plot to compare measured O2 saturation with the calculated 
# O2 saturation based on the back-calculated 2,3-DPG values.
data$Mean <- rowMeans(data[,c('Sat_measured', 'SHbO2kin')], na.rm=TRUE)
data$Difference <- data$Sat_measured - data$SHbO2kin

bland_altman <- ggplot(data, aes(x = Mean, y = Difference)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + # Zero line
  geom_hline(yintercept = mean(data$Difference, na.rm = TRUE), col = "red") +
  geom_hline(
    yintercept = mean(data$Difference, na.rm = TRUE) + 1.96 * sd(data$Difference, na.rm = TRUE),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = mean(data$Difference, na.rm = TRUE) - 1.96 * sd(data$Difference, na.rm = TRUE),
    linetype = "dashed",
    color = "black"
  ) +
  labs(x = expression("Mean of Fractional O"[2] * "Hb Saturations"), y = "Difference") +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))

print(bland_altman)


#             2,3-DPG PLAUSIBILITY CHECK: PLOT DPGrbc VS %VO2MAX
################################################################################
# Plot DPGrbc versus %VO2max, to see if this also makes physiologically sense.
DPGrbc_vs_pcVO2max <-
  ggplot(data = data, aes(x = pcVO2max, y = DPGrbc)) +
  geom_path(aes(group = 1), color = "black", linewidth = 0.5) +
  geom_point(aes(fill="Original Data"), shape = 21, color = "black", size = 1.5, fill = "black") +
  labs(x = expression("%VO"[2] ~ max), y = "2,3-DPG (M)") +
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001)) +
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001), limits = c(0, 0.008)) +  # Setting explicit y-axis limits
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))

print(DPGrbc_vs_pcVO2max)

# Creating multi-panel figure
DPG_validation <- plot_grid(bland_altman, DPGrbc_vs_pcVO2max, nrow = 1, labels = c('A', 'B'))


#            CALCULATE DASH's OUTPUTS FOR STRINGER's 12 DATA POINTS
################################################################################
# Using the successfully back-calculated 2,3-DPG values, all necessary variables 
# (PO2, PCO2, pHrbc, Temp, Hbrbc, Hct, and 2,3-DPG) are now available for each of 
# Stringer's 12 data points, allowing the application of Dash's model to calculate outputs.

# Adjusted function call within the apply loop to use DPGrbc values from the data dataframe
results <- apply(data, 1, function(row) {
  SHbO2CO2(
    PO2 = row['PO2'],        # measured by Stringer 
    PCO2 = row['PCO2'],      # measured by Stringer
    pHrbc = row['pHrbc'],    # plasma pH measured by Stringer, pHrbc calculated by rearranging Gibbs Donnan condition
    DPGrbc = row['DPGrbc'],  # Now using back calculated DPGrbc from the data frame
    Temp = Temp_global,      # increases <1 degree celcius and standardized to 37 degree celcius in blood gas analysis     
    Hbrbc = Hbrbc_global,    # does not affect saturation calculation
    Hct = Hct_global         # does not affect saturation calculation
  )
})

# Convert the list of results into a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Check the names of the results to set the column names correctly
if (!is.null(colnames(results_df))) {
  colnames(results_df) <-
    c(
      'aO2',
      'aCO2',
      'nH',
      'P50',
      'K4p',
      'KHbO2',
      'KHbCO2',
      'SHbO2',
      'SHbCO2',
      'O2tot',
      'O2cont',
      'CO2tot',
      'CO2cont',
      'CO2totfb',
      'CO2contfb',
      'HbNH2',
      'HbNH3p',
      'O2HbNH2',
      'O2HbNH3p',
      'HbNHCOOH',
      'HbNHCOOm',
      'O2HbNHCOOH',
      'O2HbNHCOOm',
      'SHbO2kin',
      'SHbCO2kin',
      "CO2bound", 
      "O2bound"
    )
}

# Check the first few rows of the results to confirm
head(results_df)


#                   RECREATING FIGURE 4 of STRINGER et al.
################################################################################
# In this section, Figure 4 of Stringer et al. is recreated. This needs 4 steps:
# 1st, calculating Dash's ODCs for Stringers 12 data points (to plot them)
# 2nd, get PO2 and O2 saturation of Stringers 12 data points (to plot them)
# 3rd, get inflections for each ODC (to plot them)
# 4th, get the P50 values to plot them as well 

#               1st, CALCULATE DASH ODCs FOR STRINGERS 12 DATA POINTS    
# ------------------------------------------------------------------------------
# Initialize a list to store the results for each set of conditions from the data
# Each element in this list will be a data frame containing a dissociation curve
dissociation_curves <- list()

# Loop over each row in the 'data' data frame to generate dissociation curves for each condition
for (i in 1:nrow(data)) {
  # Create a sequence of PO2 values from 0 to 100 with small steps (0.01) for high resolution
  po2_values <- seq(0, 100, by = 0.01)
  
  # Initialize a data frame for storing the current dissociation curve
  # This includes PO2 values, and placeholders for O2 saturation and other values
  curve_data <- data.frame(
    PO2 = po2_values,  # PO2 values for the curve
    SHbO2 = numeric(length(po2_values)),  # Placeholder for O2 saturation values
    O2bound = numeric(length(po2_values)),
    CO2bound = numeric(length(po2_values)),
    HbNH3p = numeric(length(po2_values)),
    CurveID = rep(i, length(po2_values))  # Tagging each point with the current row number from 'data'
  )
  
  res_tmp <- SHbO2CO2(
    PO2 = curve_data$PO2,
    PCO2 = data$PCO2[i],               
    pHrbc = data$pHrbc[i],             
    DPGrbc = data$DPGrbc[i],           
    Temp = Temp_global,                         
    Hbrbc = Hbrbc_global,                      
    Hct = Hct_global
  )
  
  curve_data$SHbO2 <- unlist(res_tmp["SHbO2"])
  curve_data$O2bound <- unlist(res_tmp["O2bound"])
  curve_data$CO2bound <- unlist(res_tmp["CO2bound"])
  curve_data$HbNH3p <- unlist(res_tmp["HbNH3p"])
  
  # Append the curve data frame to the list of dissociation curves
  dissociation_curves[[i]] <- curve_data
}

# Combine all the individual dissociation curve data frames into a single data frame
# This results in 'combined_curves' containing all generated curves with distinguishable CurveIDs
combined_curves <- do.call(rbind, dissociation_curves)


#             2nd, GET PO2 AND O2 SATURATION OF STRINGERS 12 DATA POINTS
# --------------------------------------------------------------------------------
# To plot the original 12 data points from Stringer, we need their PO2 and calculated
# O2 saturation (SHbO2). SO these are extracted from the results_df and assign each point
# a CurveID based on its row in the original data
original_points <-
  data.frame(
    PO2 = data$PO2,
    SHbO2 = results_df$SHbO2,
    O2bound = results_df$O2bound,
    CurveID = 1:nrow(data)
  )


#                      3rd, GET INFLECTIONS FOR EACH ODC 
# ------------------------------------------------------------------------------
# Create function that numerically calculates the second derivative
deriv_fun <- function(x,
                      i = 1,
                      returnVar = "SHbO2",
                      n = 2) {
  pracma::fderiv(
    \(x)(
      SHbO2CO2(
        PO2 = x,
        PCO2 = data$PCO2[i],
        pHrbc = data$pHrbc[i],
        DPGrbc = data$DPGrbc[i],
        Temp = Temp_global,
        Hbrbc = Hbrbc_global,
        Hct = Hct_global
      )[[returnVar]]
    ),
    x = x,
    n = n,
    h = 0,
    method = "central"
  )
}

# Plot original function, first and second derivatives
par(mfrow = c(3, 1))
# Plot original output
i <- 1
curve(
  SHbO2CO2(
    x,
    PCO2 = data$PCO2[i],
    pHrbc = data$pHrbc[i],
    DPGrbc = data$DPGrbc[i],
    Temp = Temp_global,
    Hbrbc = Hbrbc_global,
    Hct = Hct_global
  )[["O2bound"]] / Hct_global,
  from = 0,
  to = 100,
  n = 1e4,
  main = "Original"
)
curve((\(x) (
  SHbO2CO2(
    x,
    PCO2 = data$PCO2[i],
    pHrbc = data$pHrbc[i],
    DPGrbc = data$DPGrbc[i],
    Temp = Temp_global,
    Hbrbc = Hbrbc_global,
    Hct = Hct_global
  )[["CO2bound"]] / Hct_global
))(x),
from = 0,
to = 100,
n = 1e4,
col = "red",
add = TRUE
)
curve(
  SHbO2CO2(
    x,
    PCO2 = data$PCO2[i],
    pHrbc = data$pHrbc[i],
    DPGrbc = data$DPGrbc[i],
    Temp = Temp_global,
    Hbrbc = Hbrbc_global,
    Hct = Hct_global
  )[["HbNH3p"]],
  from = 0,
  to = 100,
  n = 1e4,
  col = "steelblue",
  add = TRUE
)
abline(v = 23, lty = 2)
legend(
  "topright",
  legend = c("O2bound", "CO2bound", "HbNH3p"),
  col = c("black", "red", "steelblue"),
  lwd = c(1, 1, 1),
  bty = "n"
)


# Plot first derivatives
curve(
  deriv_fun(x, i = 1, returnVar = "O2bound", n = 1) / Hct_global,
  from = 0,
  to = 100,
  n = 1e4,
  ylim = c(-1e-4, 5.5e-4),
  main = "First derivative"
)
curve((\(x) (
  deriv_fun(x, i = 1, returnVar = "CO2bound", n = 1) / Hct_global
))(x),
from = 0,
to = 100,
n = 1e4,
add = TRUE,
col = "red"
)
curve(
  deriv_fun(x, i = 1, returnVar = "HbNH3p", n = 1),
  from = 0,
  to = 100,
  n = 1e4,
  add = TRUE,
  col = "steelblue"
)
abline(v = 23, lty = 2)
legend(
  "topright",
  legend = c("O2bound", "CO2bound", "HbNH3p"),
  col = c("black", "red", "steelblue"),
  lwd = c(1, 1, 1),
  bty = "n"
)

# Plot second derivatives
curve(
  deriv_fun(x, i = 1, returnVar = "O2bound", n = 2) / Hct_global,
  from = 0,
  to = 100,
  n = 1e4,
  ylim =  c(-2e-5, 2.5e-5),
  main = "Second derivative"
)
curve((\(x) (
  deriv_fun(x, i = 1, returnVar = "CO2bound", n = 2) / Hct_global
))(x),
from = 0,
to = 100,
n = 1e4,
add = TRUE,
col = "red"
)
curve(
  deriv_fun(x, i = 1, returnVar = "HbNH3p", n = 2),
  from = 0,
  to = 100,
  n = 1e4,
  add = TRUE,
  col = "steelblue"
)
abline(h = 0, v = 23, lty = 2)
legend(
  "topright",
  legend = c("O2bound", "CO2bound", "HbNH3p"),
  col = c("black", "red", "steelblue"),
  lwd = c(1, 1, 1),
  bty = "n"
)
par(mfrow = c(1, 1))

# Numerically find the zero of the second derivative to find the inflection point
data$inflection_point_SHbO2 <-
  sapply(seq_len(nrow(data)), FUN = \(i)(
    uniroot(
      deriv_fun,
      interval = c(15, 30),
      tol = 1e-20,
      i = i,
      returnVar = "O2bound"
    )$root
  ))

# Recalculate SHbO2 for the found inflection points
inflection_points <- data.frame(PO2 = data$inflection_point_SHbO2,
                                O2bound = unlist(
                                  mapply(
                                    SHbO2CO2,
                                    PO2 = data$inflection_point_SHbO2,
                                    PCO2 = data$PCO2,
                                    pHrbc = data$pHrbc,
                                    DPGrbc = data$DPGrbc,
                                    Temp = Temp_global,
                                    Hbrbc = Hbrbc_global,
                                    Hct = Hct_global
                                  )["O2bound",]
                                ))


#                     4th, CALCULATE P50 VALUES FOR EACH ODC
#-------------------------------------------------------------------------------
# Calculate P50 values for each CurveID
p50_data <- combined_curves %>%
  group_by(CurveID) %>%
  summarize(P50 = PO2[which.min(abs(SHbO2 - 0.5))]) %>%
  ungroup()

# Ensure combined_curves has a uniquely identifying column for merging, such as CurveID
# Merge based on CurveID and the P50 value matched to the PO2 column in combined_curves
p50_data <- p50_data %>%
  left_join(
    combined_curves %>% select(CurveID, PO2, SHbO2, O2bound) %>% distinct(),
    by = c("CurveID", "P50" = "PO2")
  )

#                  PUT PREVIOUS 4 STEPS TOGETHER TO PLOT FIGURE 4 
# ------------------------------------------------------------------------------
# Define a palette for the curves 
my_curves_colors <- c(
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C",  
  "#5C5C5C"   
)

# Define colors for data points, inflection points, and P50 points
data_points_edge <- "#000000"  
data_points_fill <- "#1874CD" 

inflection_points_edge <- "#5C5C5C"
inflection_points_fill <- "#000000"  

p50_points_edge <- "#5C5C5C" 
p50_points_fill <- "#000000"  

Figure_4 <- 
  ggplot() +
  geom_line( # Plotting the O2 dissociation curves 
    data = combined_curves,
    aes(x = PO2, y = O2bound/Hct_global, group = CurveID, color = as.factor(CurveID)),
    linewidth = 0.2,  
    show.legend = FALSE
  ) +
  geom_point( # Stringers original data points
    data = original_points,
    aes(x = PO2, y = O2bound/Hct_global, fill = "Original Data"), 
    shape = 21, color = data_points_edge, size = 1.5, fill = data_points_fill
  ) +
  annotate(
    "text",
    x = 2,
    y = 0.006,
    label = "VT1 amid sample 7 & 8 -->",
    hjust = 0,
    vjust = 1,
    size = 4,  
    color = "black"
  ) +
  annotate(
    "text",
    x = 2,
    y = 0.0077,
    label = "Sample 4 passing inflection -->",
    hjust = 0,
    vjust = 1,
    size = 4,  
    color = "black"
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.0076,
    label = "Inflection points",
    hjust = 0,
    vjust = 1,
    size = 4,  
    color = "black"
  ) +
  annotate(
    "text",
    x = 37,
    y = 0.0106,
    label = "P50s",
    hjust = 0,
    vjust = 1,
    size = 4,  
    color = "black"
  ) +
  annotate(
    "text",
    x = 24,
    y = 0.0025,
    label = "Inflections and Intersection",
    hjust = 0,
    vjust = 1,
    size = 4,  
    color = "black"
  ) +
  annotate("point",   x = 1,   y = 0.0110, shape = 21, size = 1.5, color = "#000000", fill = "#1874CD") +
  annotate("text",    x = 2.6, y = 0.0110, label = "Blood Samples", hjust = 0, size = 4) +
  annotate("segment", x = 1,   y = 0.0106, xend = 2, yend = 0.0106, color = "#5C5C5C", size = 1) +
  annotate("text",    x = 2.6, y = 0.0106, label = expression(O[2]~Dissociation~Curves), hjust = 0, size = 4) +
  annotate("segment", x = 1,   y = 0.0102, xend = 2, yend = 0.0102, color = "#D55E00", size = 1) +
  annotate("text",    x = 2.6, y = 0.0102, label = expression(CO[2]~Dissociation~Curve~of~Sample~4), hjust = 0, size = 4) +
  annotate("segment", x = 1,   y = 0.0098, xend = 2, yend = 0.0098, color = "#009E73", size = 1) +
  annotate("text",    x = 2.6, y = 0.0098, label = expression(H^"+"~Dissociation~Curve~of~Sample~4), hjust = 0, size = 4) +
  geom_point( # Inflection points of the ODCs
    data = inflection_points,
    aes(x = PO2, y = O2bound/Hct_global),
    size = 1.5, shape = 21, fill = inflection_points_fill, color = inflection_points_edge
  ) +
  geom_point( # P50 values of the ODCs
    data = p50_data,
    aes(x = P50, y = O2bound/Hct_global),
    size = 1.5, shape = 21, color = p50_points_edge, fill = p50_points_fill
  ) +
  scale_color_manual(values = my_curves_colors, name = "Curves") +
  scale_fill_manual(values = c("Original Data" = data_points_fill)) +
  guides(fill = "none") +
  scale_y_continuous(breaks = seq(0, 0.0205, 0.001)) +
  scale_x_continuous(breaks = seq(10, 40, 10)) +
  labs(
    x = expression(PO[2]~"(mmHg)"),
    y = expression(Fractional~O[2]~Saturation~("%"))
  ) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.011)) +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(angle = 90, size = 10),  # Increased axis title font size
    axis.title.x = element_text(size = 10),  # Increased axis title font size
    axis.title.y = element_text(size = 10)   # Increased axis title font size
  )

Figure_4

#                      Adding CO2bound, HbNH3p TO FIGURE 4 
################################################################################
# Filter the data for CurveID 4
data_for_curve_4 <- combined_curves[combined_curves$CurveID == 4, ]

Figure_4 <- Figure_4 +
  geom_line(data = data_for_curve_4, aes(x = PO2, y = CO2bound/Hct_global, group = 1), color = "#D55E00") +
  geom_line(data = data_for_curve_4, aes(x = PO2, y = HbNH3p, group = 1), color = "#009E73") +
  scale_y_continuous(
    name = expression(Hb-bound~to~O[2] * "," ~ CO[2] * "," ~ H^"+"~(mM~"in"~RBC~space)),
    labels = function(x) x * 1000,  # Scale the y-axis labels by 1000 (convert to mmol/L)
    sec.axis = sec_axis(~ (((100 - 0)*(. - 0))/(max_bounded - 0)) + 0, name = expression(O[2]*Hb~Saturation~("%")))
  ) +
  theme_minimal()+
  theme(
    axis.title.y.right = element_text(angle = 90, size = 10),  # Increased axis title font size
    axis.title.x = element_text(size = 10),  # Increased axis title font size
    axis.title.y = element_text(size = 10)   # Increased axis title font size
  )

# Print the updated plot
print(Figure_4)

#             CALCULATE THE INFLECTION POINTS OF CO2bound & HbNH3p 
################################################################################
# Numerically find the zero of the second derivative to find the inflection point
data$inflection_point_CO2bound <-
  sapply(seq_len(nrow(data)), FUN = \(i)(
    uniroot(
      deriv_fun,
      interval = c(15, 30),
      tol = 1e-20,
      i = i,
      returnVar = "CO2bound"
    )$root
  ))
data$inflection_point_HbNH3p <-
  sapply(seq_len(nrow(data)), FUN = \(i)(
    uniroot(
      deriv_fun,
      interval = c(15, 30),
      tol = 1e-20,
      i = i,
      returnVar = "HbNH3p"
    )$root
  ))

# Recalculate SHbO2 for the found inflection points
inflection_points <- data.frame(PO2 = data$inflection_point_CO2bound[4]
                                , as.data.frame(
                                  mapply(
                                    SHbO2CO2,
                                    PO2 = data$inflection_point_CO2bound,
                                    PCO2 = data$PCO2,
                                    pHrbc = data$pHrbc,
                                    DPGrbc = data$DPGrbc,
                                    Temp = Temp_global,
                                    Hbrbc = Hbrbc_global,
                                    Hct = Hct_global
                                  )[, 4]
                                ))


# Add inflection points to the plot 
Figure_4 <- Figure_4 +
  geom_point(
    data = inflection_points,
    aes(x = PO2, y = CO2bound/Hct_global, color = "Inflection CO2"),
    size = 3,
    shape = 19
  ) +
  geom_point(
    data = inflection_points,
    aes(x = PO2, y = HbNH3p, color = "Inflection HbNH3p"),
    size = 2,
    shape = 19
  ) +
  geom_vline(xintercept = 23.37, linetype = "longdash", color = "darkgrey", linewidth = 0.4) +
  scale_color_manual(values = c("Inflection CO2" = "#D55E00", "Inflection HbNH3p" = "#009E73"))+ 
  theme(legend.position = "none")  # This line removes the legend

# Print the updated plot
print(Figure_4)

Figure_4 <- Figure_4 +
  theme(axis.title.y.right = element_text(color = "black"))+
  theme(
    axis.title.x = element_text(size = 10),      # Increased axis title font size
    axis.title.y = element_text(size = 10),      # Increased axis title font size
    axis.text.x = element_text(size = 10),       # Increase x-axis ticks font size
    axis.text.y = element_text(size = 10),       # Increase left y-axis ticks font size
    axis.text.y.right = element_text(size = 10)  # Increase right y-axis ticks font size
  )

print(Figure_4)


#                FIND INTERSECTION OF HbNH3p AND CO2bound CURVES
################################################################################
# Assuming 'combined_curves' contains CO2bound and HbNH3p for CurveID 4
data_for_curve_4 <- combined_curves[combined_curves$CurveID == 4, ]

# Interpolate the curves
interp_CO2bound <- approxfun(data_for_curve_4$PO2, data_for_curve_4$CO2bound / Hct_global)
interp_HbNH3p <- approxfun(data_for_curve_4$PO2, data_for_curve_4$HbNH3p)

# Define a function that calculates the difference between the two interpolated functions
curve_diff <- function(x) {interp_CO2bound(x) - interp_HbNH3p(x)}

# Visualize the curves to identify a good interval
plot(data_for_curve_4$PO2, data_for_curve_4$CO2bound / Hct_global, type = "l", col = "red",
     ylim = range(c(data_for_curve_4$CO2bound / Hct_global, data_for_curve_4$HbNH3p)), ylab = "Values", xlab = "PO2")
lines(data_for_curve_4$PO2, data_for_curve_4$HbNH3p, col = "blue")
legend("topright", legend = c("CO2bound", "HbNH3p"), col = c("red", "blue"), lty = 1)

intersection <- uniroot(curve_diff, lower = 15, upper = 30)  # Adjusted interval based on the visualization

print(intersection$root) # Print the PO2 where the intersection occurs


#                             Figure 5 (CREATE 12 PANEL PLOT)
################################################################################
# Add an index to results_df as CurveID
results_df$CurveID <- as.factor(seq_len(nrow(results_df)))

# Define colors for the binding types
normal_colors <- c("O2bound" = "#1874CD", "CO2bound" = "#D55E00", "HbNH3p" = "#009E73")

# Prepare the results data points for plotting
results_long <- results_df %>%
  mutate(PO2 = data$PO2) %>%
  pivot_longer(cols = c(O2bound, CO2bound, HbNH3p), 
               names_to = "BindingType", 
               values_to = "Concentration") %>%
  select(CurveID, PO2, BindingType, Concentration)

# Convert combined_curves to long format
combined_long <- combined_curves %>%
  pivot_longer(cols = c(O2bound, CO2bound, HbNH3p), 
               names_to = "BindingType", 
               values_to = "Concentration")

# Prepare the inflection points data
inflection_points <- data.frame(
  CurveID = as.factor(seq_len(nrow(data))),
  PO2 = data$inflection_point_SHbO2
)

# Create a list to store individual plots
plot_list <- list()
num_panels <- length(unique(combined_long$CurveID)) # Total number of panels
n_cols <- 4 # Number of columns in the facet grid
n_rows <- ceiling(num_panels / n_cols) # Number of rows calculated based on total panels and columns

# Loop through each CurveID (panel)
for (i in 1:num_panels) {
  current_id <- i
  col_position <- ((i - 1) %% n_cols) + 1  # Column position (1 to n_cols)
  row_position <- ceiling(i / n_cols)      # Row position (1 to n_rows)
  
  # Separate the data into previous curves and current curves
  previous_curves <- combined_long %>%
    filter(as.numeric(CurveID) < current_id)
  
  current_curves <- combined_long %>%
    filter(as.numeric(CurveID) == current_id)
  
  # Separate the data points
  current_points <- results_long %>%
    filter(as.numeric(CurveID) == current_id)
  
  # Extract previous data points
  previous_points <- results_long %>%
    filter(as.numeric(CurveID) < current_id)
  
  # Create the plot
  p <- ggplot() +
    # Plot previous curves in grey
    geom_line(data = previous_curves, aes(
      x = PO2,
      y = ifelse(
        BindingType == "HbNH3p",
        Concentration,
        Concentration / Hct_global
      ),
      group = interaction(CurveID, BindingType)
    ), linewidth = 0.1, color = "grey", alpha = 0.4) +
    
    # Plot previous data points in grey
    geom_point(data = previous_points, aes(
      x = PO2,
      y = ifelse(
        BindingType == "HbNH3p",
        Concentration,
        Concentration / Hct_global
      ),
      color = BindingType
    ), size = 1.0) +
    
    # Plot current curves in color with solid lines
    geom_line(data = current_curves, aes(
      x = PO2,
      y = ifelse(
        BindingType == "HbNH3p",
        Concentration,
        Concentration / Hct_global
      ),
      color = BindingType
    ), linewidth = 0.6) +
    
    # Plot current data points in color
    geom_point(data = current_points, aes(
      x = PO2,
      y = ifelse(
        BindingType == "HbNH3p",
        Concentration,
        Concentration / Hct_global
      ),
      color = BindingType
    ), size = 1.5) +
    
    # Add vertical line for inflection point (optional)
    geom_vline(data = inflection_points %>% filter(as.numeric(CurveID) == current_id),
               aes(xintercept = PO2),
               linetype = "dashed",
               color = "#999999",
               size = 0.3) +
    
    # Set x and y-axis limits
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.011)) +
    
    # Set colors and labels
    scale_color_manual(values = normal_colors) +
    
    # Adjust y-axis scaling to display mmol/L without modifying the data
    scale_y_continuous(
      labels = scales::label_number(scale = 1000),
      name = if (col_position == 1) expression(Hb - O[2] * "," ~ CO[2] * "," ~ H^{"+ "} ~ plain("(mM in RBC space)")) else NULL
    ) +
    
    labs(
      title = paste("Sample", current_id),
      x = if (row_position == n_rows) expression(PO[2]~"(mmHg)") else NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title.x = if (row_position == n_rows) element_text(size = 9) else element_blank(),
      axis.title.y = if (col_position == 1) element_text(size = 9) else element_blank(),
      axis.text.x = if (row_position == n_rows) element_text() else element_blank(),
      axis.text.y = if (col_position == 1) element_text() else element_blank(),
      axis.ticks.x = if (row_position == n_rows) element_line() else element_blank(),
      axis.ticks.y = if (col_position == 1) element_line() else element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Add annotations to specific samples
  # Annotation for Sample 4
  if (current_id == 4) {
    p <- p +
      annotate(
        "text",
        x = 0,
        y = 0.010,
        label = "Passing",
        color = "black",
        size = 3,
        hjust = 0
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.009,
        label = "inflections &",
        color = "black",
        size = 3,
        hjust = 0
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.008,
        label = "Intersection",
        color = "black",
        size = 3,
        hjust = 0
      )
  }
  
  # Annotation for Sample 6
  if (current_id == 6) {
    p <- p +
      annotate(
        "text",
        x = 0,
        y = 0.010,
        label = "1st appreciable",
        color = "black",
        size = 3,
        hjust = 0
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.009,
        label = "right shift",
        color = "black",
        size = 3,
        hjust = 0
      )
  }
  
  # Annotation for Sample 7
  if (current_id == 7) {
    p <- p +
      annotate(
        "text",
        x = 0,
        y = 0.010,
        label = "VT1 mid 7 to 8",
        color = "black",
        size = 3,
        hjust = 0
      )
  }
  
  # Annotation for Sample 8
  if (current_id == 8) {
    p <- p +
      annotate(
        "text",
        x = 0,
        y = 0.010,
        label = "Curves start",
        color = "black",
        size = 3,
        hjust = 0
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.009,
        label = "drifting apart",
        color = "black",
        size = 3,
        hjust = 0
      )
  }
  
  
  # Add the plot to the list
  plot_list[[i]] <- p
}

# Combine all plots into one figure
Figure_5 <- wrap_plots(plot_list, ncol = n_cols)

# Remove the overall title by omitting the 'title' argument
Figure_5 <- Figure_5 +
  plot_annotation(
    theme = theme(
      plot.title = element_blank(),  # Ensure no title is displayed
      plot.caption = element_text(size = 10, hjust = 0.5)
    )
  ) &
  theme(
    plot.margin = unit(c(5, 5, 5, 5), "pt")  # Adjust margins if necessary
  )

# Display the multi-panel plot
print(Figure_5)

ggsave("Figure_5.pdf", Figure_5, width = 18, height = 18, units = "cm", dpi = 300, device = "pdf")


#        FITTING STRINGERS DATA TO RUN FITTED VALUES THROUGH DASH's MODEL
#################################################################################
# Fitting a 3rd order polynomial model for PO2
PO2_poly <- lm(PO2 ~ poly(pcVO2max, 3), data = data)
PO2_poly_pi <- predict(PO2_poly, interval = "prediction")
data_poly_pi_PO2 <- cbind(data, PO2_poly_pi)

# Plotting PO2 as a function of %VO2max with 3rd order polynomial fit
Figure_2A <- 
  ggplot(data = data_poly_pi_PO2, aes(x = pcVO2max, y = PO2)) +
  geom_point(shape = 21, color = "grey21", fill = "black", stroke = 0.1, size = 1.5) + 
  geom_smooth(method = lm, formula = y ~ poly(x, 3), se = TRUE, linewidth = 0.6, colour = "black") + 
  geom_line(aes(y = lwr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_line(aes(y = upr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_segment(x = 64, xend = 64, y = 18, yend = 21,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 64, y = 17.3, label = "VT1", size = 2.7, hjust = 0.5) +  # Add annotation text below the arrow
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")), 
    formula = y ~ poly(x, 3), 
    size = 3,
    label.x = 25,  # Adjust the x-position of the equation
    label.y = 29.8  # Adjust the y-position of the equation
  ) +
  xlab(expression("%" * dot(V) * O[2] ~ max)) +  
  ylab(expression(PO[2]~"(mmHg)")) + 
  theme_minimal() +
  theme(text = element_text(size = 10), axis.text = element_text(size = 8))

print(Figure_2A)

# Predict PO2 for new %VO2max values (1 to 100%)
fitted_values_PO2 <- data.frame(pcVO2max = seq(10, 100, by = 0.01))
fitted_values_PO2$PO2_fitted <- predict(PO2_poly, fitted_values_PO2)


# Fitting a linear model for PCO2 as function of %VO2max
PCO2_poly <- lm(PCO2 ~ pcVO2max, data = data)
PCO2_poly_pi <- predict(PCO2_poly, interval = "prediction")
data_poly_pi_PCO2 <- cbind(data, PCO2_poly_pi)

# Plotting PCO2 as a function of %VO2max with linear fit
Figure_2B <- 
  ggplot(data = data_poly_pi_PCO2, aes(x = pcVO2max, y = PCO2)) +
  geom_point(shape = 21, color = "grey21", fill = "black", stroke = 0.1, size = 1.5) + 
  geom_smooth(method = lm, formula = y ~ x, se = TRUE, linewidth = 0.6, colour = "black") + 
  geom_line(aes(y = lwr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_line(aes(y = upr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_segment(x = 64, xend = 64, y = 55, yend = 65,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 64, y = 53, label = "VT1", size = 2.7, hjust = 0.5) +  # Add annotation text below the arrow
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")), 
    formula = y ~ x, 
    size = 3,
    label.x = 30,  # Adjust the x-position of the equation
    label.y = 83  # Adjust the y-position of the equation
  ) +
  xlab(expression("%" * dot(V) * O[2] ~ max)) +  
  ylab(expression(PCO[2]~"(mmHg)")) + 
  theme_minimal() +
  theme(text = element_text(size = 10), axis.text = element_text(size = 8))

print(Figure_2B)

# Predict PCO2 for new %VO2max values (1 to 100%)
fitted_values_PCO2 <- data.frame(pcVO2max = seq(10, 100, by = 0.01))
fitted_values_PCO2$PCO2_fitted <- predict(PCO2_poly, fitted_values_PCO2)


# Fitting a 2nd order polynomial model for pH as function of %VO2max
pHrbc_poly <- lm(pHrbc ~ poly(pcVO2max, 2), data = data)
pHrbc_poly_pi <- predict(pHrbc_poly, interval = "prediction")
data_poly_pi <- cbind(data, pHrbc_poly_pi)

Figure_2C <- 
  ggplot(data = data_poly_pi, aes(x = pcVO2max, y = pHrbc)) +
  geom_point(shape = 21, color = "grey21", fill = "black", stroke = 0.1, size = 1.5) + 
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE, linewidth = 0.6, colour = "black") + 
  geom_line(aes(y = lwr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_line(aes(y = upr), color = "black", linetype = "dashed", linewidth = 0.4) + 
  geom_segment(x = 64, xend = 64, y = 7.0, yend = 7.05,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 64, y = 6.98, label = "VT1", size = 2.7, hjust = 0.5) +  # Add annotation text below the arrow
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")), 
    formula = y ~ poly(x, 2), 
    size = 3,
    label.x = 25,  # Adjust the x-position of the equation
    label.y = 7.22  # Adjust the y-position of the equation
  ) +
  xlab(expression("%" * dot(V) * O[2] ~ max)) +  
  ylab(expression(paste("pH"[RBC]))) + 
  theme_minimal() +
  theme(text = element_text(size = 10), axis.text = element_text(size = 8))

print(Figure_2C)


fitted_values <- data.frame(pcVO2max = seq(10, 100, by = 0.01))
fitted_values$pHrbc_fitted <- predict(pHrbc_poly, fitted_values)



# Fitting a 3rd order polynomial model for pH as function of %VO2max
Sat_measured_poly <- lm(Sat_measured ~ poly(pcVO2max, 3), data = data)
Sat_measured_poly_pi <- predict(Sat_measured_poly, interval = "prediction")
data_poly_pi_Sat_measured <- cbind(data, Sat_measured_poly_pi)

# Plotting Sat_measured as a function of %VO2max with 3rd order polynomial fit
Figure_2D <-
  ggplot(data = data_poly_pi_Sat_measured, aes(x = pcVO2max, y = Sat_measured)) +
  geom_point(shape = 21, color = "grey21", fill = "black", stroke = 0.1, size = 1.5) +
  geom_smooth(method = lm, formula = y ~ poly(x, 3), se = TRUE, linewidth = 0.6, colour = "black") +
  geom_line(aes(y = lwr), color = "black", linetype = "dashed", linewidth = 0.4) +
  geom_line(aes(y = upr), color = "black", linetype = "dashed", linewidth = 0.4) +
  geom_segment(x = 64, xend = 64, y = 0.4, yend = 0.3,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 64, y = 0.43, label = "VT1", size = 2.7, hjust = 0.5) +  # Add annotation text below the arrow
  stat_regline_equation(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ poly(x, 3),
    size = 3,
    label.x = 10,  # Adjust the x-position of the equation
    label.y = 0.55  # Adjust the y-position of the equation
  ) +
  xlab(expression("%" * dot(V) * O[2] ~ max)) +  
  ylab(expression("Fract. O"[2]*Hb ~ "Saturation")) +
  theme_minimal() +
  theme(text = element_text(size = 10), axis.text = element_text(size = 8))

print(Figure_2D)


# Combine all updated subplots into a grid
Figure_2 <- plot_grid(Figure_2A, Figure_2B, Figure_2C, Figure_2D, nrow = 2, labels = c('A', 'B', 'C', 'D'))

ggsave("Figure_2.pdf", Figure_2, width = 18, height = 12, units = "cm")


# Predict Sat_measured for new %VO2max values (1 to 100%)
fitted_values_Sat_measured <- data.frame(pcVO2max = seq(10, 100, by = 0.01))
fitted_values_Sat_measured$Sat_measured_fitted <- predict(Sat_measured_poly, fitted_values_Sat_measured)

# Print the new data frame with fitted values for Sat_measured
print(fitted_values_Sat_measured)


# Create a new data frame for the predicted values
fitted_values <- data.frame(pcVO2max = seq(10, 100, by = 0.01))
fitted_values$pHrbc_fitted <- predict(pHrbc_poly, fitted_values)
fitted_values$PO2_fitted <- predict(PO2_poly, fitted_values)
fitted_values$PCO2_fitted <- predict(PCO2_poly, fitted_values)
fitted_values$Sat_measured_fitted <- predict(Sat_measured_poly, fitted_values)
print(fitted_values)


# Modified function that calculates the difference using the fitted values
find_dpg_fitted <- function(dpg, i, fitted_values) {
  PO2 <- fitted_values$PO2_fitted[i]
  PCO2 <- fitted_values$PCO2_fitted[i]
  pHrbc <- fitted_values$pHrbc_fitted[i]
  Sat_measured_fraction <- fitted_values$Sat_measured_fitted[i]  # Fitted O2 saturation
  SHbO2CO2(PO2,
           PCO2,
           pHrbc,
           dpg,
           Temp = Temp_global,
           Hbrbc = Hbrbc_global,
           Hct = Hct_global)$SHbO2kin - Sat_measured_fraction
}


# Use the fitted values to back-calculate 2,3-DPG
fitted_values$DPGrbc_fitted <- sapply(seq_len(nrow(fitted_values)), FUN = \(i) (
  uniroot(
    find_dpg_fitted, 
    interval = c(1e-10, 0.01), 
    i = i, 
    fitted_values = fitted_values, 
    tol = 1e-20
  )$root
))

# Print the new data frame with back-calculated 2,3-DPG values
print(fitted_values)


# Assuming your new dataset is also called `fitted_values`
fitted_values$SHbO2kin <-
  unlist(
    mapply(
      SHbO2CO2,
      PO2 = fitted_values$PO2,
      PCO2 = fitted_values$PCO2,
      pHrbc = fitted_values$pHrbc,
      DPGrbc = fitted_values$DPGrbc,
      Temp = Temp_global,
      Hbrbc = Hbrbc_global,
      Hct = Hct_global
    )["SHbO2kin",]
  )


# Select the relevant columns from the `data` dataframe
data_selected <- data[, c('Sat_measured', 'SHbO2kin', 'DPGrbc', 'pcVO2max')]

# Rename columns in `fitted_values` to match `data`
fitted_values_selected <- fitted_values[, c('Sat_measured_fitted', 'SHbO2kin', 'DPGrbc_fitted', 'pcVO2max')]
colnames(fitted_values_selected) <- c('Sat_measured', 'SHbO2kin', 'DPGrbc', 'pcVO2max')

# Combine the two data frames by stacking rows
combined_data <- rbind(data_selected, fitted_values_selected)

# Check the new combined data frame
print(combined_data)


# First, ensure you have the original and fitted data separated
data_selected <- data[, c('pcVO2max', 'DPGrbc')] # Original data
fitted_values_selected <- fitted_values[, c('pcVO2max', 'DPGrbc_fitted')] # Fitted data

# Plotting 2,3-DPG vs. %VO2max
Figure_3A <- ggplot() +
  # Plot original data as points
  geom_point(data = data_selected, aes(x = pcVO2max, y = DPGrbc), color = "black", size = 2) +
  # Plot fitted data as a line
  geom_line(data = fitted_values_selected, aes(x = pcVO2max, y = DPGrbc_fitted), color = "black", linewidth = 0.5) +
  labs(x = expression("%" * dot(V) * O[2] ~ max), y = "2,3-DPG (mM)") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(
    breaks = seq(0, 0.006, by = 0.001),  # Specify the breaks in original units
    limits = c(0, 0.006),                # Keep the limits in the original scale
    labels = scales::label_number(scale = 1000)  # Multiply the axis labels by 1000 to convert to mM
  ) +  
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))

# Print the combined plot
print(Figure_3A)


# Calculate means and differences
combined_data$Mean <- rowMeans(combined_data[,c('Sat_measured', 'SHbO2kin')], na.rm=TRUE)
combined_data$Difference <- combined_data$Sat_measured - combined_data$SHbO2kin

# Plot Bland-Altman Plot
Figure_3B <- ggplot(combined_data, aes(x = Mean, y = Difference)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + # Zero line
  geom_hline(yintercept = mean(combined_data$Difference, na.rm = TRUE), col = "red") +
  geom_hline(
    yintercept = mean(combined_data$Difference, na.rm = TRUE) + 1.96 * sd(combined_data$Difference, na.rm = TRUE),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = mean(combined_data$Difference, na.rm = TRUE) - 1.96 * sd(combined_data$Difference, na.rm = TRUE),
    linetype = "dashed",
    color = "black"
  ) +
  labs(x = expression("Mean of Fractional O"[2] * "Hb Saturations"), y = "Difference") +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))

print(Figure_3B)


# Combine the updated Figure 6 panels into one figure
Figure_3 <- plot_grid(
  Figure_3A,
  Figure_3B,
  nrow = 1,
  labels = c('A', 'B')
)

# Print the updated Figure 6
print(Figure_3)

ggsave("Figure_3.pdf", Figure_3, width = 18, height = 5.5, units = "cm")



#            CALCULATE DASH's OUTPUT FOR fitted_values POINTS
################################################################################
# Now, for the fitted_values points, PO2_fitted, PCO2_fitted, pHrbc_fitted, Temp, Hbrbc, Hct, and 2,3-DPG_fitted 
# have been made available/assumed so that Dash's model can be applied to them.

# Adjusted function call within the apply loop to use DPGrbc_fitted values from the fitted_values dataframe
results_fitvalues <- apply(fitted_values, 1, function(row) {
  SHbO2CO2(
    PO2 = row['PO2_fitted'],      # Fitted PO2 values
    PCO2 = row['PCO2_fitted'],    # Fitted PCO2 values
    pHrbc = row['pHrbc_fitted'],  # Fitted pHrbc values
    DPGrbc = row['DPGrbc_fitted'],# Back-calculated DPGrbc_fitted values
    Temp = Temp_global,           # Temperature (standardized to 37°C)
    Hbrbc = Hbrbc_global,         # Standard Hbrbc
    Hct = Hct_global              # Standard Hct
  )
})

# Convert the list of results into a data frame
results_fitvalues_df <- do.call(rbind, lapply(results_fitvalues, as.data.frame))

# Check the names of the results to set the column names correctly
if (!is.null(colnames(results_fitvalues_df))) {
  colnames(results_fitvalues_df) <-
    c(
      'aO2',
      'aCO2',
      'nH',
      'P50',
      'K4p',
      'KHbO2',
      'KHbCO2',
      'SHbO2',
      'SHbCO2',
      'O2tot',
      'O2cont',
      'CO2tot',
      'CO2cont',
      'CO2totfb',
      'CO2contfb',
      'HbNH2',
      'HbNH3p',
      'O2HbNH2',
      'O2HbNH3p',
      'HbNHCOOH',
      'HbNHCOOm',
      'O2HbNHCOOH',
      'O2HbNHCOOm',
      'SHbO2kin',
      'SHbCO2kin',
      "CO2bound", 
      "O2bound"
    )
}

# Check the first few rows of the results to confirm
head(results_fitvalues_df)


plot(results_fitvalues_df$O2bound/ Hct_global ~ fitted_values$PO2_fitted , xlim = c(0,40))


ggplot(data = fitted_values, aes(x = pcVO2max, y = results_fitvalues_df$KHbO2)) +
  geom_path(aes(group = 1), color = "black", size = 0.5) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  labs(x = expression("%VO"[2] ~ max), y = expression('KHbO'[2] ~ (M ^ -1))) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  geom_segment(x = 32.5, xend = 32.5, y = 12800, yend = 15000,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 32.5, y = 11000, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5)+
  geom_segment(x = 56, xend = 56, y = 15000, yend = 12800,
               arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = 56, y = 16500, label = "VT1", size = 2.7, hjust = 0.5) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 10)) 


# Adding the new dataset points to Figure 1 without title, legend, and only showing the black line
Figure_4_updated <- Figure_4 + 
  # Plot the O2bound/Hct_global from the newly fitted data (`results_fitvalues_df`) against `PO2_fitted`
  geom_line(
    data = fitted_values,
    aes(x = PO2_fitted, y = results_fitvalues_df$O2bound / Hct_global),
    color = "#1874CD",  # Change the line color to black
    linewidth = 0.8,       # Adjust line width for visibility
    show.legend = FALSE  # Remove legend
  ) +
  # Remove the overall legend
  theme(
    legend.position = "none",
    text = element_text(size = 10), 
    axis.text = element_text(size = 10)
  )

# Print the updated figure without title, points, and legend
print(Figure_4_updated)

# ggsave("Figure_4_updated.pdf", Figure_4_updated, width = 20, height = 20, units = "cm")


plot(results_fitvalues_df$O2bound / Hct_global ~ fitted_values$PO2_fitted, type = "l")


predicted_fun <- Vectorize(function(x, var = "O2bound") {
  
  assumed_pcVO2max <- uniroot(
    \(z)(predict(PO2_poly, newdata = data.frame(pcVO2max = z)) - x)
    , lower = 0
    , upper = 100
    , tol = 1e-20
  )$root
  
  tmp_frame <- data.frame(
    PO2_fitted = x
    , PCO2_fitted = predict(PCO2_poly, newdata = data.frame(pcVO2max = assumed_pcVO2max))
    , pHrbc_fitted = predict(pHrbc_poly, newdata = data.frame(pcVO2max = assumed_pcVO2max))
    , Sat_measured_fraction = predict(Sat_measured_poly, newdata = data.frame(pcVO2max = assumed_pcVO2max))
  )
  
  find_dpg <- function(dpg, data) {
    SHbO2CO2(
      PO2 = x,     
      PCO2 = data$PCO2_fitted,   
      pHrbc = data$pHrbc_fitted, 
      DPGrbc = dpg,
      Temp = Temp_global,          
      Hbrbc = Hbrbc_global,        
      Hct = Hct_global             
    )$SHbO2kin - data$Sat_measured_fraction
  }
  
  tmp_frame$DPGrbc_fitted <- uniroot(
    find_dpg, 
    interval = c(1e-10, 0.01), 
    data = tmp_frame, 
    tol = 1e-20
  )$root
  
  SHbO2CO2(
    PO2 = x,
    PCO2 = tmp_frame$PCO2_fitted,
    pHrbc = tmp_frame$pHrbc_fitted,
    DPGrbc = tmp_frame$DPGrbc_fitted,
    Temp = Temp_global,
    Hbrbc = Hbrbc_global,
    Hct = Hct_global
  )[[var]]
  
}, "x")

deriv_fun <- function(x,
                      returnVar = "O2bound",
                      n = 2) {
  pracma::fderiv(
    \(x)(
      predicted_fun(
        x
        , var = returnVar
      )
    ),
    x = x,
    n = n,
    h = 0,
    method = "central"
  )
}

inflection_point_PO2 <- uniroot(
  deriv_fun,
  interval = c(19, 28),
  tol = 1e-20,
  returnVar = "O2bound"
)$root

plot(results_fitvalues_df$O2bound / Hct_global ~ fitted_values$PO2_fitted, type = "l")
abline(v = inflection_point_PO2, col = "red")


# Retrieve the y-coordinate for the inflection point from the results_fitvalues_df
y_value_at_inflection <- approx(
  x = fitted_values$PO2_fitted,
  y = results_fitvalues_df$O2bound / Hct_global,
  xout = inflection_point_PO2
)$y


# Add red dot at the inflection point
Figure_4_updated <-
  Figure_4_updated +
  geom_point(aes(x = inflection_point_PO2, y = y_value_at_inflection),
             colour = "red",
             size = 2) # Adjust size if necessary

# Save the updated figure
ggsave("Figure_4.pdf", Figure_4_updated, width = 15, height = 17, units = "cm")


#                                   FIGURE 6 
################################################################################
Figure_6A <- 
  ggplot(data = results_df, aes(x = HbNH3p, y = P50)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values,
    aes(x = results_fitvalues_df$HbNH3p, y = results_fitvalues_df$P50), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(x = expression(H^"+"~"bound (mM)"),
       y = "P50 (mmHg)") +
  scale_x_continuous(labels = function(x) x * 1000) +  # Multiply x-axis values by 1000
  annotate("text", x = 0.002100000, y = 31.8, label = "Sample 4 \n passing inflection", size = 2.7, hjust = 0.5) +
  annotate("segment", x = 0.002167339, xend = 0.002167339, y = 30.8, yend = 29.8,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8)) 


Figure_6B <- 
  ggplot(data = results_df, aes(x = CO2bound/Hct_global, y = P50)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_path(
    data = fitted_values,
    aes(x = results_fitvalues_df$CO2bound/Hct_global, y = results_fitvalues_df$P50), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(
    x = expression(CO[2]~"bound (mM)"),
    y = "P50 (mmHg)"
  ) +
  scale_x_continuous(labels = function(x) x * 1000) +  # Multiply x-axis values by 1000
  annotate("text", x = 0.00195, y = 31.2, label = "Sample 4 \n passing inflection", size = 2.7, hjust = 0.5) +
  annotate("segment", x = 0.00207, xend = 0.00213, y = 30.5, yend = 29.5,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_6C <- 
  ggplot(data = results_df, aes(x = SHbO2, y = P50)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values,
    aes(x = results_fitvalues_df$SHbO2, y = results_fitvalues_df$P50), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(
    x = expression(Fractional~O[2]*Hb~saturation~"(% Saturation/100)"),
    y = "P50 (mmHg)"
  ) +
  annotate("text", x = 0.36, y = 32, label = "Sample 4 \n passing inflection", size = 2.7, hjust = 0.5) +
  annotate("segment", x = 0.36, xend = 0.36, y = 31, yend = 30,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_6D <- 
  ggplot(data = results_df, aes(x = CO2bound/Hct_global, y = HbNH3p)) +
  geom_point(shape = 21, color = data_points_edge, size = 1.5, fill = "black") +
  geom_path(
    data = fitted_values,
    aes(x = results_fitvalues_df$CO2bound/Hct_global, y = results_fitvalues_df$HbNH3p),
    color = "black",
    linewidth = 0.5
  ) +
  labs(x = expression(CO[2]~"bound (mM)"), y = expression(H^"+"~"bound (mM)")) +
  scale_x_continuous(labels = function(x) x * 1000) +  # Multiply x-axis values by 1000
  scale_y_continuous(
    labels = function(y)
      y * 1000,
    limits = c(0.0016, 0.0034),
    breaks = seq(0.0016, 0.0035, by = 0.0002)
  ) +  # Multiply y-axis values by 1000
  annotate(
    "text",
    x = results_df$CO2bound[4] / Hct_global - 0.00003,
    y = 0.00172,
    label = "Sample 4 \n passing inflection",
    size = 2.7,
    hjust = 0.5
  ) +
  annotate(
    "segment",
    x = results_df$CO2bound[4] / Hct_global,
    xend = results_df$CO2bound[4] / Hct_global,
    y = 0.0019,
    yend = 0.0021,
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


# Combine the updated Figure 6 panels into one figure
Figure_6 <- plot_grid(Figure_6A, Figure_6B, Figure_6C, Figure_6D, nrow = 2, labels = c('A', 'B', 'C', 'D'))

ggsave("Figure_6.pdf", Figure_6, width = 18, height = 11, units = "cm")


#                                   FIGURE 7
#                        CALCULATE HCO3- from PCO2 and pH
################################################################################
# HCO3m calculation based on Beaver equation derived from Siggard-Andersen Normogram
# Stringer, W.W., R. Casaburi, and K. Wasserman. Acid base regulation during exercise and recovery in humans. 
# J. Appl. Physiol. 72: 954-961, 1992. 
HCO3m_SA <- function(pCO2, pH, Hb = 15.4) {
  HCO3m_mmol_per_l <- 10^(pH-7.376*(1-0.00305*Hb)+0.848*(1-0.0162*Hb)*log10(pCO2))
  return(HCO3m_mmol_per_l)  
}

data$pHpl <- data$pHrbc - log10(0.69)
data$HCO3m <- HCO3m_SA(data$PCO2, data$pHpl)

# Adjust the fitted pH values back to plasma pH values by rearranging Gibbs Donnan equilibrium
fitted_values$pHpl_fitted <- fitted_values$pHrbc_fitted - log10(0.69)

# Calculate the HCO3m for the fitted values using the plasma pH
fitted_values$HCO3m_fitted <- HCO3m_SA(fitted_values$PCO2_fitted, fitted_values$pHpl_fitted)

Figure_7A <-
  ggplot(data = data, aes(x = PO2, y = results_df$HbNH2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$HbNH2), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(
    x = expression("PO"[2] ~ "(mmHg)"), 
    y =  expression(HbNH[2]~"(mM)"))+
  scale_y_continuous(
    labels = function(y) y * 1000) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))

Figure_7B <- 
  ggplot(data = data, aes(x = PO2, y = results_df$HbNHCOOm)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$HbNHCOOm), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(
    x = expression("PO"[2] ~ "(mmHg)"), 
    y = expression("HbNHCOO"^"-"~ "(mM)")) +
  scale_y_continuous(labels = function(y) y * 1000) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)  # Adjust title size and alignment
  )


Figure_7C <-
  ggplot(data = data, aes(x = PO2, y = results_df$O2HbNHCOOm)) +
  geom_point(
    aes(fill = "Original Data"), shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values,
    aes(x = PO2_fitted, y = results_fitvalues_df$O2HbNHCOOm),
    color = "black",
    linewidth = 0.5
  ) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression("O"[2] * "HbNHCOO"^"-"~"(mM)")) +
  scale_y_continuous(labels = function(y) y * 1000) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +  # Set x-axis breaks from 15 to 30 in increments of 2
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.0)  # Adjust title size and alignment
  )  


Figure_7D<-
  ggplot(data = data, aes(x = PO2, y = results_df$HbNH3p)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$HbNH3p), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression("HbNH"[3]^"+"~"(mM)")) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +  # Set x-axis breaks from 15 to 30 in increments of 2
  scale_y_continuous(labels = function(y) y * 1000) +
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)  # Adjust title size and alignment
  )      


Figure_7E <-
  ggplot(data = data, aes(x = PO2, y = results_df$O2HbNH2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$O2HbNH2), color = "black", linewidth = 0.5) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression("O"[2] * "HbNH" [2]~"(mM)")) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +  # Set x-axis breaks from 15 to 30 in increments of 2
  scale_y_continuous(labels = function(y) y * 1000) +
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)  # Adjust title size and alignment
  )    


Figure_7F <-
  ggplot(data = data, aes(x = PO2, y = results_df$O2HbNH3p)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$O2HbNH3p), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression("O"[2] * "HbNH"[3]^"+"~"(mM)")) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +  # Set x-axis breaks from 15 to 30 in increments of 2
  scale_y_continuous(labels = function(y) y * 1000) +
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)  # Adjust title size and alignment
  )     


Figure_7G <-
  ggplot(data = data, aes(x = PO2, y = HCO3m)) +
  geom_point(
    aes(fill = "Original Data"), shape = 21, color = data_points_edge, size = 1.5, fill = "black") +
  geom_line(data = fitted_values, aes(x = PO2_fitted, y = results_fitvalues_df$HbNHCOOm), color = "black", linewidth = 0.5) +
  geom_line(data = fitted_values, aes(x = PO2_fitted, y = HCO3m_fitted), color = "black", linewidth = 0.5) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression("HCO"[3]^"-"~"(mmol/L blood)")) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +  # Set x-axis breaks from 15 to 30 in increments of 2
  ylim(18, 26) +  # Set y-axis limits from 18 to 26
  theme_minimal() +
  theme(
    text = element_text(size = 9), 
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)  # Adjust title size and alignment
  )     
# Print the updated Figure 7A
print(Figure_7G)


Figure_7 <- plot_grid(Figure_7A, Figure_7B, Figure_7C, Figure_7D, Figure_7E, Figure_7F, Figure_7G,
                      nrow = 3, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G')
)

ggsave("Figure_7.pdf", Figure_7, width = 18, height = 16, units = "cm")


#                                    Figure 8
################################################################################
Figure_8A <-
  ggplot(data = data, aes(x = pcVO2max, y = results_df$KHbO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = pcVO2max, y = results_fitvalues_df$KHbO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 55, xend = 40, y = 19000, yend = 17200, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 70, y = 21000, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression("%" * dot(V) * O[2] ~ max), y = expression('KHbO'[2] ~ (M^-1))) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8B <- 
  ggplot(data = data, aes(x = pHrbc, y = results_df$KHbO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = pHrbc_fitted, y = results_fitvalues_df$KHbO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 7.1, xend = 7.15, y = results_df$KHbO2[4], yend = results_df$KHbO2[4], 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 7.01, y = 17000, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression('pH'[RBC]), y = expression('KHbO'[2] ~ (M^-1))) +
  #  scale_x_continuous(breaks = seq(15, 30, by = 2)) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8C <- 
  ggplot(data = data, aes(x = PCO2, y = results_df$KHbO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PCO2_fitted, y = results_fitvalues_df$KHbO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 60, xend = 55, y = 19000, yend = results_df$KHbO2[4] + 500, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 70, y = 21000, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression("PCO"[2] ~ "(mmHg)"), y = expression('KHbO'[2] ~ (M^-1))) +
  #  scale_x_continuous(breaks = seq(15, 30, by = 2)) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8D <-
  ggplot(data = data, aes(x = pcVO2max, y = results_df$KHbCO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = pcVO2max, y = results_fitvalues_df$KHbCO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 32.5, xend = 32.5, y = 55, yend = 60, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 32.5, y = 51, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression("%" * dot(V) * O[2] ~ max), y = expression('KHbCO'[2] ~ (M^-1))) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8E <- 
  ggplot(data = data, aes(x = pHrbc, y = results_df$KHbCO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = pHrbc_fitted, y = results_fitvalues_df$KHbCO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 7.05, xend = 7.15, y = results_df$KHbCO2[4], yend = results_df$KHbCO2[4], 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 7.0, y = 63, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression('pH'[RBC]), y = expression('KHbCO'[2] ~ (M^-1))) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8F <- 
  ggplot(data = data, aes(x = PO2, y = results_df$KHbCO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$KHbCO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 25, xend = 24, y = 55, yend = 60.5, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 25.5, y = 50, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression('KHbCO'[2] ~ (M^-1))) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8G <- 
  ggplot(data = data, aes(x = PO2, y = results_df$KHbO2)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$KHbO2), 
    color = "black", 
    linewidth = 0.5
  ) +
  geom_segment(x = 22, xend = 23, y = 19000, yend = results_df$KHbO2[4] + 500, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 21.5, y = 21000, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = expression('KHbO'[2] ~ (M^-1))) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8H <-
  ggplot(data = data, aes(x = PO2, y = results_df$nH)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = PO2_fitted, y = results_fitvalues_df$nH), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(x = expression("PO"[2] ~ "(mmHg)"), y = "nH") +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +
  theme_minimal() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8I <-
  ggplot(data = data, aes(x = pcVO2max, y = results_df$nH)) +
  geom_point(
    aes(fill = "Original Data"),
    shape = 21,
    color = data_points_edge,
    size = 1.5,
    fill = "black"
  ) +
  geom_line(
    data = fitted_values, 
    aes(x = pcVO2max, y = results_fitvalues_df$nH), 
    color = "black", 
    linewidth = 0.5
  ) +
  labs(x = expression("%" * dot(V) * O[2] ~ max), y = "nH") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  geom_segment(x = 32.5, xend = 32.5, y = 2.58, yend = 2.60, 
               arrow = arrow(length = unit(0.15, "cm")), color = "black") + # Add the upward arrow
  annotate("text", x = 32.5, y = 2.565, label = "Sample 4\npassing inflection", size = 2.7, hjust = 0.5) +
  theme(text = element_text(size = 9), axis.text = element_text(size = 8))


Figure_8 <- plot_grid(Figure_8A, Figure_8B, Figure_8C, Figure_8D, Figure_8E, Figure_8F, Figure_8G, Figure_8H, Figure_8I,
                      nrow = 3,
                      labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
)

ggsave("Figure_8.pdf", Figure_8, width = 18, height = 16, units = "cm")


