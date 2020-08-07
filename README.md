# Compound-RF-summation-data-and-model-functions
This repository contains data and model functions utilized in:   Schmidtmann, G., Kingdom, F. A. A., &amp; Loffler, G. (2019). The processing of compound radial frequency patterns. Vision Research, 161, 63–74.

Create three additional folders (Data, Figures, Models).
The key script is called "Analyse_Data.m". Please read the documentation in the script. The function "RF_compound_data.m" contains the raw data. The “Models” folder contains all model simulations, superimposed on individual data for all subjects, as Matlab figures (.fig). The models are divided into fixed - and matched attention window scenarios (FAW, MAW). A separate folder shows the model simulations for the condition where the transducer exponents are fixed (Fixed Transducers). We also provide model simulations for a single channel model. Please see Schmidtmann, Kingdom and Loffler (2019) for details.

The models described in Schmidtmann, Kingdom & Loffler (2019) employ additional functions from the Palamedes Toolbox  (Prins & Kingdom, 2018) to determine whether the data from a 5-PF (psychometric function) summation square experiment, for the detection of two stimuli in the target interval, accords more with probability summation (PS) or additive summation (AS) under the assumptions of signal-detection-theory (SDT) and assuming that the observer is monitoring both channels sensitive to the two stimuli. 

Palamedes function (PF) fitting routines:

PAL_SDT_Summ_MultiplePFML_Fit

PAL_SDT_Summ_MultiplePFML_BootstrapParametric

PAL_SDT_Summ_MultiplePFML_GoodnessOfFit

Palamedes SDT PS (probability) and AS (additive)

PAL_SDT_PS_uneqSLtoPC

PAL_SDT_AS_uneqSLtoPC

Palamedes PF fitting routines:

PAL_PFML_Fit

PAL_Logistic

Prins, N. & Kingdom, F. A. A. (2018) Applying the Model-Comparison Approach to Test Specific Research Hypotheses in Psychophysical Research Using the Palamedes Toolbox. Frontiers in Psychology, 9:1250.
