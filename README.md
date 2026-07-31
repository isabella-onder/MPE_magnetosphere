### Summer Project at the Max Planck Institute for Extraterrestrial Physics, Garching ###

"Correlating high-energy background between SRG/eROSITA and XMM-Newton/EPIC in the magnetosphere"

**Supervised by Dr M. Freyberg -- part of the High Energy Group.**


Data reduction and analysis project correlating background events detected by the X-ray telescope eROSITA onboard the satellite SRG and EPIC onboard XMM-Newton, to quantitatively contrast the magnetosphere in their environments (near-Earth and L2 respectively). I coded a pipeline algorithm to flag and sort background events from eROSITA's raw data, wherein >600 flaring events were detected - time-series were then convoluted with EPIC data to match peaks and hence estimate the velocity of a hypothesised source travelling from the Sun radially outwards. 

This project was important in view of launching ESA's new high-energy WFI NewAthena satellite, to contribute to the determination of the most optimal orbit for observations.

My results were presented at the SPIE conference on Astronomical Telescopes + Instrumentation, July 5-10 2026, and the corresponding paper is in the proceedings: 14146-293.



The repository is organised as follows:
- *ds52/* contains the code I ran remotely on a server, directly interacting with the raw eROsita data and processed XMM data. *original_codes/* are simple scripts to visualise data I wrote at the start of my project, whereas *pipeline/* contains the files ran in succession to respectively:
    0. flag and sort flare detections in eROsita data
    1. use Bayesian Blocks to delimitate flares and extract timestamps
    2. retrieve the corresponding data at XMM as well as spatial coordinates
    3. convolute the time-series to find statistically significant time delays and correlation
- *data_analysis/* contains the code used written and used locally, including time translators to go from XMM to eROSITA to UTC systems for instance, as well as drafts for the final code
- *plotting_code* and *figures* contain the essential visualisation of the data used in my paper. 


