# fNIRS-Motion-Artefact

### Motion Artefact Removal in Functional Near-Infrared Spectroscopy Signals based on Robust Estimation

*Abstract: Functional Near-InfraRed Spectroscopy (fNIRS) has gained widespread acceptance as a non-invasive neuroimaging modality for monitoring functional brain activities. FNIRS uses light in the near infra-red spectrum (600-900 nm) to penetrate human brain tissues and estimates the oxygenation conditions based on the proportion of light absorbed. In order to get reliable results, artefacts and noise need to be separated from fNIRS physiological signals. This paper focuses on removing motion-related artefacts. A new motion artefact removal algorithm based on robust parameter estimation is proposed. Results illustrate that the proposed algorithm can outperform the state-of-art algorithms in removing motion artefacts. Moreover, the proposed algorithm is robust in estimating the parameters under different interference conditions.*

- M. Wang and A. Seghouane, "Motion Artefact Removal in Functional Near-infrared Spectroscopy Signals Based on Robust Estimation," *ICASSP 2019 - 2019 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)*, Brighton, United Kingdom, 2019, pp. 1145-1149. [[Full paper](https://github.com/mengmwang/fNIRS-Motion-Artefact/blob/main/Wang2019.pdf)]

### MATLAB Implementation

This repository contains the MATLAB code for the motion artefact removal algorithm used in functional near-infrared spectroscopy signals.

**Data**

`simulated_signal.mat` - An example of a simulated signal with simulated artefacts

`experimental_signal.mat` - An example of an experimental fNIRS signal with simulated artefacts

**Code**

`mar_algo.m` - The algorithm implementation

**Demo**

`TestSimulated.m` - Test the motion artefact removal algorithm using simulated signals and simulated artefacts

`TestExperimental.m` - Test the motion artefact removal algorithm using expperimental fNIRS signals and simulated artefacts


*This is part of my PhD research. Please get in touch if you are interested to know more! 😉*
