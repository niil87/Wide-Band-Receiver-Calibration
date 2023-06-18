# Wide-Band-Receiver-Calibration
This repo contains files as part of Project work (Advanced Course in Electrical and Information Technology - EITN35) at Lund University in collaboration with Michiel Sandra. We utilize a new (as of 2023) open-source hardware, USRP X410, and perform calibration to estimate the phase contributions from the hardware. Utilizing this information, it is possible to estimate the wireless channel coefficients, thereby making it possible to utilize MIMO-specific optimization processes such as Maximum Ratio Combining (MRC), Water pouring, and other such methods.

The Presentation and Report document are located at [PPT_and_Report](https://github.com/niil87/Wide-Band-Receiver-Calibration/tree/main/PPT_and_Report)

The summary of the source codes (provided in report document Appendix I)

Step1 : Setup Requirement of Tx and Rx antennas (section 3.1 ) <br>
Purpose : Determine the Tx/Rx layout to ensure we get good LoS signal <br>
Source Code : Most of the effort is manual calculations with one useful code [ArrangementCalculator.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/ArrangementCalculator.m)


Step2 : Data Collection ( section 3.2 ) <br>
Purpose : Collect both physical measurements as well as OTA and B2B data required for calibration <br>
Source Code : [USRP Channel Sounder](https://github.com/michielsandra/openucs) (Measurement details are saved in DATA SHEET)


Step3 : Frequency Error Detection and Correction ( section 3.3 ) <br>
Purpose : Resolve frequency errors prior to the calibration process <br>
Source Code : [C1 FrequencyErrorAndCorrection.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C1_FrequencyErrorAndCorrection.m)


Step4 : Stable stationary point selection for Channel Estimation ( section : 3.4 ) <br>
Purpose : Select stable OTA samples to avoid phase jitters due to physical antenna movements <br>
Source Code : [C2 KNN Clustering.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C2_KNN_Clustering.m)


Step5 : Unwrapping error removal using outlier detection ( section : 3.5 ) <br>
Purpose : Catch and rectify phase unwrapping errors <br>
Source Code : [C3 OutlierDetectorAndPhaseCorrector.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C3_OutlierDetectorAndPhaseCorrector.m)


Step6 : LoS Separation ( section : 3.6 ) <br>
Purpose : Use SAGE algorithm to separate out LoS <br>
Source Code : [C4 Sage Implementation.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C4_Sage_Implementation.m)


Step7 : Calibration ( section : 3.8 ) <br>
Purpose : Perform calibration to determine the hardware response <br>
Source Code : [C5 FinalCalibration.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C5_FinalCalibration.m)


Step8 : MUSIC Algorithm for verifying calibration data reliability ( section : 3.9 ) <br>
Purpose : Estimate DoA on estimated channel response using calibration data; use ground truth <br>
Source Code : [C6 MusicImplementation.m](https://github.com/niil87/Wide-Band-Receiver-Calibration/blob/main/C6_MusicImplementation.m)
