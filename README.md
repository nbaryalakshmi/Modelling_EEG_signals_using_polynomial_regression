# Modelling_EEG_signals_using_polynomial_regression
Modelling EEG signals using polynomial regression
----------------------------
The objective of this project is to identify the best regression model (from a potential set of nonlinear regression models) that can adequately describe the brain activity elicited during guided meditation. The information is presumably gathered during a neuroscience experiment in which a person is encouraged to practise guided meditation while receiving instructions over speakers. The modulation of neural activity in two different brain regions during the meditation is of interest to the researchers. 

Particularly, electroencephalography (EEG) is used to gauge the brain activity in the right auditory cortex and prefrontal cortex. Whereas area (2) is linked to executive function, planning, and consciousness, area (1) is responsible for processing auditory experiences. The prefrontal cortex is projected to have a nonlinear association, in contrast to the auditory cortex, which the researchers believe is linearly related to the audio signal (i.e., the voice of the mediation guide). Using nonlinear regression modelling, this report examines these correlations.

The two distinct Excel files contain the "simulated" EEG time-series data and the sound signal. The sound signal y is contained in the y.csv file, while the X.csv file comprises the EEG signals x1 and x2 that were measured from the prefrontal and auditory cortices. The sampling times for all three signals are listed in seconds in the file time.csv. A total of 2 minutes' worth of signal data were gathered at a sampling rate of 20 Hz. Due to distortions during recording, all signals are subject to additive noise with unknown variance (assumed to be independent and identically distributed ("i.i.d") Gaussian with zero-mean).

For coding and resolving this problem, R programming language is utilised.

Full report is provided in Report.pdf
