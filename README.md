# BINF6999-Summer-Project---Modeling-CIN-Cancer-Evolution

This project aims to model chromosome counts among evolutionarily distinct cancer populations using a Markov chain and grid search. Attached are demonstration files for astrocytomas, pancreatic adenocarcinomas, and prostrate adenocarcinomas. As running a grid search takes some time, loading the supplied RData environments alongside their matching "Mark_Model" file is recommended. 

Scripting was performed in R 4.5.1, with the dplyr 1.2.1, tidyr 1.3.2, stringr 1.6.0, ggplot2 4.0.3, forreach 1.5.2, and doparallel 1.0.17 packages. R was used due to its broad support of statistical tools and analysis software, on top of prior familiarity. dplyr was used for statistical analysis. tidyr was used to facilitate piping and make scripting more legible. stringr was used for text extraction of chromosome numbers from database files. ggplot2 was used for plot creation and  data display. forreach and doparallel were used to allow multithreading in the grid search.


Each "Mark_Model" R file contains information on how the model can be adapted to explore other cancer populations. If doing so, exploring chromosome count data using one of the "Kary_Comp" scripts as a template is recommended before setting grid search parameters. 

When parameter testing, the "Maximum Log-Likelihood" is a useful metric for comparing models using the same dataset. The more positive the Maximum Log-Likelihood, the better the model fits supplied data.
