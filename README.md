# CLIApet

This package is the CLIA Performance Evaluation Toolkit: An R Library For Supporting Gene Expression Based CLIA Validation Studies. Specifically, it contains R functions that were developed in support of the CLIA validation of the Radiation Sensitivity Index assay. This assay is derived from an Affymetrix (Thermo-Fisher) gene expression array (the HG-U133 Plus GeneChip). 

To install the library in R, use

```
devtools::install_github("steveneschrich/CLIApet")
```

The RSI CLIA validation consisted of a number of different experiments performed over the course of a year. During this time, various visualizations become more useful to have for each experiment. Things like scatter plots with density, PCA plots, and others. None of them are particularly complicated but all of them took time to develop. And of course, with many different experiments it was clear that functions had to be developed so the figures would be consist across the project. Finally, it was time to develop the final report summarizing all experiments as a bookdown project. That is a whole other saga. However, it became clear that a stand-alone package would provide the needed functionality and, importantly, provide a reusable set of functions for future projects. Thus, the CLIApet library was developed.

CLIApet has undergone one or two code reorganizations (refactoring) but still remains a bit too complicated for my taste. I believe that tools like ggplot2 are excellent but are multi-tools, whereas when I'm doing a particular task I'd prefer some "nice-looking" defaults. For instance, ggpubr provides this type of activity. Interestingly, tex/latex provides this heavy-handed stylistic intervention which I quite appreciate due to my lack of innate stylistic ability. This library essentially provides a number of these defaults for various activities, inspired by that thought.
