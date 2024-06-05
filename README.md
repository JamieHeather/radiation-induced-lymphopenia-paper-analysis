# Radiotherapy induced lymphopenia:
## Photons versus protons paper analysis

### A collaboration between the Cobbold and Hong labs at MGH, 2018-2020*

This repository contains all of the scripts used to generate the data for the *\[TODO manuscript in preparation\]*. 

The input TCR data is stored in one of two places, in different formats. Users have the choice of downloading [the raw data from the immuneACCESS resource]([https://clients.adaptivebiotech.com/pub/heather-2020](https://clients.adaptivebiotech.com/pub/heather-2020)), or the AIRR-seq Community compliant format as produced by running the data through [```immunoseq2airr```](https://github.com/JamieHeather/immunoseq2airr), available from (TODO LINK TO GO HERE).

To conveniently execute all of the required analyses, and thus generate all the plots and intermediate files, follow these steps:

1) Obtain the AIRR-seq format files, either by:
    * Directly downloading them (LINK), placing them in a directory named 'Converted-Data' in the top level of this folder.
    * Download the raw data in `v2' format from Adaptive ([available at this temporary pre-publication immuneACCESS link](https://clients.adaptivebiotech.com/pub/heather-2020)), placing it in the 'Raw-Data' directory (alongside the three CSV files that should be there), and then running the ```get-data.py``` script in the  'Scripts' directory.
      * Note that this might require you to alter the path to your Python installation in the line beginning ```cmd = ...```.
   * If users which to run the final antigenic clustering script, they also need to [download the latest VDJdb release from their GitHub repo](https://github.com/antigenomics/vdjdb-db/), extract the `vdjdb.slim.txt` file and placing it in the `Raw-Data/VDJdb` directory.
   
2) Run the ```run-analysis.py``` script in the Scripts directory, which is a wrapper to run the following scripts in order. Users may run them individually if they wish, but bear in some of them require the output of preceding scripts:
    1) ```lymphopenia-incidence.py```
    2) ```analyse-metadata.py```
    3) ```analyse-tcrs.py```
    4) ```iterative-subsample-tcrs.py```
    5) ```analyse-flow.py```
    6) ```analyse-survival.py```
    7) ```cluster-specificities.py``` 
    
Note that everything should be run inside the Scripts directory; many of the scripts use relative paths, and that's where they expect to be. Most of the scripts also make use of a shared ```functions.py``` script, which must therefore also be present.

Intermediate files produced are stored in `Converted-Data/` along with the AIRR-seq formatted repertoire files. Images produced end up in the Plots directory, with each Script producing its own folder (prefixed with the date). Note that some of the analyses involve randomly subsampling the data, or work from averaged sub-sampled data, and thus will produce plots that differ slightly each time.

#### Dependencies

The following Python modules are required:

* `scipy`
* `numpy`
* `pandas`
* `seaborn`
* `Levenshtein`
* `lifelines`
* `graph_tool`


<sub><sup>&ast; Manuscript delayed due to pandemic, PIs leaving, lead authors moving labs, and life generally getting in the way.</sup></sub>