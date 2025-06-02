# ROCK1-docking-analysis
Repository for analysis performed in a bachelor thesis called "Large-Scale Analysis of Ligand-Target Interactions of Rho-associated proteinkinase 1".

## Steps:
* Structure of folders is described below
* With downloaded materials from Zenodo, each step can be ran on its own

1. Download materials from Zenodo: https://doi.org/10.5281/zenodo.15458756
    * Copy folder _materials_, so that the folder placement is the same as defined bellow

1. Create Conda Virtual Environment (for both Linux and Windows)
    * Instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
    * Linux: ```environment_linenv.yml```, further on as "linenv"
    * Windows: ```environment_venv.yml```, further on as "venv"

1. ```align.ipynb``` (linenv)
    * Align ligands to a ROCK1 structure with PDB-ID 2ETR

1. ```vina_docking.ipynb``` (linenv)
    * Docking done using Dockstring [^1], with docking algorithm being AutoDock Vina (distribution 1.2.7 [^2])
    [^1]: Miguel García-Ortegón, Gregor N. C. Simm, Austin J. Tripp, José Miguel Hernández-Lobato, Andreas Bender, and Sergio Bacallado. Dockstring: Easy molecular docking yields better benchmarks for ligand design. _Journal of Chemical Information and Modeling_, 62(15):3486–3502, Aug 2022.
    [^2]: Oleg Trott and Arthur J. Olson. Autodock vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. _Journal of Computational Chemistry_, 31(2):455–461, Jan 2010.
    * Change Vina distribution:
        * save the desired distribution to: linenv\lib\python3.10\site-packages\dockstring\resources\bin
        * change the version: linenv\lib\python3.10\site-packages\dockstring\utils.py, line 70

1. ```vina.ipynb``` (linenv)
    * Read poses docked by AutoDock Vina 1.2.7, calculate RMS between them and the co-crystallized ligands

1. ```moe.ipynb``` (linenv)
    * Read poses docked by MOE's Dock [^3], calculate RMS between them and the co-crystallized ligands
    [^3]: 2024.0601 Chemical Computing Group ULC. _Molecular Operating Environment (MOE)_. 910-1010 Sherbrooke St. W., Montreal, QC H3A 2R7, 2025.
    * \+ show graph of best RMSDs by Vina and MOE

1. ```glide.ipynb``` (linenv)
    * Read poses docked by Glide HTVS, SP and XP [^4], calculate RMS between them and the co-crystallized ligands
    [^4]:  Richard A. Friesner, Jay L. Banks, Robert B. Murphy, Thomas A. Halgren, Jasna J. Klicic, Daniel T. Mainz, Matthew P. Repasky, Eric H. Knoll, Mee Shelley, Jason K. Perry, David E. Shaw, Perry Francis, and Peter S. Shenkin. Glide: a new approach for rapid, accurate docking and scoring. 1. method and assessment of docking accuracy. _Journal of Medicinal Chemistry_, 47(7):1739–1749, Mar 2004.
    * \+ show graph of best RMSDs by all of them and MOE

1. ```ifs_plip_aligned.ipynb``` (venv), ```ifs_prolif_aligned.ipynb``` (venv)
    * Calculate interaction fingerprints (IFs) on the co-crystallized dataset using PLIP [^5]/ProLIF [^6]
    [^5]:  Sebastian Salentin, Sven Schreiber, V. Joachim Haupt, Melissa F. Adasme, and Michael Schroeder. Plip: fully automated protein–ligand interaction profiler. _Nucleic Acids Research_, 43(W1):W443-W447, 04 2015.
    [^6]: Cédric Bouysset and Sébastien Fiorucci. Prolif: a library to encode molecular interactions as fingerprints. _Journal of Cheminformatics_, 13(1):72, Sep 2021.

1. ```rescore_plip.ipynb``` (venv), ```rescore_prolif.ipynb``` (venv)
    * Optimize weights of interactions found in co-crystallized and docked poses using PLIP/ProLIF

1. ```get_fingerprints_glide.ipynb``` (venv), ```get_fingerprints_vina.ipynb``` (venv)
    * Calculate and save IFs of docked poses by Glide/Vina (both actives and decoys)

1. ```auc_and_histograms_glide.ipynb```, ```auc_and_histograms_vina.ipynb``` (venv)
    * Find the best ratio between docking (original) score and new (IFs) score

1. ```get_graphs.ipynb``` (venv)
    * Plot evaluation graphs:
        * Vina vs Glide SP: Scoring evaluation
        * Accuracies: original, PLIP and ProLIF scores
        * IFs of co-crystallized dataset (PLIP and ProLIF)

1. ```roc_all.ipynb``` (venv)
    * Plot ROC curve, calculate AUC score and Early Enrichment Factor for Vina and Glide (original and rescored)


## Structure (Folders and Files):
### GIT - subfolders containing Jupyter Notebooks for analyses
* _materials_ ... folder containing materials from Zenodo
* _env\_files_: ```environment_linenv.yml```, ```environment_venv.yml```
* _align_: ```align.ipynb```
* _docking_: ```vina_docking.ipynb```, ```vina.ipynb```, ```moe.ipynb``` , ```glide.ipynb```
* _rescore_: ```ifs_plip_aligned.ipynb```, ```ifs_prolif_aligned.ipynb```, ```rescore_plip.ipynb```, ```rescore_prolif.ipynb``` \+ ```descriptors.py```, ```pdb_reading.py``` (two helper scripts)
* _actives\_decoys_: ```get_fingerprints_glide.ipynb```, ```get_fingerprints_vina.ipynb```,  ```auc_and_histograms_glide.ipynb```, ```auc_and_histograms_vina.ipynb```
* _graphs_: ```get_graphs.ipynb```, ```roc_all.ipynb```

### Materials - subfolders with data (downloaded from Zenodo)
* ```aligned_molecules.sdf```, ```2etr.pdb```
* _ligands\_pdb_, _proteins\_pdb_ ... structures downloaded from the Protein Data Bank
* _docking_
    * _vina_ 
        * _vina\_versions_ ... distributions of AutoDock Vina 1.1.2 [^7], 1.2.7 and QuickVina 2.1 [^8]
        [^7]:  Jerome Eberhardt, Diogo Santos-Martins, Andreas F. Tillack, and Stefano Forli. Autodock vina 1.2.0: New docking methods, expanded force field, and python bindings. _Journal of Chemical Information and Modeling_, 61(8):3891–3898, Aug 2021.
        [^8]: Stephanus Daniel Alhossary, Amrand Handoko, Yuguang Mu, and Chee-Keong Kwoh. Fast, accurate, and reliable molecular docking with quickvina 2. _Bioinformatics_, 31(13):2214–2216, Jul 2015.
        * ```vina127_docked.sdf``` ... docked molecules by AutoDock Vina 1.2.7
    * _moe_
        * ```moe_docked.sdf``` ... docked molecules by MOE's Dock
    * _glide_
        * ```glide_htvs_docked.sdf```, ```glide_sp_docked.sdf```, ```glide_xp_docked.sdf``` ... docked molecules by Glide HTVS, SP and XP
* _rescore_
    * ```ifs_aligned_plip.csv```, ```ifs_aligned_prolif.csv``` ... interaction fingerprints (IFs) on co-crystallized ligands
* _actives\_decoys_
    * ```glide_actives_docked.sdf```, ```glide_decoys_docked.sdf``` ... active and decoy molecules docked using Glide SP
    * ```vina_actives_docked.sdf```, ```vina_decoys_docked.sdf``` ... active and decoy molecules docked using AutoDock Vina 1.2.7
    * _ifs\_glide_, _ifs\_vina_ ... folders with subfolders (_actives_, _decoys_) containing IFs of docked poses
* _graphs_
    * ```vina.csv```, ```moe.csv```, ```glide_xp.csv```, ```glide_sp.csv```, ```glide_htvs.csv``` ... docking results (RMSDs)
    * ```ifs_aligned_plip.csv```, ```ifs_aligned_prolif.csv``` ... interaction fingerprints (IFs) on co-crystallized ligands
    * ```accuracy_original.txt```, ```accuracy_plip.txt```, ```accuracy_prolif.txt``` ... accuracies of docking scores and IFs scores (has to be copied into a .txt file manually from outputs of ```rescore_plip.ipynb``` and ```rescore_prolif.ipynb```)
