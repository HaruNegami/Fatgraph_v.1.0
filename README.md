# Fatgraph_v.1.0
#################################
##### First, be sure to write the path of the working directory to /configuration/classFilePath.py
#################################

First settings:

You can change parameters at configuration files in configuration folder

To change parameters for thresholds of potentials, see configuration.py

To change directories to save results, see classFilePath.py

Usage:

    python main.py model_option test_option update_option

Usage example:

    python main.py chain test difference


model_option:

    chain : calculate fatgraph invariants for each chain
    
    whole : calculate fatgraph invariants of whole structure

test_option:

    test       : calculate fatgraph invariants of file at ./test_dssp/
    
    production : calculate fatgraph invariants of data at ./dssp/

update_option:

    overwrite  : calculate all files
    
    difference : calculate invariants iff not yet calculated

