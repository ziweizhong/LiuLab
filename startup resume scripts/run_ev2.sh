#!/bin/bash

# Enables interrupts
trap '' INT

# Saves the current working directory to file
savecwd () (
    trap - INT
    pwd >> /home/liusynevolab/curr_ev2_exps
)

# Runs the python script for eVOLVER 1.0
runpy () (
    trap - INT
    python /home/liusynevolab/scripts/ev2.0/eVOLVER.py
)

# If the python script for eVOLVER terminated correctly (via user input or other method)
#  removes the current working directory from running experiments file
removecwd () (
    trap - INT
    python /home/liusynevolab/scripts/remove_cwd.py 2
)

savecwd

runpy

removecwd
