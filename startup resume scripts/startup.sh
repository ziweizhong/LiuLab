#!/bin/bash

# Initializes the miniconda command. Without this, unable to activate conda environments upon reboot
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false '/home/liusynevolab/miniconda3/bin/conda' shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/liusynevolab/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/liusynevolab/miniconda3/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        export PATH="/home/liusynevolab/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<

# evolver 1.0 experiments to restart
# experiments that were running prior to unexpected shutdown/restart are stored in the input file
input='/home/liusynevolab/curr_ev1_exps'
while IFS= read -r line
do
	echo $line
	echo "cd $line;eval $(conda shell.bash hook);conda activate py27;python ~/scripts/ev1.0/main_eVOLVER.py --emergency-restart;exec bash" > run.sh
	chmod +x run.sh
	gnome-terminal -- ~/run.sh
	rm -f run.sh
done < $input


# evolver 2.0 experiments to restart
input='/home/liusynevolab/curr_ev2_exps'
while IFS= read -r line
do
	echo $line
	echo "cd $line;eval $(conda shell.bash hook);conda activate py365;python ~/scripts/ev2.0/eVOLVER.py --emergency-restart;exec bash" > run.sh
	chmod +x run.sh
	gnome-terminal -- ~/run.sh
	rm -f run.sh
done < $input
