import os,sys

cwd = os.getcwd()

with open("/home/liusynevolab/curr_ev%s_exps" % sys.argv[1], "r") as f:
    lines = f.readlines()
with open("/home/liusynevolab/curr_ev%s_exps" % sys.argv[1], "w") as f:
    for line in lines:
        if line.strip("\n") != cwd:
            f.write(line)
