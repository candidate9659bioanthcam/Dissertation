from mpi4py import MPI
import numpy as np
import os
import sys
import subprocess

def create_file(name, content):
    f = open(name+".txt", "w")
    f.write(content)
    f.close()

def convert_to_txt(c):
    time1 = [-6123,-5535,-5051,-4507,-4001,-3517,-3044,-2527,-1947,-1504,-1094,-400,0]
    time2 = [-5500, -4906, -4240, -3935, -3461, -3013, -2471, -2045, -1550, -971,-492,-400,0]
    res = ""
    for i, (t1, t2) in zip(range(1, 13,1), zip(time1, time2)):
        res += '{}\t{}\t{}\t{}\n'.format(c[2*i],c[2*i+1],t1,t2)#s f"{c[2*i]}\t{c[2*i+1]}\t{t1}\t{t2}\n" 
    return res

def get_dict(file_name):
    d = dict()
    f = open(file_name)
    for line in f:
        line = line. strip('\n')
        val = line. split(",")
        d[val[1]] = convert_to_txt(val)
    del d["SNP"]
    return d


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

if rank == 0:

    data = get_dict(sys.argv[1])#{i: i for i in range(16)}#

    # determine the size of each sub-task
    ave, res = divmod(len(data), nprocs)
    counts = [ave + 1 if p < res else ave for p in range(nprocs)]

    # determine the starting and ending indices of each sub-task
    starts = [sum(counts[:p]) for p in range(nprocs)]
    ends = [sum(counts[:p+1]) for p in range(nprocs)]

    key_list = list(data.keys())
    data = [{key_list[i]: data[key_list[i]] for i in range(s, e, 1)} for s, e in zip(starts, ends)]
else:
    data = None

data = comm.scatter(data, root=0)

for k, v in data.items():
    create_file(k, v)
    cli_com = "./sr -D {file_name}.txt -G 25 -N 10000 -n 500000 -d 0.001 -F 20 -f 1000 -s 100 -P constant.pop -e 8067 -a -o {file_name}".format(file_name = k)
    cli_com = cli_com.split(" ")
    cli_com[-1] = cli_com[-1].replace("\n", "")
    result = subprocess.check_output(cli_com)
    if result.split(",")[0]==k:
        print(result.replace("\n", ""))
    os.remove(k+".txt")
