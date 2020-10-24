import os
import subprocess

def load_data():
    nodes = ["1n", "2n", "4n", "6n", "8n", "10n", "12n"]
    out = ["10", "100", "1000", "1500", "2500", "5000"]

    files = [i+"1c_"+j+"_nlogn.out" for i in nodes for j in out]

    files = list(filter(lambda x: os.path.isfile(x), files))

    timing = []

    for f in files:
        process = subprocess.Popen(['tail', '-n', '1', f],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)

        stdout, _ = process.communicate()
        timing.append(str(stdout).split(" ")[-2])

    group_res = []
    cores = []
    
    p = len(out)

    for i in range(int(len(files) / p)):
        group_res.append(list(map(lambda x: float(x),timing[i*p:(i+1)*p])))
        if 'n' in files[i*p][:2]:
            cores.append(files[i*p][:1])
        else:
            cores.append(files[i*p][:2])
    
    print(group_res)

    return group_res, cores


if __name__ == "__main__":
    load_data()