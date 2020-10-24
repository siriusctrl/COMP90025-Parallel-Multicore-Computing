import os
import subprocess

if __name__ == "__main__":
    nodes = ["1n", "2n", "4n", "6n", "8n", "10n", "12n"]
    out = ["10", "100", "1000", "5000"]

    files = [i+"1c_"+j+"_nlogn.out" for i in nodes for j in out]

    files = list(filter(lambda x: os.path.isfile(x), files))

    timing = []

    for f in files:
        process = subprocess.Popen(['tail', '-n', '1', f],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)

        stdout, _ = process.communicate()
        timing.append(str(stdout).split(" ")[-2])
    
    for i in range(len(timing)):
        print(files[i], timing[i])

    group_res = []
    cores = []
    
    for i in range(int(len(files) / 4)):
        group_res.append(timing[i*4:(i+1)*4])
        if 'n' in files[i*4][:2]:
            print(files[i*4][:1])
        else:
            print(files[i*4][:2])

        