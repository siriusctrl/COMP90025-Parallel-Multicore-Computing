def test(filename="easy"):
    p = "parallel_" + filename

    pxy = 1
    pgap = 1

    with open(p, "r") as f:
        for i in range(5):
            f.readline()
        p_string = f.readline()
        o_string = f.readline()
    
    penalty = 0

    for i in range(len(p_string)):
        if p_string[i] == o_string[i]:
            pass
        elif p_string[i] == "_" or o_string[i] == "_":
            penalty += pgap
        else:
            penalty += pxy

    print(penalty)

test()