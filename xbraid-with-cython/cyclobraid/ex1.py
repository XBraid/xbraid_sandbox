def step(u, tstart, tstop):
    print("py step!")
    return 1./(1. + tstop-tstart)*u

def init():
    return 0.6
