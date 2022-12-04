from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import atexit, time

# A simple timing decorator adapted from TensorMol.
PRTIMER = {}
PRSTARTTIME = time.time()
def PrintPRTIMER():
    print("=======    Accumulated Time Information    =======")
    print("Category ||| Time Per Call ||| Total Elapsed     ")
    for key in PRTIMER.keys():
        if (PRTIMER[key][1]>0):
            print(key+" ||| {:.5f} ||| {:.5f} ".format(PRTIMER[key][0]/(PRTIMER[key][1]),PRTIMER[key][0]))
def PRTiming(nm_="Obs"):
    if (not nm_ in PRTIMER.keys()):
        PRTIMER[nm_] = [0.,0]
    def wrap(f):
        def wf(*args,**kwargs):
            t0 = time.time()
            output = f(*args,**kwargs)
            PRTIMER[nm_][0] += time.time()-t0
            PRTIMER[nm_][1] += 1
            return output
        print("PRTimed "+nm_+str(PRTIMER[nm_]))
        return wf
    return wrap
@atexit.register
def exitProgram():
    PrintPRTIMER()
    print("Total Time in seconds: {:.5f}".format(time.time()-PRSTARTTIME))
