# jobwatcher  - submit x number of jobs at once, when any of the jobs finishes, submit a new one

# todo: check for error , resubmit?

    
from subprocess import Popen, check_call, call
import multiprocessing as mp

samples = subprocess.check_output("cd ~/projects/fap/expands/; ls *combined | cut -d. -f1", shell=True)
samples = samples.split("\n")
cmd_list = ["Rscript expands.R " [sample] for sample in samples]


def worker(cmd):
    #cmd = "sh ./test.sh"
    print "a;sldkjfal;sdjkf"
    p = Popen(cmd, shell=True)
    p.wait()
    
    #check_call("exit 1", shell=True)
    #print "return code " + str(p.returncode)
def main():

    # start processes
    pool = mp.Pool( processes = 2 );
    results =[pool.apply_async(worker, [cmd]) for cmd in cmd_list];
    ans = [res.get() for res in results];
    type(ans)
    for an in ans:
        print type(an)
        print an
    
if __name__=="__main__":
    main()

        