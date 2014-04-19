# TODO: Move the management of timer to different module.
import time;

def create_timers():
    return [];

def tic(timers):
    timers.append(time.time());
    
def toc(timers):
    timers[-1] = time.time() - timers[-1];

def print_timers(header, timers):
    timers_with_sum = timers;
    timers_with_sum.append(sum(timers));
    print "%s:" % (header) + ",".join(map(lambda f: "%.2f" % (f), timers_with_sum));
    
