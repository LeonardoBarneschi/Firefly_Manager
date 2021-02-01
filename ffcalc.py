#!/usr/bin/env python3

import errno
import os
import schedule
import subprocess
import sys
import time

def check_done(out, check):

    if os.path.exists(out):
        with open(out, 'r') as o:
            for line in o:
                if check in line:
                    return True
    else:
        pass

def check_done_wrap():

    done = check_done(ffout, check)
    if done:
        print(f'The file {ffout} has been checked succesfully')
        return schedule.CancelJob 


def fflaunch(ffinp):
    launcher = os.environ.get('FFLAUNCHER')
    nprocs = 16
    queue = 'normal'
    cmd = f'{launcher} {ffinp} {nprocs} {queue}'
    return cmd 

if __name__ == "__main__":
    
    check = 'TERMINATED NORMALLY' 

    # PSEUDO-CODE 
    # List of directories

    paths = ['/home/Leo/FIREFLY/auto-test/rs1/',
             '/home/Leo/FIREFLY/auto-test/rs2/',
             '/home/Leo/FIREFLY/auto-test/rs3/']

    for path in paths:

        os.chdir(path)
        print(f'current path is: {os.getcwd()}')
        time.sleep(2)

        if len([f.endswith('inp') for f in os.listdir()]) == 1:
            
            ffinp = [f for f in os.listdir() if f.endswith('inp')][0]
            ffout = ffinp.split('.')[0]+'.out'
            print(f'FF input file is {ffinp}')
            time.sleep(2)
            
            cmd = fflaunch(ffinp)
            print(f'FF launcher command is: {cmd}')
            time.sleep(2)

            subprocess.run([cmd], shell=True)
            schedule.every(5).seconds.do(check_done_wrap)

            while True:
                schedule.run_pending()
                if not schedule.jobs:
                    break
                time.sleep(2)

        else:
            continue

        print()

