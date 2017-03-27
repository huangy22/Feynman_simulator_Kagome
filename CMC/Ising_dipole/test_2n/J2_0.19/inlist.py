# monte carlo job defintion
MonteCarlo={
"Control": {
    "__Execute" : "./simulator.exe",
    "__Duplicate" :  1,
    "__IsCluster" : False, 
    "__AutoRun" : True,
    "__PBSCommand": "#PBS -l walltime=200:00:00"
    },

"Job": {
    "Beta": 1.0,
    "ExternalField": 0.0,

    #2D lattice
    "Dim": 2,
    "LatticeName": "Kagome", "NSublat": 3,
    "L": [16,16],

    "Sample" : 100000000,  ##0.8 min for 1000000(*1000) Samples in MC
    "Sweep" : 100, "Toss" : 100000,
    "Timer":{
        "PrinterTimer": 100,
        "DiskWriterTimer": 100,},
    }  
}

import job_class as job
'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
TO_DO = []
TO_DO.append(job.JobMonteCarlo(MonteCarlo))
CPU = 4
SLEEP = 1    #check job status for every SLEEP seconds
