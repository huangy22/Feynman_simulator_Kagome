''' This is a file to define the class of all jobs,
    You have to modify this file if you want to add new type of jobs.'''
import sys, os, random, re, copy

def get_current_PID(KeyWord):
    workspace=os.path.abspath(".")
    filelist=sorted([int(e.split('_')[0]) for e in os.listdir(workspace) if (KeyWord in e) and e[0] is not '_'])
    if len(filelist)==0:
        NextPID=0
    else:
        NextPID=filelist[-1]+1
    return filelist, NextPID

#base class of all jobs
class Job:
    '''Class Job is the base class of all job objects.
       You have to add new subclass of Job if you want to
       create new type of job.'''
    def __init__(self, para):
        # interesting point here, python will actually call the 
        #__check_parameters__ of subclass,
        # So you don't have to call __check_parameters__ 
        #in subclass initialization again!
        if self.__check_parameters__(para) is False:
            print "Something is wrong with the inlist! Abandon!"
            sys.exit()
        self.control=para.pop("Control")
        self.para=para.pop("Job")
        self.pid = []

    def to_dict(self):
        '''output the corresponding string of the job class'''
        pid=self.pid.pop(0)
        self.para["PID"] = pid
        self.__set_model_specific__()
        return pid, {"Job": copy.deepcopy(self.para)}

    def __check_parameters__(self, para):
        if para["Control"]["__Execute"] is "":
            print "Please specify the executive file name!"
            return False

    def __set_model_specific__(self):
        pass

class JobMonteCarlo(Job):
    '''job subclass for monte carlo jobs'''
    def __init__(self, para):
        Job.__init__(self, para)
        self.control["__KeepCPUBusy"]=True
        self.control["__Type"]="MC"

        #search folder for old jobs, the new pid=largest old pid+1
        PIDList, NextPID=get_current_PID("statis")
        if len(PIDList) is not 0:
            self.pid=PIDList[:self.control["__Duplicate"]]
        else:
            self.pid=range(NextPID, NextPID+self.control["__Duplicate"])

    def __check_parameters__(self, para):
        if Job.__check_parameters__(self, para) is False:
            return False

    def to_dict(self):
        pid, Dict=Job.to_dict(self)
        #set Seed here so that each job has it own rng seed
        Dict["Job"]["Seed"] = int(random.random()*2**30)
        Timer=Dict["Job"]["Timer"]
        #don't let MC process output at the same time
        for e in Timer.keys():
            Timer[e]=int(Timer[e]*random.uniform(0.8, 1.2))
        return pid, Dict
