#!/usr/bin/env python

#from __future__ import with_statement
import pypar
import corrcperm
import sys
import time
import cPickle
#import cProfile
import random
import getopt
import math


#this code uses mpi to calculate the pvalues using the nstatistic.
#a vector of nstatistics for each permutation of the genetic data is generated
# by distributing the work for each permutation to different nodes.
#the nstatistics are then collected at the server and combined into a vector of pvalues

mu = 0.0
std = 2.0
numnodes = pypar.size()

class task:
    #these are the default run parameters.
    # 'Hyperdip_n.pck' and 'Tel_n.pck': The Hyperdip and Tel childhood leukemia data provided by St. Jude Children's research hospital.
    simulatedata = False
    file1 = 'Hyperdip_n.pck'
    file2 = 'Tel_n.pck'
    outfilename = 'pvals'
    numthreads = 1
    kern = 2
    seed = 12345
    logprefix = ''
    genes = 1000
    columns = 80
    perms = 1000
    groupsize = columns/8

    fisher = True                       # Fisher's transformation
    pvaluethreshold = -1

    def __init__(self, argv):
        lastincluded = False
        for n, arg in enumerate(argv):#this is a hack so that we don't get the mpi arguments
            if arg == 'last': 
                lastincluded = True
                argv = argv[1:n]
        if not lastincluded:
            abnormalexit("please end the arguments to the script with 'last'")

        try:
            opts, args = getopt.getopt(argv, 't' ,['genes=', 'columns=', 'groupsize=', 'permutations=', 'kernel=', 'numthreads=','seed=', 'file1=', 'file2=', 'outfile=', 'logprefix=', 'threshold=', 'fisher', 'nofisher'])
        except getopt.GetoptError, err:
            abnormalexit(str(err))
        
        for o, a in opts:
            if o == '-t':
                self.simulatedata = True
            elif o == '--genes':
                self.genes = eval(a)
            elif o == '--columns':
                self.columns = eval(a)
            elif o == '--groupsize':
                self.groupsize = eval(a)
            elif o == '--permutations':
                self.perms = eval(a)
            elif o == '--kernel':
                self.kern == eval(a)
            elif o == '--numthreads':
                self.numthreads = eval(a)
            elif o == '--seed':
                self.seed = eval(a)
            elif o == '--file1':
                self.file1 = a
            elif o == '--file2':
                self.file2 = a
            elif o == '--outfile':
                self.outfilename = a
            elif o == '--logprefix':
                self.logprefix = a
            elif o == '--fisher':
                self.fisher = True
            elif o == '--nofisher':
                self.fisher = False 
            elif o == '--threshold':
                self.pvaluethreshold = eval(a)
            else:
                abnormalexit('unrecognized argument: '+o)

        if self.pvaluethreshold == -1:
            self.pvaluethreshold = self.perms
        

#        q, r = divmod(self.columns,self.groupsize)
#        if r == 0:
#            self.groups = q
#        else:
#            self.groups = q+1

        if self.groupsize > self.columns:
            abnormalexit('groupsize greater than number of columns')

        self.initdata()

    def initdata(self):
        random.seed(self.seed)#if this isn't done, python uses the current system time for the seed, or the os random source.  Since each mpi node generates simulated data separately, this must be done so that they generate the same data.

        if self.simulatedata == True:
            #set up the program with generated data
            self.data = corrcperm.corrcperm(self.genes, self.columns, self.groupsize, self.fisher)
            self.data.makerandomdata(self.seed, mu, std)
        else:
            #use real data
            #with open(self.file1, 'rb') as f1:
            #    condition1 = cPickle.load(f1)
            f1 = open(self.file1, 'rb')
            condition1 = cPickle.load(f1)
            #with open(self.file2, 'rb') as f2:
            #    condition2 = cPickle.load(f2)
            f2 = open(self.file2, 'rb')
            condition2 = cPickle.load(f2)
            if condition1.shape != condition2.shape:
                sys.exit('Error: conditions are not the same size')
            self.genes = condition1.shape[0]
            self.columns = condition1.shape[1]
            self.data = corrcperm.corrcperm(self.genes, self.columns, self.groupsize, self.fisher)
            for rownum, row in enumerate(condition1):
                for col, value in enumerate(row):
                    self.data.data_set(rownum, col, value)
            for rownum, row in enumerate(condition2):
                for col, value in enumerate(row):
                    self.data.data_set(rownum, col+self.columns, value)

    def dojob(self, job):
        """this is where the actual computation in the mpi nodes occurs"""
        vec, invalidated = job
        
        for i in invalidated:
            self.data.ignoregene(i)

        setpermutevec(self.data, vec)
        self.data.rearrangeall()
        self.data.threadedallNstats(self.kern, self.numthreads)
        #return the data to its original arrangement.
        #this is not strictly necessary, but it makes the program deterministic
        #(and otherwise, the data will be completely scrambled at the end)
        undorearrange(self.data,vec)

        return [self.data.getNstat(gene) for gene in xrange(self.data.genes)]

    def printruninfo(self):
        print 'genes',self.genes
        print 'columns', self.columns
        print 'perms', self.perms
        print 'kern', self.kern
        print 'numthreads', self.numthreads
        print 'groupsize', self.groupsize
#        print 'number of groups', self.groups
        print 'numnodes', numnodes
        print 'threshold', self.pvaluethreshold
        print 'fisher', self.fisher


def setpermutevec(data, vec):
    for i, v in enumerate(vec):
        data.permutevecset(i, v)

def undorearrange(data,vec):
    for i, v in enumerate(vec):
        data.permutevecset(v,i)
    data.rearrangeall()

def sendtoall(msg):
    for proc in xrange(1, numnodes):
        pypar.send(msg, proc, tag=OUT)

def abnormalexit(reason):
    """this tells each worker node to exit, then kills the server process.
       this should only be called by the server node"""
    print 'abnormal exit'
    print reason
    sendtoall(('Die', 0))
    pypar.barrier()
    pypar.finalize()
    sys.exit(2)



############

rank = pypar.rank()

procname = pypar.get_processor_name()

#these are message tags
#"OUT" messages go from the server to workers
#"RETURN" messages go from the workers to the server
#the value of these variables should remain constant!
OUT = 0
RETURN = 1


 

class invalidatedgenes:
    def __init__(self, numworkers, numgenes):
        self.workerindices = []
        self.invalidatedlist = []
        self.invalidated = {}
        for i in range(0, numgenes):
            self.invalidated[i] = False
        for i in range(0, numworkers+1):
            self.workerindices.append(0)  #every worker's index starts at zero
    
    def invalidategene(self, genenum):
        if self.invalidated[genenum] != True:
            self.invalidated[genenum] = True
            self.invalidatedlist.append(genenum)

    def getnewinvalidatedgenes(self, workernum):
        previousindex = self.workerindices[workernum]
        self.workerindices[workernum] = len(self.invalidatedlist)
        return self.invalidatedlist[previousindex:] #return every invalidated gene from the index to the end

class serverdata:
    def __init__(self, t):#t is of type task
        self.opts = t
        self.invalidatedgeneset = invalidatedgenes(numnodes-1, self.opts.genes)

class jobgenerator:
    """generates jobs for dojob"""
    jobindex = 0
    def __init__(self, server):
        self.numjobs = server.opts.perms
        self.server = server
        self.opts = server.opts
    def __iter__(self):
        return self
    def hasnext(self, workernum):
        """are there new jobs available?"""
        return self.jobindex < self.numjobs
    def next(self, workernum):
        """get the next job
           in this case a job is a permutation vector
        """
        if self.jobindex < self.numjobs:
            self.jobindex += 1
            job = range(self.opts.columns*2)
            random.shuffle(job)
            print self.jobindex
            invalidated = self.server.invalidatedgeneset.getnewinvalidatedgenes(workernum)
            #print invalidated
            return (self.jobindex, (job, invalidated))
        else:
            raise StopIteration

class resultcollector:
    """collects the responses from the worker nodes and combines the results from them at the server"""
    def __init__(self, server):#server is of type serverdata
        self.server = server
        self.nstats = server.opts.dojob((range(server.opts.columns*2), []))#calculate the unpermuted nstatistics
        self.t = server.opts
        self.pvals = [0.0 for i in xrange(server.opts.genes)]
    def collect(self, response):
        jobnum, result = response
        for i, v in enumerate(result):
            if v >= self.nstats[i]:
                self.pvals[i] += 1
            if self.pvals[i] > self.t.pvaluethreshold:
                self.server.invalidatedgeneset.invalidategene(i)

    def finish(self):
        """call this when all the results are in"""
        #print self.pvals
        self.pvals = [x/self.t.perms for x in self.pvals]

        #write the results out to a file
        outfile = open(self.t.outfilename, 'w')
        cPickle.dump(self.pvals, outfile)
        outfile.close()
    
        
    
######################################################
def main():
    #--------------------#
    # server code
    #--------------------#
    if rank == 0:
        print 'server running on ', procname

        opts = task(sys.argv)

        opts.printruninfo()

        sendtoall(('Start', sys.argv))
        server = serverdata(opts)

        #set up the collector and generator
        start = time.time()

        collector = resultcollector(server)
        end = time.time()
        print end-start
        
        jobs = jobgenerator(server)

        numjobsreceived = 0
        #begin distributing work
        for proc in xrange(1, min(numnodes, jobs.numjobs+1)):
            job = jobs.next(proc)
            pypar.send(('job',job), proc, tag=OUT)
        while numjobsreceived < jobs.jobindex:#while any job is still running
            #wait for any node to send a result
            msg, status = pypar.receive(pypar.any_source, return_status=True, tag=RETURN)
            numjobsreceived += 1
            proc, response = msg

            if jobs.hasnext(proc):#see if there is more work to be done
                job = jobs.next(proc)
                pypar.send(('job',job), proc, tag=OUT)#send it to the node that just completed

            #combine the results *after* sending the new job
            #(this way the worker can proceed while the results are being combined)
            collector.collect(response)


        #all jobs collected, kill the workers
        sendtoall(('Done', 0))

        #finish up the computation
        collector.finish()
        
    #--------------------#    
    # worker code
    #--------------------#
    else:
        while True:
            start = time.time()
            (code, msg), status = pypar.receive(0, return_status=True, tag=OUT)
            end = time.time()
            print 'waiting', end-start
            if code == 'Done':#all work is done
                opts.printruninfo()
                break
            elif code == 'Die':#abnormal exit
                break
            elif code == 'Start':
                opts = task(msg)
                sys.stdout = open(opts.logprefix+'%02d.log'%rank, 'w') #logfile
                print 'client', rank, 'running on', procname                
            else:
                start = time.time()
                jobnum, job = msg
                print jobnum
                result = opts.dojob(job)#do the job
                end = time.time()
                print 'working',msg[0], end-start

                start = time.time()
                pypar.send((rank, (jobnum, result)), 0, tag=RETURN)#return the result to the server
                end = time.time()
                print 'sending', end-start

    #------------------#
    #end of parallel code
    pypar.barrier()
    pypar.finalize()


start = time.time()
#cProfile.run('main()', 'nstatprof')
main()
end = time.time()
print end - start
