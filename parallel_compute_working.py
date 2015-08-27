__author__ = 'jlu96'


import multiprocessing as mup
import time
import sys
import os

# class partitionThread (threading.Thread):
#     def __init__(self, function, args, input, input_index, threadID):
#         threading.Thread.__init__(self)
#         self.threadID = threadID
#         self.function = function
#         self.args = args
#         self.input = input
#         self.input_index = input_index
#         self.threadID = threadID
#     def run(self):
#         args = self.args
#         args[self.input_index] = self.input
#         return (function)
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
        return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


# def chunks(l, n):
#     """Yield successive n-sized chunks from l."""
#     for i in xrange(0, len(l), n):
#         yield l[i:i+n]

def partition_inputs(input, number):
    num_inputs = len(input)
    return [input[num_inputs * i/number:num_inputs * (i+1)/number] for i in range(number)]

def partition_dictionary(inputdict, number):
    # print "Length of old inputdict is", len(inputdict)
    newkeys_list = partition_inputs(inputdict.keys(), number)
    newdict_list = []
    for i in range(number):
        newdict = inputdict.copy()
        newkeys = newkeys_list[i]
        for key in inputdict:
            if key not in newkeys:
                newdict.pop(key)
        newdict_list.append(newdict)
    # print "Length of new inputdicts is ", [len(newdict) for newdict in newdict_list]
    print len(set(newdict_list[0].keys()).difference(set(newdict_list[1].keys())))
    print set(inputdict.keys()).difference(set(newdict_list[0].keys())) == set(newdict_list[1].keys())
    print set(inputdict.keys()).difference(set(newdict_list[1].keys())) == set(newdict_list[0].keys())
    # print "hello"
    # print "hello"
    # print "help"

    return newdict_list



# Left off here 7/22/15 -jlu
# Make the worker function create a return_dict
# Then use the manager to put everything together.
# def worker_function(function, args, procnum, return_dict):
#     # print "Process ", procnum, " running function ", function, " with args ", args
#     return_dict[procnum] = function(args)
#
# def parallel_compute(function, args, input, input_index, join_functions_dict, number=32):
#     # Partition the data
#     # Only return those indices with a corresponding join function.
#
#     new_inputs = partition_inputs(input, number)
#
#     manager = mup.Manager()
#     return_dict = manager.dict()
#
#     jobs = []
#
#     for i in range(number):
#         new_input = new_inputs[i]
#         new_args = args[:]
#         new_args[input_index] = new_input
#         p = mup.Process(target=worker_function, args=(function, new_args, i, return_dict))
#         jobs.append(p)
#         p.start()
#
#     for proc in jobs:
#         proc.join()
#
#     return_values = return_dict[0]
#
#     for index in join_functions_dict:
#         for proc in range(1, number):
#             return_values[index] = join_functions_dict[index](return_values[index], return_dict[proc][index])
#
#     return return_values

def worker_function_new(args):
    function = args[0]
    function_args = args[1]
    procnum = args[2]

    print "Process ", os.getpid(), "number ", procnum, "with function ", function, " begun at ", time.strftime("%H:%M:%S")
    returnvalue = function(*function_args)
    print "Process ", os.getpid(), "number ", procnum, "with function ", function, " finished at ", time.strftime("%H:%M:%S")
    return returnvalue

# def parallel_compute_new(function, args, input, input_index, join_functions_dict, number=32, procnumber=32):
#     # Partition the data
#     # Only return those indices with a corresponding join function.
#     t = time.time()
#
#     new_inputs = partition_inputs(input, number)
#     t1 = time.time()
#     print "time to partition inputs ", t1- t
#
#     worker_function_args_list = []
#
#     #print "Args were ", args
#
#     for i in range(number):
#         new_args = args[:]
#         new_args[input_index] = new_inputs[i]
#         #print "New_args are ", new_args
#         worker_function_args_list.append((function, new_args, i))
#
#     t2 = time.time()
#     print "Time to prepare function args ", t2-t1
#
#     pool = mup.Pool(processes=procnumber)
#
#     t2_5 = time.time()
#     print "Time to prepare processes ", t2_5 - t2
#
#     #print worker_function_args_list
#
#     returnlist = pool.map(worker_function_new, worker_function_args_list)
#     t3 = time.time()
#     print "Time to get return values ", t3-t2_5
#
#     return_values = list(returnlist[0])
#
#     for index in join_functions_dict:
#         for proc in range(1, number):
#             return_values[index] = join_functions_dict[index](return_values[index], returnlist[proc][index])
#
#     t4 = time.time()
#     print "Time to integrate return values ", t4-t3
#
#     return return_values

def initialize_pool(procnumber=32):
    pool = mup.Pool(processes=procnumber)
    print "Pool initialized at ", time.strftime("%H:%M:%S")
    return pool

# def parallel_compute_new(function, args, input, input_index, partition_input_function, join_functions_dict, pool,
#                          number=32, procnumber=32):
#     # Partition the data
#     # Only return those indices with a corresponding join function.
#     t = time.time()
#
#     new_inputs = partition_input_function(input, number)
#     t1 = time.time()
#     print "time to partition inputs ", t1- t
#
#     worker_function_args_list = []
#
#     #print "Args were ", args
#
#     for i in range(number):
#         new_args = args[:]
#         new_args[input_index] = new_inputs[i]
#         #print "New_args are ", new_args
#         worker_function_args_list.append((function, new_args, i))
#
#     t2 = time.time()
#     print "Time to prepare function args ", t2-t1
#
#
#     t2_5 = time.time()
#     # print "Time to prepare processes ", t2_5 - t2
#
#     #print worker_function_args_list
#
#     print "Mapping begun at ", time.strftime("%H:%M:%S")
#     returnlist = pool.map(worker_function_new, worker_function_args_list)
#     # pool.close()
#     # pool.join()
#     t3 = time.time()
#     print "Mapping finished at ", time.strftime("%H:%M:%S")
#     print "Time to get return values ", t3-t2_5
#
#     return_values = list(returnlist[0])
#
#     for index in join_functions_dict:
#         for proc in range(1, number):
#             return_values[index] = join_functions_dict[index](return_values[index], returnlist[proc][index])
#
#     t4 = time.time()
#     print "Time to integrate return values ", t4-t3
#
#     return return_values

# Old working parallel compute new with pool passed as argument above.

def parallel_compute_new(function, args, input, input_index, partition_input_function, join_functions_dict,
                         number=32, procnumber=32):
    # Partition the data
    # Only return those indices with a corresponding join function.
    t = time.time()

    new_inputs = partition_input_function(input, number)
    t1 = time.time()
    print "time to partition inputs ", t1- t

    worker_function_args_list = []

    #print "Args were ", args

    for i in range(number):
        new_args = args[:]
        new_args[input_index] = new_inputs[i]
        #print "New_args are ", new_args
        worker_function_args_list.append((function, new_args, i))

    t2 = time.time()
    print "Time to prepare function args ", t2-t1



    pool = initialize_pool(number)

    t2_5 = time.time()
    # print "Time to prepare processes ", t2_5 - t2

    #print worker_function_args_list

    print "Mapping begun at ", time.strftime("%H:%M:%S")
    returnlist = pool.map(worker_function_new, worker_function_args_list)
    pool.close()
    pool.join()
    t3 = time.time()
    print "Mapping finished at ", time.strftime("%H:%M:%S")
    print "Time to get return values ", t3-t2_5

    return_values = list(returnlist[0])

    for index in join_functions_dict:
        for proc in range(1, number):
            return_values[index] = join_functions_dict[index](return_values[index], returnlist[proc][index])

    t4 = time.time()
    print "Time to integrate return values ", t4-t3

    return return_values





#
# def basic_worker_function(i):
#     a = []
#
#     for b in i:
#         a.append(b*3.0)
#     return b

# def basic_parallel_compute(function, args, input, input_index, join_functions_dict, number=32, procnumber=32):
#     # Partition the data
#     # Only return those indices with a corresponding join function.
#     t = time.time()
#
#     new_inputs = partition_inputs(input, number)
#     t1 = time.time()
#     print "time to partition inputs ", t1- t
#
#     worker_function_args_list = []
#
#     #print "Args were ", args
#
#     for i in range(number):
#         new_args = args[:]
#         new_args[input_index] = new_inputs[i]
#         #print "New_args are ", new_args
#         worker_function_args_list.append((function, new_args, i))
#
#     t2 = time.time()
#     print "Time to prepare function args ", t2-t1
#
#     pool = mup.Pool(processes=procnumber)
#
#     t2_5 = time.time()
#     print "Time to prepare processes ", t2_5 - t2
#
#     #print worker_function_args_list
#
#     returnlist = pool.map(worker_function_new, worker_function_args_list)
#     t3 = time.time()
#     print "Time to get return values ", t3-t2_5
#
#     return_values = list(returnlist[0])
#
#     for index in join_functions_dict:
#         for proc in range(1, number):
#             return_values[index] = join_functions_dict[index](return_values[index], returnlist[proc][index])
#
#     t4 = time.time()
#     print "Time to integrate return values ", t4-t3
#
#     return return_values
#




    # jobs = []
    #
    # for i in range(number):
    #     new_input = new_inputs[i]
    #     new_args = args
    #     new_args[input_index] = new_input
    #     p = mup.Process(target=worker_function, args=(function, new_args, i, return_dict))
    #     jobs.append(p)
    #     p.start()
    #
    # for proc in jobs:
    #     proc.join()
    #
    # return_values = return_dict[0]
    # #
    # for index in join_functions_dict:
    #     for proc in range(1, number):
    #         return_values[index] = join_functions_dict[index](return_values[index], return_dict[proc][index])
    #
    # return return_values



def test_function(a, b, c):
    d = []
    e = 0
    f = {}
    for i in b:
        d.append(a*i*c)
        f[i] = set([a, i, c])
    return d, e, f

def test_parallel_compute(a, b, c, number=32, procnumber=32):

    join_functions_dict = {}
    join_functions_dict[0] = lambda a, b: a + b
    join_functions_dict[2] = combine_dictionaries

    # a = 100.
    # b = range(1000000000)
    # c = 15.


    t0 = time.time()
    d, e, f = test_function(a, b, c)
    t1 = time.time()
    print "Time to run with no parallel compute: ", t1 - t0


    d1, e1, f1 = parallel_compute_new(test_function, [a, b, c], b, 1, join_functions_dict, number=number, procnumber=procnumber)
    t2 = time.time()
    print "Time to run with parallel compute: ", t2 - t1, " with number ", number

    print "d and d1 equality: ", d == d1
    print "e and e1 equality: ", e == e1
    print "f and f1 equality: ", f == f1


def combine_dictionaries(a,b):
    a.update(b)
    return a

def combine_lists(a, b):
    a.extend(b)
    return a

def main():
    # test_parallel_compute(100., range(1000000), 15., number=2)
    args = get_parser().parse_args(sys.argv[1:])
    b = range(args.b)
    print
    print "*********************************************************************************"
    print "Parallel_compute called with arguments ", args.a, args.b, args.c, args.number
    test_parallel_compute(args.a, b, args.c, number=args.number, procnumber=args.procnumber)
    print "*********************************************************************************"
    print

def get_parser():
    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--a', type=float, default = 100.)

    parser.add_argument('-b', '--b', type=int, default=1000000, help="length of b list.")

    parser.add_argument('-c', '--c', type=int, default=15.)

    parser.add_argument('-n', '--number', type=int, default=32)

    parser.add_argument('-p', '--procnumber', type=int, default=32)
    return parser

if __name__ == '__main__':
    main()


    # For each partition, start a new thread that implements that function

    # Integrate the return values
