"""
Wrapper on multiprocessing Pool.

Utility that will will take care of mapping job into pool of workers and allow
for user interrupts.
"""
import logging
import pathos.multiprocessing as mp
import signal


def how_many_cores():
    return mp.cpu_count()

def curr_process():
    """
    Get current process id.

    Returns:
        pid of multiprocessing process.
    """
    return mp.current_process().pid


def submit_jobs(function, inputs, num_threads):
    """
    Submit the inputs to the function, splitting the inputs across cores.

    Args:
        function (function):
            function we wish to apply in parallel.
        inputs (list of tuples):
            list of arguments we wish to pass in, each entry a different set.
        num_threads (int):
            number of threads to run with.

    Returns list of corresponding return values of the called function. None
    if timeout or user interrupt occurs.
    """
    if len(inputs) == 0:
        return

    logging.info('Processing %d inputs.', len(inputs))

    if num_threads > 1:
        logging.info('Parallel Mode: Using %d threads for writing.',
                     num_threads)
        pool = mp.ProcessPool(num_threads, maxtasksperchild=50)

        inputs = [(function,) + x for x in inputs]

        res = pool.amap(
            submit_helper,
            inputs,
            chunksize=1)

        out = res.get(0xFFFF)
    else:
        logging.info('Sequential Mode.')
        out = [function(*args) for args in inputs]

    return out


def submit_helper(args):
    """
    Helper function that allows us to use map_async.

    map_async only allows one argument to be passed in, so we pass in a
    tuple to this helper and then unpack the arguments here to pass to final
    function.  The function is the first argument of tuple.
    """
    function = args[0]
    inputs = args[1:]
    return function(*inputs)
