__author__ = 'marcos'
import os
import time
import sys


class Simulator:
    def __init__(self):
        pass

    @staticmethod
    def execute(commands, number_of_processes=2, delay_between=2, output_dir=None):
        if not output_dir:
            output_dir = "."
        if number_of_processes > len(commands):
            number_of_processes = len(commands)
        epoch_time = str(int(time.time()))
        pairs_and_files = dict([(i, (commands[i],
                                     "%s/%s_touch_%d" % (output_dir, epoch_time, i)))
                                for i in range(len(commands))])
        pairs_indices = range(0, len(commands))
        # remove old files
        for index in pairs_indices:
            _, touch = pairs_and_files[index]
            os.system("rm " + touch + " 2> /dev/null")
        threads_started = []
        # pairs_indices can change through the iterations because we can try to retrieve data again
        for index in pairs_indices:
            command, touch = pairs_and_files[index]
            os.system('( ' + command + ' ; touch ' + touch + ' ) &')
            print "> started [%d]: '%s'" % (index, command)
            threads_started.append(index)
            if len(threads_started) < number_of_processes:
                continue
            # if we already started the maximum numbers of processes, we join
            thread_done = None
            while thread_done is None:
                for thread in threads_started:
                    _, touch = pairs_and_files[thread]
                    if os.path.exists(touch):
                        thread_done = thread
                        break
                time.sleep(delay_between)
            print "> finished [%d]" % thread_done
            threads_started.remove(thread_done)
            # total_executed += len(threads_started)
            # sys.stdout.write("%d/%d/%d " % (len(threads_started), total_executed, len(commands)))
            # sys.stdout.flush()
            # threads_started = []
        # for _, touch in pairs_and_files
        while threads_started:
            threads_done = []
            for thread in threads_started:
                _, touch = pairs_and_files[thread]
                if os.path.exists(touch):
                    threads_done.append(thread)
            for thread in threads_done:
                threads_started.remove(thread)
                print "> finished [%d]" % thread
            time.sleep(delay_between)
        sys.stdout.write("\n")
        # remove old files
        for index in pairs_indices:
            _, touch = pairs_and_files[index]
            os.system("rm " + touch + " 2> /dev/null")

    @staticmethod
    def pso():
        """
        PSO.jar runs particles evaluations dimensions function topology mechanism [mechanism_parameter]
        > TOPOLOGY
         0: GLOBAL 1: RING 2: RANDOM 3: VON_NEUMANN 4: THREESOME_PARTNERS 5: NSOME_PARTNERS
        > MECHANISM
         0: NONE 1: DYNAMIC_2011
        """
        commands = []
        topologies = [("global", 0), ("ring", 1), ("vonneumann", 4)]
        runs = 30
        runs = 2
        functions = range(1, 21)
        functions = range(1, 3)
        evaluations = 100000
        evaluations = 1000
        dimensions = 1000
        particles = 100
        for topology in topologies:
            for function in functions:
                commands += ["java -jar PSO/jar/pso.jar 1 %d %d %d %d %d 0 > %s_F%02d_%02d.teste " %
                             (particles, evaluations, dimensions, function, topology[1],
                             topology[0], function, r) for r in range(runs)]
        # dynamic topology
        for topology in topologies:
            for function in functions:
                    commands += ["java -jar PSO/jar/pso.jar 1 %d %d %d %d %d 1 > dynamic%s_F%02d_%02d.teste " %
                                 (particles, evaluations, dimensions, function, topology[1],
                                 topology[0], function, r) for r in range(runs)]
        return commands
"""
execfile("simulator.py")
commands = ["java -jar PSO/jar/pso.jar > t%02d.teste " % i for i in range(50)]
commands = ["sleep %d" % ((i+3)%4) for i in range(10)]
Simulator.execute(Simulator.pso(), number_of_processes=4, delay_between=0)
Simulator.pso()

"""
