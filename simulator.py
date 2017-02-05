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
        evaluations = 600000
        dimensions = 1000
        particles = 100
        functions = range(21, 29)
        functions = [23]
        runs = 30
        # for topology in [("kregular", 6)]:
        #     for function in [23]:
        #         for k in [12, 14, 16, 18, 22, 24, 26, 28]:
        #             commands += ["java -jar PSO/jar/pso.jar 1 %d %d %d %d %d 0 %d True > %s%d_F%02d_%02d.with_positions " %
        #                          (particles, evaluations, dimensions, function, topology[1], k,
        #                          topology[0], k, function, r) for r in range(runs)]

        # runs = 2
        # # for topology in [("global", 0, 0), ("ring", 1, 0), ("regular30", 6, 30)]:
        # for topology in [("noc2", 0, 0)]:
        #     for function in [23]:
        #         commands += ["java -jar PSO/jar/pso_no_c2.jar 1 %d %d %d %d %d 0 %d True > %s_F%02d_%02d.with_positions " %
        #                      (particles, evaluations, dimensions, function, topology[1], topology[2],
        #                      topology[0], function, r) for r in range(runs)]
        # runs = 30
        # for topology in [("global", 0, 0), ("ring", 1, 0), ("noc2", 0, 0)]:
        #     for function in functions:
        #         commands += ["java -jar PSO/jar/pso_no_c2.jar 1 %d %d %d %d %d 0 %d False > %s_F%02d_%02d.normal " %
        #                      (particles, evaluations, dimensions, function, topology[1], topology[2],
        #                      topology[0], function, r) for r in range(runs)]
        for topology in [("kregular", 6)]:
            for function in [23]:
                # for k in [12, 14, 16, 18, 22, 24, 26, 28] + range(10, 100, 10):
                for k in range(3, 10, 1):
                    commands += ["java -jar PSO/jar/pso_no_info.jar 1 %d %d %d %d %d 0 %d False > %s%d_F%02d_%02d.no_info " %
                                 (particles, evaluations, dimensions, function, topology[1], k,
                                 topology[0], k, function, r) for r in range(runs)]
        return commands

"""
execfile("simulator.py")
Simulator.execute(Simulator.pso(), number_of_processes=6, delay_between=10)
Simulator.pso()[0]
range(10, 21)
"""
