import argparse
from unittest import result
from olsq_ori.device import qcdevice
from olsq_ori import OLSQ
import json
import timeit


def get_nnGrid(n: int, swap_duration):
    my_coupling = []
    for i in range(n):
        for j in range(n):
            if j < n-1:
                my_coupling.append((i*n +j,i*n +j+1))
            if i < n-1:
                my_coupling.append((i*n +j,i*n +j+n))
   
    return qcdevice(name="grid", nqubits=n*n,
        connection=my_coupling, swap_duration=swap_duration)

def get_device_by_name(name, swap_duration):
    device_set_edge = { "qx" : [(0,2), (0,1), (1,2), (2,3), (2,4), (3,4)],
                        "ourense": [(0, 1), (1, 2), (1, 3), (3, 4)],
                       "sycamore": [(0, 6), (1, 6), (1, 7), (2, 7), (2, 8), (3, 8), (3, 9), (4, 9), (4, 10), (5, 10), (5, 11),
                                    (6, 12), (6, 13), (7, 13), (7, 14), (8, 14), (8, 15), (9, 15), (9, 16), (10, 16), (10, 17), (11, 17),
                                    (12, 18), (13, 18), (13, 19), (14, 19), (14, 20), (15, 20), (15, 21), (16, 21), (16, 22), (17, 22), (17, 23),
                                    (18, 24), (18, 25), (19, 25), (19, 26), (20, 26), (20, 27), (21, 27), (21, 28), (22, 28), (22, 29), (23, 29),
                                    (24, 30), (25, 30), (25, 31), (26, 31), (26, 32), (27, 32), (27, 33), (28, 33), (28, 34), (29, 34), (29, 35),
                                    (30, 36), (30, 37), (31, 37), (31, 38), (32, 38), (32, 39), (33, 39), (33, 40), (34, 40), (34, 41), (35, 41),
                                    (36, 42), (37, 42), (37, 43), (38, 43), (38, 44), (39, 44), (39, 45), (40, 45), (40, 46), (41, 46), (41, 47),
                                    (42, 48), (42, 49), (43, 49), (43, 50), (44, 50), (44, 51), (45, 51), (45, 52), (46, 52), (46, 53), (47, 53)],
                       "rochester": [(0, 1), (1, 2), (2, 3), (3, 4),
                                     (0, 5), (4, 6), (5, 9), (6, 13),
                                     (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15),
                                     (7, 16), (11, 17), (15, 18), (16, 19), (17, 23), (18, 27),
                                     (19, 20), (20, 21), (21, 22), (22, 23), (23, 24), (24, 25), (25, 26), (26, 27),
                                     (21, 28), (25, 29), (28, 32), (29, 36),
                                     (30, 31), (31, 32), (32, 33), (33, 34), (34, 35), (35, 36), (36, 37), (37, 38),
                                     (30, 39), (34, 40), (38, 41), (39, 42), (40, 46), (41, 50),
                                     (42, 43), (43, 44), (44, 45), (45, 46), (46, 47), (47, 48), (48, 49), (49, 50),
                                     (44, 51), (48, 52)],
                       "tokyo": [(0, 1), (1, 2), (2, 3), (3, 4),
                                 (0, 5), (1, 6), (1, 7), (2, 6), (2, 7), (3, 8), (3, 9), (4, 8), (4, 9),
                                 (5, 6), (6, 7), (7, 8), (8, 9),
                                 (5, 10), (5, 11), (6, 10), (6, 11), (7, 12), (7, 13), (8, 12), (8, 13), (9, 14),
                                 (10, 11), (11, 12), (12, 13), (13, 14),
                                 (10, 15), (11, 16), (11, 17), (12, 16), (12, 17), (13, 18), (13, 19), (14, 18), (14, 19),
                                 (15, 16), (16, 17), (17, 18), (18, 19)],
                       "aspen-4": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7),
                                   (0, 8), (3, 11), (4, 12), (7, 15),
                                   (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15)],
                        "eagle": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (9, 10), (10, 11), (11, 12), (12, 13),
                                   (0, 14), (14, 18), (4, 15), (15, 22), (8, 16), (16, 26), (12, 17), (17, 30),
                                   (18, 19), (19, 20), (20, 21), (21, 22), (22, 23), (23, 24), (24, 25), (25, 26), (26, 27), (27, 28), (28, 29), (29, 30), (30, 31), (31, 32),
                                   (20, 33), (33, 39), (24, 34), (34, 43), (28, 35), (35, 47), (32, 36), (36, 51),
                                   (37, 38), (38, 39), (39, 40), (40, 41), (41, 42), (42, 43), (43, 44), (44, 45), (45, 46), (46, 47), (47, 48), (48, 49), (49, 50), (50, 51),
                                   (37, 52), (52, 56), (41, 53), (53, 60), (45, 54), (54, 64), (49, 55), (55, 68),
                                   (56, 57), (57, 58), (58, 59), (59, 60), (60, 61), (61, 62), (62, 63), (63, 64), (64, 65), (65, 66), (66, 67), (67, 68), (68, 69), (69, 70),
                                   (58, 71), (71, 77), (62, 72), (72, 81), (66, 73), (73, 85), (70, 74), (74, 89),
                                   (75, 76), (76, 77), (77, 78), (78, 79), (79, 80), (80, 81), (81, 82), (82, 83), (83, 84), (84, 85), (85, 86), (86, 87), (87, 88), (88, 89),
                                   (75, 90), (90, 94), (79, 91), (91, 98), (83, 92), (92, 102), (87, 93), (93, 106),
                                   (94, 95), (95, 96), (96, 97), (97, 98), (98, 99), (99, 100), (100, 101), (101, 102), (102, 103), (103, 104), (104, 105), (105, 106), (106, 107), (107, 108),
                                   (96, 109), (100, 110), (110, 118), (104, 111), (111, 112), (108, 112), (112, 126),
                                   (113, 114), (114, 115), (115, 116), (116, 117), (117, 118), (118, 119), (119, 120), (120, 121), (121, 122), (122, 123), (123, 124), (124, 125), (125, 126)]
                       }
    
    device_set_qubit_num = {"qx": 5,
                        "ourense": 5,
                       "sycamore": 54,
                       "rochester": 53,
                       "tokyo": 20,
                       "aspen-4": 16,
                       "eagle": 127}
    
    device = qcdevice(name=name, nqubits=device_set_qubit_num[name],
                        connection=device_set_edge[name], swap_duration=swap_duration)
    return device

def run_olsq_tbolsq(obj_is_swap, circuit_info, mode, device, use_sabre, encoding, swap_bound = -1, thread = 1):
    lsqc_solver = OLSQ(obj_is_swap = obj_is_swap, mode=mode, encoding = encoding, swap_up_bound=swap_bound, thread = thread)
    lsqc_solver.setprogram(circuit_info)
    lsqc_solver.setdevice(device)
    start = timeit.default_timer()
    result = lsqc_solver.solve(use_sabre, output_mode="IR", memory_max_size=0, verbose=0)
    stop = timeit.default_timer()
    print('Time: ', stop - start)  
    return result

def dump_olsq(obj_is_swap, circuit_info, device, folder, encoding):
    lsqc_solver = OLSQ(objective = "depth", mode="transition", encoding = encoding)
    lsqc_solver.setprogram(circuit_info)
    lsqc_solver.setdevice(device)
    lsqc_solver.dump(folder)
    return

if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser()
    # Adding optional argument
    parser.add_argument("--dt", dest='device_type', type=str,
        help="grid or heavy-hexagon or ring")
    parser.add_argument("--d", dest='device', type=int,
        help="device (x-by-x grid)")
    parser.add_argument("--f", dest='folder', type=str,
        help="the folder to store results")
    parser.add_argument("--b", dest='benchmark', type=str,
        help="qaoa or others")
    parser.add_argument("--qf", dest="qasm", type=str,
        help="Input file name")
    parser.add_argument("--encoding", dest="encoding", type=int, default=1,
        help="seqcounter  = 1, sortnetwrk  = 2, cardnetwrk  = 3, totalizer   = 6, mtotalizer  = 7. kmtotalizer = 8, native = 9")
    parser.add_argument("--sabre", action='store_true', default=False,
        help="Use sabre to get SWAP upper bound")
    parser.add_argument("--tran", action='store_true', default=False,
        help="Use TB-OLSQ")
    parser.add_argument("--swap", action='store_true', default=False,
        help="Optimize SWAP")
    parser.add_argument("--dump", action='store_true', default=False,
        help="Only dump constraint")
    parser.add_argument("--test_sabre", action='store_true', default=False,
        help="test sabre scalability")
    parser.add_argument("--swap_bound", dest="swap_bound", type=int, default=-1,
        help="user define swap bound")
    parser.add_argument("--thread", dest="thread", type=int, default=-1,
        help="number of thread used for SWAP optimization")
    # Read arguments from command line
    args = parser.parse_args()

    if args.benchmark == "qaoa":
        swap_duration = 1
    else:
        swap_duration = 3
    circuit_info = open(args.qasm, "r").read()
    if args.device_type == "grid":
        device = get_nnGrid(args.device, swap_duration)
    else:
        device = get_device_by_name(args.device_type, swap_duration)

    if args.dump:
        dump_olsq(args.swap, circuit_info, device, args.folder, args.encoding)
    elif args.test_sabre:
        import datetime
        from olsq.run_h_compiler import run_sabre
        start_time = datetime.datetime.now()
        swap_num, depth, _ = run_sabre(args.benchmark, circuit_info, device.list_qubit_edge, device.count_physical_qubit)
        print("Run heuristic compiler sabre to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
        print("Heuristic optimization time = {}".format(datetime.datetime.now() - start_time))
    else:
        data = dict()
        b_file = args.qasm.split('.')
        b_file = b_file[-2]
        b_file = b_file.split('/')
        b_file = b_file[-1]
        file_name = args.folder+"/"+str(args.device_type)+"_"+b_file+".json"

        mode = "normal"
        if args.tran:
            mode = "transition"
        result = run_olsq_tbolsq(args.swap, circuit_info, mode, device, args.sabre, args.encoding, args.thread)
        data["device"] = str(args.device)
        data["mode"] = mode
        data["depth"] = result[0]
        data["gate_spec"] = result[1]
        data["gate"] = result[2]
        data["final_mapping"] = result[3]
        data["initial_mapping"] = result[4]
        
        with open(file_name, 'w') as file_object:
            json.dump(data, file_object, default=int)
    
