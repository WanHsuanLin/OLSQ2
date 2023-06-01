import os
import argparse
from run_olsq2 import get_output

if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser()
    # Adding optional argument
    # parser.add_argument("--dt", dest='device_type', type=str,
    #     help="grid, ourense, sycamore, rochester, tokyo, aspen-4, or eagle")
    # parser.add_argument("--d", dest='device', type=int,
    #     help="device (x-by-x grid)")
    parser.add_argument("--o", dest='output_folder', type=str, default='.',
        help="the folder to store results")
    parser.add_argument("--i", dest="input_folder", type=str,
        help="Input file folder")
    # parser.add_argument("--encoding", dest="encoding", type=int, default=1,
    #     help="seqcounter = 1, sortnetwrk  = 2, cardnetwrk  = 3, totalizer   = 6, mtotalizer  = 7. kmtotalizer = 8, native = 9")
    # parser.add_argument("--sabre", action='store_true', default=False,
    #     help="Use sabre to get SWAP upper bound")
    # parser.add_argument("--tran", action='store_true', default=False,
    #     help="Use TB-OLSQ")
    # parser.add_argument("--swap", action='store_true', default=False,
    #     help="Optimize SWAP")
    # parser.add_argument("--all_commute", action='store_true', default=False,
    #     help="All gates  are commute. e.g., qaoa")
    # parser.add_argument("--swap_bound", dest="swap_bound", type=int, default=-1,
    #     help="user define swap bound")
    # parser.add_argument("--swap_duration", dest="swap_duration", type=int, default=1,
    #     help="swap duration")
    # parser.add_argument("--IR_output", action='store_true', default=False,
    #     help="user define output type")
    # # add the argument to verify the result
    # parser.add_argument("--check", action='store_true', default=False,
    #     help="Verify the result")
    # Read arguments from command line

    args = parser.parse_args()
    # swap = args.swap
    # mode_type = False
    # if args.tran:
    # #     mode_type = True
    # device_type = args.device_type
    # use_sabre = args.sabre
    # encoding = args.encoding
    # swap_duration = args.swap_duration
    # swap_bound = args.swap_bound
    output_directory = args.output_folder
    # output_mode = args.IR_output
    # device_size = args.device
    # verify_result = args.check
    
    input_directory = args.input_folder
    for filename in os.listdir(input_directory):
        file_name = os.path.join(input_directory, filename)
        get_output(file_name, output_directory)
        # get_output(file_name, swap, mode_type, device_type, use_sabre, encoding, swap_duration, swap_bound, store, output_mode, device_size, verify_result)


