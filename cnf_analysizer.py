import argparse

if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser()
    # Adding optional argument
    parser.add_argument("--f", dest='file', type=str,
        help="cnf file name")
    # Read arguments from command line
    args = parser.parse_args()
    file = open(args.file, 'r')
    max_length = 0
    min_length = 1E10
    total_clause_length = 0
    clause_count = 0
    line = file.readline()
    line = line.strip().split()
    clause_count = int(line[-1])
    line = file.readline()
    while line:
        line = line.strip().split()
        max_length = max(max_length, len(line))
        min_length = min(max_length, len(line))
        total_clause_length += len(line)
        line = file.readline()
    file.close
    
    print(f"File: {args.file}")
    print(f"#clause: {clause_count}")
    print(f"Max clause length: {max_length}")
    print(f"Min clause length: {min_length}")
    print(f"Average clause length: {total_clause_length/clause_count}")