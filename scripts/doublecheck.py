def lines_exist(file1, file2, output_file):
    lines_found = []
    with open(file1, 'r') as f1:
        lines_file1 = [line.strip().lower() for line in f1]
        print(lines_file1)
    with open(file2, 'r') as f2:
        for line in f2:
            if line.strip().lower() in lines_file1:
                lines_found.append(line.strip())

    with open(output_file, 'w') as output:
        for line in lines_found:
            output.write(line + '\n')

if __name__ == "__main__":
    file1_path = "/home/kerney/wolb-recomb/simulations/testing_5_genome/possible_recomb.txt"
    file2_path = "/home/kerney/wolb-recomb/output.txt"
    output_file_path = "doublecheck.txt"  # Output file path

    lines_exist(file1_path, file2_path, output_file_path)
    print(f"Lines from file 2 found in file 1 have been written to {output_file_path}.")