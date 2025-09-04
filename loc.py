import os

emptiness_criterion_level = 3

def is_empty_line(line):
    stripped = line.strip()
    c = False
    if emptiness_criterion_level >= 1:
        c = c or stripped in ['', '\n']
    if emptiness_criterion_level >= 2:
        c = c or stripped in ['}', '{']
    if emptiness_criterion_level >= 3:
        c = c or stripped.startswith('//')
    return c

def count_nonempty_lines(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        # Read all lines and filter out empty ones (after stripping whitespace)
        return sum(1 for line in f if not is_empty_line(line))

def count_all_c_files(directory):
    total_c_lines = 0
    total_h_lines = 0
    
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if "build" in file_path:
                continue
            try:
                if file.endswith('.c') or file.endswith('.cpp'):
                    lines = count_nonempty_lines(file_path)
                    print(f"{file_path}: {lines} lines")
                    total_c_lines += lines
                elif file.endswith('.h') or file.endswith('.hpp'):
                    lines = count_nonempty_lines(file_path)
                    print(f"{file_path}: {lines} lines")
                    total_h_lines += lines
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
    
    return total_c_lines, total_h_lines

if __name__ == "__main__":
    # Use current directory as default
    directory = "."
    total_c, total_h = count_all_c_files(directory)
    print(f"\nTotal non-empty lines in object files: {total_c}")
    print(f"Total non-empty lines in header files: {total_h}")
    print(f"Total non-empty lines in all files: {total_c + total_h}")
