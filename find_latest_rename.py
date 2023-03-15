import glob
import os
import sys
import shutil

if len(sys.argv) < 3:
    print("Please provide a path to the folder and destination folder")
    sys.exit()

folder_path = sys.argv[1]
dest_folder = sys.argv[2]
new_header = '>' + sys.argv[1]


if not os.path.exists(dest_folder):
    os.makedirs(dest_folder)

list_of_files = glob.glob(folder_path + '/*/pilon/*.fasta') # * means all if need specific format then *.csv
latest_file = max(list_of_files, key=os.path.getmtime)
filename = os.path.basename(latest_file)
new_filename = os.path.join(dest_folder, sys.argv[1] + '.fasta')


with open(latest_file, 'r') as infile, open(new_filename, 'w') as outfile:
    first_line = True
    for line in infile:
        if first_line:
            outfile.write(new_header + '\n')
            first_line = False
        else:
            outfile.write(line)

shutil.copystat(latest_file, new_filename)

print(f"{filename} copied to {new_filename} with new header: {new_header}")



