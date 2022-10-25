#!/usr/bin/env python

import pandas as pd
import shutil
import os
import argparse


def copy_rename_raw_files(input_dir, file_name_description, output_dir):
    '''
    '''
    data_dir = "data"

    # Extract the correspondance between sample name and current file names
    file_names = pd.read_csv(file_name_description)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_nb = 0
    renamed_sample_files = {}
    for path, dirs, files in os.walk(input_dir):
        for f in files:
            if not f.endswith("fastq.gz"):
                continue
            # Retrieve the name from the path
            split_path = os.path.split(path)
            name = split_path[1]
            if not any(file_names["Name in project"] == name):
                raise ValueError("A name in a project not found")
            # Retrieve the project id from the path
            # (to know if there is several directory for the same sample)
            project_id = os.path.split(split_path[0])[1]
            num = "a"
            if project_id.endswith("a"):
                num = "a"
            if project_id.endswith("b"):
                num = "b"
            # Extract the sample name based on the project name
            line = file_names[file_names["Name in project"] == name]
            sample_name = str(line["Sample name"].tolist()[0])
            # Increment the num if the same sample was already found
            if sample_name in renamed_sample_files:
                previous_num = renamed_sample_files[sample_name]
                num = chr(ord(previous_num) + 1)
            # Build the filename
            filename = "%s%s.fastq.gz" % (sample_name, num)
            # Copy the file
            shutil.copy2(
                os.path.join(path, f),
                os.path.join(output_dir, filename))
            # Just some stats
            file_nb += 1
            renamed_sample_files.setdefault(sample_name, "")
            renamed_sample_files[sample_name] = num
            print(filename)
    print("Number of renamed files: %s" %(file_nb))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--file_name_description', required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()
    copy_rename_raw_files(
        args.input_dir,
        args.file_name_description,
        args.output_dir)
    