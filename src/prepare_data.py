from interact_galaxy import *
from bioblend.galaxy import GalaxyInstance
import pandas as pd


if __name__ == '__main__':
    '''
    Import the files from the data library, merge the files sequenced on 2
    different lanes (for Project_S178 and Project_S225) and move the input files
    into collections
    '''
    # Get config
    with open('config.yaml', 'r') as config_f:
        config = yaml.safe_load(config_f)
    # Connect to Galaxy and retrieve the history
    gi = GalaxyInstance(config["galaxy_url"], config["api_key"])
    hist = get_hist(gi, config["hist_name"])
    # Get tools in the Galaxy instance
    tools = gi.tools.get_tools()
    # Extract the sample names
    finename_desc_df = pd.read_csv("data/file_description.csv", index_col = 0)
    sample_names = list(finename_desc_df.index)
    out = open("logs/prepare_data.log", "w")
    # Find the data library
    lib = gi.libraries.get_libraries(name=config["library_names"]["input_data"])
    assert len(lib) > 0, "No library found for Prinz lab"
    lib_id = lib[0]["id"]
    # Parse the data library datasets
    for ds in gi.libraries.show_library(lib_id, contents=True):
        # Eliminate the folder
        if ds['type'] != 'file':
            continue
        # Eliminate the files from other folders
        if ds["name"].find(config["folder_names"]["input_data"]) == -1:
            continue
        # Add the files to the history
        gi.histories.upload_dataset_from_library(hist, ds["id"])
    out.write("Added the files to the history\n")
    # Retrieve the name of samples to merge
    to_merge = {}
    for sample in sample_names:
        project_id = finename_desc_df["Project id"][sample]
        if project_id not in ["Project_S178", "Project_S198", "Project_S225"]:
            continue
        to_merge.setdefault(sample, [])
    # Parse the dataset in history to extract the ids of dataset to merge,
    # rename the other files and add them to a collection
    raw_dataset_ids = []
    input_data = gi.histories.show_matching_datasets(hist)
    for dataset in input_data:
        name = dataset['name']
        if not name.endswith("fastq"):
            continue
        sample_name = os.path.splitext(name)[0][:-1]
        if sample_name in to_merge:
            to_merge[sample_name].append(dataset["id"])
            # Hide the file
            gi.histories.update_dataset(hist, dataset["id"], visible = False)
        else:
            # Rename the file and hide it
            gi.histories.update_dataset(hist,
                dataset["id"],
                name="%s%s" % (config["collection_names"]["raw_data"], sample_name),
                visible = False)
            # Add the file to the collection
            raw_dataset_ids.append({'id': dataset["id"], 'name': sample_name, 'src': 'hda'})
    out.write("Got info about the files to merge\n")
    # Get concatenate tool
    tool_id = get_working_tool_id(config["tool_ids"]["merging"], tools)
    # Merge datasets
    unmerged_dataset_ids = []
    for dataset in to_merge:
        if len(to_merge[dataset]) != 2:
            print("Issue with %s" %(dataset))
            continue
        # Create the input datamap
        datamap = dict()
        datamap["inputs"] = [
            {'src':'hda', 'id': to_merge[dataset][0]},
            {'src':'hda', 'id': to_merge[dataset][1]}]
        # Run the tool
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Rename the dataset and hide it
        gi.histories.update_dataset(
            hist,
            info['outputs'][0]['id'],
            name="%s%s" % (config["collection_names"]["raw_data"], dataset),
            visible = False)
        # Add the merge file to the collection of raw dataset
        raw_dataset_ids.append({'id': info['outputs'][0]['id'], 'name': dataset, 'src': 'hda'})
        # Add the input files (before merging) to the collection of unmerged
        # datasets
        unmerged_dataset_ids.append({'id': to_merge[dataset][0], 'name': "%s 0" % dataset, 'src': 'hda'})
        unmerged_dataset_ids.append({'id': to_merge[dataset][1], 'name': "%s 1" % dataset, 'src': 'hda'})
    out.write("Merged the files\n")
    # Prepare and create the collections
    raw_data_collection = {
        'collection_type': 'list',
        'element_identifiers': raw_dataset_ids,
        'name': config["collection_names"]["raw_data"]
    }
    gi.histories.create_dataset_collection(hist, raw_data_collection)
    unmerged_collection = {
        'collection_type': 'list',
        'element_identifiers': unmerged_dataset_ids,
        'name': "Unmerged files"
    }
    gi.histories.create_dataset_collection(hist, unmerged_collection)
    out.write("Organized files into collections\n")
    # Find the data library with the annotations
    lib = gi.libraries.get_libraries(name=config["library_names"]["genome_annotations"])
    assert len(lib) > 0, "No library found for %s lab" % config["library_names"]["genome_annotations"]
    lib_id = lib[0]["id"]
    # Add the annotation file to the history
    for ds in gi.libraries.show_library(lib_id, contents=True):
        # Eliminate the folder
        if ds['type'] != 'file':
            continue
        # Eliminate the files from other folders
        if ds["name"].find(config["folder_names"]["annotation"]) == -1:
            continue
        # Find only the "Mus_musculus.GRCm38.87.gtf (mm10)"
        if ds["name"].find(config["annotation_name"]) == -1:
            continue
        # Add the files to the history
        gi.histories.upload_dataset_from_library(hist, ds["id"])
        annotation_id = ds["id"]
    assert annotation_id != '', "No annotation file for %s" % config["annotation_name"]
    out.write("Imported the annotation file\n")
    # Close log file
    out.close()
