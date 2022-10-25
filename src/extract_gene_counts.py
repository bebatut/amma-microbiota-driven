#!/usr/bin/env python

import yaml
from bioblend.galaxy import GalaxyInstance

from interact_galaxy import *

if __name__ == '__main__':
    '''
    Apply the worklow to extract the gene counts on the data collection
    '''
    # Get config
    with open('config.yaml', 'r') as config_f:
        config = yaml.safe_load(config_f)
    # Connect to Galaxy and retrieve the history
    gi = GalaxyInstance(config["galaxy_url"], config["api_key"])
    hist = get_hist(gi, config["hist_name"])
    out = open(str("logs/extract_gene_counts.log"), "w")
    # Search for the collection id with the raw data
    raw_data_coll_id = get_working_collection_id(config["collection_names"]["raw_data"], hist, gi)
    # Get the workflow id or import it if it is not there
    wf_id = get_wf_id(config["workflow_names"]["extract_gene_counts"], gi)
    if wf_id != '':
        wf_id = gi.workflows.import_workflow_from_local_path("src/from_sequences_to_counts.ga")['id']
    out.write("Got workflow id\n")
    # Launch the workflow
    annotation_id = get_annotation_id(config["annotation_name"], hist, gi)
    inputs = {'0': {'src':'hdca', 'id': raw_data_coll_id}, '1': {'id': annotation_id, 'src': 'hda'}}
    invoc = gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=hist)
    out.write("Invoked the workflow\n")
    # Close log file
    out.close()
