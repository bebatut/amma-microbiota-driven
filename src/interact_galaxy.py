from bioblend.galaxy import GalaxyInstance


def get_tool_id(tool_ids, tools):
    '''
    Retrieve the id of a tool
    '''
    tool_id = ''
    for tool in tools:
        if tool["id"] == tool_ids:
            tool_id = tool["id"]
    return tool_id


def get_working_tool_id(tool_name, tools):
    '''
    Retrieve the id of a tool and test if it exists
    '''
    tool_id = get_tool_id(tool_name, tools)
    assert tool_id != '', "No %s tool" % tool_name
    return tool_id


def get_collection_id(collection_name, hist_id, gi):
    '''
    Retrieve the id of a collection
    '''
    coll_id = ''
    for ds in gi.histories.show_history(hist_id, contents=True, visible = True):
        if ds["history_content_type"] != 'dataset_collection':
            continue
        if ds["name"] == collection_name:
            coll_id = ds["id"]
    return coll_id

def get_working_collection_id(coll_name, hist, gi):
    '''
    Retrieve the id of a collection and test if it exists
    '''
    coll_id = get_collection_id(coll_name, hist, gi)
    assert coll_id != '', "No collection for %s" % coll_name
    return coll_id


def fill_multiqc_inputs(coll_id,hist, gi):
    '''
    Extract the list of element in the collection and format it as input for tool
    '''
    inputs = []
    for el in gi.histories.show_dataset_collection(hist, coll_id)["elements"]:
        inputs.append({'src':'hda','id': el["object"]["id"]})
    return inputs

def get_annotation_id(annot_name, hist, gi):
    '''
    Extract the dataset id in the history
    '''
    annotation_id = ''
    for ds in gi.histories.show_history(hist, contents=True, visible=True):
        if ds["name"].find(annot_name) != -1:
            annotation_id = ds["id"]
    assert annotation_id != '', "No annotation file for %s in history" % annot_name
    return annotation_id


def rename_generated_collection(info, new_name, visible):
    '''
    Rename a generated collection
    '''
    col_id = get_output_collection_ids(info)
    assert type(col_id) == str
    gi.histories.update_dataset_collection(hist, col_id, name=new_name, visible=visible)
    # Rename the generated files too
    for out in info['outputs']:
        out_id = out['id']
        prov = gi.histories.show_dataset_provenance(hist, out_id)
        print(prov)
        if not 'input|__identifier__' in prov['parameters']:
            print("Issue to find a correct name")
            print(prov)
            continue
        new_name = prov['parameters']['input|__identifier__'].replace('"','')
        gi.histories.update_dataset(hist, out_id, name=new_name)
    return col_id


def get_hist(gi, hist_name):
    '''
    Get history id based on a name
    '''
    histories = gi.histories.get_histories()
    hist = ''
    for history in histories:
        if history["name"] == hist_name:
            hist = history['id']
    if hist == '':
        hist = gi.histories.create_history(hist_name)["id"]
    return hist


def get_wf_id(wf_name, gi):
    '''
    Get worklow id based on its name
    '''
    workflows = gi.workflows.get_workflows()
    wf_id = ""
    for wf in workflows:
        if wf["name"] == wf_name:
            wf_id = wf["id"]
    return wf_id
