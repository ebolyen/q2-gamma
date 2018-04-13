import os
import pkg_resources
import json
from distutils.dir_util import copy_tree
import collections
import math
import base64
import bz2

import biom
import skbio
import numpy as np
import qiime2 as q2

PLOT = pkg_resources.resource_filename('q2_gamma', 'visualizers/plot')


def plot(output_dir, table: biom.Table, metadata: q2.Metadata,
         case_where: str, control_where: str,
         feature_tree: skbio.TreeNode=None):

    with open('/tmp/tree.nwk', 'w') as fh:
        feature_tree.write(fh)

    copy_tree(os.path.join(PLOT, 'assets', 'dist'), output_dir)
    data_dir = os.path.join(output_dir, 'data')
    os.mkdir(data_dir)

    metadata = metadata.filter_ids(table.ids(axis='sample'))
    case_samples = sorted(list(metadata.get_ids(case_where)))
    control_samples = sorted(list(metadata.get_ids(control_where)))

    table.filter(case_samples + control_samples)
    table.remove_empty('observation')
    features = list(table.ids(axis='observation'))

    if feature_tree is not None:
        feature_tree = shear_no_prune(feature_tree, features)
    else:
        feature_tree = TreeNode()

    tree_data = tree_to_array(feature_tree)
    idx, = np.where(np.asarray(tree_data['children']) == 0)
    tree_data['lookup'] = dict(zip(map(str, idx), range(len(idx))))

    tip_order = np.asarray(tree_data['names'])[idx]
    table = table.sort_order(tip_order, axis='observation')
    table = table.sort_order(case_samples + control_samples, axis='sample')

    with open(os.path.join(data_dir, 'packed_table.jsonp'), 'w') as fh:
        fh.write('LOAD_PACKED_TABLE(')
        fh.write(json.dumps(table_to_b64pa(table)))
        fh.write(');')

    with open(os.path.join(data_dir, 'tree.jsonp'), 'w') as fh:
        fh.write('LOAD_TREE(')
        fh.write(json.dumps(tree_data))
        fh.write(');')


def shear_no_prune(tree, names):
    tcopy = tree.deepcopy()
    all_tips = {n.name for n in tcopy.tips()}
    ids = set(names)

    if not ids.issubset(all_tips):
        raise ValueError("ids are not a subset of the tree.")

    marked = set()
    for tip in tcopy.tips():
        if tip.name in ids:
            marked.add(tip)
            for anc in tip.ancestors():
                if anc in marked:
                    break
                else:
                    marked.add(anc)

    for node in list(tcopy.traverse()):
        if node not in marked:
            node.parent.remove(node)

    return tcopy


def table_to_b64pa(table):
    array = table.pa().matrix_data.astype(np.uint8).toarray()

    row_uint8_len = array.shape[0]
    col_uint8_len = math.ceil(array.shape[1]/8)  # np.packbits pads to nearest byte boundary

    row_pad = (32 - row_uint8_len % 32) % 32
    col_pad = (4 - col_uint8_len % 4) % 4

    row_uint8_len += row_pad
    col_uint8_len += col_pad

    row_uint32_len = row_uint8_len
    col_uint32_len = col_uint8_len // 4

    packed_array = np.packbits(array, axis=1)
    padded_array = np.pad(packed_array, [(0, row_pad), (0, col_pad)], mode='constant')
    bytes_ = padded_array.tobytes()
    data = base64.b64encode(bz2.compress(bytes_)).decode('ascii')
    return {
        'shape': {
            'bitmatrix': list(table.shape),
            'uint8': [row_uint8_len, col_uint8_len],
            'uint32': [row_uint32_len, col_uint32_len]
        },
        'encoding': 'base64+bzip2',
        'data': data
    }


def tree_to_array(tree):
    children = []
    parents = []
    names = []
    distances = []

    tree.parent_idx = 0
    idx = 0
    offset = 0
    queue = collections.deque([tree])
    while queue:
        node = queue.pop()

        distance_to_longest_tip = 0
        for tip in node.tips():
            distance_to_longest_tip = max(distance_to_longest_tip,
                                          node.distance(tip))
        for child in node.children:
            child.parent_idx = idx

        children.append(1 + offset if node.children else 0)
        parents.append(node.parent_idx)
        if not node.name:
            node.name = ""
        names.append(node.name.strip())
        distances.append(distance_to_longest_tip)

        offset += len(node.children)
        idx += 1
        queue.extendleft(node.children)

    return {
        "children": children,
        "parent": parents,
        "names": names,
        "distances": distances
    }

def collapse(x, depth):
    children = x['children']
    distances = x['distances']
    mapping = {}

    search = collections.deque([None])
    start, end, g = 0, 0, []
    pending_start = 0
    pending_g = []
    for i, child in enumerate(children):
        # finalize last iteration
        if pending_start and child:
            search.appendleft([pending_start, child, pending_g])
            pending_start = 0
        if i == end:
            search.pop()
        if search:
            start, end, g = search[-1]

        # current iteration
        if child:
            if start <= i < end:
                pending_g = g
                pending_start = child
            elif distances[i] <= depth:
                pending_g = []
                mapping[i] = pending_g
                pending_start = child
            else:
                pass
        else:
            if start <= i < end:
                g.append(i)
            elif pending_start and pending_start <= i:
                pending_g.append(i)
            else:
                mapping[i] = [i]

    return mapping
