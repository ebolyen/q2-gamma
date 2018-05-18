import os
import collections
import itertools

import biom
import skbio
import ternary
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import qiime2 as q2

from q2_gamma.geommed import geometric_median
from q2_feature_table import group


def simple_plot(output_dir, table: biom.Table, feature_tree: skbio.TreeNode,
                metadata: q2.Metadata, case_where: str, control_where: str,
                n_transects: int=10, stratify_by: str=None):
    print("Data extracted")

    metadata = metadata.filter_ids(table.ids(axis='sample'))
    case_samples = sorted(list(metadata.get_ids(case_where)))
    control_samples = sorted(list(metadata.get_ids(control_where)))
    get_pairs = comparisons(metadata, control_samples, case_samples, stratify_by)

    table.filter(case_samples + control_samples)
    table.remove_empty('observation')
    features = list(table.ids(axis='observation'))
    feature_tree = shear_no_prune(feature_tree, features)
    print("Extraneous features removed")

    for n in feature_tree.traverse():
        if not n.length:
            n.length = 0
    tree = tree_to_array(feature_tree)
    print("Tree index created")

    possible_transects = len(np.unique(np.asarray(tree['distances'])))
    tree_length = tree['distances'][0]  # the root
    if n_transects > possible_transects:
        n_transects = possible_transects
        print("Only %d transects exist, using that instead" % n_transects)

    transects = list(np.linspace(0, tree_length, num=n_transects))
    print("Will transect at: %s" % ", ".join(map(str, transects)))

    figure_gen = prepare_plot(tree_length)
    figure_gen.send(None)  # initialize co-routine

    feature_ids = table.ids(axis='observation')
    points, ranks = pairwise_components(table, get_pairs())
    write_ranks(output_dir, pd.Series(feature_ids, index=feature_ids), ranks, None)
    figure_gen.send(points)

    collapsed_groups = pd.DataFrame()
    rank_files = []
    for distance in transects:
        collapsed_table, collapsed_counts, groups = group_by_transect(
            table, tree, distance)
        collapsed_groups[groups.name] = groups
        print("Table collapsed at transect %s" % distance)

        points, ranks = pairwise_components(collapsed_table, get_pairs())

        filename = write_ranks(output_dir, collapsed_counts, ranks, distance)
        rank_files.append(filename)

        figure_gen.send((points, distance))

    print("Finalizing visualization")
    figure = figure_gen.send((None, None))
    figure.savefig(os.path.join(output_dir, 'figure.png'))

    with open(os.path.join(output_dir, 'collapsed_groups.tsv'), 'w') as fh:
        collapsed_groups.to_csv(fh, sep='\t')

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write("""<html><body><img src="figure.png" />""")
        fh.write("<p><a href='collapsed_groups.tsv'>collapsed_groups.tsv</a></p>")
        fh.write("""<ul>""")
        for filename in rank_files:
            fh.write("<li><a href='%s'>%s</a></li>" % (filename, filename))
        fh.write("</ul></body></html>")


def write_ranks(output_dir, collapsed_counts, ranks, distance):
    name = "feature_ranks_no_transect.tsv"
    if distance is not None:
        name = "feature_ranks_transect_%s.tsv" % distance

    path = os.path.join(output_dir, name)

    rank_table = collapsed_counts.to_frame()
    rank_table.index.name = "id"
    # These vectors are aligned to the table, which is aligned to this series
    # from the group call
    rank_table['Case'] = ranks[0]
    rank_table['Control'] = ranks[1]
    rank_table['Shared'] = ranks[2]
    rank_table['Missing'] = ranks[3]
    with open(path, 'w') as fh:
        rank_table.to_csv(fh, sep='\t')

    return name


def prepare_plot(tree_length):
    figure, tax = ternary.figure(scale=1.0)
    figure.set_size_inches(9, 9)
    tax.set_title("", fontsize=20)
    tax.bottom_axis_label("Case", fontsize=14)
    tax.left_axis_label("Control", fontsize=14)
    tax.right_axis_label("Shared", fontsize=14)
    tax.ticks(axis='lbr', linewidth=1, multiple=0.1, tick_formats="%.1f")
    tax.clear_matplotlib_ticks()

    tax.gridlines(multiple=0.05, color="grey")
    tax.line((0, 1, 0), (0.5, 0, 0.5), linewidth=.5, color='black', linestyle="-")
    line = []

    t_null = yield
    tax.scatter(t_null, marker='.', s=10, color='red')
    t_null_center = [tuple(geometric_median(np.asarray(t_null)))]
    tax.scatter(t_null_center, marker='X', color='red', label="Original", s=50)
    line += t_null_center

    points, distance = yield
    while points is not None:
        color = cm.viridis(distance/tree_length)
        tax.scatter(points, marker='.', color=color, s=10)
        c = [tuple(geometric_median(np.asarray(points)))]
        line += c
        tax.scatter(c, marker='X', color=color, label="T_%s" % distance, s=50)

        points, distance = yield

    tax.plot(line, linewidth=2, color='red', label='trajectory')
    tax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figure.set_tight_layout(True)
    yield figure


def group_by_transect(table, tree, distance):
    otu_map = collapse_tree(tree, distance)
    use_names = not (not tree['names'][1] or is_float(tree['names'][1]))
    def name(otu):
        if use_names:
            i = otu
            names = []
            while i != 0:
                names.append(tree['names'][i])
                i = tree['parent'][i]
            if names:
                return '; '.join(reversed(names))

        return "T_%s:%d" % (distance, otu)

    groups = pd.Series(
        {tree['names'][tip]: name(otu)
            for otu, tips in otu_map.items() for tip in tips},
        name='T_%s' % distance).sort_index()
    groups.index.name = 'id'
    collapsed_counts = pd.Series(
        {name(otu): len(tips) for otu, tips in otu_map.items()},
        name='Cluster Size').sort_index()

    column = q2.CategoricalMetadataColumn(groups)
    collapsed_table = group(table, axis='feature', metadata=column, mode='sum')

    return collapsed_table, collapsed_counts, groups

def is_float(n):
    try:
        float(n)
    except:
        return False
    else:
        return True

def pairwise_components(table, comparisons):
    points = []
    feature_count = table.length('observation')
    r_g1_and_g2 = np.zeros(feature_count, dtype=int)
    r_g2_sans_g1 = np.zeros(feature_count, dtype=int)
    r_g1_sans_g2 = np.zeros(feature_count, dtype=int)
    r_sans_g1_sans_g2 = np.zeros(feature_count, dtype=int)

    for s1, s2 in comparisons:
        v1 = table.data(s1).astype(bool)
        v2 = table.data(s2).astype(bool)

        cF = (v1 | v2).sum()

        v1_and_v2 = v1 & v2
        r_g1_and_g2 += v1_and_v2
        c1_and_c2 = v1_and_v2.sum() / cF

        v2_sans_v1 = v2 & ~v1
        r_g2_sans_g1 += v2_sans_v1
        c2_sans_c1 = v2_sans_v1.sum() / cF

        v1_sans_v2 = v1 & ~v2
        r_g1_sans_g2 += v1_sans_v2
        c1_sans_c2 = v1_sans_v2.sum() / cF

        r_sans_g1_sans_g2 += ~v1 & ~v2

        points.append((c2_sans_c1, c1_and_c2, c1_sans_c2))

    return points, (r_g2_sans_g1, r_g1_sans_g2, r_g1_and_g2, r_sans_g1_sans_g2)


def comparisons(metadata, control_samples, case_samples, stratify_by):
    if stratify_by is None:
        def iterator():
            yield from itertools.product(control_samples, case_samples)
    else:
        control_samples = set(control_samples)
        case_samples = set(case_samples)
        groups = metadata.to_dataframe().groupby(stratify_by)
        def iterator():
            for _, df in groups:
                in_group = set(df.index)
                yield from itertools.product(control_samples & in_group,
                                             case_samples & in_group)
    return iterator

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

def collapse_tree(x, depth):
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
