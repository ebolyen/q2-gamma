import csv

from q2_gamma.plugin_setup import plugin

import skbio
from q2_types.feature_data import TSVTaxonomyFormat


@plugin.register_transformer
def _0(ff: TSVTaxonomyFormat) -> skbio.TreeNode:
    root = skbio.TreeNode('root', length=0)
    with ff.open() as fh:
        reader = iter(csv.reader(fh, delimiter='\t'))
        next(reader)  # skip header
        for row in reader:
            id_, taxonomy = row[:2]
            taxonomy = taxonomy.split(';')
            node = root
            for taxon in taxonomy:
                for child in node.children:
                    if child.name == taxon:
                        node = child
                        break
                else:
                    child = skbio.TreeNode(taxon, length=1)
                    node.append(child)
                    node = child

            node.append(skbio.TreeNode(id_, length=1))

    return root



