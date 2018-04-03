import os
import pkg_resources
from distutils.dir_util import copy_tree

import biom
import skbio
import qiime2 as q2

PLOT = pkg_resources.resource_filename('q2_gamma', 'visualizers/plot')


def plot(output_dir, table: biom.Table, metadata: q2.Metadata,
         case_where: str, control_where: str,
         feature_tree: skbio.TreeNode=None):

    copy_tree(os.path.join(PLOT, 'assets', 'dist'), output_dir)

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
