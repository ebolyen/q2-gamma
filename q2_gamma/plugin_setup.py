import importlib

import qiime2.plugin as plg
from q2_types.feature_table import FeatureTable, Frequency, PresenceAbsence
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.tree import Phylogeny, Rooted, Hierarchy

import q2_gamma
import q2_gamma.visualizers

plugin = plg.Plugin(
    name='gamma',
    version=q2_gamma.__version__,
    website='https://github.com/ebolyen/q2-gamma',
    package='q2_gamma',
    description='TODO',
    short_description='TODO'
)

plugin.visualizers.register_function(
    function=q2_gamma.visualizers.simple_plot,
    inputs={
        'table': FeatureTable[Frequency | PresenceAbsence],
        'feature_tree': Hierarchy | Phylogeny[Rooted] | FeatureData[Taxonomy]
    },
    parameters={
        'metadata': plg.Metadata,
        'case_where': plg.Str,
        'control_where': plg.Str,
        'n_transects': plg.Int % plg.Range(1, 50),
        'stratify_by': plg.Str
    },
    input_descriptions={
        'table': '',
        'feature_tree': ''
    },
    parameter_descriptions={
        'metadata': '',
        'case_where': '',
        'control_where': ''
    },
    name='TODO',
    description='TODO'
)

importlib.import_module('q2_gamma.transformers');
