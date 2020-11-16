import logging
log = logging.getLogger(__name__)
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from bshap import meth5py
from bshap.core import combinemeths

import seaborn as sns
sns.set(style="white", color_codes=True)

import palettable
context_colors = palettable.colorbrewer.qualitative.Set2_3.hex_colors



sra_names = ["stem_positive_7days", "stem_negative_7days", "stem_positive_14days", "stem_negative_14days", "stem_positive_35days", "stem_negative_35days"]
sra_files = [ "allc_" + i + ".hdf5"  for i in sra_names]
i = (2, 3) ## 14days stem vs nonstem


meths_file_list = [sra_files[i[0]], sra_files[i[1]]]
cmeths = combinemeths.CombinedMethsTable(meths_file_list)


read_threshold = 10
## Plotting for all the chromosomes.
cmeths.derive_methylated_identical_pos_ix( read_depth=read_threshold )
