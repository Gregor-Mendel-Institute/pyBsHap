## Functions used for plotting
import matplotlib.pylab as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from scipy import stats
import numpy as np
import logging

log = logging.getLogger(__name__)

## A joint plot using seaborn

def meths_jointplot(x, y, reqcond, kde=True, hexplt=False):
    ## reqcond is a dictionry with values given below.
    # 1) reqcond['color']
    # 2) reqcond['xlab']
    # 3) reqcond['ylab']
    # 4) reqcond['plt_limits']
    # 4) reqcond['size'] -- not yet implemented
    import seaborn as sns
    sns.set(style="white", color_codes=True)
    if kde:
        p = sns.jointplot(x = x, y = y, kind = "kde", joint_kws={'gridsize':30}, color=reqcond['color'])
        _ = p.ax_marg_x.hist(p.x, bins = np.arange(0, 1.1, 0.05), color= reqcond['color'])
        _ = p.ax_marg_y.hist(p.y, bins = np.arange(0, 1.1, 0.05), orientation="horizontal", color = reqcond['color'])
    else:
        if hexplt:
            p = sns.jointplot(x = x, y = y, kind="hex", color = reqcond['color'], joint_kws={'gridsize':15})
        else:
            p = sns.jointplot(x = x, y = y, size = 5, kind = "scatter", joint_kws={"s": 4}, alpha = 0.1, color=reqcond['color'])
    p.ax_joint.plot(reqcond['plt_limits'], reqcond['plt_limits'], ':k')
    p.set_axis_labels(reqcond['xlab'], reqcond['ylab'])
    p = p.annotate(stats.pearsonr, template="{stat}: {val:.2f}; nCs: %s" % len(p.x), fontsize=12)
    p.ax_marg_x.axis([reqcond['plt_limits'][0], reqcond['plt_limits'][1], 0, len(p.x)/5])
    p.ax_marg_y.axis([0, len(p.x)/5, reqcond['plt_limits'][0], reqcond['plt_limits'][1]])
    return(p)
