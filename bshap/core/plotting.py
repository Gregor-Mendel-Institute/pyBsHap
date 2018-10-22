## Functions used for plotting
import matplotlib.pylab as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd
import logging

log = logging.getLogger(__name__)

## A joint plot using seaborn

def meths_jointplot(x, y, reqcond, filter_pos=False, kde=True, hexplt=False):
    ## reqcond is a dictionry with values given below.
    # 1) reqcond['color']
    # 2) reqcond['xlab']
    # 3) reqcond['ylab']
    # 4) reqcond['plt_limits']
    # 4) reqcond['size']
    if filter_pos:
        filter_inds = np.where(x + y > 0)[0]
        x = x[filter_inds]
        y = y[filter_inds]
    if 'color' not in reqcond:
        reqcond['color'] = "#43a2ca"
    if 'xlab' not in reqcond:
        reqcond['xlab'] = ''
    if 'ylab' not in reqcond:
            reqcond['ylab'] = ''
    if 'plt_limits' not in reqcond:
        reqcond['plt_limits'] = (-0.05, 1.05)
    sns.set(style="white", color_codes=True)
    if kde:
        if 'size' in reqcond:
            p = sns.jointplot(x = x, y = y, kind = "kde", joint_kws={'gridsize':reqcond['size']}, color=reqcond['color'])
        else:
            p = sns.jointplot(x = x, y = y, kind = "kde", color=reqcond['color'])
        _ = p.ax_marg_x.hist(p.x, bins = np.linspace(reqcond['plt_limits'][0], reqcond['plt_limits'][1], 20), color= reqcond['color'])
        _ = p.ax_marg_y.hist(p.y, bins = np.linspace(reqcond['plt_limits'][0], reqcond['plt_limits'][1], 20), orientation="horizontal", color = reqcond['color'])
    else:
        if hexplt:
            if 'size' in reqcond:
                p = sns.jointplot(x = x, y = y, kind="hex", color = reqcond['color'], joint_kws={'gridsize':reqcond['size']})
            else:
                p = sns.jointplot(x = x, y = y, kind="hex", color = reqcond['color'])
        else:
            if 'size' in reqcond:
                p = sns.jointplot(x = x, y = y, kind = "scatter", joint_kws={"s": reqcond['size']}, alpha = 0.1, color=reqcond['color'])
            else:
                p = sns.jointplot(x = x, y = y, kind = "scatter", alpha = 0.1, color=reqcond['color'])
    p.set_axis_labels(reqcond['xlab'], reqcond['ylab'])
    if 'annotate' in reqcond:
        t_anno = "%s; npts: %s" % (reqcond['annotate'], len(p.x))
        p = p.annotate(stats.pearsonr, template="%s" % t_anno, fontsize=12)
    else:
        p = p.annotate(stats.pearsonr, template="{stat}: {val:.2f}; npts: %s" % len(p.x), fontsize=12)
    p.ax_joint.plot(reqcond['plt_limits'], reqcond['plt_limits'], ':k')
    p.ax_marg_x.axis([reqcond['plt_limits'][0], reqcond['plt_limits'][1], 0, len(p.x)/5])
    p.ax_marg_y.axis([0, len(p.x)/5, reqcond['plt_limits'][0], reqcond['plt_limits'][1]])
    return(p)

class PlotMethylationContexts(object):

    def __init__(self, t_req_gene, cg_thres = 0.08, chg_thres = 0.012, chh_thres = 0.015):
        assert type(t_req_gene) is pd.DataFrame, 'please provide pandas dataframe'
        self.meths = t_req_gene
        self._cg_thres = cg_thres
        self._chg_thres = chg_thres
        self._chh_thres = chh_thres

    def plot_cg_chg_chh(self, color = "#74a9cf"):
        chh_plt_limits = get_context_limits(max( self.meths.iloc[:,2] ), "chh", self._cg_thres, self._chg_thres, self._chh_thres)
        chg_plt_limits = get_context_limits(max( self.meths.iloc[:,1] ), "chg", self._cg_thres, self._chg_thres, self._chh_thres)
        p1 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,1], stat_func=stats.spearmanr, color = color)
        p1.ax_joint.plot((self._cg_thres,self._cg_thres), chg_plt_limits, ':k')
        p1.ax_joint.plot( (-0.01,max( self.meths.iloc[:,0] ) ), (self._chg_thres,self._chg_thres), ':k')
        p2 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,2], stat_func=stats.spearmanr, color = color)
        p2.ax_joint.plot((self._cg_thres,self._cg_thres), chh_plt_limits, ':k')
        p2.ax_joint.plot( (-0.01,max( self.meths.iloc[:,0] ) ), (self._chh_thres,self._chh_thres), ':k')
        self.p1 = p1
        self.p2 = p2

    def add_points(self, plt_fig, x_pts, y_pts, pts_label, pts_color):
        assert hasattr(self, plt_fig)
        x_pts = np.array(x_pts)
        y_pts = np.array(y_pts)
        pts_label = np.array(pts_label)
        pts_color = np.array(pts_color)
        for eind in range(len(x_pts)):
            self.__getattribute__(plt_fig).ax_joint.plot(x_pts[eind], y_pts[eind], 'o', label = pts_label[eind], color = pts_color[eind])
        legend = self.__getattribute__(plt_fig).ax_joint.legend(loc="upper left")

def get_context_limits(max_value, context, cg_thres, chg_thres, chh_thres):
    if max_value < eval(context + "_thres"):
        #return((0, eval(context + "_thres")))
        return((0, 0.2))
    return((-0.01, max_value ))


class SeabornFig2Grid(object):

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())
