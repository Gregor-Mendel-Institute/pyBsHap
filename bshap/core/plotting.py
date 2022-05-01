
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
import pygenome
from . import the1001g, meth5py
import pybedtools as pybed


log = logging.getLogger(__name__)

def scale_colors(minval, maxval, val, safe_colors = None):
    import palettable
    if safe_colors is None:
        safe_colors = palettable.colorbrewer.sequential.BuGn_7.colors
    import sys
    EPSILON = sys.float_info.epsilon  # Smallest possible difference.
    i_f = float(val-minval) / float(maxval-minval) * (len(safe_colors)-1)
    i, f = int(i_f // 1), i_f % 1  # Split into whole & fractional parts.
    if f < EPSILON:
        ret_col = safe_colors[i]
    else:  # Otherwise return a color within the range between them.
        (r1, g1, b1), (r2, g2, b2) = safe_colors[i], safe_colors[i+1]
        ret_col = int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))
    return('#%02x%02x%02x' % (ret_col[0], ret_col[1], ret_col[2]))

np_scale_colors = np.vectorize(scale_colors)


class metaPlot(object):
    """
    Class function to do metaplots on Genes/TEs
        * first divide the genome into windows
        * Denote each window as either +100 or -100 based on the GFF file for genes
        * Then calculate the methylation averages for each of those windows 
        * 

    # gene list.
    """
    def __init__(self, genome_obj, meth_winds_obj):
        assert type(genome_obj) is pygenome.genome.GenomeClass, "please provide a pygenome.genome object"
        # assert type(meth_winds_obj) is methh5py.HDF5MethTable, "please provide methh5py.HDF5MethTable object"
        assert type(meth_winds_obj) is the1001g.ContextsHDF51001gTable, "please provide the1001g.ContextsHDF51001gTable object"
        ## make sure genes and TEs are added to genome object
        self.genome = genome_obj
        self.meth_windows = meth_winds_obj

    def get_windows(self):
        return( pybed.BedTool.from_dataframe(self.meth_windows.cg.get_bed_df(None)) )

    def get_tss_tes_positions(self):
        tss_df = self.genome.genes.iloc[:, [0, 1, 2, 3]].values
        tss_df = pd.DataFrame(tss_df, index = self.genome.genes.index, columns = ["chr", "start", "end", "gene.id"])
        for ef_gene in self.genome.genes.iterrows():
            if ef_gene[1].iloc[5] == "-":
                tss_df.loc[ef_gene[0], 'start'] = ef_gene[1].iloc[2]
        # tss_df = tss_df.iloc[:,[0, 1, 1, 3]]
        tss_df['end'] = tss_df['start'].values
        tss_df = tss_df.sort_values(["chr", "start"])
        # tss_df.columns = np.array(("chr", "start", "end", "gene.id"))
        return( tss_df ) #pybed.BedTool.from_dataframe(tss_df) )

    def get_closest_TSS_and_TES(self):
        #https://dalzer.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.closest.html
        tasdasd


    def then_calculate_methylation_averages(self):
        asdash

def plot_pie(x, ax, r=0.07, colors=['#3182bd', '#bdbdbd']): 
    # radius for pieplot size on a scatterplot
    ax.pie([x['prop'], 1 - x['prop']], center=(x['x'],x['y']), radius=r, colors=colors, startangle = 90, counterclock = False, frame = True)

def classify_array(np_arr, break_point = 0):
    np_arr_str = np.repeat("nan", np_arr.shape[0])
    np_arr_str[ np_arr < break_point ] = "neg"
    np_arr_str[ np_arr > break_point ] = "pos"
    ef_zeros = np.where(np_arr == 0)[0]
    ### Randomly put zeros in either neg or positive bins
    np_arr_str[ ef_zeros ] = np.array(("neg", "pos"))[np.random.binomial(1, 0.5, ef_zeros.shape[0])]
    return(np_arr_str)

def sort_df(df, columns_to_sort = None):
    if columns_to_sort is not None:
        return( df.sort_values( columns_to_sort ) )
    df_sorted_index = df.loc[:,columns_to_sort].sum(axis = 1).argsort().sort_values()
    return( df.loc[df_sorted_index.index,:] )

def effect_plot_with_ci(effects_df):

    ef_str = "eff_temp_"

    ef_data = sort_df(effects_df.loc[effects_df[ef_str + 'diff.CG'] < 0], [ef_str + 'upr.CG']).loc[:,[ef_str + "lwr.CG", ef_str + "diff.CG", ef_str + "upr.CG"]].dropna()
    ef_data = ef_data.append( sort_df(effects_df.loc[effects_df[ef_str + 'diff.CG'] > 0], [ef_str + "lwr.CG"] ).loc[:,[ef_str + "lwr.CG", ef_str + "diff.CG", ef_str + "upr.CG"]].dropna() )

    p = sns.scatterplot( x = np.arange(ef_data.shape[0]), y = 1, marker = "" )
    p.axes.fill_between( np.arange(ef_data.shape[0]), ef_data.iloc[:,0].values, ef_data.iloc[:,2].values, interpolate=False, color = cb.qualitative.Dark2_5.hex_colors[1], alpha = 0.3 )
    p.plot( np.arange(ef_data.shape[0]), ef_data.iloc[:,1].values, '-', color = cb.qualitative.Dark2_5.hex_colors[1])

    p.set_ylim( (-1.5, 1.5) )
    p.set_xlim( (0, ef_data.shape[0]) )
    plt.plot((0, ef_data.shape[0]), (0, 0), '--', c = cb.sequential.Greys_5.hex_colors[3])

    plt.xlabel("#genes")
    plt.ylabel("Temp. effect on CG gbM in F2s")

    plt.show()
    return(p)


def quadrant_plot(x, y, plt_limits = {}, plt_options = {}, **kwargs):
    if 'xlim' not in plt_limits:
        plt_limits['xlim'] = 0.2
    if 'ylim' not in plt_limits:
        plt_limits['ylim'] = 0.2
    if 'xpie' not in plt_limits:
        plt_limits['xpie'] = 0.1
    if 'ypie' not in plt_limits:
        plt_limits['ypie'] = 0.1
    if 'pie_rad' not in plt_limits:
        plt_limits['pie_rad'] = 0.05
    
    if 'color' not in plt_options:
        plt_options['color'] = "#fe9929"
    if 'plt_center' not in plt_options:
        plt_options['plt_center'] = 0
    if 'marginal_bins' not in plt_options:
        plt_options['marginal_bins'] = 200
    if 'annotation' not in plt_options:
        plt_options['annotation'] = "#pts: %s"

    p = sns.jointplot( x=x, y = y, marginal_kws={"bins": plt_options['marginal_bins']}, xlim=(-plt_limits['xlim'], plt_limits['xlim']), ylim=(-plt_limits['ylim'], plt_limits['ylim']), color=plt_options['color'], **kwargs)

    p.ax_joint.plot((plt_options['plt_center'],plt_options['plt_center']), (-plt_limits['ylim'], plt_limits['ylim']), '--', color = "#636363")
    p.ax_joint.plot((-plt_limits['xlim'],plt_limits['xlim']), (plt_options['plt_center'], plt_options['plt_center']), '--', color = "#636363")

    t_effects = pd.Series(classify_array(x, plt_options['plt_center'])) + "_" + pd.Series(classify_array(y, plt_options['plt_center']))
    t_effects = t_effects.value_counts()

    plot_pie( {"prop": t_effects['neg_pos']/float(len(x)), 'x': -plt_limits['xpie'], 'y': plt_limits['ypie'] }, r = plt_limits['pie_rad'], ax = p.ax_joint, colors = [plt_options['color'], "#cccccc"] )
    plot_pie( {"prop": t_effects['pos_pos']/float(len(x)), 'x': plt_limits['xpie'], 'y': plt_limits['ypie'] }, r = plt_limits['pie_rad'], ax = p.ax_joint, colors = [plt_options['color'], "#cccccc"] )
    plot_pie( {"prop": t_effects['neg_neg']/float(len(x)), 'x': -plt_limits['xpie'], 'y': -plt_limits['ypie'] }, r = plt_limits['pie_rad'], ax = p.ax_joint, colors = [plt_options['color'], "#cccccc"] )
    plot_pie( {"prop": t_effects['pos_neg']/float(len(x)), 'x': plt_limits['xpie'], 'y': -plt_limits['ypie'] }, r = plt_limits['pie_rad'], ax = p.ax_joint, colors = [plt_options['color'], "#cccccc"] )

    t_anno = plt_options['annotation'] % len(x)
    # \nchi-square pval: %s sp.stats.chisquare(t_effects.values, ef_data.shape[0] * np.repeat(0.25, 4))[1]) # ; pearsonr: %.2f, sp.stats.pearsonr(p.x, p.y)[0] )
    plt.text( 0.1, 0.9, t_anno, transform=p.ax_joint.transAxes, size = 15 )
    # p.ax_joint.text(0, 1, t_anno, ha='center', va='center', )
    # text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)
    return(p)

## A joint plot using seaborn

def meths_jointplot(x, y, reqcond = {}, axs = None, min_sum=None, **kwargs):
    ## reqcond is a dictionry with values given below.
    # 1) reqcond['color']
    # 2) reqcond['xlab']
    # 3) reqcond['ylab']
    # 4) reqcond['plt_limits']
    # 4) reqcond['size']
    if min_sum is not None:
        filter_inds = np.where(x + y > min_sum)[0]
        x = pd.Series(x).iloc[filter_inds]
        y = pd.Series(y).iloc[filter_inds]
    if 'color' not in reqcond:
        reqcond['color'] = "#43a2ca"
    if 'plt_limits' not in reqcond:
        reqcond['plt_limits'] = (-0.05, 1.05)
    if 'annotate' not in reqcond: 
        reqcond['annotate'] = False
    if "alpha" not in reqcond:
        reqcond['alpha'] = 0.1
    marginal_kws = { "bins": np.linspace(reqcond['plt_limits'][0], reqcond['plt_limits'][1], 20)  }
    if "kind" in reqcond: 
        if "size" in reqcond:
            p = sns.jointplot(x = x, y = y, kind = reqcond['kind'], joint_kws={'gridsize':reqcond['size']}, color=reqcond['color'], marginal_ticks=True, **kwargs)
        else:
            p = sns.jointplot(x = x, y = y, kind = reqcond['kind'], color=reqcond['color'], marginal_ticks=True, **kwargs)
    else:
        if 'size' in reqcond:
            p = sns.jointplot(x = x, y = y, kind = "scatter", joint_kws={'s':reqcond['size']}, color=reqcond['color'], alpha = reqcond['alpha'], marginal_ticks=True, **kwargs)
        else:
            p = sns.jointplot(x = x, y = y, kind = "reg", joint_kws = {'scatter_kws':dict(alpha=reqcond['alpha'])}, color=reqcond['color'], marginal_ticks=True, **kwargs)
    # p.set_axis_labels(reqcond['xlab'], reqcond['ylab'])
    # _ = p.ax_marg_x.hist(p.x, bins = np.linspace(reqcond['plt_limits'][0], reqcond['plt_limits'][1], 20), color= reqcond['color'])
    # _ = p.ax_marg_y.hist(p.y, bins = np.linspace(reqcond['plt_limits'][0], reqcond['plt_limits'][1], 20), orientation="horizontal", color = reqcond['color'])
    if reqcond['annotate']:
        t_anno = "npts: %s; pearsonr: %.2f" % (len(p.x), stats.pearsonr( x, y )[0])
        plt.text( 0.1, 0.9, t_anno, transform=p.ax_joint.transAxes, size = 12 )
        # p.ax_joint.text(0.1, 0.9, t_anno, ha='center', va='center', size = 12 )
        # p = p.annotate(stats.pearsonr, template="{stat}: {val:.2f}; npts: %s" % len(p.x), fontsize=12)
    p.ax_joint.plot(reqcond['plt_limits'], reqcond['plt_limits'], ':k')
    p.ax_joint.set_xlim(reqcond['plt_limits'])
    p.ax_joint.set_ylim(reqcond['plt_limits'])
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
        if isinstance(color, basestring):
            p1 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,1], stat_func=stats.spearmanr, color = color)
            p2 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,2], stat_func=stats.spearmanr, color = color)
        else:
            p1 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,1], stat_func=stats.spearmanr, scatter = True, color = None)
            p1.ax_joint.scatter(x = self.meths.iloc[:,0], y = self.meths.iloc[:,1], c = pd.Series(color))
            p2 = sns.jointplot(x = self.meths.iloc[:,0], y = self.meths.iloc[:,2], stat_func=stats.spearmanr, scatter = True, color = None)
            p2.ax_joint.scatter(x = self.meths.iloc[:,0], y = self.meths.iloc[:,2], c = pd.Series(color))
        p1.ax_joint.plot((self._cg_thres,self._cg_thres), chg_plt_limits, ':k')
        p1.ax_joint.plot( (-0.01,max( self.meths.iloc[:,0] ) ), (self._chg_thres,self._chg_thres), ':k')
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
