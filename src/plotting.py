import os,sys
import numpy as np
import warnings
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib import colors
from matplotlib import ticker as t
from matplotlib import patches
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

# Turn off warnings
warnings.filterwarnings("ignore", module="matplotlib")

# A discrete colormap
cdict1 = {'red':   ((0.0, 0.0, 0.0),
                    #((0.124, 255/255., 255/255.),
                    (0.125, 50/255., 50/255.),
                    (0.25, 102/255., 102/255.),
                    (0.375, 171/255., 171/255.),
                    (0.5, 230/255., 230/255.),
                    (0.625, 254/255., 254/255.),
                    (0.75, 253/255., 253/255.),
                    (0.875, 244/255., 244/255.),
                    (1.0, 213/255., 213/255.)),

          'green': ((0.0, 0.0, 0.0),
                    #((0.124, 255/255., 255/255.),
                    (0.125, 136/255., 136/255.),
                    (0.25, 194/255., 194/255.),
                    (0.375, 221/255., 221/255.),
                (0.5, 245/255., 245/255.),
                    (0.625, 224/255., 224/255.),
                    (0.75, 174/255., 174/255.),
                    (0.875, 109/255., 109/255.),
                    (1.0, 62/255., 62/255.)),

          'blue':   ((0.0, 0.0, 0.0),
                    #((0.124, 255/255., 255/255.),
                     (0.125, 189/255., 189/255.),
                     (0.25, 165/255., 165/255.),
                     (0.375, 164/255., 164/255.),
                     (0.5, 152/255., 152/255.),
                     (0.625, 139/255., 139/255.),
                    (0.75, 97/255., 97/255.),
                    (0.875, 67/255., 67/255.),
                     (1.0, 79/255., 79/255.)),
          }
discrete = colors.LinearSegmentedColormap('BlueRed1', cdict1) 

def plot_fitting_parameter_histograms(datalist, titlelist, bins=100):
    """
    Plotting a figure with 5 histograms.
    """
    fig = plt.figure()
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(234)
    ax4 = fig.add_subplot(235)
    ax5 = fig.add_subplot(133)

    ax1.hist(datalist[0].flatten(), bins=bins)
    ax1.set_title(titlelist[0])
    ax2.hist(datalist[1].flatten(), bins=bins)
    ax2.set_title(titlelist[1])
    ax3.hist(datalist[2].flatten(), bins=bins)
    ax3.set_title(titlelist[2])
    ax4.hist(datalist[3].flatten(), bins=bins)
    ax4.set_title(titlelist[3])
    ax5.hist(datalist[4].flatten(), bins=bins)
    ax5.set_title(titlelist[4])

    plt.show()

class Plot(object):
    """
    A plotting class.
    """
    def __init__(self, save_pdf=False, save_png=False, save_ps=False, savename=None, text_width=16, window_width=800,rows=1, cols=1, max_nr=None, exclude=[], aspect=1., header=False, colorbar=False, colorbar_left=False, legend=False, legend_location=1, legend_frameon=True, fontsize=11, border_in=0.05, border_out=0.1, axes_visible=True, usetex=False):

        # SAVE PLOT AS PDF or PNG?
        self.save_pdf  = save_pdf
        self.save_png  = save_png
        self.save_ps   = save_ps
        self.save_plot = save_pdf | save_png | save_ps

        # SWITCH BACKEND
        #if self.save_ps: plt.switch_backend('PS')
        #else:
        #    if plt.get_backend() != 'TkAgg':
        #        plt.switch_backend('TkAgg')

        # MATPLOTLIB MODULES
        self.patches = patches

        # SET PROPERTIES OF CANVAS
        self.dpi = 100 # fix the dpi (keeps also all labels with same size)
        self.text_width = text_width # the desired with of the figure in cm
        self.window_width = window_width # the desired width of the figure in px
        self.aspect = aspect # aspect ratio (height/width) of plot(s)
        self.rows = rows # nr. of rows
        self.cols = cols # nr. of columns
        self.max_nr = max_nr # maximum nr. of axes
        self.exclude = exclude # A list of axes to exclude
        self.border_inside = border_in # border between subfigures
        self.border_outside = border_out # border around the figure(s)

        # SET PROPERTIES OF AXES
        self.title_label = [r'self.title\_label'] * cols * rows
        self.xlabel = [r'self.xlabel'] * cols * rows
        self.ylabel = [r'self.ylabel'] * cols * rows
        self.axes_visible = [axes_visible] * cols * rows

        # SET PROPERTIES OF HEADER
        self.header = header # show header on top
        self.header_height = 0.1 # height of header (as fraction of total height)

        # SET PROPERTIES OF COLORBAR
        self.colorbar = colorbar # show colorbar
        self.colorbar_left = colorbar_left # show colorbar on the left
        self.colorbar_width = 0.1 # width of colorbar (as fraction of total width)
        self.colorbar_label = r'self.colorbar\_label'

        # SET MORE PROPERTIES
        self.fontsize = fontsize
        self.legend = legend
        self.legend_location = legend_location
        self.legend_frameon  = legend_frameon

        # USE LATEX FONT
        if usetex:
            rc('text', usetex=True)
            rc('text.latex', preamble = '\usepackage{amsmath}\usepackage{amssymb}\usepackage{units}\usepackage[utf8x]{inputenc}\usepackage{kerkis}\usepackage{upgreek}')
            rc('font',**{'family':'sans-serif', 'serif':['Kerkis'], 'sans-serif':['Kerkis']})

        # CREATE FIGURE CANVAS
        self.figsize = self.get_proper_figsize()
        self.fig = plt.figure(figsize=self.figsize, dpi=self.dpi)

        # CREATE FIGURE(S)
        self.set_layout()

    def get_proper_figsize(self):

        # adjust aspect ratio due to multiple subfigures
        self.hfract = self.rows + (self.rows-1)*self.border_inside
        self.wfract = self.cols + (self.cols-1)*self.border_inside
        aspect = self.aspect * (self.hfract  / self.wfract )

        # adjust aspect ratio due to additional canvas elements
        self.width = 1. - (2*self.border_outside + 0.025 + self.colorbar*self.colorbar_width)
        self.height = 1. - (2*self.border_outside + self.header*self.header_height)
        aspect *= (self.width  / self.height)

        # use aspect ratio to define proper figure size
        if self.save_plot: figure_width = self.text_width*0.393701 # fixed figure width in cm
        else: figure_width = self.window_width / self.dpi # fixed figure width in px
        figure_height = figure_width * aspect
        return (figure_width, figure_height)


    # LAYOUT FUNCTIONS #
    ####################
    def set_layout(self):

        x0 = self.border_outside + 0.025 + self.colorbar_left*self.colorbar_width
        y0 = self.border_outside

        w = self.width / self.wfract
        h = self.height / self.hfract
        b = self.border_inside

        if self.max_nr is None: self.max_nr = self.rows*self.cols
        self.axes = [self.fig.add_axes([x0+i*w*(b+1),y0+j*h*(b+1),w,h])  for j in range(self.rows-1, -1, -1) for i in range(0, self.cols) if (((self.rows-j)*self.cols - self.cols + i) < self.max_nr) & (((self.rows-j)*self.cols - self.cols + i) not in self.exclude) ]
        self.nr_axes = len(self.axes)

        if self.header: self.hax = self.fig.add_axes([x0, 1. - (self.header_height+self.border_outside) + b*h, self.width, self.header_height])
        if self.colorbar: 
            if self.colorbar_left:
                self.cax = self.fig.add_axes([self.border_outside + 0.025, y0, 0.25*self.colorbar_width , self.height])
            else:
                self.cax = self.fig.add_axes([1. - (self.colorbar_width + self.border_outside) + 0.1*w, y0, 0.25*self.colorbar_width , self.height])

        # different layout depending on purpose
        self.axes = [self.set_axes_layout(ax) for ax in self.axes]

    def add_axes(self, (i,j), nw, nh, padx=0, pady=0, wfrac=1, hfrac=1):
        x0 = self.border_outside + 0.025 + self.colorbar_left*self.colorbar_width + padx
        y0 = self.border_outside + pady

        w = wfrac * self.width / self.wfract
        h = hfrac * self.height / self.hfract
        b = self.border_inside

        self.axes.append(self.fig.add_axes([x0+i*w*(b+1), y0+j*h*(b+1),w*nw+(nw-1)*w*b, h*nh+(nh-1)*h*b]))
        self.nr_axes = len(self.axes)
        self.axes = [self.set_axes_layout(ax) for ax in self.axes]

    def add_twiny_axes(self, axid):
        self.axes.append(self.axes[axid].twiny())
        self.axes = [self.set_axes_layout(ax) for ax in self.axes]
        
    def set_axes_layout(self, ax, visible=True):
        if not visible:
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_xticks([])
        ax.tick_params(labelsize=self.fontsize)
        return ax

    def set_header_layout(self, title=None, color='k', border=1, fontcolor='k', fontsize=None, titlepos='center'):
        self.hax.xaxis.set_visible(False)
        self.hax.yaxis.set_visible(False)
        if fontsize is None: fontsize = 1.5*self.fontsize
        plt.setp(self.hax.spines.values(), color=color, lw=border)

        if title is not None:
            if titlepos=='left': pos=0.05
            elif titlepos=='right': pos=0.95
            else: pos=0.5
            self.hax.text(pos,0.5, title, color=fontcolor, fontsize=fontsize, va='center', ha=titlepos, transform=self.hax.transAxes)

    def add_overview_panel_to_header(self, pos='left', facecolor='k', fontcolor='w'):
        aspect = (self.figsize[1] / self.figsize[0]) * (self.header_height / self.width)
        hborder = 0.1
        wborder = hborder * aspect
        width = 1. / self.cols - 2*wborder
        height = 1. - 2*hborder
        hb = 0.06
        wb = hb * aspect
        h = (height  - (self.rows - 1)*hb ) / self.rows
        w = h * aspect
        vertices = [(wborder + i*(w+wb), hborder + j*(h+hb)) for j in range(self.rows-1, -1, -1) for i in range(0, self.cols)][:self.max_nr]
        [self.hax.add_patch(plt.Rectangle(vertex,w,h, facecolor=facecolor, edgecolor=facecolor, transform=self.hax.transAxes)) for vertex in vertices]

        labels_all = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        labels = labels_all[:self.nr_axes]
        [self.hax.text(v[0]+w/2, v[1]+h/2, labels[i], color=fontcolor, va='center', ha='center', transform=self.hax.transAxes) for i,v in enumerate(vertices) ]

    def set_colorbar_layout(self, tickpos='right', direction='in', pad=5, labelpad=5, ticks=None, color='black', labelcolor='black', ticklabels=None):
        self.cb.ax.yaxis.set_ticks_position(tickpos)
        if ticks is not None: self.cb.set_ticks(ticks)
        if ticklabels is not None: self.cb.set_ticklabels(ticklabels)
        self.cb.ax.set_ylabel(self.colorbar_label, labelpad=labelpad, fontsize=self.fontsize, color=labelcolor)
        self.cb.ax.tick_params(direction=direction, pad=pad, colors=color, labelsize=self.fontsize)

    def optimize_colorbar_for_mask(self, ticks, labels, colors, pad=-15):
        self.cb.set_ticks(ticks)
        self.cb.ax.set_yticklabels(labels, rotation='vertical', fontsize=self.fontsize)
        self.cb.ax.tick_params(width=0, pad=pad)
        self.cb.ax.yaxis.get_ticklabels()[0].set_color(colors[0])
        self.cb.ax.yaxis.get_ticklabels()[1].set_color(colors[1])

    def cmap_discretize(self, cmap, N):
        """Return a discrete colormap from the continuous colormap cmap.
        
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
        Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
        """
        import matplotlib.colors
        cdict = cmap._segmentdata.copy()
        # N colors
        colors_i = np.linspace(0,1.,N)
        # N+1 indices
        indices = np.linspace(0,1.,N+1)
        for key in ('red','green','blue'):
            # Find the N colors
            D = np.array(cdict[key])
            I = interp1d(D[:,0], D[:,1])
            colors = I(colors_i)
            # Place these colors at the correct indices.
            A = np.zeros((N+1,3), float)
            A[:,0] = indices
            A[1:,1] = colors
            A[:-1,2] = colors
            # Create a tuple for the dictionary.
            L = []
            for l in A:
                L.append(tuple(l))
            cdict[key] = tuple(L)
            # Return colormap object.
        return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

    # PLOTTING FUNCTIONS #
    ######################
    def plotting_a_mask(self,mask):

        # check for single/multiple mask(s)
        if not isinstance(mask, list): mask = [mask]

        # use imshow to display mask(s)
        ims = [self.axes[i].imshow(mask[i], interpolation='none', cmap=cm.get_cmap('gray', 2)) for i in range(0, self.nr_axes)]

        # header
        if self.header:
            self.set_header_layout(title='This is a title', titlepos='right')
            self.add_overview_panel_to_header(facecolor='red')

        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(ims[0], cax = self.cax)
            self.set_colorbar_layout(tickpos='right', labelpad=13)
            self.optimize_colorbar_for_mask([0.25,0.75], ['False', 'True'], ['white', 'black'], pad=-13)

    def plotting_patterns(self, map, vmin=None, vmax=None, cmap='jet', type=None, bad='k', under=None, over=None, log=False, mask=None):

        # check for single/multiple map(s)
        if not isinstance(map, list): map = [map]

        # colormap
        cmap = cm.get_cmap(cmap)
        if under is not None: cmap.set_under(under,1.)
        if over is not None: cmap.set_over(over,1.)
            
        # mask out stuff
        if mask is not None:
            map = [np.ma.array(map[i], mask=~mask) for i in range(0, self.nr_axes)]
            cmap.set_bad(bad, 1.)

        # use imshow to display map(s)
        if log:
            ims = [self.axes[i].imshow(map[i], interpolation='bicubic', cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax, clip=False)) for i in range(0, self.nr_axes)]
        else:
            ims = [self.axes[i].imshow(map[i], interpolation='nearest', vmin=vmin, vmax=vmax, cmap=cmap) for i in range(0, self.nr_axes)]

        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(ims[0], cax = self.cax)
            self.set_colorbar_layout(tickpos='right')

        [self.set_axes_layout(self.axes[i], visible=self.axes_visible[i]) for i in range(self.nr_axes)]

        # put a title for single figure
        [self.axes[i].set_title(self.title_label[i], fontsize=self.fontsize) for i in range(self.nr_axes)]

        # x/y labels
        [self.axes[i].set_xlabel(self.xlabel[i], fontsize=self.fontsize) for i in range(self.nr_axes)]
        [self.axes[i].set_ylabel(self.ylabel[i], fontsize=self.fontsize) for i in range(self.nr_axes)]
        
    def plotting_a_map(self, axid, map, vmin=None, vmax=None, cmap='jet', type=None, bad='k', under=None, over=None, log=False, logx=False, logy=False, mask=None, extent=None, aspect=None, origin='upper', interpolation='nearest', colorbar=False, cax=None, colorbar_orientation='vertical', discrete_colors=None):
        # Colormap
        if cmap == 'discrete':
            cmap = discrete
        cmap = cm.get_cmap(cmap, discrete_colors)
        if discrete_colors is not None:
            cmap = self.cmap_discretize(cmap, discrete_colors)
        if under is not None: cmap.set_under(under,1.)
        if over is not None: cmap.set_over(over,1.)
            
        # mask out stuff
        if mask is not None:
            map = np.ma.array(map, mask=~mask)
            cmap.set_bad(bad, 1.)

        # use imshow to display map(s)
        if log:
            ims = self.axes[axid].imshow(map, interpolation=interpolation, cmap=cmap, extent=extent, aspect=aspect, origin=origin, norm=colors.LogNorm(vmin=vmin, vmax=vmax, clip=False))
        else:
            ims = self.axes[axid].imshow(map, interpolation=interpolation, vmin=vmin, vmax=vmax, cmap=cmap, extent=extent, aspect=aspect, origin=origin)

        # colorbar
        if self.colorbar or colorbar:
            if cax is None: cax = self.cax
            self.cb = self.fig.colorbar(ims, cax = cax, orientation=colorbar_orientation)
            self.set_colorbar_layout(tickpos='right')

        self.set_axes_layout(self.axes[axid], visible=self.axes_visible[axid])

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()


    def circles(self, axid, x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        if isinstance(c,basestring):
            color = c     # ie. use colors.colorConverter.to_rgba_array(c)
        else:
            color = None  # use cmap, norm after collection is created
        kwargs.update(color=color)
            
        if np.isscalar(x):
            patches = [Circle((x, y), s),]
        elif np.isscalar(s):
            patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
        else:
            patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
        collection = PatchCollection(patches, **kwargs)
                
        if color is None:
            collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)
                        
        self.axes[axid].add_collection(collection)
        self.axes[axid].autoscale_view()
        return collection

    def plotting_circles(self, circles, center, color='r', alpha=1., width=1., linestyle='dashed'):
        print self.nr_axes
        for i in range(self.nr_axes):
            x0, y0 = center[i]
            circle = circles[i]
            circle_patches = [patches.Circle((x0,y0), circle_patch, facecolor='none', edgecolor=color, linewidth=width, linestyle=linestyle,alpha=alpha) for circle_patch in circle]
            [self.axes[i].add_patch(circle_patch) for circle_patch in circle_patches]
        
    def plotting_a_coordinate_map(self, shape, mask=None, option=2):

        map = self.get_color_grid(shape, mask=mask, option=option)
        self.axes[0].imshow(map)

    def get_color_grid(self, shape, mask=None, flat=False, option=1):

        if mask is None: mask = np.zeros(shape).astype(np.bool)

        rmap = m.rgrid(shape)
        rmap /= rmap.max()

        tmap = m.tgrid(shape)
        tmap -= tmap.min()
        tmap /= tmap.max()

        if option == 1:
            H = tmap
            V = rmap
            S = np.ones_like(V)

        elif option == 2:
            H = rmap
            V = np.ones_like(H)
            S = np.ones_like(V)

        if flat:
            HSV = np.dstack((H.flatten(), S.flatten(), V.flatten()))
            RGB = c.hsv_to_rgb(HSV)[0]
            RGB[mask.flatten()] = np.array([1.,1.,1.])
        else:
            HSV = np.dstack((H,S,V))
            RGB = c.hsv_to_rgb(HSV)
            RGB[mask] = np.array([1.,1.,1.])

        return RGB


    def plotting_a_gaussian_kde(self, axid, data, xlim=(0,1), bins=100, bw=0.5, weights=None, color='k', label='', logx=False, logy=False):
        density = gaussian_kde(data, bw_method=True)
        density.set_bandwidth(bw)
        if isinstance(bins, np.ndarray):
            xs = bins
        else:
            xs = np.linspace(xlim[0],xlim[1],bins)
        self.axes[axid].plot(xs, density(xs), color=color, label=label)

        # logarithmic axis
        if logy: self.axes[axid].semilogy(nonposy='clip')
        if logx: self.axes[axid].semilogx()

        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location)
    
    def plotting_a_histogram(self, axid, hist, edges, type='step', cmap='jet', color=None, edgewidth=0, linestyle='-', label=None, xlim=None, ylim=None, logy=False, logx=False, zorder=0):

        if not isinstance(hist, list): hist = [hist]
        if not isinstance(edges, list): edges = [edges]
        nr_hist = len(hist)
        if label is None: label = range(nr_hist)
        if color is None: color = [cm.get_cmap(cmap, nr_hist)(i) for i in range(nr_hist)]

        # draw histograms in multiple axes
        #if self.nr_axes != 1:
        #   [self.axes[i].step(0.5*(edges[i][1:] + edges[i][:-1]), hist[i]) for i in range(0, self.nr_axes)]

        # histograms in one axes
        #else:
        if type == 'step': [self.axes[axid].step(0.5*(edges[i][1:] + edges[i][:-1]), hist[i], color=color[i], label=str(label[i]), lw=edgewidth, zorder=zorder) for i in range(0,nr_hist)]
        elif type == 'line': [self.axes[axid].plot(0.5*(edges[i][1:] + edges[i][:-1]), hist[i], color=color[i], label=str(label[i]), lw=edgewidth, ls=linestyle, zorder=zorder) for i in range(0,nr_hist)]
        elif type == 'bar':
            [self.axes[axid].bar(edges[i][:-1], hist[i], width=((edges[i][-1] - edges[i][0])/(edges[i].shape[0]-1)), edgecolor='w', lw=edgewidth, label=str(label[i]), color=color[i], zorder=zorder) for i in range(0,nr_hist)]

        # Legend
        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location, frameon=self.legend_frameon)

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # set x limits
        if xlim is None:
            self.axes[axid].set_xlim([edges[0][0], edges[0][-1]])
        else:
            self.axes[axid].set_xlim(xlim)
        self.axes[axid].set_ylim(ylim)
        
        # Set tick params
        self.axes[axid].tick_params(axis='x', which='both', bottom=1, top=0)
        self.axes[axid].tick_params(axis='y', which='both', right=0, left=1)

        # logarithmic axis
        if logy: self.axes[axid].semilogy(nonposy='clip')
        if logx: self.axes[axid].semilogx()
    
    def plotting_correlation(self, axid, x, y, color=None, cmap='jet', alpha=1., xlim=(None, None), ylim=(None, None), logx=False, logy=False, markersize=10, marker='o', markers=None, labels=None, zorder=0):

        if not isinstance(x, list): x = [x]
        if not isinstance(y, list): y = [y]
        nr = len(x)
        if markers is None:
            markers = [marker, marker, marker, '<', '>', '^', 'v'][:nr]
        if labels is None:
            labels = nr*['']
        if color is None:
            color = cm.get_cmap(cmap, nr)
            [self.axes[axid].scatter(x[i].flatten(),y[i].flatten(), color=color(i), marker=markers[i], s=markersize, alpha=alpha, edgecolors='None', label=labels[i], zorder=zorder) for i in range(0,nr)]
        else:
            [self.axes[axid].scatter(x[i].flatten(),y[i].flatten(), color=color[i], marker=markers[i], s=markersize, alpha=alpha, edgecolors='none', label=labels[i], zorder=zorder) for i in range(0,nr)]

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)
        
        # limits
        if xlim is not None: self.axes[axid].set_xlim([xlim[0], xlim[1]])
        if ylim is not None: self.axes[axid].set_ylim([ylim[0], ylim[1]])

        #Legend
        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location, scatterpoints=1, markerscale=2, frameon=self.legend_frameon)
        
        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()

    def plotting_correlation3d(self, axid, x, y, z, cmap='jet', alpha=1., ms=20, marker='o', vmin=None, vmax=None, xlim=[None, None], ylim=[None, None],scale_size=False, logx=False, logy=False, log=False, cticks=None, cticklabels=None):
        if scale_size: ms = 100 * (z - np.nanmin(z)) / np.nanmax(z - np.nanmin(z))
            
        if log:
            sc = self.axes[axid].scatter(x.flatten(),y.flatten(), c=z, marker=marker, s=ms,cmap=cmap, alpha=alpha, linewidths=1, edgecolors='none', norm=colors.LogNorm(vmin=vmin, vmax=vmax, clip=False))
        else:
            sc = self.axes[axid].scatter(x.flatten(),y.flatten(), c=z, marker=marker, s=ms,cmap=cmap, alpha=alpha, linewidths=1, edgecolors='none', vmin=vmin, vmax=vmax)
            
        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # limits
        #if xlim[0] is None: xlim[0] = x.min()
        #if xlim[1] is None: xlim[1] = x.max()
        #if ylim[0] is None: ylim[0] = y.min()
        #if ylim[1] is None: ylim[1] = y.max()
        self.axes[axid].set_xlim([xlim[0], xlim[1]])
        self.axes[axid].set_ylim([ylim[0], ylim[1]])
        
        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()
        
        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(sc, cax = self.cax)
            self.set_colorbar_layout(tickpos='right', ticks=cticks, ticklabels=cticklabels)

    def plotting_griddata(self, axid, x, y, z, xyrange, bins=[50,50], N=15, cmap='jet', alpha=1., ms=5, lw=0.5, ls='dashed', c='k', vmin=None, vmax=None, logx=False, logy=False, log=False, interpolation='linear', fill=False, bad='k', under=None, over=None):

        # Range and bins
        [[xmin, xmax],[ymin,ymax]] = xyrange
        [xbins, ybins] = bins

        # colormap
        cmap = cm.get_cmap(cmap)
        if under is not None: cmap.set_under(under,1.)
        if over is not None: cmap.set_over(over,1.)

        # Fill border
        if fill:
            xr = np.linspace(xmin, xmax, xbins)
            yr = np.linspace(ymin, ymax, ybins)
            x0 = np.ones(yr.shape)*xmin
            x1 = np.ones(yr.shape)*xmax
            y0 = np.ones(xr.shape)*ymin
            y1 = np.ones(xr.shape)*ymax
            zxb = np.zeros(xr.shape)
            zyb = np.zeros(yr.shape)
            xf = np.hstack([x, xr, x1, xr, x0])
            yf = np.hstack([y, y0, yr, y1, yr])
            zf = np.hstack([z, zxb, zyb, zxb, zyb])
        else:
            xf = x
            yf = y
            zf = z


        # Discretize data
        xi = np.linspace(xmin, xmax, xbins)
        yi = np.linspace(ymin, ymax, ybins)
        xgrid, ygrid = np.meshgrid(xi,yi)
        xy = np.vstack([xf,yf]).T
        zi = griddata(xy, zf, (xgrid, ygrid), method=interpolation)

        print zi.shape
        
        zi[zi<vmin] = vmin
        zi[zi>vmax] = vmax

        #if fill:
        #    zi = np.ma.array(zi, mask=(zi<0))
        #    cmap.set_bad(bad, 1.)

        if log:
            levs = np.logspace(np.log10(z.min()), np.log10(z.max()), 20)
            lev_exp = np.arange(np.floor(np.log10(z.min())-1), np.ceil(np.log10(z.max())+1))
            levs10 = np.power(10, lev_exp)
            levs = np.hstack([l*np.arange(1,10,1) for l in levs10])
        #else:
        #    levs = np.linspace(vmin, vmax, 10)
            
        # Contour plot
        #self.axes[axid].contour(xi, yi, zi, N, linewidths=lw, linestyles=ls, colors='k')
        if log:
            cs = self.axes[axid].contourf(xi, yi, zi, levs, cmap=cmap, norm=colors.LogNorm())
        else:
            cs = self.axes[axid].contourf(xi, yi, zi, N, vmin=vmin, vmax=vmax, cmap=cmap)

        # plot data points
        self.axes[axid].scatter(x, y, marker='o', c=c, s=ms, alpha=alpha, linewidths=0, zorder=10)

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()
        
        # limits
        self.axes[axid].set_xlim([xmin, xmax])
        self.axes[axid].set_ylim([ymin, ymax])

        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(cs, cax = self.cax)
            self.set_colorbar_layout(tickpos='left')
            if log: self.cb.set_ticks([levs10])


    def plotting_contourf(self, axid, x, y, z, N=15, cmap='jet', alpha=1., vmin=None, vmax=None, logx=False, logy=False, log=False):

        if vmin is None: vmin = z.min()
        if vmax is None: vmax = z.max()
        z[z>vmax] = vmax
        
        if log:
            lev_exp = np.arange(np.floor(np.log10(vmin)-1), np.ceil(np.log10(vmax)+1))
            levs10 = np.power(10, lev_exp)
            print levs10
            levs = np.hstack([l*np.arange(1,10,1) for l in levs10])
        else:
            levs = np.linspace(vmin, vmax, N)
            
        # Contour plot
        #self.axes[axid].contour(xi, yi, zi, N, linewidths=lw, linestyles=ls, colors='k')
        if log:
            cs = self.axes[axid].contourf(x, y, z, levs, cmap=cmap, norm=colors.LogNorm())
        else:
            cs = self.axes[axid].contourf(x, y, z, levs, cmap=cmap)

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()
        
        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(cs, cax = self.cax)
            self.set_colorbar_layout(tickpos='right')
            if log: self.cb.set_ticks([levs10])

    def plotting_a_contour(self, axid, x, y, z, levels, colors='k', linewidths=1, linestyles='-', logx=False, logy=False, label=False, format='%.2f', manual=None):
        cs = self.axes[axid].contour(x, y, z, levels, colors = colors, linewidths=linewidths, linestyles=linestyles,hold='on')

        if label:
            self.axes[axid].clabel(cs, inline=1, fontsize=self.fontsize, fmt=format, manual=manual)
        
        # logarithmic axis
        if logx: self.axes[axid].semilogx()
        if logy: self.axes[axid].semilogy()
        
    def plotting_correlation_color(self, x, y, markersize=10, option=1):

        [self.axes[i].scatter(x[i].flatten(),y[i].flatten(), s=markersize, c=self.get_color_grid(x[i].shape, flat=True, option=option), marker='x') for i in range(0, self.nr_axes)]

    def plotting_a_heatmap(self, axid, x, y, hrange, bins=100, cmaplist=None, alpha=1., vmin=1, vmax=None, log=False, under='k', bad='k', colorbar=False, cax=None, colorbar_orientation='vertical', zorder=0):

        if not isinstance(x, list): x = [x]
        if not isinstance(y, list): y = [y]
        nr = len(x)
        if cmaplist is None: cmaplist = ['Blues', 'Reds', 'Greens', 'Greys']
        
        # colormap
        cmaplist = [cm.get_cmap(cmap) for cmap in cmaplist]
        for cmap in cmaplist: cmap.set_under(under, 1.)
        for cmap in cmaplist: cmap.set_bad(bad, 1.)
        
        H = [np.histogram2d(x[i].flatten(), y[i].flatten(), range=hrange, bins=bins) for i in range(0,nr)]
        img = H[i][0].T
        if log:
            ims = [self.axes[axid].imshow(img, extent=[H[i][1][0], H[i][1][-1], H[i][2][0], H[i][2][-1]],interpolation = 'nearest', aspect='auto', cmap=cmaplist[i], alpha = alpha, norm=colors.LogNorm(vmin=vmin, vmax=vmax), origin = 'lower', zorder=zorder) for i in range(0,nr)]
        else:
            ims = [self.axes[axid].imshow(img, extent=[H[i][1][0], H[i][1][-1], H[i][2][0], H[i][2][-1]],interpolation = 'nearest', aspect='auto', cmap=cmaplist[i], alpha = alpha, vmin=vmin, vmax=vmax,origin = 'lower', zorder=zorder) for i in range(0,nr)]
            
        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # Colorbar axis
        if cax is not None:
            self.cax = cax
        
        # colorbar
        if self.colorbar or colorbar:
            self.cb = self.fig.colorbar(ims[0], cax = self.cax, orientation=colorbar_orientation)
            self.set_colorbar_layout(tickpos='right')

    def plotting_a_meanmap(self, axid, x, y, z, hrange, bins=100, cmaplist=None, alpha=1., vmin=0.0001, vmax=None, log=False):
        
        if not isinstance(x, list): x = [x]
        if not isinstance(y, list): y = [y]
        if not isinstance(z, list): z = [z]
        nr = len(x)
        if cmaplist is None: cmaplist = ['Greens', 'Reds', 'Blues', 'Greys']
        
        # colormap
        cmaplist = [cm.get_cmap(cmap) for cmap in cmaplist]
        for cmap in cmaplist: cmap.set_under('w', 0.)

        S = [np.histogram2d(x[i].flatten(), y[i].flatten(), range=hrange, weights=z[i], bins=bins) for i in range(0,nr)]
        H = [np.histogram2d(x[i].flatten(), y[i].flatten(), range=hrange, bins=bins) for i in range(0,nr)]
        for i in range(0,nr): H[i][0][S[i][0]==0] = 1. 

        if log:
            ims = [self.axes[axid].imshow(S[i][0].T/H[i][0].T, extent=[H[i][1][0], H[i][1][-1], H[i][2][0], H[i][2][-1]],interpolation = 'nearest', aspect='auto', cmap=cmaplist[i], alpha = alpha, norm=colors.LogNorm(vmin=vmin, vmax=vmax, clip=False), origin = 'lower') for i in range(0,nr)]
        else:
            ims = [self.axes[axid].imshow(S[i][0].T/H[i][0].T, extent=[H[i][1][0], H[i][1][-1], H[i][2][0], H[i][2][-1]],interpolation = 'nearest', aspect='auto', cmap=cmaplist[i], alpha = alpha, vmin=vmin, vmax=vmax, origin = 'lower') for i in range(0,nr)]

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # colorbar
        if self.colorbar:
            self.cb = self.fig.colorbar(ims[0], cax = self.cax)
            self.set_colorbar_layout(tickpos='right')
            
            
    def plotting_traces(self, axid, x,y, labels, linestyles=None, linewidths=None, xlim=None, ylim=None, markersize=5, markers=None, mecs=None, mfcs=None, colors=None, alphas=None, logx=False, logy=False, zorder=0):
        if not isinstance(x,list): x = list(x)
        if not isinstance(y,list): y = list(y)
        nr = len(x)
        if colors is None: colors=['b', 'r', 'g', 'k', 'm', 'y'][:nr]
        if markers is None: markers = nr*[None]
        if mecs is None: mecs=nr*['k']
        if mfcs is None: mfcs=nr*['k']
        if linestyles is None: linestyles = nr*['-']
        if linewidths is None: linewidths = nr*[1]
        if alphas is None: alphas = nr*[1]
        [self.axes[axid].plot(x[i],y[i], label=labels[i], ls=linestyles[i], lw=linewidths[i], color=colors[i], marker=markers[i], markersize=markersize, mec=mecs[i], mfc=mfcs[i], alpha=alphas[i], zorder=zorder) for i in range(nr)]

        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx(nonposx='clip')
        if logy: self.axes[axid].semilogy(nonposy='clip')
        
        # legend
        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location, numpoints=1, markerscale=2, frameon=self.legend_frameon)

        # limits
        if xlim is not None: self.axes[axid].set_xlim([xlim[0], xlim[1]])
        if ylim is not None: self.axes[axid].set_ylim([ylim[0], ylim[1]])

    def plotting_points_with_error(self, axid, x, y, l, xerr=None, yerr=None, marker='o', markersize=5, color='b', ecolor='k', xlim=None, ylim=None, logx=False, logy=False):

        # Plotting points with errorbars
        self.axes[axid].errorbar(x, y, xerr=xerr, yerr=yerr, fmt=marker, color=color, ecolor=ecolor, label=l, markersize=markersize)
        
        # put a title for single figure
        self.axes[axid].set_title(self.title_label[axid], fontsize=self.fontsize)

        # x/y labels
        self.axes[axid].set_xlabel(self.xlabel[axid], fontsize=self.fontsize)
        self.axes[axid].set_ylabel(self.ylabel[axid], fontsize=self.fontsize)

        # logarithmic axis
        if logx: self.axes[axid].semilogx(nonposx='clip')
        if logy: self.axes[axid].semilogy(nonposy='clip')
        
        # legend
        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location, numpoints=1, markerscale=2, frameon=self.legend_frameon)

        # limits
        if xlim is not None: self.axes[axid].set_xlim([xlim[0], xlim[1]])
        if ylim is not None: self.axes[axid].set_ylim([ylim[0], ylim[1]])

        
    def plotting_a_circle(self, axid, r, x,y, color='w', alpha=1., width=1., linestyle='dashed', label=''):
        circle = patches.Circle((x,y), r, facecolor=color, fill=False, edgecolor=color, alpha=alpha, linewidth=width, linestyle=linestyle, label=label) 
        self.axes[axid].add_patch(circle)
            
        # legend (default)
        if self.legend: self.axes[axid].legend(prop={'size':self.fontsize-1}, loc=self.legend_location, numpoints=1, markerscale=2, frameon=self.legend_frameon)

        
    def plotting_a_rectangle(self, axid, x0, y0, w, h, alpha=1., linestyle='-', edgecolor='r', facecolor='none', linewidth=1):
        self.axes[axid].add_patch(plt.Rectangle((x0,y0),w,h, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, alpha=alpha))

    # I/O FUNCTIONS #
    #################
    def show(self): plt.show()
    def save(self, filename, facecolor='w'):
        if self.save_pdf:
            self.fig.savefig(filename, dpi=300, format='pdf', facecolor=facecolor, bbox_inches='tight')
        elif self.save_png:
            self.fig.savefig(filename, dpi=300, format='png', facecolor=facecolor, bbox_inches='tight')
        elif self.save_ps:
            self.fig.savefig(filename, dpi=300, format='eps', facecolor=facecolor, bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close(self.fig)
