# coding: utf-8

from __future__ import division
from __future__ import print_function

from __future__ import division

import numpy as np
import sys

import matplotlib.widgets as mwidgets
from matplotlib import lines
import pylab as pl


class RectangleSelector(mwidgets.RectangleSelector):
    """Widget widget for selecting a rectangular region in a plot.

    Unlike :class:`matplotlib.widgets.RectangleSelector`, this widget remains
    visible after selection and can be resized using corner and edge handles.

    After making the desired selection, press "Enter" to accept the selection
    and call the `onselect` callback function.

    Note the difference from :class:`matplotlib.widgets.RectangleSelector`:
    The callback method `onselect` is called on "Enter", *not* after release of
    the mouse button.  In addition, the `onselect` function now takes a single
    argument `extents`, which is a tuple specifying (xmin, xmax, ymin, ymax) of
    the rectangle.

    TODO: Allow rectangle to flip over itself when resizing.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes

    onselect : function
        Function accepting rectangle extents as the only argument. If None,
        print extents of rectangle.

    mindist : float
        Minimum distance in pixels for selection of a control (i.e. corner or
        edge) handle.

    rectprops : dict
        Properties for :class:`matplotlib.patches.Rectangle`. This class
        redefines defaults in :class:`matplotlib.widgets.RectangleSelector`.

    kwargs : see :class:`matplotlib.widgets.RectangleSelector`.

    Attributes
    ----------
    extents : tuple
        Rectangle extents: (xmin, xmax, ymin, ymax).
    """

    def __init__(self, ax, mouserelease, extent_init, onselect=None, rectprops=None,
                 mindist=10, **kwargs):

        if 'drawtype' in kwargs and not kwargs['drawtype'] == 'box':
            raise ValueError('"drawtype" must be "box"')

        self.mindist = mindist
        self.active_handle = None
        self._on_mouse_release = mouserelease

        rectprops_defaults = dict(edgecolor='k', facecolor='r', alpha=0.2)
        
        if rectprops is None:
            rectprops = {}
            
        rectprops.update(rectprops_defaults)

        mwidgets.RectangleSelector.__init__(self, ax,
                                            self._on_mouse_release,
                                            rectprops=rectprops,
                                            **kwargs)
        # Alias rectangle attribute.
        self._rect = self.to_draw

        if onselect is None:
            def onselect(extents):
                print("(xmin=%.3g, xmax=%.3g, ymin=%.3g, ymax=%.3g)" % extents)
                
        self.onenter = onselect

        handle_props = dict(mfc='none', mec='k', ls='none', alpha=0.7,
                            visible=False)

        self._corner_order = ['NW', 'NE', 'SE', 'SW']
        xc, yc = self.corner_coords
        self._corner_handles = lines.Line2D(xc, yc, marker='o', **handle_props)
        # replace with widget method for clean up
        self.ax.add_line(self._corner_handles)

        self._edge_order = ['W', 'N', 'E', 'S']
        xe, ye = self.edge_coords
        self._edge_handles = lines.Line2D(xe, ye, marker='s', **handle_props)
        self.ax.add_line(self._edge_handles)

        self.connect_event('key_press_event', self.onkeypress)
        [self.x0, self.x1, self.y0, self.y1] = extent_init

    @property
    def _rect_bbox(self):
        x0 = self._rect.get_x()
        y0 = self._rect.get_y()
        width = self._rect.get_width()
        height = self._rect.get_height()
        return x0, y0, width, height

    @property
    def corner_coords(self):
        """Corners of rectangle from lower left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        xc = x0, x0 + width, x0 + width, x0
        yc = y0, y0, y0 + height, y0 + height
        return xc, yc

    @property
    def edge_coords(self):
        """Midpoint of rectangle edges from lower left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        w = width / 2.
        h = height / 2.
        xe = x0, x0 + w, x0 + width, x0 + w
        ye = y0 + h, y0, y0 + h, y0 + height
        return xe, ye

    @property
    def extents(self):
        xmin = min(self.x0, self.x1)
        xmax = max(self.x0, self.x1)
        ymin = min(self.y0, self.y1)
        ymax = max(self.y0, self.y1)
        return xmin, xmax, ymin, ymax


    def release(self, event):
        mwidgets.RectangleSelector.release(self, event)
        # Undo hiding of rectangle and redraw.
        self.set_visible(True)
        self.update()
        self.set_animated(False)

    def press(self, event):
        dist = []
        for h in (self._corner_handles, self._edge_handles):
            pts = np.transpose((h.get_xdata(), h.get_ydata()))
            pts = self.ax.transData.transform(pts)
            diff = pts - ((event.x, event.y))
            dist.append(np.sqrt(np.sum(diff**2, axis=1)))

        dist = np.asarray(dist)
        idx = np.argmin(dist)

        close_to_handle = dist.flat[idx] < self.mindist
        if idx < 4:
            handle = 'corner'
        else:
            handle = 'edge'
            idx -= 4

        if close_to_handle and handle == 'corner':
            self.active_handle = self._corner_order[idx]
        elif close_to_handle and handle == 'edge':
            self.active_handle = self._edge_order[idx]
        else:
            self.active_handle = None

        # Clear previous rectangle before drawing new rectangle.
        self.set_animated(True)
        if not close_to_handle:
            self.set_visible(False)
            self.update()
            self.set_visible(True)

        mwidgets.RectangleSelector.press(self, event)

    def onkeypress(self, event):
        if event.key == 'enter':
            self.onenter(self.extents)
            self.set_visible(False) 
            self.update()

    def onmove(self, event):

        if self.eventpress is None or self.ignore(event):
            return

        if self.active_handle is None:
            xmin, ymin = event.xdata, event.ydata
            xmax = self.eventpress.xdata
            ymax = self.eventpress.ydata
        else:
            xmin, ymin, width, height = self._rect_bbox
            xmax = xmin + width
            ymax = ymin + height

            if self.active_handle in ('W', 'SW', 'NW'):
                xmin = event.xdata
            if self.active_handle in ('E', 'SE', 'NE'):
                xmax = event.xdata
            if self.active_handle in ('N', 'NW', 'NE'):
                ymin = event.ydata
            if self.active_handle in ('S', 'SW', 'SE'):
                ymax = event.ydata

        # Order by value instead of time.
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        self._rect.set_x(xmin)
        self._rect.set_y(ymin)
        self._rect.set_width(xmax-xmin)
        self._rect.set_height(ymax-ymin)

        xc, yc = self.corner_coords
        self._corner_handles.set_xdata(xc)
        self._corner_handles.set_ydata(yc)
        xe, ye = self.edge_coords
        self._edge_handles.set_xdata(xe)
        self._edge_handles.set_ydata(ye)
        self.update()
        return False

    def set_visible(self, val):
        self._rect.set_visible(val)
        self._corner_handles.set_visible(val)
        self._edge_handles.set_visible(val)

    def set_animated(self, val):
        self._rect.set_animated(val)
        self._corner_handles.set_animated(val)
        self._edge_handles.set_animated(val)

    def update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self._rect)
            self.ax.draw_artist(self._edge_handles)
            self.ax.draw_artist(self._corner_handles)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()
        return False



def smooth(x,window_len=11,window='hanning'):
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]
    
    
def movingSmooth(array, axis, windowLen, windowType):
    
    [m, n] = np.shape(array)
    smoothedArray = np.zeros_like(array)
    
    if axis == 0:
        for idx in range(n):       
            smoothedArray[:, idx] = smooth(array[:, idx], windowLen, windowType)
    elif axis == 1:
        for idx in range(m):       
            smoothedArray[idx, :] = smooth(array[idx, :], windowLen, windowType)
    else:
        raise ValueError, "You didn't enter a valid value for axis number (0, 1)"
        
    return smoothedArray
        
    


def calcLMS(vec):
	
    L = pl.ceil(vec.size/2) - 1
    alpha = 1.
    M = pl.rand(L,vec.size) + alpha
    	
    for k in range(1,int(L)):
    	for i in range(int(k+2), int(vec.size-k+1)):			
    		if (vec[i-1] > vec[i-k-1]) and (vec[i-1]>vec[i+k-1]):
    			M[k,i] = 0.
     
    return M
	 
def getPeaks(x,y):
	
    linfit = pl.polyfit(x, y, 1)
    yfit = pl.polyval(linfit, x)
    
    ynew = y - yfit
    
    M = calcLMS(ynew)

    gamma_k = M.sum(1)
    
    lambda_m = gamma_k.argmin()
    	# rescale the LMS

    Mr = M[1:lambda_m, :]
    
    sigma = 1./(lambda_m-1) * pl.sum( pl.sqrt(( Mr - 1./lambda_m* (pl.ones((lambda_m-1,1))* Mr.sum(0)) )**2), 0)
    
    p = pl.find(sigma == 0) - 1
    	
    return p
    


