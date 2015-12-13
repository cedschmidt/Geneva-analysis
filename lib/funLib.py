# coding: utf-8

from __future__ import division
from __future__ import print_function

import numpy as np
from matplotlib.widgets import _SelectorWidget, ToolHandles
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.sparse import csc_matrix, spdiags
from scipy.sparse.linalg import spsolve





############################################################ Rectangle selector


class RectangleSelector(_SelectorWidget):
    """
    Select a rectangular region of an axes.

    For the cursor to remain responsive you much keep a reference to
    it.

    Example usage::

        from matplotlib.widgets import  RectangleSelector
        from pylab import *

        def onselect(eclick, erelease):
          'eclick and erelease are matplotlib events at press and release'
          print(' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata))
          print(' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata))
          print(' used button   : ', eclick.button)

        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print(' RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print(' RectangleSelector activated.')
                toggle_selector.RS.set_active(True)

        x = arange(100)/(99.0)
        y = sin(x)
        fig = figure
        ax = subplot(111)
        ax.plot(x,y)

        toggle_selector.RS = RectangleSelector(ax, onselect, drawtype='line')
        connect('key_press_event', toggle_selector)
        show()
    """

    _shape_klass = Rectangle
    
    def __init__(self, ax, onselect, drawtype='box',
                     minspanx=None, minspany=None, useblit=False,
                     lineprops=None, rectprops=None, spancoords='data',
                     button=None, maxdist=10, marker_props=None,
                     interactive=False, state_modifier_keys=None):
    
            """
            Create a selector in *ax*.  When a selection is made, clear
            the span and call onselect with::
    
              onselect(pos_1, pos_2)
    
            and clear the drawn box/line. The ``pos_1`` and ``pos_2`` are
            arrays of length 2 containing the x- and y-coordinate.
    
            If *minspanx* is not *None* then events smaller than *minspanx*
            in x direction are ignored (it's the same for y).
    
            The rectangle is drawn with *rectprops*; default::
    
              rectprops = dict(facecolor='red', edgecolor = 'black',
                               alpha=0.2, fill=True)
    
            The line is drawn with *lineprops*; default::
    
              lineprops = dict(color='black', linestyle='-',
                               linewidth = 2, alpha=0.5)
    
            Use *drawtype* if you want the mouse to draw a line,
            a box or nothing between click and actual position by setting
    
            ``drawtype = 'line'``, ``drawtype='box'`` or ``drawtype = 'none'``.
    
            *spancoords* is one of 'data' or 'pixels'.  If 'data', *minspanx*
            and *minspanx* will be interpreted in the same coordinates as
            the x and y axis. If 'pixels', they are in pixels.
    
            *button* is a list of integers indicating which mouse buttons should
            be used for rectangle selection.  You can also specify a single
            integer if only a single button is desired.  Default is *None*,
            which does not limit which button can be used.
    
            Note, typically:
             1 = left mouse button
             2 = center mouse button (scroll wheel)
             3 = right mouse button
    
            *interactive* will draw a set of handles and allow you interact
            with the widget after it is drawn.
    
            *state_modifier_keys* are keyboard modifiers that affect the behavior
            of the widget.
    
            The defaults are:
            dict(move=' ', clear='escape', square='shift', center='ctrl')
    
            Keyboard modifiers, which:
            'move': Move the existing shape.
            'clear': Clear the current shape.
            'square': Makes the shape square.
            'center': Make the initial point the center of the shape.
            'square' and 'center' can be combined.
            """
            _SelectorWidget.__init__(self, ax, onselect, useblit=useblit,
                                     button=button,
                                     state_modifier_keys=state_modifier_keys)
    
            self.to_draw = None
            self.visible = True
            self.interactive = interactive
    
            if drawtype == 'none':
                drawtype = 'line'                        # draw a line but make it
                self.visible = False                     # invisible
    
            if drawtype == 'box':
                if rectprops is None:
                    rectprops = dict(facecolor='red', edgecolor='black',
                                     alpha=0.2, fill=True)
                rectprops['animated'] = self.useblit
                self.rectprops = rectprops
                self.to_draw = self._shape_klass((0, 0),
                                         0, 1, visible=False, **self.rectprops)
                self.ax.add_patch(self.to_draw)
            if drawtype == 'line':
                if lineprops is None:
                    lineprops = dict(color='black', linestyle='-',
                                     linewidth=2, alpha=0.5)
                lineprops['animated'] = self.useblit
                self.lineprops = lineprops
                self.to_draw = Line2D([0, 0], [0, 0], visible=False,
                                      **self.lineprops)
                self.ax.add_line(self.to_draw)
    
            self.minspanx = minspanx
            self.minspany = minspany
    
            if spancoords not in ('data', 'pixels'):
                msg = "'spancoords' must be one of [ 'data' | 'pixels' ]"
                raise ValueError(msg)
    
            self.spancoords = spancoords
            self.drawtype = drawtype
    
            self.maxdist = maxdist
    
            if rectprops is None:
                props = dict(mec='r')
            else:
                props = dict(mec=rectprops.get('edgecolor', 'r'))
            self._corner_order = ['NW', 'NE', 'SE', 'SW']
            xc, yc = self.corners
            self._corner_handles = ToolHandles(self.ax, xc, yc, marker_props=props,
                                               useblit=self.useblit)
    
            self._edge_order = ['W', 'N', 'E', 'S']
            xe, ye = self.edge_centers
            self._edge_handles = ToolHandles(self.ax, xe, ye, marker='s',
                                             marker_props=props, useblit=self.useblit)
    
            xc, yc = self.center
            self._center_handle = ToolHandles(self.ax, [xc], [yc], marker='s',
                                              marker_props=props, useblit=self.useblit)
    
            self.active_handle = None
    
            self.artists = [self.to_draw, self._center_handle.artist,
                            self._corner_handles.artist,
                            self._edge_handles.artist]
    
            if not self.interactive:
                self.artists = [self.to_draw]
    
            self._extents_on_press = None
    
    def _press(self, event):
        """on button press event"""
        # make the drawed box/line visible get the click-coordinates,
        # button, ...
        if self.interactive and self.to_draw.get_visible():
            self._set_active_handle(event)
        else:
            self.active_handle = None
        
        if self.active_handle is None or not self.interactive:
            # Clear previous rectangle before drawing new rectangle.
            self.update()
        
        self.set_visible(self.visible)
    
    def _release(self, event):
            """on button release event"""
            if not self.interactive:
                self.to_draw.set_visible(True)
    
            if self.spancoords == 'data':
                xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
                xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
                # calculate dimensions of box or line get values in the right
                # order
            elif self.spancoords == 'pixels':
                xmin, ymin = self.eventpress.x, self.eventpress.y
                xmax, ymax = self.eventrelease.x, self.eventrelease.y
            else:
                raise ValueError('spancoords must be "data" or "pixels"')
    
            if xmin > xmax:
                xmin, xmax = xmax, xmin
            if ymin > ymax:
                ymin, ymax = ymax, ymin
    
            spanx = xmax - xmin
            spany = ymax - ymin
            xproblems = self.minspanx is not None and spanx < self.minspanx
            yproblems = self.minspany is not None and spany < self.minspany
    
            if (((self.drawtype == 'box') or (self.drawtype == 'line')) and
                    (xproblems or yproblems)):
                # check if drawn distance (if it exists) is not too small in
                # neither x nor y-direction
                self.extents = [0, 0, 0, 0]
                return
    
            # update the eventpress and eventrelease with the resulting extents
            x1, x2, y1, y2 = self.extents
            self.eventpress.xdata = x1
            self.eventpress.ydata = y1
            xy1 = self.ax.transData.transform_point([x1, y1])
            self.eventpress.x, self.eventpress.y = xy1
    
            self.eventrelease.xdata = x2
            self.eventrelease.ydata = y2
            xy2 = self.ax.transData.transform_point([x2, y2])
            self.eventrelease.x, self.eventrelease.y = xy2
    
            self.onselect(self.eventpress, self.eventrelease)
                                                  # call desired function
            self.update()
    
            return False
    
    def _onmove(self, event):
        """on motion notify event if box/line is wanted"""
        # resize an existing shape
        if self.active_handle and not self.active_handle == 'C':
            x1, x2, y1, y2 = self._extents_on_press
            if self.active_handle in ['E', 'W'] + self._corner_order:
                x2 = event.xdata
            if self.active_handle in ['N', 'S'] + self._corner_order:
                y2 = event.ydata

        # move existing shape
        elif (('move' in self.state or self.active_handle == 'C')
              and self._extents_on_press is not None):
            x1, x2, y1, y2 = self._extents_on_press
            dx = event.xdata - self.eventpress.xdata
            dy = event.ydata - self.eventpress.ydata
            x1 += dx
            x2 += dx
            y1 += dy
            y2 += dy

        # new shape
        else:
            center = [self.eventpress.xdata, self.eventpress.ydata]
            center_pix = [self.eventpress.x, self.eventpress.y]
            dx = (event.xdata - center[0]) / 2.
            dy = (event.ydata - center[1]) / 2.

            # square shape
            if 'square' in self.state:
                dx_pix = abs(event.x - center_pix[0])
                dy_pix = abs(event.y - center_pix[1])
                if not dx_pix:
                    return
                maxd = max(abs(dx_pix), abs(dy_pix))
                if abs(dx_pix) < maxd:
                    dx *= maxd / (abs(dx_pix) + 1e-6)
                if abs(dy_pix) < maxd:
                    dy *= maxd / (abs(dy_pix) + 1e-6)

            # from center
            if 'center' in self.state:
                dx *= 2
                dy *= 2

            # from corner
            else:
                center[0] += dx
                center[1] += dy

            x1, x2, y1, y2 = (center[0] - dx, center[0] + dx,
                              center[1] - dy, center[1] + dy)

        self.extents = x1, x2, y1, y2

    @property
    def _rect_bbox(self):
        if self.drawtype == 'box':
            x0 = self.to_draw.get_x()
            y0 = self.to_draw.get_y()
            width = self.to_draw.get_width()
            height = self.to_draw.get_height()
            return x0, y0, width, height
        else:
            x, y = self.to_draw.get_data()
            x0, x1 = min(x), max(x)
            y0, y1 = min(y), max(y)
            return x0, y0, x1 - x0, y1 - y0

    @property
    def corners(self):
        """Corners of rectangle from lower left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        xc = x0, x0 + width, x0 + width, x0
        yc = y0, y0, y0 + height, y0 + height
        return xc, yc

    @property
    def edge_centers(self):
        """Midpoint of rectangle edges from left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        w = width / 2.
        h = height / 2.
        xe = x0, x0 + w, x0 + width, x0 + w
        ye = y0 + h, y0, y0 + h, y0 + height
        return xe, ye

    @property
    def center(self):
        """Center of rectangle"""
        x0, y0, width, height = self._rect_bbox
        return x0 + width / 2., y0 + height / 2.

    @property
    def extents(self):
        """Return (xmin, xmax, ymin, ymax)."""
        x0, y0, width, height = self._rect_bbox
        xmin, xmax = sorted([x0, x0 + width])
        ymin, ymax = sorted([y0, y0 + height])
        return xmin, xmax, ymin, ymax

    @extents.setter
    def extents(self, extents):
        # Update displayed shape
        self.draw_shape(extents)
        # Update displayed handles
        self._corner_handles.set_data(*self.corners)
        self._edge_handles.set_data(*self.edge_centers)
        self._center_handle.set_data(*self.center)
        self.set_visible(self.visible)
        self.update()

    def draw_shape(self, extents):
        x0, x1, y0, y1 = extents
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])
        xlim = sorted(self.ax.get_xlim())
        ylim = sorted(self.ax.get_ylim())

        xmin = max(xlim[0], xmin)
        ymin = max(ylim[0], ymin)
        xmax = min(xmax, xlim[1])
        ymax = min(ymax, ylim[1])

        if self.drawtype == 'box':
            self.to_draw.set_x(xmin)
            self.to_draw.set_y(ymin)
            self.to_draw.set_width(xmax - xmin)
            self.to_draw.set_height(ymax - ymin)

        elif self.drawtype == 'line':
            self.to_draw.set_data([xmin, xmax], [ymin, ymax])

    def _set_active_handle(self, event):
        """Set active handle based on the location of the mouse event"""
        # Note: event.xdata/ydata in data coordinates, event.x/y in pixels
        c_idx, c_dist = self._corner_handles.closest(event.x, event.y)
        e_idx, e_dist = self._edge_handles.closest(event.x, event.y)
        m_idx, m_dist = self._center_handle.closest(event.x, event.y)

        if 'move' in self.state:
            self.active_handle = 'C'
            self._extents_on_press = self.extents

        # Set active handle as closest handle, if mouse click is close enough.
        elif m_dist < self.maxdist * 2:
            self.active_handle = 'C'
        elif c_dist > self.maxdist and e_dist > self.maxdist:
            self.active_handle = None
            return
        elif c_dist < e_dist:
            self.active_handle = self._corner_order[c_idx]
        else:
            self.active_handle = self._edge_order[e_idx]

        # Save coordinates of rectangle at the start of handle movement.
        x1, x2, y1, y2 = self.extents
        # Switch variables so that only x2 and/or y2 are updated on move.
        if self.active_handle in ['W', 'SW', 'NW']:
            x1, x2 = x2, event.xdata
        if self.active_handle in ['N', 'NW', 'NE']:
            y1, y2 = y2, event.ydata
        self._extents_on_press = x1, x2, y1, y2

    @property
    def geometry(self):
        if hasattr(self.to_draw, 'get_verts'):
            xfm = self.ax.transData.inverted()
            y, x = xfm.transform(self.to_draw.get_verts()).T
            return np.array([x[:-1], y[:-1]])
        else:
            return np.array(self.to_draw.get_data())
    
###############################################################################       
 









        
        
####################################################    Draggable vertical line
        
class DraggableVLine:
    def __init__(self, onmovefun, ax, canvas, posinit, colorV = "black", lineVWidth = 2):
        self.canvas = canvas
        self.press = None
        self.ax = ax
        self.pos = posinit
        self.color = colorV
        self.linewidth = lineVWidth
        self.fun = onmovefun
        self.background = None
        self.Vline = self.ax.axvline(posinit, color = colorV, linewidth = lineVWidth)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        contains, attrd = self.Vline.contains(event)
        
        if not contains: 
            return
            
        self.pos = self.Vline.get_xdata()
        self.press = self.pos, event.xdata
        
        self.Vline.set_animated(True)
        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.Vline)
        self.canvas.blit(self.ax.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: 
            return
            
        if event.inaxes != self.Vline.axes: 
            return
            
        self.pos = event.xdata 
        self.Vline.set_xdata(self.pos)
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.Vline)
        self.canvas.blit(self.ax.bbox)



    def on_release(self, event):
        'on release we reset the press data'
        self.press = None    
        self.background = None
        self.Vline.set_animated(False)
        self.fun(self.pos)
        
        
    def update(self, position):
        self.pos = position
        self.Vline = self.ax.axvline(self.pos, color = self.color, linewidth = self.linewidth)  

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.disconnect_events()

###############################################################################
        
   




     
        
####################################################  Draggable horizontal line
        
class DraggableHLine:
    def __init__(self, onmovefun, ax, canvas, posinit, colorH = "black", lineVWidth = 2):
        self.canvas = canvas
        self.press = None
        self.ax = ax
        self.pos = posinit
        self.color = colorH
        self.linewidth = lineVWidth
        self.fun = onmovefun
        self.background = None
        self.Hline = self.ax.axhline(posinit, color = colorH, linewidth = lineVWidth)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        contains, attrd = self.Hline.contains(event)
        
        if not contains: 
            return
            
        self.pos = self.Hline.get_xdata()
        self.press = self.pos, event.xdata
        self.Hline.set_animated(True)
        self.canvas.draw()
        
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)

        # now redraw just the rectangle
        self.ax.draw_artist(self.Hline)

        # and blit just the redrawn area
        self.canvas.blit(self.ax.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: 
            return
            
        if event.inaxes != self.Hline.axes: 
            return
            
        #self.pos, xpress = self.press
        self.pos = event.ydata 
        
        
        self.Hline.set_ydata(self.pos)
        
        self.canvas.restore_region(self.background)

        # redraw just the current rectangle
        self.ax.draw_artist(self.Hline)

        # blit just the redrawn area
        self.canvas.blit(self.ax.bbox)


    def on_release(self, event):
        'on release we reset the press data'
        self.press = None    
        self.background = None
        self.Hline.set_animated(False)
        self.fun(self.pos)
        
    def update(self, position):
        self.pos = position
        self.Hline = self.ax.axhline(self.pos, color = self.color, linewidth = self.linewidth)  

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.disconnect_events()

###############################################################################
        





############################################################ Rectangle selector






###############################################################################


      

###############################################################  Moving average


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
    
    if window_len%2==0:
        return y[(window_len/2):-(window_len/2)+1]
    else:
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
    
###############################################################################
    
    
    
    
    
    
    
    
    
##############################################################  Peak extrcation
       
def baseline_als(y, lam, p, niter=10):
    '''
    Baseline C orrection with Asymmetric Least Squares SmoothingPaul 
    H. C. Eilers and Hans F.M. Boelens
    October 21, 2005
    '''
    L = len(y)
    D = csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    
    for i in xrange(niter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
        
    return z

def peakExtract(array, lam, p, niter=10):  
    [m, n] = np.shape(array)
    backArray = np.zeros_like(array)
    
    for idx in range(m):
        z = baseline_als(array[idx, :], lam, p, niter=10)
        backArray[idx, :] = z
        
    return -backArray

###############################################################################
 


