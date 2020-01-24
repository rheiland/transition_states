"""
  step_svg_substrate_grid.py - step through 2D cell plots on top of substrate plots, with the grid

"""
import sys,pathlib
import xml.etree.ElementTree as ET
import os
import math
import matplotlib.colors as mplc
from collections import deque
import scipy.io
import matplotlib
#import matplotlib.pyplot as plt  # NB! do this AFTER the TkAgg line below!
#import matplotlib.colors as mplc
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection
import numpy as np
from numpy.random import randn


try:
  # apparently we need mpl's Qt backend to do keypresses 
#  matplotlib.use("Qt5Agg")
  matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
except:
  print("\n---Error: cannot use matplotlib's TkAgg backend")
#  print("Consider installing Anaconda's Python 3 distribution.")
  raise

current_idx = 0

# if (len(sys.argv) < 3):
#   current_idx = 0
#   field_idx = 4
#   fix_cmap = 0
#   print("Usage: %s start_idx field_idx fix_cmap(0/1) vmin vmax" % sys.argv[0])
#   print("e.g. %s start_idx field_idx fix_cmap(0/1) vmin vmax\n" % sys.argv[0])
# else:
#   kdx = 1
#   current_idx = int(sys.argv[kdx]); kdx += 1
#   field_idx = int(sys.argv[kdx]); kdx += 1
#   fix_cmap = int(sys.argv[kdx]); kdx += 1
#   vmin = float(sys.argv[kdx]); kdx += 1
#   vmax = float(sys.argv[kdx]); kdx += 1


current_idx = 0
print("# args=",len(sys.argv)-1)

#for idx in range(len(sys.argv)):
use_defaults = True
show_nucleus = 0
current_idx = 0
xmin = 0.0
xmax = 1000  # but overridden by "width" attribute in .svg

vmin = 0.0
vmax = 1050
fix_cmap = 0

scale_radius = 1.0
if (len(sys.argv) == 12):
  use_defaults = False
  kdx = 1
  show_nucleus = int(sys.argv[kdx])
  kdx += 1
  current_idx = int(sys.argv[kdx])
  kdx += 1
  svg_xmin = float(sys.argv[kdx])
  kdx += 1
  svg_xmax = float(sys.argv[kdx])
  kdx += 1
  svg_ymin = float(sys.argv[kdx])
  kdx += 1
  svg_ymax = float(sys.argv[kdx])
  kdx += 1
  xmin = float(sys.argv[kdx])
  kdx += 1
  xmax = float(sys.argv[kdx])
  kdx += 1
  ymin = float(sys.argv[kdx])
  kdx += 1
  ymax = float(sys.argv[kdx])
  # kdx += 1
  # scale_radius = float(sys.argv[kdx])
  kdx += 1
  field_idx = int(sys.argv[kdx])
else:
  print("Usage:")
  # usage_str = "show_nucleus start_index svg_xmin svg_xmax svg_ymin svg_ymax xmin xmax ymin ymax scale_radius field_idx"
  usage_str = "show_nucleus start_index svg_xmin svg_xmax svg_ymin svg_ymax xmin xmax ymin ymax field_idx"
  print(usage_str)
  print("e.g.,")
  eg_str = "%s 0 0 0 200 0 200 -100 100 -100 100 0" % (sys.argv[0])
  print(eg_str)
  sys.exit(1)

#field_idx = 0
field_idx += 4

print('current_idx, field_idx = ',current_idx, field_idx)

# figure out the domain sizes (might not be square)
ifname = "initial.xml"
tree = ET.parse(ifname)
xml_root = tree.getroot()
xcoord_vals = xml_root.find(".//x_coordinates").text.split()
ycoord_vals = xml_root.find(".//y_coordinates").text.split()
xmin = float(xcoord_vals[0])
xmax = float(xcoord_vals[-1])   # should be 999.0
ymin = float(ycoord_vals[0])
ymax = float(ycoord_vals[-1])   # should be 999.0

# numx = int((xmax - xmin) / xdel)  # need to also round maybe?
# numy = int((ymax - ymin) / ydel)
numx = len(xcoord_vals)
numy = len(ycoord_vals)
print("numx, numy = ",numx,numy)


fig = plt.figure(figsize=(7,5.8))
#ax = fig.gca()

time_delay = 0.1
count = -1

cbar = None

#-----------------------------------------------------
def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83 

    Make a scatter plot of circles. 
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_)
               for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    return collection

#-----------------------------------------------------
def plot_substrate():
    global current_idx, axes_max, cbar

    # select whichever substrate index you want, e.g., for one model:
    # 4=tumor cells field, 5=blood vessel density, 6=growth substrate

    xml_file = "output%08d.xml" % current_idx
    tree = ET.parse(xml_file)
    root = tree.getroot()
#    print('time=' + root.find(".//current_time").text)
    mins = float(root.find(".//current_time").text)
    hrs = mins/60.
    days = hrs/24.
    title_str = '%d days, %d hrs, %d mins' % (int(days),(hrs%24), mins - (hrs*60))
#    print(title_str)

    fname = "output%08d_microenvironment0.mat" % current_idx
    output_dir_str = '.'
    fullname = output_dir_str + "/" + fname
    if not pathlib.Path(fullname).is_file():
        print("file not found",fullname)
        return

    info_dict = {}
    scipy.io.loadmat(fullname, info_dict)
    M = info_dict['multiscale_microenvironment']
    print('plot_substrate: field_idx=',field_idx)
    f = M[field_idx,:]   # 
    
    #N = int(math.sqrt(len(M[0,:])))
    #grid2D = M[0,:].reshape(N,N)
    xgrid = M[0, :].reshape(numy, numx)
    ygrid = M[1, :].reshape(numy, numx)

#    xvec = grid2D[0,:]
    #xvec.size
    #xvec.shape
    num_contours = 30
    num_contours = 10
#    vmin = 30.
#    vmax = 38.

    levels = MaxNLocator(nbins=30).tick_values(vmin, vmax)
#    cmap = plt.get_cmap('PiYG')
    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#    my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), num_contours, cmap='viridis') #'viridis'
    if fix_cmap > 0:
      # my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), levels=levels, cmap=cmap)
      my_plot = plt.contourf(xgrid, ygrid, M[field_idx, :].reshape(numy, numx), levels=levels, extend='both', cmap=cmap)
    else:
      # my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), cmap=cmap)
      my_plot = plt.contourf(xgrid, ygrid, M[field_idx, :].reshape(numy, numx), cmap=cmap)

    if cbar == None:
#      cbar = plt.colorbar(my_plot, boundaries=np.arange(vmin, vmax, 1.0))
      cbar = plt.colorbar(my_plot)
    else:
      cbar = plt.colorbar(my_plot, cax=cbar.ax)

#    plt.axis('equal')
    plt.title(title_str)

#    plt.show()

    png_file = "aaa%08d.png" % current_idx
#    fig.savefig(png_file)
#    plt.pause(time_delay)

#--------------------------------------------
def plot_svg():
  global current_idx, axes_max
  fname = "snapshot%08d.svg" % current_idx
  if (os.path.isfile(fname) == False):
    print("File does not exist: ",fname)
    return

  xlist = deque()
  ylist = deque()
  rlist = deque()
  rgb_list = deque()

#  print('\n---- ' + fname + ':')
  tree = ET.parse(fname)
  root = tree.getroot()
#  print('--- root.tag ---')
#  print(root.tag)
#  print('--- root.attrib ---')
#  print(root.attrib)


#  print('--- child.tag, child.attrib ---')
  numChildren = 0
  for child in root:
#    print(child.tag, child.attrib)
#    print("keys=",child.attrib.keys())
    if use_defaults and ('width' in child.attrib.keys()):
      axes_max = float(child.attrib['width'])
#      print("--- found width --> axes_max =", axes_max)
    if child.text and "Current time" in child.text:
      svals = child.text.split()
      title_str = "(" + str(current_idx) + ") Current time: " + svals[2] + "d, " + svals[4] + "h, " + svals[7] + "m"

#    print("width ",child.attrib['width'])
#    print('attrib=',child.attrib)
#    if (child.attrib['id'] == 'tissue'):
    if ('id' in child.attrib.keys()):
#      print('-------- found tissue!!')
      tissue_parent = child
      break

#  print('------ search tissue')
  cells_parent = None

  for child in tissue_parent:
#    print('attrib=',child.attrib)
    if (child.attrib['id'] == 'cells'):
#      print('-------- found cells, setting cells_parent')
      cells_parent = child
      break
    numChildren += 1


  num_cells = 0
  svg_xrange = svg_xmax - svg_xmin
  svg_yrange = svg_ymax - svg_ymin
  x_range = xmax - xmin
  y_range = ymax - ymin
#  print('------ search cells')
  for child in cells_parent:
#    print(child.tag, child.attrib)
#    print('attrib=',child.attrib)
    for circle in child:  # two circles in each child: outer + nucleus
    #  circle.attrib={'cx': '1085.59','cy': '1225.24','fill': 'rgb(159,159,96)','r': '6.67717','stroke': 'rgb(159,159,96)','stroke-width': '0.5'}
#      print('  --- cx,cy=',circle.attrib['cx'],circle.attrib['cy'])
      xval = float(circle.attrib['cx'])
      # map into desired coord sys (same as substrate mesh)
      xval = (xval-svg_xmin)/svg_xrange * x_range + xmin

      s = circle.attrib['fill']
#      print("s=",s)
#      print("type(s)=",type(s))
      if (s[0:3] == "rgb"):  # if an rgb string, e.g. "rgb(175,175,80)" 
        rgb = list(map(int, s[4:-1].split(",")))  
        rgb[:]=[x/255. for x in rgb]
      else:     # otherwise, must be a color name
        rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
        rgb = [x for x in rgb_tuple]

      # test for bogus x,y locations (rwh TODO: use max of domain?)
      too_large_val = 10000.
      if (math.fabs(xval) > too_large_val):
        print("bogus xval=",xval)
        break
      yval = float(circle.attrib['cy'])
      if (math.fabs(yval) > too_large_val):
        print("bogus xval=",xval)
        break
      # map into desired coord sys (same as substrate mesh)
      yval = (yval-svg_ymin)/svg_yrange * y_range + ymin

      rval = float(circle.attrib['r'])
#      if (rgb[0] > rgb[1]):
#        print(num_cells,rgb, rval)
      xlist.append(xval)
      ylist.append(yval)
      rlist.append(rval)
      rgb_list.append(rgb)

#     For .svg files with cells that *have* a nucleus, there will be a 2nd
      if (show_nucleus == 0):
        break

    num_cells += 1

#    if num_cells > 3:   # for debugging
#      print(fname,':  num_cells= ',num_cells," --- debug exit.")
#      sys.exit(1)
#      break

  print(fname,':  num_cells= ',num_cells)

  xvals = np.array(xlist)
  yvals = np.array(ylist)
  rvals = np.array(rlist)
  rgbs =  np.array(rgb_list)
#print("xvals[0:5]=",xvals[0:5])
#print("rvals[0:5]=",rvals[0:5])
#  print("rvals.min, max=",rvals.min(),rvals.max())

#  plt.cla()
  title_str += " (" + str(num_cells) + " agents)"
  plt.title(title_str)
  # axes range labels
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)
#  plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs)
  circles(xvals,yvals, s=rvals, color=rgbs, alpha=1.0, edgecolor='black')
#  circles(xvals,yvals, s=rvals, color=rgbs)


  xs = np.linspace(xmin,xmax,numx)
  ys = np.linspace(ymin,ymax,numy)
  hlines = np.column_stack(np.broadcast_arrays(xs[0], ys, xs[-1], ys))
  vlines = np.column_stack(np.broadcast_arrays(xs, ys[0], xs, ys[-1]))
  grid_lines = np.concatenate([hlines, vlines]).reshape(-1, 2, 2)
  line_collection = LineCollection(grid_lines, color="red", linewidths=0.5)
  ax = plt.gca()
  ax.add_collection(line_collection)
  ax.set_xlim(xs[0], xs[-1])
  ax.set_ylim(ys[0], ys[-1])

#plt.xlim(0,2000)  # TODO - get these values from width,height in .svg at top
#plt.ylim(0,2000)
  plt.pause(time_delay)

step_value = 1
def press(event):
  global current_idx, step_value
#    print('press', event.key)
  sys.stdout.flush()
  if event.key == 'escape':
    sys.exit(1)
  elif event.key == 'h':  # help
    print('esc: quit')
    print('right arrow: increment by step_value')
    print('left arrow:  decrement by step_value')
    print('up arrow:   increment step_value by 1')
    print('down arrow: decrement step_value by 1')
    print('0: reset to 0th frame')
    print('h: help')
  elif event.key == 'left':  # left arrow key
#    print('go backwards')
#    fig.canvas.draw()
    current_idx -= step_value
    if (current_idx < 0):
      current_idx = 0
    plot_substrate()
    plot_svg()
  elif event.key == 'right':  # right arrow key
#        print('go forwards')
#        fig.canvas.draw()
    current_idx += step_value
    plot_substrate()
    plot_svg()
  elif event.key == 'up':  # up arrow key
    step_value += 1
    print('step_value=',step_value)
  elif event.key == 'down':  # down arrow key
    step_value -= 1
    if (step_value <= 0):
      step_value = 1
    print('step_value=',step_value)
  elif event.key == '0':  # reset to 0th frame/file
    current_idx = 0
    plot_substrate()
    plot_svg()
  else:
    print('press', event.key)

plot_substrate()
plot_svg()
print("\nNOTE: click in plot window to give it focus before using keys.")

fig.canvas.mpl_connect('key_press_event', press)
#plot_substrate(frame_idx)
plt.show()

