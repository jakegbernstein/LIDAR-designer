#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:52:49 2018

@author: jake
"""
### Import Modules
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import tkinter as tk
from tkinter import ttk
#from matplotlib.widgets import TextBox
from LIDAR_designer import LIDARoptics, ureg, U_, Q_

### Initialize Figure
lidar = LIDARoptics()
window = tk.Tk()
mainframe = ttk.Frame(window)
in_frame = ttk.Frame(mainframe)
static_frame = ttk.Frame(mainframe)
output_frame = ttk.Frame(mainframe)
graph_frame = ttk.Frame(mainframe)

### Set Parameters
param_names = ('f','N','nu')
params_in = dict()
for name in param_names:
    params_in[name] = str(getattr(lidar,name).magnitude)

### Test code
def update_design():
    for key in param_names:
        setattr(lidar, key, float(entryvars[key].get()) * getattr(lidar,key).units)
        print(getattr(lidar,key))
    
def quitloop():
    window.quit()


mainframe.grid(column = 0, row = 0)
in_frame.grid(column = 0, row = 0)
output_frame.grid(column = 1, row = 0)
graph_frame.grid(column = 2, row = 0)

f=6
### Draw Input Parameter frame
inrows = 1+len(params_in)
in_title = ttk.Label(in_frame,text = 'Input Parameters')
in_title.grid(column = 0, row = 0, columnspan = 2)

i =0
entrywidgets = dict()
entryvars = dict()
for key in param_names:
    i = i+1
    entryvars[key] = tk.StringVar(window)
    entryvars[key].set(getattr(lidar,key).magnitude)
    templabel = ttk.Label(in_frame, text = key)
    templabel.grid(column = 0, row = i)
    templabel = ttk.Label(in_frame, text = str(getattr(lidar,key).units))
    templabel.grid(column = 2, row = i)
    #tempentry = ttk.Entry(in_frame, textvariable = getattr(lidar,params_in[i]))
    #print(getattr(lidar,params_in[i]))
    entrywidgets[key] = ttk.Entry(in_frame, textvariable = entryvars[key])
    entrywidgets[key].grid(column = 1, row = i)
    
    
updatebutton = ttk.Button(in_frame, text = 'Update Design', command = update_design)
updatebutton.grid(column = 0, row = len(params_in)+1, columnspan = 2)


quitbutton = ttk.Button(in_frame, text = 'QUIT', command = quitloop)
quitbutton.grid(column = 0, row = len(params_in)+2, columnspan = 2)



fig = plt.figure()
fig.suptitle('LIDAR System Design')
figax =fig.add_axes([0,0,1,1])
figax.plot([1,2,3,4,5])

graph_canvas = FigureCanvasTkAgg(fig,graph_frame)
graph_canvas.get_tk_widget().grid(column = 0, row = 0)
graph_canvas.draw()

#graph_canvas = FigureCanvasAgg(fig,graph_frame)

#graph_frame.grid(column=2,row=1)

#ax_inputs = fig.add_subplot(1,3,1)
#ax_outputs = fig.add_subplot(1,3,2)
#ax_graphs = fig.add_subplot(1,3,3)
#
#for ax in fig.get_axes():
#    ax.set_xticks([])
#    ax.set_yticks([])
#
#ax_inputs.set_title('Input Parameters')
#ax_outputs.set_title('Derived Parameters')
#ax_graphs.set_title('Performance Graphics')

#ax_inputs.set_ylim(0,len(params_in))
#ax_inputs.set_xlim(0,1)

#for p in params_in:
#    ax_inputs.text(0,1+params_in.index(p),p)
#    
#
#
#
#### Do stuff
#
#lidar.array_active = [100,100]
##ureg.setup_matplotlib(True)
#b = lidar.blur_defocus(Q_(1,'m'))
#res = []
#bits = []
#Drange = np.linspace(.6,10,200)
#for D in Drange:
#    res.append(lidar.res_lat(D*U_('m')).to('cm').magnitude)
#    bits.append(lidar.pixel_response(D*U_('m')).magnitude)
#    
#
#ax_graphs.plot(Drange,res,'.-')
#ax_graphs.set_ylabel('Lateral Resolution (cm)')

#subplot(2,1,2)
#semilogy(Drange,bits,'.-')
#xlabel('Distance (m)')
#ylabel('Bits per pixel per LED W ms')
window.mainloop()