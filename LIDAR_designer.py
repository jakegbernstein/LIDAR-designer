#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class defintion for LIDAR optics system

Created on Wed Mar 29
@author: Jake Bernstein
"""
### Import necessary modules    
import sys
from math import pi, atan, degrees
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import tkinter as tk
from tkinter import ttk
import json

### Have to be careful to only have one unit registry for pint
if not (('pint' in locals()) or ('pint' in globals())):
    import pint
    
#if not (('ureg' in locals()) or ('ureg' in globals())):
#    global ureg, U_, Q_
#    ureg = pint.UnitRegistry()
#    U_ = ureg.parse_expression
#    Q_ = ureg.Quantity 

defaultparams = dict()
dp = defaultparams
#Lens Parameters
dp['f'] = (6,'mm')
dp['N'] = (1,'')
dp['D_max'] = (10,'m')
dp['D_min'] = (.6,'m')
#LED Parameters
#dp['LED_label'] = 'Custom'
dp['LED_lambda'] = (473,'nm')
dp['lumfn_eff'] = (.2,'')
dp['LED_Vref'] = (3.1,'volt')
dp['LED_Iref'] = (0.350,'amp')
dp['LED_Poutref'] = (39.8,'lumen')
dp['LED_spotsize'] = (2,'degrees')
dp['albedo'] = (0.25,'')
dp['pixel_w'] = (20,'micrometer')
dp['array_size'] = ([320,240],'')
dp['array_active'] = ([320,240],'')
dp['pixel_relsens'] = (0.25,'')
dp['pixel_sens'] = (150e3,'1/(lux*s)')
dp['t_int'] = (1,'ms')
dp['nu_max'] = (0.1,'1/s')
dp['v_max'] = (0.5,'m/s')

LED_partnums = ('CXA1820','CXA1830','Custom')
datadir = './Data/'
luminosityfile = 'luminosity_spectrum.json'
sensitivityfile = 'epc660_sensitivity.json'

class LIDARoptics:

    def __init__(
            self,
            ureg = None,
            **kwargs
            # Lens Parameters
#            f = 6 * U_('mm'), # focal distance of lens [millimeters]
#            N = 1 * U_(''), # f-stop of lens aperture
#            D_max = 10 * U_('m'), # maximum imaging distance specified
#            D_min = 0.6 * U_('m'), # minimum imaging distance specified             
#            # Illumination Parameters (Cree XLAMP XB-D)
#            LED_lambda = 473 * U_('nm'), # wavelength of illumination
#            lumfn = .2 * U_(''), # luminosity function value at LED_lambda
#            LED_Vref = 3.1 * U_('volt'), # Input voltage at datasheet reference
#            LED_Iref = 0.350 * U_('amp'), # Input current at datasheet reference
#            LED_Poutref = 39.8 * U_('lumen'), # Optical power output at datasheet reference
#            LED_spotsize = 2 * U_('degrees'), # Spot size of each LED
#            albedo = 0.25 * U_(''), # reflectivity of cave wall                
#            # Image Sensor Parameters (epc660)
#            pixel_w = 20 * U_('um'), # width of pixel [microns]            
#            array_size = [320, 240], # size of image array [pixels]
#            array_active = [320,240], # size of active image area [pixels]
#            pixel_relsens = 0.25 * U_(''), # relative sensitivity of pixel at lambda_illum
#            pixel_sens = 150e3 /(U_('lux') * U_('s')), # pixel sensitivy in LSBs
#            t_int = 0.5 * U_('ms'), # amount of time light is pulsed for LIDAR image
#            # Mechanical Parameters
#            nu = 1 / U_('s'), # rotational frequency of lens               
#            v_max = 0.5 * U_('m') / U_('s'),            
            ):
#        self.f = f
#        self.N = N
#        self.D_max = D_max
#        self.D_min = D_min
#        self.nu = nu
#        self.LED_lambda = LED_lambda
#        self.lumfn = lumfn
#        self.LED_Vref = LED_Vref
#        self.LED_Iref = LED_Iref
#        self.LED_Poutref = LED_Poutref
#        self.LED_spotsize = LED_spotsize
#        self.albedo = albedo
#        self.pixel_w = pixel_w
#        self.array_size = array_size
#        self.array_active = array_active
#        self.pixel_relsens = pixel_relsens
#        self.pixel_sens = pixel_sens
#        self.t_int = t_int
        
        if (ureg == None) or (type(ureg) != type(pint.UnitRegistry())):
            self.ureg = pint.UnitRegistry()
        self.U_ = self.ureg.parse_expression
        self.Q_ = self.ureg.Quantity
        
        for par in defaultparams.keys():
            if par in kwargs.keys():
                setattr(self,par,self.Q_(kwargs[par][0],kwargs[par][1]))
            else:
                setattr(self,par,self.Q_(dp[par][0],dp[par][1]))
            
        
        ### Calculate derived values for object attributes
    def LED_Pinref(self):
        return self.LED_Vref * self.LED_Iref # Input power at datasheet reference
        
    def confusion(self):
        return self.pixel_w # max circle of confusion on image array
        
    def pixel_lambdacorrection(self):
        return self.pixel_relsens/self.lumfn # wavelength correction
        
    ######################################################################
    ### Calculations #####################################################
    ######################################################################
    ### Calculate depth of field based on lens equation
    ### Eqs from https://en.wikipedia.org/wiki/Depth_of_field
    
    #Field of View (FoV): alpha = 2*arctan(d/2f)
    ## d = sensor size in direction
    ## f = focal length of lens

    def arrayFoV(self, printans=False):
        array_FoV = 2*atan(self.array_active[0]*self.pixel_w/(2*self.f))*self.U_('radian')
        if printans:
            print("Image Sensor Field of View = %.3f degrees".format(degrees(array_FoV)))
        return array_FoV
    
    def pixelFoV(self, printans=False):
        pixel_FoV = 2*atan(self.pixel_w/(2*self.f))*self.U_('radian')
        if printans:
            print("Single Pixel Field of View = %.3f degrees".format(degrees(pixel_FoV)))
        return pixel_FoV
    
    def pixelRes(self, D, printans=False):
        pixel_res = D*self.pixelFoV()
        if printans:
            print("Lateral Resolution at {}s = {:4.3f}s".format(D,pixel_res.to('cm')))
        return pixel_res
    
    # Hyperfocal distance: H = f + f**2/(Nc)
    ## f = focal distance of lens [length]
    ## c = circle of confusion [length]
    ## N = f-number of lens aperture
    
    def hyperfocal(self, printans=False):
        H = self.f + self.f**2/(self.N*self.confusion())
        H.ito(self.U_('m'))
        if printans:
            print("Hyperfocal distance = {:4.3f}s".format(H))
        return H
        
    # Near limit DOF D_N = H*s/(H+s), s = subject distance
    ## Assume s is set to H
    ## D_N = H/2
    def dnear(self, s = None, printans=False):
        H = self.hyperfocal()
        if s is None : s = H
        D_N = H*s/(H+s)
        if printans:
            print("Near limit of DOF = {:4.3f}s".format(D_N))
            print("Lateral resolution at {:4.3f}s = {:4.3f}s ".format(D_N, D_N.to('cm')*self.pixel_FoV()))
        return D_N
    
    #Calculate blur at D_min < D_N
    # b = (f*mag/N)*(x_D/D)
    ## D = distance to foreground (e.g. D_min)
    ## x_d = difference between subject distance and D
    ## mag = f/(s-f)
    ## set s = H
    def blur_defocus(self, D, s = None, printans=False):
        if s is None : s = self.hyperfocal()
        x_d = abs(D-s)
        mag = self.f/(s-self.f)
        blur_length = (self.f*mag/self.N)*(x_d/D)
        blur_length.ito('micrometer')
        blur_pixels = blur_length/self.pixel_w
        blur_pixels.ito_base_units()
        blur_resolution = blur_pixels.magnitude*self.pixelRes(D)
        if printans:
            print("Blur at {}s = {} pixels".format(D,blur_pixels.magnitude))
            print("Blur at {}s = {}s".format(D,blur_resolution.to('cm')))
        return {'pixel':blur_pixels,'space':blur_resolution}
    
    def blur_rotation(self, D, printans=False):
        b_r = 2*pi*self.nu*self.t_int
        return b_r
    
    #Light capture from lens
    ## Effective aperture radius = f/(2*N)
    ## Lambertian radiation pattern ~ cos(theta)
    ## P_collectedbylens/P_total = Area_lens,effective*cos(theta)/(pi*Distance**2)
    ## P/P_total = (f/(2*N*Distance))**2
    def loss_geometric(self, D = None, printans=False):
        if D is None : D = self.D_max
        lg = (self.f/(2*self.N*D))**2
        lg.ito_base_units()
        if printans:
            print("Fraction of reflected light collected by lens at {}s = {:4e}".format(D, lg.magnitude))
        return lg
        
     # Illumination will spread out w/ same FOV as lens, so spots imaged by each pixel will receive constant illumination.
    def pixel_response(self, D = None, printans=False):
        if D is None : D = self.D_max
        #illum_spot = self.LED_Poutref/(self.array_active[0]*self.array_active[1])
        illum_spot = self.LED_Poutref*(pi/4)*((self.pixelFoV()/self.LED_spotsize).to_base_units())**2
        illum_spot.ito(self.U_('lumen'))
        illum_pixel = self.albedo*illum_spot*self.loss_geometric(D)        
        flux_pixel = illum_pixel/(self.pixel_w)**2
        flux_pixel.ito(self.U_('lux'))
        #bits_per_pixel_per_watt = pixel_sens*flux_pixel*t_int*pixel_lambdacorrection/LED_Pinref
        bits_per_pixel_per_watt = self.pixel_sens*flux_pixel*self.pixel_lambdacorrection()/self.LED_Pinref()
        bits_per_pixel_per_watt.ito(1/(self.U_('W*ms')))
        if printans:
            print("Illumination per spot imaged by each pixel = {:.3e}s".format(illum_spot))
            print("Illumination per pixel at {}s = {:.3e}s".format(D,illum_pixel)) 
            print("Radiant flux per pixel at {}s = {:.3e}".format(D,flux_pixel)) 
            print("Bits per pixel/LED input energy at {}s= {:.2f}".format(D, bits_per_pixel_per_watt))
        return bits_per_pixel_per_watt
    
    def res_lat(self, D, s = None, printans=False):
        b = self.blur_defocus(D, s)
        if b['pixel'] > 1:
            return b['space']
        else:
            return self.pixelRes(D)
   
class LED:
    def __init__(self, ureg,
                 partnum = None):
        if (partnum == None) or (partnum not in LED_partnums):
            self.partnum = LED_partnums[0]
        else:
            self.partnum = partnum
        with open(datadir+luminosityfile) as lumfile:
            self.lumfn = json.load(lumfile)
        self.ureg = ureg
        self.Q_ = ureg.quantity
        self.load_LED()
        
    def load_LED(self):
        #self.outspecfilename = self.partnum+'_spectrum.json'
        self.specfilename = self.partnum+'.json'
        with open(datadir+self.specfilename) as specfile:
            tempspec = json.load(specfile)
        for spec in tempspec.keys():
            if type(tempspec[spec]) is tuple:
                setattr(self,spec,self.Q_(tempspec[spec][0],tempspec[spec][1]))
            else:
                setattr(self,spec,tempspec[spec])
        with open(datadir+self.outspecfilename) as outspecfile:
            self.outspec = json.load(outspecfile)
        #Normalize output spectrum to integrate to 1
        self.outspecnorm = np.asarray(self.outspec['y'])
        self.outspecnorm = self.outspecnorm/np.sum(self.outspecnorm)
        self.outspeclambda = np.asarray(self.outspec['x'],dtype=int)
        #Match luminosity function range to LED output spectrum
        self.lumfnlambda = np.asarray(self.lumfn['x'],dtype=int)
        self.lumfnarray  = np.asarray(self.lumfn['y'])
        self.lumfnlambdaoffset = self.lumfnlambda.tolist().index(self.outspeclambda[0])
        self.lumfn_eff = np.sum(self.outspecnorm*
                                self.lumfnarray[self.lumfnlambdaoffset+np.arange(len(self.outspecnorm))])
        
     
class LIDARsummary(ttk.Frame):
    param_names = ('f','N','nu','t_int','LED_spotsize')
    def __init__(self,parent, lidar_optic):
        ttk.Frame.__init__(self,parent,borderwidth=2, relief='solid')
        self.lidar = lidar_optic
        self.ureg = lidar_optic.ureg
        self.Q_ = lidar_optic.Q_
        self.U_ = lidar_optic.U_
        
        self.initialize()
        
    def initialize(self):
        self.grid(column = 0, row = 0)
        # Set up three subframes
        self.in_frame = ttk.Frame(self,borderwidth=1, relief='solid')
        self.output_frame = ttk.Frame(self,borderwidth=1, relief='solid')
        self.graph_frame = ttk.Frame(self,borderwidth=1, relief='solid')
        self.in_frame.grid(column = 0, row = 0, sticky='n')
        self.output_frame.grid(column = 1, row = 0, sticky ='n')
        self.graph_frame.grid(column = 2, row = 0, sticky = 'nsew')
        self.rowconfigure(0,weight=1)
        self.columnconfigure(2,weight=1)
        # Set up input parameter frame
        in_title = ttk.Label(self.in_frame,text = 'Input Parameters')
        in_title.grid(column = 0, row = 0, columnspan = 3)
        i =0
        self.entrywidgets = dict()
        self.entryvars = dict()
        for key in self.param_names:
            i = i+1
            self.entryvars[key] = tk.StringVar(self)
            self.entryvars[key].set(getattr(self.lidar,key).magnitude)
            templabel = ttk.Label(self.in_frame, text = key)
            templabel.grid(column = 0, row = i)
            templabel = ttk.Label(self.in_frame, text = str(getattr(self.lidar,key).units))
            templabel.grid(column = 2, row = i)
            #tempentry = ttk.Entry(in_frame, textvariable = getattr(lidar,params_in[i]))
            #print(getattr(lidar,params_in[i]))
            self.entrywidgets[key] = ttk.Entry(self.in_frame, textvariable = self.entryvars[key])
            self.entrywidgets[key].grid(column = 1, row = i)
        # Make buttons
        self.updatebutton =  ttk.Button(self.in_frame, text = 'Update Design', command = self.update_design)
        self.updatebutton.grid(column = 0, row = len(self.param_names)+1, columnspan = 2)        
        # Set up output parameter frame
        out_title = ttk.Label(self.output_frame, text = 'Design Output')
        out_title.grid(column = 1, row = 0)
        #Set up matplotlib figure
        self.fig = plt.figure()
        self.fig.suptitle('Graphs')
        self.ax_latres = self.fig.add_subplot(2,1,1)
        self.ax_bitres = self.fig.add_subplot(2,1,2)
        #Import matplotlib figure to tk frame
        self.graph_canvas = FigureCanvasTkAgg(self.fig,self.graph_frame)
        #self.graph_canvas = FigureCanvasTkAgg(self.fig,self)
        self.graph_canvas.get_tk_widget().grid(column = 0, row = 0, sticky='nsew')
        self.graph_frame.rowconfigure(0, weight=1)
        self.graph_frame.columnconfigure(0, weight=1)

        self.plotresults()
        #self.graph_canvas.draw
        #self.resizable(True,False)
        #self.update()


    def update_design(self):
        for key in self.param_names:
            setattr(self.lidar, key, float(self.entryvars[key].get()) * getattr(self.lidar,key).units)
            print(getattr(self.lidar,key))
        self.plotresults()
        
    def plotresults(self):
        res = []
        bits =[]
        D_range = np.linspace(self.lidar.D_min.magnitude,self.lidar.D_max.magnitude,200)
        for D in D_range:
            res.append(self.lidar.res_lat(D*self.U_('m')).to('cm').magnitude)
            bits.append(self.lidar.pixel_response(D*self.U_('m')).magnitude)
        self.ax_latres.clear()
        self.ax_bitres.clear()
        self.ax_latres.plot(D_range,res)
        self.ax_latres.set_ylabel('Lateral Resolution at D_max [cm]')
        self.ax_bitres.semilogy(D_range,bits)
        # Plot best sensitivity lines for bit resolution
#        self.ax_bitres.lines([self.lidar.D_min.magnitude, self.lidar.D_max.magnitude],
#                             [200, 200])
#        self.ax_bitres.lines([self.lidar.D_min.magnitude, self.lidar.D_max.magnitude],
#                             [1500, 1500])
        self.ax_bitres.hlines([200,1500],self.lidar.D_min.magnitude, self.lidar.D_max.magnitude)
        self.ax_bitres.set_xlabel('Distance [m]')
        self.ax_bitres.set_ylabel('Bits per pixel per frame')        
        #self.fig.draw()
        self.graph_canvas.draw()
        self.update()
        
class LIDARwork(ttk.Frame):
    def __init__(self,parent,lidar_optic):
        ttk.Frame.__init__(self,parent,borderwidth=2, relief='solid')
        self.lidar = lidar_optic
        self.ureg = lidar_optic.ureg
        self.Q_ = lidar_optic.Q_
        self.U_ = lidar_optic.U_
        
        self.initialize()
        
    def initialize(self):
        self.grid(column = 0, row = 0, stickey='nsew')
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        
class LIDARGUI(tk.Tk):    
    def __init__(self, lidar_optic, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        #tk.Tk.wm_title('LIDAR Design GUI')
        
        self.lidar_optic = lidar_optic
        self.mainframe =  LIDARsummary(self,self.lidar_optic)
        self.mainframe.grid(column=0, row=0, sticky='nsew')
        self.rowconfigure(0,weight = 1)
        self.columnconfigure(0,weight = 1)
        
        self.statusframe = ttk.Frame(self,borderwidth=2,relief='solid')
        self.statusframe.grid(column=0, row=1, sticky='w')
        
        #self.showwork = tk.BooleanVar()
        self.showwork = False
        self.showworkbutton = ttk.Checkbutton(self.statusframe, text = 'Show Work?',
                                        command = work_changed, variable = self.showwork,
                                        onvalue = True, offvalue = False)
        self.stopbutton = ttk.Button(self.statusframe, text = 'STOP', command = self.quit)
        self.stopbutton.grid(column=0, row=0)
        self.quitbutton = ttk.Button(self.statusframe, text = 'EXIT', command = self.exit)
        self.quitbutton.grid(column = 1, row = 0)
        self.showworkbutton.grab_current(column = 2, row = 0)
        self.update()
        self.mainloop()
        
    def work_changed(self):
        if self.showwork:
            self.workwindow = tk.Toplevel(self)
            self.workframe = LIDARwork(self,self.lidar_optic)
            self.update()
        else:
            self.workwindow.destroy()
        
        
    def exit(self):
        plt.close()
        self.destroy()
#    def quitloop(self.):
#        window.quit()