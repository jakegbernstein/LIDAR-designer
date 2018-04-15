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
#dp['LED_lambda'] = (473,'nm')
#dp['lumfn_eff'] = (.2,'')
#dp['LED_Vref'] = (3.1,'volt')
#dp['LED_Iref'] = (0.350,'amp')
#dp['LED_Poutref'] = (39.8,'lumen')
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


outopticsmethods = ('hyperfocal','dnear','blur_defocus','blur_rotation','loss_geometric')
outpowermethods = ()
outdatamethods = ()
outmethods = (outopticsmethods,outpowermethods,outdatamethods)

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
                
        self.LED = LED(self.ureg)
            
        
        ### Calculate derived values for object attributes
    def LED_Pinref(self):
        return self.LED.Vref * self.LED.Iref # Input power at datasheet reference
        
    def confusion(self):
        return self.pixel_w # max circle of confusion on image array
        
    def pixel_lambdacorrection(self):
        return self.pixel_relsens/self.LED.lumfn_eff # wavelength correction
        
    ######################################################################
    ### Calculations #####################################################
    ######################################################################
    ### Calculate depth of field based on lens equation
    ### Eqs from https://en.wikipedia.org/wiki/Depth_of_field
    
    #Field of View (FoV): alpha = 2*arctan(d/2f)
    ## d = sensor size in direction
    ## f = focal length of lens

    def arrayFoV(self, printans=False, returnstring=False):
        array_FoV = 2*atan(self.array_active[0]*self.pixel_w/(2*self.f))*self.U_('radian')
        outstring = "Image Sensor Field of View = %.3f degrees".format(degrees(array_FoV))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return array_FoV
    
    def pixelFoV(self, printans=False, returnstring=False):
        pixel_FoV = 2*atan(self.pixel_w/(2*self.f))*self.U_('radian')
        outstring = "Single Pixel Field of View = %.3f degrees".format(degrees(pixel_FoV))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return pixel_FoV
    
    def pixelRes(self, D, printans=False, returnstring=False):
        pixel_res = D*self.pixelFoV()
        outstring = "Lateral Resolution at {}s = {:4.3f}s".format(D,pixel_res.to('cm'))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return pixel_res
    
    # Hyperfocal distance: H = f + f**2/(Nc)
    ## f = focal distance of lens [length]
    ## c = circle of confusion [length]
    ## N = f-number of lens aperture
    
    def hyperfocal(self, printans=False, returnstring=False):
        H = self.f + self.f**2/(self.N*self.confusion())
        H.ito(self.U_('m'))
        outstring = "Hyperfocal distance = {:4.3f}s".format(H)
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return H
        
    # Near limit DOF D_N = H*s/(H+s), s = subject distance
    ## Assume s is set to H
    ## D_N = H/2
    def dnear(self, s = None, printans=False, returnstring=False):
        H = self.hyperfocal()
        if s is None : s = H
        D_N = H*s/(H+s)
        outstring = "Near limit of DOF = {:4.3f}s\n".format(D_N)
        outstring = outstring + "Lateral resolution at {:4.3f}s = {:4.3f}s ".format(D_N, D_N.to('cm')*self.pixelFoV()/self.U_('radian'))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return D_N
    
    #Calculate blur at D_min < D_N
    # b = (f*mag/N)*(x_D/D)
    ## D = distance to foreground (e.g. D_min)
    ## x_d = difference between subject distance and D
    ## mag = f/(s-f)
    ## set s = H
    def blur_defocus(self, D = None, s = None, printans=False, returnstring=False):
        if s is None : s = self.hyperfocal()
        if D is None : D = self.D_max
        x_d = abs(D-s)
        mag = self.f/(s-self.f)
        blur_length = (self.f*mag/self.N)*(x_d/D)
        blur_length.ito('micrometer')
        blur_pixels = blur_length/self.pixel_w
        blur_pixels.ito_base_units()
        blur_resolution = blur_pixels.magnitude*self.pixelRes(D)
        outstring = "Blur at {}s = {:.1f} pixels\n".format(D,blur_pixels.magnitude)
        outstring = outstring + "Blur at {}s = {:.3}s".format(D,blur_resolution.to('cm'))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return {'pixel':blur_pixels,'space':blur_resolution}
    
    def blur_rotation(self, D=None, printans=False, returnstring=False):
        if D is None: D = self.D_max
        b_r = 2*pi*self.U_('radian')*self.nu_max*self.t_int
        b_r.ito_base_units()
        outstring = "Blur due to rotation at {}s = {:.3}s".format(D,b_r*D.to('cm')/self.U_('radian'))
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return b_r
    
    #Light capture from lens
    ## Effective aperture radius = f/(2*N)
    ## Lambertian radiation pattern ~ cos(theta)
    ## P_collectedbylens/P_total = Area_lens,effective*cos(theta)/(pi*Distance**2)
    ## P/P_total = (f/(2*N*Distance))**2
    def loss_geometric(self, D = None, printans=False, returnstring=False):
        if D is None : D = self.D_max
        lg = (self.f/(2*self.N*D))**2
        lg.ito_base_units()
        outstring = "Fraction of reflected light collected by lens at {}s = {:4e}".format(D, lg.magnitude)
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return lg
        
     # Illumination will spread out w/ same FOV as lens, so spots imaged by each pixel will receive constant illumination.
    def pixel_response(self, D = None, printans=False, returnstring=False):
        if D is None : D = self.D_max
        #illum_spot = self.LED_Poutref/(self.array_active[0]*self.array_active[1])
        illum_spot = self.LED.Poutref*(pi/4)*((self.pixelFoV()/self.LED_spotsize).to_base_units())**2
        illum_spot.ito(self.U_('lumen'))
        illum_pixel = self.albedo*illum_spot*self.loss_geometric(D)        
        flux_pixel = illum_pixel/(self.pixel_w)**2
        flux_pixel.ito(self.U_('lux'))
        #bits_per_pixel_per_watt = pixel_sens*flux_pixel*t_int*pixel_lambdacorrection/LED_Pinref
        bits_per_pixel_per_watt = self.pixel_sens*flux_pixel*self.pixel_lambdacorrection()/self.LED_Pinref()
        bits_per_pixel_per_watt.ito(1/(self.U_('W*ms')))
        outstring = "Illumination per spot imaged by each pixel = {:.3e}s\n".format(illum_spot)
        outstring = outstring + "Illumination per pixel at {}s = {:.3e}s\n".format(D,illum_pixel)
        outstring = outstring + "Radiant flux per pixel at {}s = {:.3e}\n".format(D,flux_pixel)
        outstring = outstring + "Bits per pixel/LED input energy at {}s= {:.2f}".format(D, bits_per_pixel_per_watt)
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return bits_per_pixel_per_watt
    
    def res_lat(self, D, s = None, printans=False, returnstring=False):
        b = self.blur_defocus(D, s)
        if b['pixel'] > 1:
            res = b['space']
        else:
            res = self.pixelRes(D)
        outstring = "Lateral resolution = {}".format(res)
        if printans:
            print(outstring)
        if returnstring:
            return outstring
        else:
            return res
   
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
        self.Q_ = self.ureg.Quantity
        self.U_ = self.ureg.parse_expression
        self.load_LED()
        
    def load_LED(self):
        #self.outspecfilename = self.partnum+'_spectrum.json'
        self.specfilename = self.partnum+'.json'
        with open(datadir+self.specfilename) as specfile:
            tempspec = json.load(specfile)
        for spec in tempspec.keys():
            if (type(tempspec[spec]) is list) or (type(tempspec[spec]) is tuple):
                setattr(self,spec,self.Q_(tempspec[spec][0],tempspec[spec][1]))
            else:
                setattr(self,spec,tempspec[spec])
        with open(datadir+self.outspecfile) as outspecfile:
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
                                self.lumfnarray[self.lumfnlambdaoffset+np.arange(len(self.outspecnorm))])*self.U_('')
        
     
class LIDARsummary(ttk.Frame):
    param_names = ('f','N','nu_max','t_int','LED_spotsize')
    param_LED = ('Iref','Vref','Poutref','lumfn_eff')
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
        self.in_frame = ttk.Frame(self,borderwidth=2, relief='solid',padding=5)
        self.output_frame = ttk.Frame(self,borderwidth=2, relief='solid',padding=5)
        self.graph_frame = ttk.Frame(self,borderwidth=2, relief='solid')
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
        
        # Make LED info input box        
        self.LED_frame = ttk.Labelframe(self.in_frame,text='LED Parameters',padding=5)
        self.LEDvar = tk.StringVar(self)
        self.LEDvar.set(self.lidar.LED.partnum)
        templabel = ttk.Label(self.LED_frame, text='LED partnum:')
        templabel.grid(column=0, row=0)
        self.LED_choice = ttk.Combobox(self.LED_frame, textvariable=self.LEDvar)
        self.LED_choice['values'] = LED_partnums
        self.LED_choice.state = 'readonly'
        self.LED_choice.bind('<<ComboboxSelected>>',self.changeLED)
        self.LED_choice.grid(column=1, row=0, columnspan = 2)
        j=0
        self.LEDwidgets = dict()
        self.LEDvars = dict()
        for key in self.param_LED:
            j=j+1
            self.LEDvars[key] = tk.StringVar(self)
            self.LEDvars[key].set(getattr(self.lidar.LED,key).magnitude)
            templabel = ttk.Label(self.LED_frame, text=key)
            templabel.grid(column = 0, row = j)
            templabel = ttk.Label(self.LED_frame, text=str(getattr(self.lidar.LED,key).units))
            templabel.grid(column=2, row=j)
            self.LEDwidgets[key] = ttk.Entry(self.LED_frame, textvariable = self.LEDvars[key])
            self.LEDwidgets[key].grid(column =1, row=j)
        self.LED_frame.grid(column=0, row = i, columnspan=3)
        self.changeLED('')
        # Make buttons
        self.updatebutton =  ttk.Button(self.in_frame, text = 'Update Design', command = self.update_design)
        self.updatebutton.grid(column = 0, row = len(self.param_names)+2, columnspan = 2)        
        # Set up output parameter frame
        out_title = ttk.Label(self.output_frame, text = 'Design Output')
        out_title.grid(column = 0, row = 0)
        #out_optical = ttk.Frame(self.output_frame,borderwidth=1, relief='solid')        
        #out_power = ttk.Frame(self.output_frame,borderwidth=1, relief='solid')
        #out_data = ttk.Frame(self.output_frame,borderwidth=1, relief='solid')
        
        ### Usinge Outbox class to keep code manageable
        #The Outbox is a frame with a number of message boxes inside, which are
        #linked to StringVars kept in a dictionary keyed to the method names
        #that produce the desired output strings
        #Outbox is passed in a list of method names and an empty dictionary
        outopticsvars = dict()
        outpowervars = dict()
        outdatavars = dict()
        self.outvars = (outopticsvars,outpowervars,outdatavars)
        self.out_optical = Outbox(self.output_frame,'Optical Specs',outopticsmethods,outopticsvars)
        self.out_power = Outbox(self.output_frame,'Power Specs',outpowermethods,outpowervars)
        self.out_data = Outbox(self.output_frame,'Data Specs',outdatamethods,outdatavars)
        self.update_outputs()
        self.out_optical.grid(column=0, row=1, sticky='ew')
        self.out_power.grid(column=0, row=2, sticky='ew')
        self.out_data.grid(column=0, row=3, sticky='ew')

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
        
    def changeLED(self,event):
        #print(event)
        self.LED_choice.selection_clear()
        if self.LEDvar.get() != 'Custom':
            self.lidar.LED.partnum = self.LEDvar.get()
            self.lidar.LED.load_LED()
        for par in self.param_LED:
            if self.LEDvar.get() == 'Custom':
                #print('Custom')
                self.LEDwidgets[par]['state'] = 'normal'
            else:
                #print('Read only')
                #print(self.LEDwidgets[par]['state'])
                self.LEDwidgets[par]['state'] = 'normal'
                self.LEDvars[par].set(getattr(self.lidar.LED,par).magnitude)
                self.LEDwidgets[par]['state'] = 'readonly'
        
    def update_outputs(self):
        for i in range(len(outmethods)):
            #map(lambda x: outvars[i][x] = getattr(self,x)(returnstring=True), outmethods[i])
            #map(lambda x: getattr(self,x)(returnstring=True), outmethods[i])
            for meth in outmethods[i]:
                #print(type(self.outvars[i][meth]))
                self.outvars[i][meth].set(getattr(self.lidar,meth)(returnstring=True))
            
        
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
        self.grid(column = 0, row = 0, sticky='nsew')
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
        
        self.showwork = tk.BooleanVar(self)
        self.showwork.set(False)
        #print(self.showwork.get())
        self.showworkbutton = ttk.Checkbutton(self.statusframe, text = 'Show Work?',
                                        command = self.work_changed, variable = self.showwork,
                                        onvalue = True, offvalue = False)
        self.stopbutton = ttk.Button(self.statusframe, text = 'STOP', command = self.quit)
        self.stopbutton.grid(column=0, row=0)
        self.quitbutton = ttk.Button(self.statusframe, text = 'EXIT', command = self.exit)
        self.quitbutton.grid(column = 1, row = 0)
        self.showworkbutton.grid(column = 2, row = 0)
        self.update()
        self.mainloop()
        
    def work_changed(self):
        #print(self.showwork.get())
        if self.showwork.get():
            self.workwindow = tk.Toplevel(self)
            self.workframe = LIDARwork(self.workwindow,self.lidar_optic)
            self.update()
        else:            
            self.workwindow.destroy()
        
        
    def exit(self):
        plt.close()
        self.destroy()
#    def quitloop(self.):
#        window.quit()
        
class Outbox(ttk.Labelframe):
    def __init__(self,parent,label,methodnames,vardict):
        ttk.Labelframe.__init__(self,parent,text=label)
        self.widgetdict = dict()
        i = 0
        for var in methodnames:
            vardict[var] = tk.StringVar()
            #print(type(vardict[var]))
            self.widgetdict[var] = tk.Message(self,textvariable=vardict[var], width=500)
            self.widgetdict[var].grid(column = 1, row = i, sticky='nw')
            templabel = ttk.Label(self,text=var)
            templabel.grid(column=0, row = i, sticky='nw')
            i = i+1
            
            