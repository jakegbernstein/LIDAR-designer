#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 20:52:49 2018

@author: jake
"""
from matplotlib.pyplot import *
from LIDAR import LIDARoptics, ureg, U_, Q_
import numpy as np


lidar = LIDARoptics()
lidar.array_active = [100,100]
#ureg.setup_matplotlib(True)
fig = figure()
b = lidar.blur_defocus(Q_(1,'m'))
res = []
bits = []
Drange = np.linspace(.6,10,200)
for D in Drange:
    res.append(lidar.res_lat(D*U_('m')).to('cm').magnitude)
    bits.append(lidar.pixel_response(D*U_('m')).magnitude)
    
subplot(2,1,1)
plot(Drange,res,'.-')
ylabel('Lateral Resolution (cm)')

subplot(2,1,2)
semilogy(Drange,bits,'.-')
xlabel('Distance (m)')
ylabel('Bits per pixel per LED W ms')