# -*- coding: utf-8 -*-
import sunpy
import sunpy.map
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy import ndimage
import matplotlib.colors as colors
import matplotlib.cm as cm
import sunpy.physics.differential_rotation as diffrot
from scipy import signal
from astropy.utils.data import download_file
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


import sunpy.data.sample

map1=sunpy.map.Map('/home/talafha/Shared_Virtual/SunPy_Project/AIA20150405_0000_0094.fits') #1st Map 
map2=sunpy.map.Map('/home/talafha/Shared_Virtual/SunPy_Project/AIA20150406_0000_0094.fits') #2nd Map


'''
print(map1.meta['cdelt1'])
print(map1.meta['cdelt2'])
print(map1.meta['CUNIT1'])
print(map1.meta['CUNIT2'])
print(map1.meta['CDELT1'])
print(map1.meta['CDELT2'])

print(map1.meta['CRPIX1'])
print(map1.meta['CRPIX2'])
print(map1.meta['CRVAL1'])
print(map1.meta['CRVAL2'])
#print(map1.meta['TDIM1'])
#print(map1.meta['TDIM2'])
#print(map1.meta['PSCAL1'])
#print(map1.meta['PSCAL2'])

'''


#rotate the map 30 degrees  #### Obliquity :7.25°(to the ecliptic) , 67.23° (to the galactic plane)
map1_rotated = map1.rotate(angle = 7.25 * u.deg)
#map1_rotated.peek(draw_limb=True, draw_grid=True)


map2_rotated = map2.rotate(angle = 7.25 * u.deg)
#map2_rotated.peek(draw_limb=True, draw_grid=True)

# reduce the resolution of the image by combining the number of pixels (in each dimension) in the dimensions argument into one single pixel.
dimensions = u.Quantity(map1_rotated.dimensions) / 18
map1_rotated_superpixel = map1_rotated.superpixel(dimensions)
#map1_rotated_superpixel.peek(draw_limb=True, draw_grid=True)


map2_rotated_superpixel = map2_rotated.superpixel(dimensions)
#map2_rotated_superpixel.peek(draw_limb=True, draw_grid=True)

#headers = {'cdelt1': 2.4, 'cdelt2': 2.4, 'CUNIT1' : 'degree', 'CUNIT2' : 'degree', 'CDELT1' : 2.4, 'CDELT2' : 2.4, 'CRPIX1' : 512.5, 'CRPIX2' : 512.5, 'CRVAL1' : 0, 'CRVAL2' : 0 , 'telescop':'sunpy'}
#map22 = sunpy.map.Map(map2_rotated_superpixel.data, headers)
#map22.peek(draw_limb=True, draw_grid=True)

#calculate the corss correlation between the two maps
corr = signal.correlate2d(map1_rotated_superpixel.data, map2_rotated_superpixel.data, boundary='symm', mode='same')
#print corr
#print corr.max()
header = {'cdelt1': 152.8, 'cdelt2': 152.8, 'CROTA1' : 7.25, 'CROTA2' : -7.25, 'CRPIX1' : 9.5, 'CRPIX2' : 9.5, 'CRVAL1' : 1.3418447, 'CRVAL2' : 1.03896718, 'TDIM1' : 18, 'TDIM2' : 18 , 'telescop':'sunpy'}

corrmap = sunpy.map.Map(corr, header)
#print corrmap.coordinate_frame

c = SkyCoord(100 * u.arcsec, 10*u.arcsec, frame=corrmap.coordinate_frame)
ax = plt.subplot(projection=corrmap)
print corrmap.max()
corrmap.plot()
ax.plot_coord(c, 'o')
plt.show()
'''
corrmap.peek(draw_limb=True)#, draw_grid=True)
fig=plt.figure()
ax=fig.add_subplot(111)#, projection='3d')
corrmap.plot()
plt.colorbar()
#plt.show()



# In case of running difference, we loop through all the maps in the mapcube and 
#differentially rotate each map in the mapcube with respect to the previous map while
# in case of base difference we differentially rotate each map in the mapcube to the 
#time of the base map
mc = sunpy.map.Map([map1_rotated_superpixel, map2_rotated_superpixel], cube=True)   #Creat Mapcube
base_diffmap = []
running_diffmap = []
for i, map_i in enumerate(mc[1:]):
    aiamap_rot = diffrot.diffrot_map(map_i, time=mc[0].date)
    mc_rot = diffrot.diffrot_map(mc[i+1], time=mc[i].date)
    diffdata = map_i.data - mc_rot.data
    smap_base = sunpy.map.Map(diffdata, map_i.meta)
    diffdata = mc_rot.data - map_i.data
    smap_run = sunpy.map.Map(diffdata, map_i.meta)
    smap_base.plot_settings['cmap'] = plt.get_cmap('Greys_r')
    smap_base.plot_settings['norm'] = colors.LogNorm(100, smap_base.max())
    smap_run.plot_settings['cmap'] = plt.get_cmap('Greys_r')
    smap_run.plot_settings['norm'] = colors.LogNorm(100, smap_run.max())
    base_diffmap.append(smap_base)
    running_diffmap.append(smap_run)


result_mapcube = sunpy.map.MapCube(running_diffmap)
result_mapcube.peek()
plt.show()
'''
raw_input("Press [enter] to continue.")
