from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from siphon.catalog import TDSCatalog
import numpy as np
import os
from scipy.interpolate import interp1d

# get zone for temperature
def get_zone(zone_inds, mint_zone, maxt_zone, tmp):
    nzone = len(mint_zone)
    zi_vals = np.arange(nzone)+1
    zi = zi_vals[(tmp>=np.array(mint_zone))&(tmp<np.array(maxt_zone))][0]
    #print(zi)
    zone_ind = zone_inds[zi]
    return zone_ind

# define crystal habit zones
zone_inds = {1:{'min_tmp':-70.,'max_tmp':-25.,'name':'Polycrystalline','color':'gray'},
             2:{'min_tmp':-25.,'max_tmp':-18.,'name':'Planar','color':'red'},
             3:{'min_tmp':-18.,'max_tmp':-12.,'name':'Dendritic','color':'yellow'},
             4:{'min_tmp':-12.,'max_tmp':-8.,'name':'Planar','color':'red'},
             5:{'min_tmp':-8.,'max_tmp':-5.,'name':'Columnar','color':'blue'},
             6:{'min_tmp':-5.,'max_tmp':0.,'name':'Aggregate/Planar','color':'m'},
             7:{'min_tmp':0.,'max_tmp':60,'name':'Rain','color':'green'}}

mint_zone = []
maxt_zone = []
for zk in zone_inds.keys():
    mint_zone.append(zone_inds[zk]['min_tmp'])
    maxt_zone.append(zone_inds[zk]['max_tmp'])

'''    
# read temperature profile
data = np.genfromtxt('profile.txt', skip_header=1)
t = data[:,1]
gph = data[:,0]
nz = len(t)
'''

# set date and time to download
#dt_request = datetime.utcnow()
year = 2022
month = 1
day = 19
hour = 12
dt_request = datetime(year, month, day, hour)

# read in latest gfs data
#best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/'
#                      'Global_0p5deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p5deg/Best')
best_gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/'
                      'Global_0p25deg/latest.xml?dataset=grib/NCEP/GFS/Global_0p25deg/'
                      f'GFS_Global_0p25deg_{year:04d}{month:02d}{day:02d}_{hour:02d}00.grib2')
best_ds = best_gfs.datasets[0]
print(best_ds)
ncss = best_ds.subset()
query = ncss.query()
#print(ncss.variables)

# get subset
print('downloading data...')


query.lonlat_point(-76.11,35.95).time(dt_request+timedelta(hours=72))
#query.lonlat_point(-89.30,48.76)
query.variables('Temperature_isobaric',
                'Geopotential_height_isobaric',
                'Vertical_velocity_pressure_isobaric').accept('netcdf')
query.add_lonlat(True)
#query.vertical_level(50000.)
#query.all_times()
data = ncss.get_data(query)

# get variable fields
print(data.variables)
t = np.squeeze(data.variables['Temperature_isobaric'][:])-273.15
gph = np.squeeze(data.variables['Geopotential_height_isobaric'][:])
omeg = np.squeeze(data.variables['Vertical_velocity_pressure_isobaric'][:])
pres = np.squeeze(data.variables['isobaric1'][:])
nz = len(t)

lon = data.variables['longitude'][0]
lat = data.variables['latitude'][0]
time_hrs = data.variables['time'][0][0]
time_units = data.variables['time'].units
print(time_units, time_hrs)

dt_units = datetime.strptime(time_units, 'Hour since %Y-%m-%dT%H:%M:%SZ')
dt_forecast = dt_units+timedelta(hours=time_hrs)
dstr_forecast = dt_forecast.strftime('%d %b %Y at %H UTC')
print(dt_units, dt_forecast)

'''
plt.plot(omeg, pres/100., 'k-')
ax = plt.gca()
ax.set_ylim([200.,1000.])
ax.invert_yaxis()

plt.savefig('test.png')
'''

# subsample (interpolate) gph-t data
nss = 501
gph_ss = np.linspace(np.min(gph), np.max(gph), nss)
fzt = interp1d(gph, t)
t_ss = fzt(gph_ss)
fzo = interp1d(gph, omeg)
omeg_ss = fzo(gph_ss)

# plot profile

# change fonts
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text', usetex=True)

fig = plt.figure(figsize=(8,5))
#ax1 = fig.add_subplot()
#ax2 = fig.add_subplot(spec[0,3])
ax1 = plt.axes([0.1,0.1,0.7,0.8])
ax2 = plt.axes([0.85,0.1,0.1,0.8])

ax1.plot(t, gph*3.28e-3, 'k--', lw=1.)
#plt.plot([-45.,0.], [z_mint_kft[-1],z_mint_kft[-1]], 'm--', lw=2.)
ax1.set_xlim([-45.,0.])
ax1.set_ylim([0.,25.])

#ax2.plot(omeg, gph*3.28e-3, 'k-', lw=1.)
#ax2.plot([0.,0.], [0.,25.], 'b--', lw=1.)

max_z = gph[np.argmin(np.abs(gph-25./3.28e-3))]

# loop over sounding and color by zone
cmap = plt.get_cmap('magma_r')

for i in range(nss-1):
    if gph_ss[i]<max_z+50.:
        #print(i, gph[i], t[i])
        
        # get zone by mean temp
        mn_tmp = 0.5*(t_ss[i]+t_ss[i+1]) 
        mn_omeg = 0.5*(omeg_ss[i]+omeg_ss[i+1])
        print(gph_ss[i], mn_tmp, mn_omeg)    
        zone_ind = get_zone(zone_inds, mint_zone, maxt_zone, mn_tmp)
        #print(zone_ind)
        
        # add habit region
        rect = patches.Rectangle((-45.,gph_ss[i+1]*3.28e-3), 45.,(gph_ss[i]-gph_ss[i+1])*3.28e-3,
                                 facecolor=zone_ind['color'], edgecolor='none', alpha=0.3)
        ax1.add_patch(rect)
        
        # add lift region
        dz = (gph_ss[i]-gph_ss[i+1])*3.28e-3
        
        if mn_omeg<0:
            rect = patches.Rectangle((0.,gph_ss[i+1]*3.28e-3-dz*0.25), mn_omeg, dz/2,
                                     facecolor=cmap(np.clip(np.abs(mn_omeg)/3., 0., 1.)), linewidth=0)
        else:
            rect = patches.Rectangle((0.,gph_ss[i+1]*3.28e-3-dz*0.25), mn_omeg, dz/2,
                                     facecolor='gray', linewidth=0)
        ax2.add_patch(rect)
        
        
# add dummy plot for legend
names = []
for zk in zone_inds.keys():
    color = zone_inds[zk]['color']
    name = zone_inds[zk]['name']
    
    if not (name in names):
        names.append(name)
        ax1.scatter(t, gph+1.e6, c=color, s=100, alpha=0.3, marker='s', label=name)
'''
# plot omega values
tvals = np.linspace(-45., 0., 46)
t2, gp2 = np.meshgrid(tvals, gph, indexing='ij')
t2, omeg2 = np.meshgrid(tvals, omeg, indexing='ij')

plt.scatter(t2[omeg2<0.], gp2[omeg2<0.]*3.28e-3, c=omeg2[omeg2<0.], s=20, cmap='inferno_r', alpha=1., marker='o', label='Omega')
'''
#plt.legend(loc=1)
h, l = ax1.get_legend_handles_labels()
ax1.legend(h, l, loc=1)

ax1.set_xlabel('Temperature ($^{\circ}$ C)')
ax1.set_ylabel('Height ASL (kft)')

ax1.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))

title = f'Crystal habit zones ({lat:.2f}N, {abs(lon):.2f}W) - GFS valid {dstr_forecast}'
ax1.set_title(title, ha='left', x=0.)

ax2.set_title('Omega')
ax2.set_xlabel('(Pa s$^{\sf{-1}}$)')
ax2.set_xlim([-3.,1.])
ax2.set_ylim([0.,25.])
ax2.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))
ax2.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
#ax2.invert_xaxis()
plt.savefig('habits.png', dpi=150)

