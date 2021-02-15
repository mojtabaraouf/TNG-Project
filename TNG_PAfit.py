import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
# import pafit
from pafit.fit_kinematic_pa import fit_kinematic_pa
# from plotbin.symmetrize_velfield import symmetrize_velfield
# from plotbin.plot_velfield import plot_velfield

# start up your interface of choice and define a helper function, whose purpose is to make a HTTP GET request to a specified URL ("endpoint"), and verify that the response is successful.
#---------------------------------------------------------------------------------------------
import requests

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"f0f388f8e939bf305fa9346ce99cecd1"}

def get(path, params=None):
    # make HTTP GET request to path
    headers = {"api-key":"f0f388f8e939bf305fa9346ce99cecd1"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r
'''   
def fit_kinematic_pa(x, y, vel, debug=False, nsteps=361,
                     quiet=False, plot=True, dvel=10):

    x, y, vel = map(np.ravel, [x, y, vel])

    assert x.size == y.size == vel.size, 'Input vectors (x, y, vel) must have the same size'

    nbins = x.size
    n = nsteps
    angles = np.linspace(0, 180, n) # 0.5 degrees steps by default
    chi2 = np.empty_like(angles)
    for j, ang in enumerate(angles):
        velSym = symmetrize_velfield(x, y, vel, sym=1, pa=ang)
        chi2[j] = np.sum(((vel - velSym)/dvel)**2)
        if debug:
            print('Ang: %5.1f, chi2/DOF: %.4g' % (ang, chi2[j]/nbins))
            plt.cla()
            plot_velfield(x, y, velSym)
            plt.pause(0.01)
    k = np.argmin(chi2)
    angBest = angles[k]

    # Compute fit at the best position
    #
    velSym = symmetrize_velfield(x, y, vel, sym=1, pa=angBest)
    if angBest < 0:
        angBest += 180

    # 3sigma confidence limit, including error on chi^2
    #
    f = chi2 - chi2[k] <= 9 + 3*np.sqrt(2*nbins)
    minErr = max(0.5, (angles[1] - angles[0])/2.0)
    if f.sum() > 1:
        angErr = (np.max(angles[f]) - np.min(angles[f]))/2.0
        if angErr >= 45:
            good = np.degrees(np.arctan(np.tan(np.radians(angles[f]))))
            angErr = (np.max(good) - np.min(good))/2.0
        angErr = max(angErr, minErr)
    else:
        angErr = minErr

    vSyst = np.median(vel - velSym)

    if not quiet:
        print('  Kin PA: %5.1f' % angBest, ' +/- %5.1f' % angErr, ' (3*sigma error)')
        print('Velocity Offset: %.2f' % vSyst)

    # Plot results
    #
    if plot:

        mn, mx = np.percentile(velSym, [2.5, 97.5])
        mx = min(mx, -mn)
        plt.subplot(121)
        plot_velfield(x, y, velSym, vmin=-mx, vmax=mx)
        plt.title('Symmetrized')

        plt.subplot(122)
        plot_velfield(x, y, vel - vSyst, vmin=-mx, vmax=mx)
        plt.title('Data and best PA')
        rad = np.sqrt(np.max(x**2 + y**2))
        ang = [0,np.pi] + np.radians(angBest)
        plt.plot(rad*np.cos(ang), rad*np.sin(ang), 'k--', linewidth=3) # Zero-velocity line
        plt.plot(-rad*np.sin(ang), rad*np.cos(ang), color="limegreen", linewidth=3) # Major axis PA

    return angBest, angErr, vSyst
'''
# cmap_reversed = matplotlib.cm.get_cmap('jet')
#2Dplot (star):   ( his2d + contour)
#---------------------------------------------------------------------------------------------
ids   = [468590]
# level = [0.1 , 0.2 , 0.5 , 0.75 , 0.95]
level = [5,10,25,50,70,100]

i = 0.0

redshift = 0.0
scale_factor = 1.0 / (1+redshift)
numRows = 1 # increase to trace back further
# plt.figure(figsize=[7*3,numRows*3])
plt.figure(figsize= (6,5))
count = 1
'''
for id in ids:
    baseUrl = "http://www.tng-project.org/api/TNG100-1/"
    start_url = baseUrl + "snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
    sub = get(start_url)

    # get cutout of gas positions and masses
    # we downloaded cutout of all ids
    cutout_request = {'gas':'Coordinates,Masses,Velocities','stars':'Coordinates,Masses,Velocities'}

    cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
    # cutout = "cutout_1256402.hdf5"
    # make 2d histogram visualization of gas distribution
    with h5py.File(cutout,'r') as f:
    # with h5py.File('cutout_'+str(id) + '_' + str(i)+'.hdf5','r') as f:

        dx_gas          = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
        dy_gas          = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
        dz_gas          = f['PartType0']['Coordinates'][:,2] - sub['pos_z']

        dvx_gas          = f['PartType0']['Velocities'][:,0] - sub['vel_x']
        dvy_gas          = f['PartType0']['Velocities'][:,1] - sub['vel_y']
        dvz_gas          = f['PartType0']['Velocities'][:,2] - sub['vel_z']
        mass_gas        = np.log10(f['PartType0']['Masses'][:] *1e10/0.704)

        dx_star          = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
        dy_star          = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
        dz_star          = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
        dvx_star         = f['PartType4']['Velocities'][:,0] - sub['vel_x']
        dvy_star         = f['PartType4']['Velocities'][:,1] - sub['vel_y']
        dvz_star         = f['PartType4']['Velocities'][:,2] - sub['vel_z']
        mass_star        = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)

        X   = dz_gas
        Y   = dy_gas
        Z   = mass_gas

    xrange = [-200,200]
    yrange = [-200,200]
    # add to plot
    # plt.subplot(numRows,5,count)

    counts_star,ybins_s,xbins_s,image_s = plt.hist2d(dz_star,dy_star,bins=[100,100],norm=LogNorm(),cmin = 1, range=[xrange,yrange])
    plt.clf()
    counts,ybins,xbins,image = plt.hist2d(X,Y,bins=[40,40],norm=LogNorm(),cmin = 5,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
    # plt.scatter(X,Y, s= 1, c = Z,cmap = 'jet')
    cb = plt.colorbar()
    # cb.ax.set_yticklabels(['4','5','6','7','8','9','10'])
    # cb.set_label(r'$\log_{10}~(\mathrm{M_{BH}}\ [M_{\odot}])$', rotation=270,labelpad=20)

    plt.contour(np.transpose(counts_star),extent=[xbins_s.min(),xbins_s.max(),ybins_s.min(),ybins_s.max()],linewidths=1,colors='black',levels = level)
    #plt.clabel(CS, inline=1, fontsize=7, fmt = '%1.1f')


    plt.xlabel('$\Delta z$ [ckpc/h]')
    plt.ylabel('$\Delta y$ [ckpc/h]')
    #plt.text(-2,-2+5,"snap="+str(sub['snap'])+'  /  id='+str(id), color = 'white')
    plt.xlim(-200, 200)
    plt.ylim(-200, 200)
    #plt.gca().axes.get_xaxis().set_ticks([])
    #plt.gca().axes.get_yaxis().set_ticks([])
    count += 1

    plt.savefig('P7_2Dplot_star_on_gas.png')

    if count > numRows*5:
        break

plt.close()


############################################
#Trace ID Illustris-1
output_I1 = "./plots/IL-1/"
id   = 104798 # non fossil 1256402
snap = 135
#baseUrl = "http://www.tng-project.org/api/TNG100-1/"
baseUrl = "http://www.illustris-project.org/api/Illustris-1/"
start_url = baseUrl + "snapshots/" + str(snap) + "/subhalos/" + str(id)
sub = get(start_url)

numRows = 4 # increase to trace back further
# plt.figure(figsize=[5*3,numRows*3])
count = 1
#list = [99, 91,84,78,72,67,59,50]

#for sub['snap'] in list:
while sub['prog_sfid'] >= -1:
# while sub['snap'] >= 135 :
    # request the full subhalo details of the progenitor by following the sublink URL
    sub = get(sub['related']['sublink_progenitor'])

    print sub['id'], sub['snap']

    # cutout_request = {'gas':'Coordinates,Masses,Velocities','stars':'Coordinates,Masses,Velocities'}
    # cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)

    # cutout = "cutout_1256402.hdf5"
    cutout = "./IL-1_cutout_104798/cutout_"+str(sub['id'])+".hdf5"
    # make 2d histogram visualization of gas distribution
    with h5py.File(cutout,'r') as f:
    # with h5py.File('cutout_'+str(id) + '_' + str(i)+'.hdf5','r') as f:

        dx_gas          = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
        dy_gas          = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
        dz_gas          = f['PartType0']['Coordinates'][:,2] - sub['pos_z']

        dvx_gas          = f['PartType0']['Velocities'][:,0] - sub['vel_x']
        dvy_gas          = f['PartType0']['Velocities'][:,1] - sub['vel_y']
        dvz_gas          = f['PartType0']['Velocities'][:,2] - sub['vel_z']
        mass_gas        = np.log10(f['PartType0']['Masses'][:] *1e10/0.704)

        dx_star          = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
        dy_star          = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
        dz_star          = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
        dvx_star         = f['PartType4']['Velocities'][:,0] - sub['vel_x']
        dvy_star         = f['PartType4']['Velocities'][:,1] - sub['vel_y']
        dvz_star         = f['PartType4']['Velocities'][:,2] - sub['vel_z']
        mass_star        = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)

        X   = dz_gas
        Y   = dy_gas
        Z   = mass_gas

    # add to plot
    # plt.subplot(numRows,5,count)
    xrange = [-200,200]
    yrange = [-200,200]

    counts_star,ybins_s,xbins_s,image_s = plt.hist2d(dz_star,dy_star,bins=[100,100],norm=LogNorm(),cmin = 1, range=[xrange,yrange])
    plt.clf()
    counts,ybins,xbins,image = plt.hist2d(X,Y,bins=[40,40],norm=LogNorm(),cmin = 2,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
    # plt.scatter(X,Y, s= 1, c = Z,cmap = 'jet')
    cb = plt.colorbar()
    # cb.ax.set_yticklabels(['4','5','6','7','8','9','10'])
    # cb.set_label(r'$\log_{10}~(\mathrm{M_{BH}}\ [M_{\odot}])$', rotation=270,labelpad=20)

    plt.contour(np.transpose(counts_star),extent=[xbins_s.min(),xbins_s.max(),ybins_s.min(),ybins_s.max()],linewidths=2,colors='black',levels = level)
    #plt.clabel(CS, inline=1, fontsize=7, fmt = '%1.1f')

    print 'finish'
    plt.xlabel('$\Delta z$ [ckpc/h]')
    plt.ylabel('$\Delta y$ [ckpc/h]')
    #plt.text(-2,-2+5,"snap="+str(sub['snap'])+'  /  id='+str(id), color = 'white')
    plt.text(-180,180,"snap="+str(sub['snap'])+'(id='+str(id)+')', color = 'red')

    plt.xlim(-200, 200)
    plt.ylim(-200, 200)
    #plt.gca().axes.get_xaxis().set_ticks([])
    #plt.gca().axes.get_yaxis().set_ticks([])
    count += 1
    plt.savefig(output_I1 + 'Illustris_gas-star' + str(sub['snap']) +'-'+ str(id) +'.pdf')
    plt.close()
    if count > numRows*5:
        break
# plt.savefig('Illustris_gas-star.pdf')
plt.close()
'''

############################################
#Trace ID TNG-300
output_TNG300 = "./plots/TNG300/"
id   =  1256402  #468590  sample # TNG 100  288014
snap = 99
baseUrl = "http://www.tng-project.org/api/TNG300-1/"
# baseUrl = "http://www.illustris-project.org/api/Illustris-1/"
start_url = baseUrl + "snapshots/" + str(snap) + "/subhalos/" + str(id)
sub = get(start_url)
# cutout_request = {'gas':'Coordinates,Masses,Velocities','stars':'Coordinates,Masses,Velocities'}
cutout_request = {'gas':'Coordinates,Masses,Velocities,InternalEnergy,SubfindVelDisp,ElectronAbundance,GFM_Metallicity,Density','stars':'Coordinates,Masses,Velocities'}
cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)

# cutout = "./TNG300_cutout_1256402/cutout_"+str(sub['id'])+".hdf5"
print sub['id'], sub['snap']
# For the last snapshot
# make 2d histogram visualization of gas distribution
with h5py.File(cutout,'r') as f:

            dx_gas    = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
            dy_gas    = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
            dz_gas    = f['PartType0']['Coordinates'][:,2] - sub['pos_z']
            dvx_gas   = f['PartType0']['Velocities'][:,0] #- sub['vel_x']
            dvy_gas   = f['PartType0']['Velocities'][:,1] #- sub['vel_y']
            dvz_gas   = f['PartType0']['Velocities'][:,2] #- sub['vel_z']
            mass_gas  = np.log10(f['PartType0']['Masses'][:] *1e10/0.704)

            dx_star   = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
            dy_star   = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
            dz_star   = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
            dvx_star  = f['PartType4']['Velocities'][:,0]  #- sub['vel_x']
            dvy_star  = f['PartType4']['Velocities'][:,1]  #- sub['vel_y']
            dvz_star  = f['PartType4']['Velocities'][:,2]  #- sub['vel_z']
            mass_star = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)

            mgas      = f['PartType0']['Masses'][:] *1e10
            u         = f['PartType0']['InternalEnergy'][:]    #  the Internal Energy
            #VelDisp   = f['PartType0']['SubfindVelDisp'][:]
            Xe        = f['PartType0']['ElectronAbundance'][:]  # xe (=ne/nH)  the electron abundance
            XH        = 0.76             # the hydrogen mass fraction
            gama      = 5.0/3.0          # the adiabatic index
            KB        = 1.3807e-16       # the Boltzmann constant in CGS units  [cm^2 g s^-2 K^-1]
            KB_kev    = 8.6173324e-8 
            mp        = 1.6726e-24       # the proton mass  [g]
            mu        = (4*mp)/(1+3*XH+4*XH*Xe)

            temperature = (gama-1)* (u/KB)* mu* 1e10

            #gas_profile = temperature
            Msun      = 1.989e+33
            kpc_in_cm = 3.08568e+21
            radius    = np.sqrt(np.power(dx_gas-0,2) +  np.power(dy_gas-0,2) +  np.power(dz_gas-0,2))
            radius_s  = np.sqrt(dx_star**2 + dy_star**2 + dz_star**2)
            rsize     = sub['halfmassrad_stars']


            #cold& hot gas
            #mgas_c    = np.sum(mgas[((radius <= rsize)& (temperature < 1e6))])
            #mgas_h    = np.sum(mgas[((radius <= rsize)& (temperature > 1e6))])


            #Cold gas accretion
            mask_c    = ((radius <= 2*rsize) ) #& (temperature < 1e6))
            vx_gas    = dvx_gas[mask_c] - sub['vel_x']
            vy_gas    = dvy_gas[mask_c] - sub['vel_y']
            vz_gas    = dvz_gas[mask_c] - sub['vel_z']
            #rx_gas    = f['PartType0']['Coordinates'][mask_c,0] - sub['pos_x']
            #ry_gas    = f['PartType0']['Coordinates'][mask_c,1] - sub['pos_y']
            #rz_gas    = f['PartType0']['Coordinates'][mask_c,2] - sub['pos_z']
            rx_gas    = dx_gas [mask_c]
            ry_gas    = dy_gas [mask_c]
            rz_gas    = dz_gas [mask_c]
            #Macc_c    = np.sum((mgas[mask_c]/rsize) * (vx_gas*rx +vy_gas*ry + vz_gas*rz)/np.sqrt((rx)**2 + (ry)**2 + (rz)**2)) # Msun  /kpc(km 3.086e16) km/ s(yr 3.17e-8) : 1e-9 Msun/yr 


            #stars
            limstar      = np.where( radius_s <= 2*rsize )
            vx_star      = dvx_star[limstar] - sub['vel_x']
            vy_star      = dvy_star[limstar] - sub['vel_y']
            vz_star      = dvz_star[limstar] - sub['vel_z']
            rx_star      = dx_star [limstar]
            ry_star      = dy_star [limstar]
            rz_star      = dz_star [limstar] 

xrange = [-200,200]
yrange = [-200,200]

# counts_star,ybins_s,xbins_s,image_s = plt.hist2d(ry_star,rz_star,bins=[100,100],norm=LogNorm(),cmin = 1, range=[xrange,yrange])
# plt.clf()
# counts,ybins,xbins,image = plt.hist2d(dx_gas,dy_gas ,bins=[100,100],norm=LogNorm(),cmin = 1,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
# plt.clf()
# plt.imshow(counts,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],cmap = 'jet_r')
#x,y,z= np.meshgrid(X,Y,Z,indexing = 'ij')
#z = np.array([dvz_star[i][j] for i in X for J in Y])

#z =  dvz_star.reshape(len(X),len(Y))

#plt.pcolormesh(x,y,z)
# plt.scatter(dx_gas,dy_gas, s= 300, c =dvz_gas ,cmap = 'jet_r', edgecolor = '', alpha = 0.2)
# plt.hexbin(dy_star,dz_star, C=dvx_star, gridsize=100, cmap='jet_r', bins=None,edgecolors=None,mincnt = 0.1)
# cb = plt.colorbar()
angBest, angErr, vSyst = fit_kinematic_pa(ry_star,rz_star, vx_star, debug=False, nsteps=361,quiet=False, plot=True, dvel=10)
print (angBest, angErr, vSyst)

# cb.ax.set_yticklabels(['4','5','6','7','8','9','10'])
# cb.set_label(r'$\log_{10}~(\mathrm{M_{BH}}\ [M_{\odot}])$', rotation=270,labelpad=20)

# plt.contour(np.transpose(counts_star),extent=[xbins_s.min(),xbins_s.max(),ybins_s.min(),ybins_s.max()],linewidths=2,colors='black',levels = level)
#plt.clabel(CS, inline=1, fontsize=7, fmt = '%1.1f')

print 'finish'
plt.xlabel('$\Delta y$ [ckpc/h]')
plt.ylabel('$\Delta z$ [ckpc/h]')
plt.text(-180,180,"snap="+str(sub['snap'])+'(id='+str(id)+')', color = 'red')
plt.xlim(-150, 150)
plt.ylim(-150, 150)
#plt.gca().axes.get_xaxis().set_ticks([])
#plt.gca().axes.get_yaxis().set_ticks([])
count += 1
plt.savefig(output_TNG300 + 'Illustris_gas-star' + str(sub['snap']) +'-'+ str(id) +'.pdf')
plt.close()

'''
numRows = 4 # increase to trace back further
# plt.figure(figsize=[5*3,numRows*3])
count = 1
#list = [99-79,72,67,59,50]

#for trace the evolution
while sub['prog_sfid'] != -1:
# while sub['snap'] >= 50 :
    # request the full subhalo details of the progenitor by following the sublink URL
    sub = get(sub['related']['sublink_progenitor'])
    if ((sub['snap'] >= 79)|((sub['snap'] == 72)|(sub['snap'] == 67)|(sub['snap'] == 59)|(sub['snap'] == 50))):
    # if ((sub['snap'] == 98)|(sub['snap'] == 91)|(sub['snap'] == 84)|(sub['snap'] == 72)|(sub['snap'] == 67)|(sub['snap'] == 59)|(sub['snap'] == 50)):
     print sub['id'], sub['snap']
     # cutout_request = {'gas':'Coordinates,Masses,Velocities','stars':'Coordinates,Masses,Velocities'}
     # cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
     cutout = "./TNG300_cutout_1256402/cutout_"+str(sub['id'])+".hdf5"

    # make 2d histogram visualization of gas distribution
     with h5py.File(cutout,'r') as f:

        dx_gas          = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
        dy_gas          = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
        dz_gas          = f['PartType0']['Coordinates'][:,2] - sub['pos_z']

        dvx_gas          = f['PartType0']['Velocities'][:,0] - sub['vel_x']
        dvy_gas          = f['PartType0']['Velocities'][:,1] - sub['vel_y']
        dvz_gas          = f['PartType0']['Velocities'][:,2] - sub['vel_z']
        mass_gas        = np.log10(f['PartType0']['Masses'][:] *1e10/0.704)

        dx_star          = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
        dy_star          = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
        dz_star          = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
        dvx_star         = f['PartType4']['Velocities'][:,0] - sub['vel_x']
        dvy_star         = f['PartType4']['Velocities'][:,1] - sub['vel_y']
        dvz_star         = f['PartType4']['Velocities'][:,2] - sub['vel_z']
        mass_star        = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)

        X   = dx_gas
        Y   = dy_gas
        Z   = dz_gas

    # add to plot
     xrange = [-200,200]
     yrange = [-200,200]

     counts_star,ybins_s,xbins_s,image_s = plt.hist2d(dy_star,dz_star,bins=[100,100],norm=LogNorm(),cmin = 1, range=[xrange,yrange])
     plt.clf()
     # counts,ybins,xbins,image = plt.hist2d(Z,Y,bins=[40,40],norm=LogNorm(),cmin = 2,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
     # plt.scatter(X,Y, s= 1, c = Z,cmap = 'jet')

     plt.hexbin(dy_gas,dz_gas, C=dvx_gas, gridsize=100, cmap='jet_r', bins=None,edgecolors=None,mincnt = 0.1)

     cb = plt.colorbar()
     # cb.ax.set_yticklabels(['4','5','6','7','8','9','10'])
     # cb.set_label(r'$\log_{10}~(\mathrm{M_{BH}}\ [M_{\odot}])$', rotation=270,labelpad=20)

     plt.contour(np.transpose(counts_star),extent=[xbins_s.min(),xbins_s.max(),ybins_s.min(),ybins_s.max()],linewidths=2,colors='black',levels = level)
     #plt.clabel(CS, inline=1, fontsize=7, fmt = '%1.1f')

     print 'finish'
     plt.xlabel('$\Delta z$ [ckpc/h]')
     plt.ylabel('$\Delta y$ [ckpc/h]')
     plt.text(-180,180,"snap="+str(sub['snap'])+'(id='+str(id)+')', color = 'red')
     plt.xlim(-200, 200)
     plt.ylim(-200, 200)
     #plt.gca().axes.get_xaxis().set_ticks([])
     #plt.gca().axes.get_yaxis().set_ticks([])
     count += 1
     plt.savefig(output_TNG300 + 'Illustris_gas-star' + str(sub['snap']) +'-'+ str(id) +'.pdf')
     plt.close()

'''
