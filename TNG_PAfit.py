#Estimate the kinematic position angle for gas, stars and dark matter at 2Re Illustris-TNG simualtion
#----------------------------------------------------------------------------------------------------
#--------------------------------------------mojtabaraouf@gmail.com----------------------------------
#TNG project-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from   pafit.fit_kinematic_pa import fit_kinematic_pa
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

plt.figure(figsize= (6,5))
TNG_dm_mass = 5.9 * 1e7  # TNG 300  (5.9 * 1e7) & TNG100(7.5 * 1e6) & TNG 50 (4.5 * 1e5)
temp_lim = 1e6
#Trace ID TNG-300 

output = './output_data/Gals_TNG300_Pafit.asc'
plot = './plots/'
id   =  1256402    # TNG 100  288014
snap = 99
baseUrl = "http://www.tng-project.org/api/TNG300-1/"
# baseUrl = "http://www.illustris-project.org/api/Illustris-1/"
start_url = baseUrl + "snapshots/" + str(snap) + "/subhalos/" + str(id)
sub = get(start_url)

# redshift = [0, 0.1, ....]
sub_halo = get(sub['related']['parent_halo'])
sub_sub_halo = get(sub_halo['meta']['info'])

#print 'Rvir , Mahalo = ',sub_sub_halo['Group_R_Crit200'], np.log10(sub_sub_halo['Group_M_Crit200']*1e10)
#----------------------------------------------------------------------------------------------
#Get the postion of galaxy in group in unit of r200
hc_Posx =  sub_sub_halo['GroupPos'][0]
hc_Posy =  sub_sub_halo['GroupPos'][1]
hc_Posz =  sub_sub_halo['GroupPos'][2]

gal_dist_center = (np.sqrt((hc_Posx - sub['pos_x'])**2 + (hc_Posy - sub['pos_y'])**2 + (hc_Posz - sub['pos_z'])**2))
hcr200 = sub_sub_halo['Group_R_Crit200']
#----------------------------------------------------------------------------------------------
#define central and satelite
sub_meta = get(sub['meta']['info'])
parentID = sub_meta['SubhaloParent']
if (parentID == id):
    cent = 1
else: 
    cent = 0
#----------------------------------------------------------------------------------------------
#redshift
sub_snap = get(sub['meta']['snapshot'])
redshift1  = sub_snap['redshift']

#----------------------------------------------------------------------------------------------
#start to get cutout    
count = 0
#for last snap 
cutout_request = {'dm':'Coordinates,Velocities' ,'gas':'Coordinates,Masses,Velocities,InternalEnergy,SubfindVelDisp,ElectronAbundance,GFM_Metallicity,Density','stars':'Coordinates,Masses,Velocities,GFM_StellarFormationTime'}

cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
#cutout = "./cutout_ramp/cutout_"+str(sub['id'])+".hdf5"
print sub['id'], sub['snap']
profile_frest = './output_data/fresor_TNG300_'+str(sub['snap']).zfill(3)+'.asc'
ff = open(profile_frest, 'w+')
ff.write('# galID r df_ramp frestor r_r200 M200'+"\n")
#print '#id  sub_id snap ramP  frest  massr  rsize  densr  dens velwrtgas gmass_cold gmass_hot gacc_cold gacc_hot' 

main_id=[];sub_ID=[];snap=[];r_size=[];hc_dist_center=[];hc_r200 =[];central=[];redshift = [];age = []
PA_gas =[];PA_gas_err =[];v_gas_sys =[];PA_s =[];PA_s_err =[];v_s_sys =[];PA_dm =[];PA_dm_err =[];v_dm_sys =[];deltaPA_GS =[];deltaPA_DS =[];deltaPA_DG =[]
    # make 2d histogram visualization of gas distribution
with h5py.File(cutout,'r') as f:
        dx_DM          = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
        dy_DM          = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
        dz_DM          = f['PartType1']['Coordinates'][:,2] - sub['pos_z']
        dvx_DM          = f['PartType1']['Velocities'][:,0] - sub['vel_x']
        dvy_DM          = f['PartType1']['Velocities'][:,1] - sub['vel_y']
        dvz_DM          = f['PartType1']['Velocities'][:,2] - sub['vel_z']

        dx_gas          = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
        dy_gas          = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
        dz_gas          = f['PartType0']['Coordinates'][:,2] - sub['pos_z']

        dvx_gas          = f['PartType0']['Velocities'][:,0] #- sub['vel_x']
        dvy_gas          = f['PartType0']['Velocities'][:,1] #- sub['vel_y']
        dvz_gas          = f['PartType0']['Velocities'][:,2] #- sub['vel_z']

        dx_star          = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
        dy_star          = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
        dz_star          = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
        dvx_star         = f['PartType4']['Velocities'][:,0] - sub['vel_x']
        dvy_star         = f['PartType4']['Velocities'][:,1] - sub['vel_y']
        dvz_star         = f['PartType4']['Velocities'][:,2] - sub['vel_z']
        mass_star        = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)
                
        mgas        = f['PartType0']['Masses'][:] *1e10
        u           = f['PartType0']['InternalEnergy'][:]    #  the Internal Energy
        VelDisp     = f['PartType0']['SubfindVelDisp'][:]
        Xe          = f['PartType0']['ElectronAbundance'][:]  # xe (=ne/nH)  the electron abundance
        XH          = 0.76             # the hydrogen mass fraction
        gama        = 5.0/3.0          # the adiabatic index
        KB          = 1.3807e-16       # the Boltzmann constant in CGS units  [cm^2 g s^-2 K^-1]
        KB_kev      = 8.6173324e-8 
        mp          = 1.6726e-24       # the proton mass  [g]
        mu          = (4*mp)/(1+3*XH+4*XH*Xe)
        Msun = 1.989e+33
        kpc_in_cm=3.08568e+21
        
        
        temperature = (gama-1)* (u/KB)* mu* 1e10 
        #CoolingRate         = f['PartType0']['GFM_CoolingRate'][:]
        Metallicity         = f['PartType0']['GFM_Metallicity'][:] / 0.0127
        Density             = f['PartType0']['Density'][:] # (10^10M/h)/(ckpc/h)3
        rsize = sub['halfmassrad_stars']

        radius_s      = np.sqrt(np.power(dx_star,2) +  np.power(dy_star,2) +  np.power(dz_star,2))
        radius_g     = np.sqrt(np.power(dx_gas,2) +  np.power(dy_gas,2) +  np.power(dz_gas,2))
        radius_dm   = np.sqrt(np.power(dx_DM,2) +  np.power(dy_DM,2) +  np.power(dz_DM,2))
        
       
#----------------------------------------------------------------------------------------------              
        #star formation history test
        SFHp = f['PartType4']['GFM_StellarFormationTime'][:]
        SFH = np.mean(SFHp[(radius_s <= rsize)])
        age1 = cosmo.age((1/SFH)+1).value

#----------------------------------------------------------------------------------------------
#PA angle for gas, star and DM
        #Cold gas 
        mask_g = ((radius_g<= 2.* rsize)) 
        vx_gas = dvx_gas[mask_g] - sub['vel_x']
        angBest_g, angErr_g, vSyst_g = fit_kinematic_pa(dy_gas[mask_g],dz_gas[mask_g], vx_gas, debug=False, nsteps=361,quiet=False, plot=False, dvel=10) 
        #stars 
        mask_s = ((radius_s<= 2.* rsize)) ; 
        angBest_s, angErr_s, vSyst_s = fit_kinematic_pa(dy_star[mask_s],dz_star[mask_s], dvx_star[mask_s], debug=False, nsteps=361,quiet=False, plot=False, dvel=10) 
        #DM 
        mask_dm = ((radius_dm<= 2.* rsize))
        angBest_dm, angErr_dm, vSyst_dm = fit_kinematic_pa(dy_DM[mask_dm],dz_DM[mask_dm], dvx_DM[mask_dm], debug=False, nsteps=361,quiet=False, plot=False, dvel=10)      
#----------------------------------------------------------------------------------------------

        print  'PA_g=',angBest_g,'PA_s=',angBest_s,'PA_dm=',angBest_dm,
        main_id.append(id)
        sub_ID.append(sub['id'])
        snap.append(sub['snap'])
        r_size.append(rsize)
        hc_dist_center.append(gal_dist_center)
        hc_r200.append(hcr200)
        central.append(cent)
        redshift.append(redshift1)
        age.append(age1)
        PA_gas.append(angBest_g); PA_gas_err.append(angErr_g); v_gas_sys.append(vSyst_g)
        PA_s.append(angBest_s); PA_s_err.append(angErr_s); v_s_sys.append(vSyst_s)
        PA_dm.append(angBest_dm); PA_dm_err.append(angErr_dm); v_dm_sys.append(vSyst_dm)
        deltaPA_GS.append(abs(angBest_g - angBest_s))
        deltaPA_DS.append(abs(angBest_dm - angBest_s))
        deltaPA_DG.append(abs(angBest_dm - angBest_g))     
        
 
#------------------------------------------------------------------------------------------------------------------------
#Plot Kinematics of galaxies
fig, axs = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=False , sharex=True, sharey=True)
        
plt.axis([-2*rsize-6,2*rsize+6, -2*rsize-6,2*rsize+7])

# star    
# add circle around virial radius
circle1 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
axs[0].add_artist(circle1)
circle2 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
axs[0].add_artist(circle2)

# plot velocity map
vel = axs[0].hexbin(dy_star[mask_s],dz_star[mask_s], C=dvx_star[mask_s], gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
#cb = fig.colorbar()
fig.colorbar( vel, ax=axs[0])

max_bin = 20
tet = np.linspace(np.pi/2.,(angBest_s)*np.pi/180.,max_bin)
x0  = 0 ; y0 = 0
x1  = 20 * np.cos(tet) + 0.0
y1  = 20 * np.sin(tet) + 0.0
#axs[0].plot(x1,y1, color = 'grey')
x2  = np.array([x0,x1[max_bin-1]])
y2  = np.array([y0,y1[max_bin-1]])
z   = np.polyfit(x2,y2,1)
x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
y3  = z[0] * x3 + z[1]
#axs[1].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
axs[0].plot(x3,y3,'b--', alpha = 0.99)

axs[0].set_title('SK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
axs[0].set_xlabel('$\Delta y$ [ckpc/h]')
axs[0].set_ylabel('$\Delta z$ [ckpc/h]')
axs[0].text( -2*rsize-5,-2*rsize-5  , 'PA_SK = '          + str(angBest_s) + (u"\u00B1") + str(angErr_s) )
#------------------------------------------------------------------------------------------------------------------------    
# gas 
# add circle around virial radius
circle3 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
axs[1].add_artist(circle3)
    
circle4 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
axs[1].add_artist(circle4)

# plot velocity map
vel = axs[1].hexbin(dy_gas[mask_g],dz_gas[mask_g], C=vx_gas, gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
#cb = fig.colorbar()
fig.colorbar( vel, ax=axs[1])

max_bin = 20
tet = np.linspace(np.pi/2.,(angBest_g)*np.pi/180.,max_bin)
x0  = 0 ; y0 = 0
x1  = 20 * np.cos(tet) + 0.0
y1  = 20 * np.sin(tet) + 0.0
#axs[0].plot(x1,y1, color = 'grey')
x2  = np.array([x0,x1[max_bin-1]])
y2  = np.array([y0,y1[max_bin-1]])
z   = np.polyfit(x2,y2,1)
x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
y3  = z[0] * x3 + z[1]
axs[1].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
#axs[0].plot(x3,y3,'r--', color = 'red', alpha = 0.99)

axs[1].set_title('GK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
axs[1].set_xlabel('$\Delta y$ [ckpc/h]')
#axs[1].set_ylabel('$\Delta z$ [ckpc/h]')
axs[1].text( -2*rsize-5,-2*rsize-5  , 'PA_GK = '          + str(angBest_g) + (u"\u00B1") + str(angErr_g) )

#------------------------------------------------------------------------------------------------------------------------
#DM
# add circle around virial radius
circle3 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
axs[2].add_artist(circle3)
    
circle4 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
axs[2].add_artist(circle4)

# plot velocity map
vel = axs[2].hexbin(dy_DM[mask_dm],dz_DM[mask_dm], C=dvx_DM[mask_dm], gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
#cb = fig.colorbar()
fig.colorbar( vel, ax=axs[2])

max_bin = 20
tet = np.linspace(np.pi/2.,(angBest_dm)*np.pi/180.,max_bin)
x0  = 0 ; y0 = 0
x1  = 20 * np.cos(tet) + 0.0
y1  = 20 * np.sin(tet) + 0.0
#axs[0].plot(x1,y1, color = 'grey')
x2  = np.array([x0,x1[max_bin-1]])
y2  = np.array([y0,y1[max_bin-1]])
z   = np.polyfit(x2,y2,1)
x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
y3  = z[0] * x3 + z[1]
axs[2].plot(x3,y3,'g--', alpha = 0.99)
#axs[0].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
axs[2].set_title('DMK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
axs[2].set_xlabel('$\Delta y$ [ckpc/h]')
#axs[2].set_ylabel('$\Delta z$ [ckpc/h]')
axs[2].text( -2*rsize-5,-2*rsize-5  , 'PA_DMK = '          + str(angBest_dm) + (u"\u00B1") + str(angErr_dm) )
plt.subplots_adjust(wspace=0.06,hspace=0.6, bottom=0.13, right=0.97,left=0.05, top=0.92)

plt.savefig(plot +'Kinematics-'+str(id)+'-'+str(snap)+'-TNG100.pdf')        

#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#for other trace list:
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
ff.close()
# while sub['prog_sfid'] >= -1:
while sub['snap'] >= 93 :    
    # request the full subhalo details of the progenitor by following the sublink URL
    sub = get(sub['related']['sublink_progenitor'])

    if ((sub['snap'] >= 79)|((sub['snap'] == 72)|(sub['snap'] == 67)|(sub['snap'] == 59)|(sub['snap'] == 50))): 
     sub_halo = get(sub['related']['parent_halo'])
     sub_sub_halo = get(sub_halo['meta']['info'])
     print sub['id'], sub['snap']  
     profile_frest = './fresor_TNG300_'+str(sub['snap']).zfill(3)+'.asc'
     ff = open(profile_frest, 'w+')
     ff.write('# galID r df_ramp frestor r_r200 M200'+"\n")

#----------------------------------------------------------------------------------------------        
 #Get the postion of galaxy in group in unit of r200
     hc_Posx =  sub_sub_halo['GroupPos'][0]
     hc_Posy =  sub_sub_halo['GroupPos'][1]
     hc_Posz =  sub_sub_halo['GroupPos'][2]

     gal_dist_center = (np.sqrt((hc_Posx - sub['pos_x'])**2 + (hc_Posy - sub['pos_y'])**2 + (hc_Posz - sub['pos_z'])**2))
     hcr200 = sub_sub_halo['Group_R_Crit200']
#define central and satelite
     sub_meta = get(sub['meta']['info'])
     parentID = sub_meta['SubhaloParent']
     if (parentID == id):
         cent = 1
     else: 
         cent = 0
#----------------------------------------------------------------------------------------------
#redshift
     sub_snap = get(sub['meta']['snapshot'])
     redshift1  = sub_snap['redshift']
#----------------------------------------------------------------------------------------------
       
     cutout_request = {'dm':'Coordinates,Velocities' ,'gas':'Coordinates,Masses,Velocities,InternalEnergy,ElectronAbundance,GFM_Metallicity,Density','stars':'Coordinates,Masses,Velocities,GFM_StellarFormationTime'}
     cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
     #cutout = "./cutout_ramp/cutout_"+str(sub['id'])+".hdf5"

    # make 2d histogram visualization of gas distribution
     with h5py.File(cutout,'r') as f:
        dx_DM          = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
        dy_DM          = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
        dz_DM          = f['PartType1']['Coordinates'][:,2] - sub['pos_z']
        dvx_DM          = f['PartType1']['Velocities'][:,0] - sub['vel_x']
        dvy_DM          = f['PartType1']['Velocities'][:,1] - sub['vel_y']
        dvz_DM          = f['PartType1']['Velocities'][:,2] - sub['vel_z']

        dx_gas          = f['PartType0']['Coordinates'][:,0] - sub['pos_x']
        dy_gas          = f['PartType0']['Coordinates'][:,1] - sub['pos_y']
        dz_gas          = f['PartType0']['Coordinates'][:,2] - sub['pos_z']

        dvx_gas          = f['PartType0']['Velocities'][:,0] #- sub['vel_x']
        dvy_gas          = f['PartType0']['Velocities'][:,1] #- sub['vel_y']
        dvz_gas          = f['PartType0']['Velocities'][:,2] #- sub['vel_z']
        mass_gas        = np.log10(f['PartType0']['Masses'][:] *1e10/0.704)

        dx_star          = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
        dy_star          = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
        dz_star          = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
        dvx_star         = f['PartType4']['Velocities'][:,0] - sub['vel_x']
        dvy_star         = f['PartType4']['Velocities'][:,1] - sub['vel_y']
        dvz_star         = f['PartType4']['Velocities'][:,2] - sub['vel_z']
        mass_star        = np.log10(f['PartType4']['Masses'][:] *1e10/0.704)
                
        mgas        = f['PartType0']['Masses'][:] *1e10
        u           = f['PartType0']['InternalEnergy'][:]    #  the Internal Energy
        # VelDisp     = f['PartType0']['SubfindVelDisp'][:]
        Xe          = f['PartType0']['ElectronAbundance'][:]  # xe (=ne/nH)  the electron abundance
        XH          = 0.76             # the hydrogen mass fraction
        gama        = 5.0/3.0          # the adiabatic index
        KB          = 1.3807e-16       # the Boltzmann constant in CGS units  [cm^2 g s^-2 K^-1]
        KB_kev      = 8.6173324e-8 
        mp          = 1.6726e-24       # the proton mass  [g]
        mu          = (4*mp)/(1+3*XH+4*XH*Xe)
        Msun = 1.989e+33
        kpc_in_cm=3.08568e+21
       
        temperature = (gama-1)* (u/KB)* mu* 1e10 
        #CoolingRate         = f['PartType0']['GFM_CoolingRate'][:]
        Metallicity         = f['PartType0']['GFM_Metallicity'][:] / 0.0127
        Density             = f['PartType0']['Density'][:] # (10^10M/h)/(ckpc/h)3
        rsize = sub['halfmassrad_stars']
        radius_s    = np.sqrt(np.power(dx_star,2) +  np.power(dy_star,2) +  np.power(dz_star,2))
        radius_g     = np.sqrt(np.power(dx_gas,2) +  np.power(dy_gas,2) +  np.power(dz_gas,2))
        radius_dm   = np.sqrt(np.power(dx_DM,2) +  np.power(dy_DM,2) +  np.power(dz_DM,2))
        
       
#----------------------------------------------------------------------------------------------              
        #star formation history test
        SFHp = f['PartType4']['GFM_StellarFormationTime'][:]
        SFH = np.mean(SFHp[(radius_s <= rsize)])
        age1 = cosmo.age((1/SFH)+1).value

#----------------------------------------------------------------------------------------------
#PA angle for gas, star and DM
        #Cold gas 
        mask_g = ((radius_g<= 2.* rsize)) 
        vx_gas = dvx_gas[mask_g] - sub['vel_x']
        angBest_g, angErr_g, vSyst_g = fit_kinematic_pa(dy_gas[mask_g],dz_gas[mask_g], vx_gas, debug=False, nsteps=361,quiet=False, plot=False, dvel=10) 
        #stars 
        mask_s = ((radius_s<= 2.* rsize)) ; 
        angBest_s, angErr_s, vSyst_s = fit_kinematic_pa(dy_star[mask_s],dz_star[mask_s], dvx_star[mask_s], debug=False, nsteps=361,quiet=False, plot=False, dvel=10) 
        #DM 
        mask_dm = ((radius_dm<= 2.* rsize))
        angBest_dm, angErr_dm, vSyst_dm = fit_kinematic_pa(dy_DM[mask_dm],dz_DM[mask_dm], dvx_DM[mask_dm], debug=False, nsteps=361,quiet=False, plot=False, dvel=10)      
#----------------------------------------------------------------------------------------------

        print  'PA_g=',angBest_g,'PA_s=',angBest_s,'PA_dm=',angBest_dm,
        main_id.append(id)
        sub_ID.append(sub['id'])
        snap.append(sub['snap'])
        r_size.append(rsize)
        hc_dist_center.append(gal_dist_center)
        hc_r200.append(hcr200)
        central.append(cent)
        redshift.append(redshift1)
        age.append(age1)
        PA_gas.append(angBest_g); PA_gas_err.append(angErr_g); v_gas_sys.append(vSyst_g)
        PA_s.append(angBest_s); PA_s_err.append(angErr_s); v_s_sys.append(vSyst_s)
        PA_dm.append(angBest_dm); PA_dm_err.append(angErr_dm); v_dm_sys.append(vSyst_dm)
        deltaPA_GS.append(abs(angBest_g - angBest_s))
        deltaPA_DS.append(abs(angBest_dm - angBest_s))
        deltaPA_DG.append(abs(angBest_dm - angBest_g))                
        
#------------------------------------------------------------------------------------------------------------------------
#Plot Kinematics of galaxies
     fig, axs = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=False , sharex=True, sharey=True)
        
     plt.axis([-2*rsize-6,2*rsize+6, -2*rsize-6,2*rsize+7])

# star    
# add circle around virial radius
     circle1 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
     axs[0].add_artist(circle1)
     circle2 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
     axs[0].add_artist(circle2)
     
     # plot velocity map
     vel = axs[0].hexbin(dy_star[mask_s],dz_star[mask_s], C=dvx_star[mask_s], gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
     #cb = fig.colorbar()
     fig.colorbar( vel, ax=axs[0])

     max_bin = 20
     tet = np.linspace(np.pi/2.,(angBest_s)*np.pi/180.,max_bin)
     x0  = 0 ; y0 = 0
     x1  = 20 * np.cos(tet) + 0.0
     y1  = 20 * np.sin(tet) + 0.0
    #axs[0].plot(x1,y1, color = 'grey')
     x2  = np.array([x0,x1[max_bin-1]])
     y2  = np.array([y0,y1[max_bin-1]])
     z   = np.polyfit(x2,y2,1)
     x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
     y3  = z[0] * x3 + z[1]
     #axs[1].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
     axs[0].plot(x3,y3,'b--', alpha = 0.99)

     axs[0].set_title('SK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
     axs[0].set_xlabel('$\Delta y$ [ckpc/h]')
     axs[0].set_ylabel('$\Delta z$ [ckpc/h]')
     axs[0].text( -2*rsize-5,-2*rsize-5  , 'PA_SK = '          + str(angBest_s) + (u"\u00B1") + str(angErr_s) )
#------------------------------------------------------------------------------------------------------------------------    
# gas 
# add circle around virial radius
     circle3 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
     axs[1].add_artist(circle3)
    
     circle4 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
     axs[1].add_artist(circle4)

# plot velocity map
     vel = axs[1].hexbin(dy_gas[mask_g],dz_gas[mask_g], C=vx_gas, gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
     #cb = fig.colorbar()
     fig.colorbar( vel, ax=axs[1])

     max_bin = 20
     tet = np.linspace(np.pi/2.,(angBest_g)*np.pi/180.,max_bin)
     x0  = 0 ; y0 = 0
     x1  = 20 * np.cos(tet) + 0.0
     y1  = 20 * np.sin(tet) + 0.0
     #axs[0].plot(x1,y1, color = 'grey')
     x2  = np.array([x0,x1[max_bin-1]])
     y2  = np.array([y0,y1[max_bin-1]])
     z   = np.polyfit(x2,y2,1)
     x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
     y3  = z[0] * x3 + z[1]
     axs[1].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
     #axs[0].plot(x3,y3,'r--', color = 'red', alpha = 0.99)

     axs[1].set_title('GK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
     axs[1].set_xlabel('$\Delta y$ [ckpc/h]')
     #axs[1].set_ylabel('$\Delta z$ [ckpc/h]')
     axs[1].text( -2*rsize-5,-2*rsize-5  , 'PA_GK = '          + str(angBest_g) + (u"\u00B1") + str(angErr_g) )

#------------------------------------------------------------------------------------------------------------------------
#DM
# add circle around virial radius
     circle3 = plt.Circle((0 , 0 ), 2*rsize, color='black',ls = '--',  fill=False)
     axs[2].add_artist(circle3)
    
     circle4 = plt.Circle((0 , 0 ), rsize, color='black',  fill=False)
     axs[2].add_artist(circle4)

# plot velocity map
     vel = axs[2].hexbin(dy_DM[mask_dm],dz_DM[mask_dm], C=dvx_DM[mask_dm], gridsize=10, cmap='jet_r', bins=None, edgecolors=None, mincnt = 0.1)
     #cb = fig.colorbar()
     fig.colorbar( vel, ax=axs[2])

     max_bin = 20
     tet = np.linspace(np.pi/2.,(angBest_dm)*np.pi/180.,max_bin)
     x0  = 0 ; y0 = 0
     x1  = 20 * np.cos(tet) + 0.0
     y1  = 20 * np.sin(tet) + 0.0
     #axs[0].plot(x1,y1, color = 'grey')
     x2  = np.array([x0,x1[max_bin-1]])
     y2  = np.array([y0,y1[max_bin-1]])
     z   = np.polyfit(x2,y2,1)
     x3  = np.linspace( -rsize-8 ,rsize+8 , 100)
     y3  = z[0] * x3 + z[1]
     axs[2].plot(x3,y3,'g--', alpha = 0.99)
     #axs[0].plot(x3,y3,'r--', color = 'red', alpha = 0.99)
     axs[2].set_title('DMK ( ID = '+str(id)+', z ='+str(np.round(redshift[0])).zfill(3)+')')
     axs[2].set_xlabel('$\Delta y$ [ckpc/h]')
     #axs[2].set_ylabel('$\Delta z$ [ckpc/h]')
     axs[2].text( -2*rsize-5,-2*rsize-5  , 'PA_DMK = '          + str(angBest_dm) + (u"\u00B1") + str(angErr_dm) )
     plt.subplots_adjust(wspace=0.06,hspace=0.6, bottom=0.13, right=0.97,left=0.05, top=0.92)

     plt.savefig(plot +'Kinematics-'+str(id)+'-'+str(snap)+'-TNG100.pdf')   

cc = np.array([main_id,sub_ID,snap,central,r_size,hc_dist_center,hc_r200,redshift,age,PA_gas,PA_gas_err,v_gas_sys,PA_s,PA_s_err,v_s_sys,PA_dm,PA_dm_err,v_dm_sys,deltaPA_GS,deltaPA_DS,deltaPA_DG]).T
np.savetxt(output, cc, fmt='%i %i %i %i %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f' , delimiter='\t',header='ID\tSubID\tsnap\tcentral\tr_size\thc_dist_center\thc_r200\tredshift\tage\tPA_gas\tPA_gas_err\tv_gas_sys\tPA_s\tPA_s_err\tv_s_sys\tPA_dm\tPA_dm_err\tv_dm_sys\tdeltaPA_GS\tdeltaPA_DS\tdeltaPA_DG')    

#----------------------------------------------------------------------------------------------
#plots
plt.figure(figsize= (6,5))
plt.plot(redshift, deltaPA_GS,deltaPA_DS,deltaPA_DG, lw= 3, color = 'blue',label = 'gas-strar')
plt.plot(redshift, deltaPA_DS, lw= 3, color = 'blue',label = 'DM-strar')
plt.plot(redshift, deltaPA_DG, lw= 3, color = 'blue',label = 'gas-DM')
#plt.plot(redshift, f_rest, lw= 3, color = 'red',label = 'Theory')
plt.xlabel('z')
plt.ylabel(r'$Position Angle$')
plt.xlim(0.01, 0.5) 
plt.ylim(0, 180)
plt.xscale('log')
plt.savefig('TNG_PAfit_trace'+str(id) +'.pdf')
plt.close()


