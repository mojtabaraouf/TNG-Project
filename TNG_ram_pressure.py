#Estimate the Ram pressure striping, cold and hot gas accretion for Illustris-TNG simualtion
#----------------------------------------------------------------------------------------------
#@mojtabaraouf TNG project---------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
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
#Trace ID TNG-300 

output = './output_data/Gals_TNG300_ramp.asc'
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
#cutout_request = {'dm':'Coordinates' ,'gas':'Coordinates,Masses,Velocities,InternalEnergy,SubfindVelDisp,ElectronAbundance,GFM_Metallicity,Density','stars':'Coordinates,Masses,Velocities,GFM_StellarFormationTime'}

#cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
cutout = "./cutout_ramp/cutout_"+str(sub['id'])+".hdf5"
print sub['id'], sub['snap']
profile_frest = './output_data/fresor_TNG300_'+str(sub['snap']).zfill(3)+'.asc'
ff = open(profile_frest, 'w+')
ff.write('# galID r df_ramp frestor r_r200 M200'+"\n")
#print '#id  sub_id snap ramP  frest  massr  rsize  densr  dens velwrtgas gmass_cold gmass_hot gacc_cold gacc_hot' 

main_id=[];sub_ID=[];snap=[];ram_p=[];f_rest=[];gmassr=[];r_size=[];dens_r=[];densg=[];velwrtgas=[];gmass_cold=[];gmass_hot=[];gacc_cold=[];gacc_hot=[];hc_dist_center=[];hc_r200 =[];central=[];redshift = [];age = [];r_trunc = []
    # make 2d histogram visualization of gas distribution
with h5py.File(cutout,'r') as f:
        dx_DM          = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
        dy_DM          = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
        dz_DM          = f['PartType1']['Coordinates'][:,2] - sub['pos_z']
        #mass_DM        = f['PartType1']['Masses'][:] *1e10

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
        
        temperature = (gama-1)* (u/KB)* mu* 1e10 
        #AGNRadiation        = f['PartType0']['GFM_AGNRadiation'][:]
        #CoolingRate         = f['PartType0']['GFM_CoolingRate'][:]
        Metallicity         = f['PartType0']['GFM_Metallicity'][:] / 0.0127
        #WindDMVelDisp       = f['PartType0']['GFM_WindDMVelDisp'][:]
        #StarFormationRate   = f['PartType0']['StarFormationRate'][:]
        Density             = f['PartType0']['Density'][:] # (10^10M/h)/(ckpc/h)3
        #NeutralHydrogenAbundance = f['PartType0']['NeutralHydrogenAbundance'][:]
        rsize = sub['halfmassrad_stars']
        #gas_profile = temperature
        radius_s      = np.sqrt(np.power(dx_star,2) +  np.power(dy_star,2) +  np.power(dz_star,2))
        radius      = np.sqrt(np.power(dx_gas,2) +  np.power(dy_gas,2) +  np.power(dz_gas,2))
        Msun = 1.989e+33
        kpc_in_cm=3.08568e+21
        temp_lim = 1e6
#----------------------------------------------------------------------------------------------              
        #star formation history test
        SFHp = f['PartType4']['GFM_StellarFormationTime'][:]
        SFH = np.mean(SFHp[(radius_s <= rsize)])
        age1 = cosmo.age((1/SFH)+1).value
    
#----------------------------------------------------------------------------------------------              
       
        #cold& hot gas
        mgas_c = np.sum(mgas[((radius <= rsize)& (temperature < temp_lim))])
        mgas_h = np.sum(mgas[((radius <= rsize)& (temperature > temp_lim))]) 
        
 #----------------------------------------------------------------------------------------------       
        #ram pressure
        # massr = sub['massinhalfrad_gas'] * 1e10 # M/h
        massr = np.sum(mgas[(radius <= rsize)])
        densr = np.nanmedian(Density[(radius <= rsize)]* 1e10 * Msun/kpc_in_cm**3)   # (M/h)/(ckpc/h)3
        max_r       = round(max(radius))
        
        Grav = 4.3 * 1e-9 * 1000 # Km2 kpc Msun-1 s-2
        Pi = 3.14

        frest = (1e10 * (Pi/2.)*Grav) * (massr) * densr / (rsize)  # [g s-2 cm]  [g s-2 kpc-3]
        
        mask = ((radius <= 2.* rsize)&(radius > rsize))
        gasVelx = np.sum(mgas[mask]*dvx_gas[mask])/np.sum(mgas[mask])
        gasVely = np.sum(mgas[mask]*dvy_gas[mask])/np.sum(mgas[mask])
        gasVelz = np.sum(mgas[mask]*dvz_gas[mask])/np.sum(mgas[mask])
        vel_wrt_gas2 = (sub['vel_x'] - gasVelx)**2 + (sub['vel_y'] - gasVely)**2 + (sub['vel_z'] - gasVelz)**2
        dens = np.nanmedian(Density[mask]* 1e10 * Msun/kpc_in_cm**3)
        ramp = (dens * vel_wrt_gas2)
#----------------------------------------------------------------------------------------------
#save as array in file FOR TRUNCATION RADII 
        bin_radii = -0.3
        rx = dx_gas
        ry = dy_gas
        rz = dz_gas
        
        rx_dm = dx_DM
        ry_dm = dy_DM
        rz_dm = dz_DM
        Grav = 4.3 * 1e-9 * 1000 # Km2 kpc Msun-1 s-2
        Pi = 3.14
        rtrunc_all = [];gap = [];ramii = [];gid = [];m200 = [];rsize1 = []
        for i in np.arange(2.0*(rsize),1,bin_radii):

            if i < rsize:
                maski = ( (np.sqrt((rx)**2 + (ry)**2 + (rz)**2)<= i))
                # maski_s = ( (np.sqrt((rx_s)**2 + (ry_s)**2 + (rz_s)**2)<= i))
                maski_dm = ( (np.sqrt((rx_dm)**2 + (ry_dm)**2 + (rz_dm)**2)<= i))
            else:
                maski = ( (np.sqrt((rx)**2 + (ry)**2 + (rz)**2)<= rsize))
                maski_dm = ( (np.sqrt((rx_dm)**2 + (ry_dm)**2 + (rz_dm)**2)<= rsize))
                
            densi =  np.nanmedian(Density[maski]* 1e10 * Msun/kpc_in_cm**3)
            massdmi = len(rx_dm[maski_dm]) * TNG_dm_mass
            massgasi = np.sum(mgas[maski])
            massi = massgasi + massdmi
            rami = ((Pi/2.)*Grav) * (massi) * densi / (i)
            rtrunc_all.append(i)
            gap.append(abs(np.log10(ramp * 1e10) - np.log10(rami * 1e10)))
            ramii.append(rami)
            gid.append(id)
            m200.append(hcr200)
            rsize1.append(rsize)

        zz = np.array([gid,rtrunc_all,gap,np.log10(ramii),rsize1,np.log10(m200)]).T
            
        np.savetxt(ff, zz, fmt='%i %1.8f %1.8f %1.8f %1.8f %1.8f',delimiter='\t')
        if len(zz[:,2]) > 0:
             w2 = np.where((zz[:,2] == np.nanmin(zz[:,2])))[0]
        if len(w2)>0:
             rtrunc =   zz[w2[0],1]
             frestor_trunc = zz[w2[0],3]
             gap1 = zz[w2[0],2]
        else:
             rtrunc = 2*int(rsize)
             frestor_trunc = np.log10(ram*1e10)
       
#----------------------------------------------------------------------------------------------
        #Cold gas accretion
        mask_c = ((radius <= 2.* rsize)&(radius > rsize)&(temperature < temp_lim))
        if (mgas_c) > 0: 
            vx_gas = dvx_gas[mask_c] - sub['vel_x']
            vy_gas = dvy_gas[mask_c] - sub['vel_y']
            vz_gas = dvz_gas[mask_c] - sub['vel_z']
            rx = dx_gas[mask_c]
            ry = dy_gas[mask_c]
            rz = dz_gas[mask_c]
            Macc_c = np.sum((mgas[mask_c]/rsize) * (vx_gas*rx +vy_gas*ry + vz_gas*rz)/np.sqrt((rx)**2 + (ry)**2 + (rz)**2))
        else: 
            Macc_c = 0.0
        #Hot gas accretion
        
        mask_h = ((radius <= 2.* rsize)&(radius > rsize)&(temperature > temp_lim))
        if (mgas_h) > 0: 
            vx_gas = dvx_gas[mask_h] - sub['vel_x']
            vy_gas = dvy_gas[mask_h] - sub['vel_y']
            vz_gas = dvz_gas[mask_h] - sub['vel_z']
            rx = dx_gas[mask_h]
            ry = dy_gas[mask_h]
            rz = dz_gas[mask_h]
            Macc_h = np.sum((mgas[mask_h]/rsize) * (vx_gas*rx +vy_gas*ry + vz_gas*rz)/np.sqrt((rx)**2 + (ry)**2 + (rz)**2))
        else: 
            Macc_h = 0.0
#----------------------------------------------------------------------------------------------
        print '#parentid id central sub_id snap ramP  frest  massr  rsize  densr  dens velwrtgas gmass_cold gmass_hot gacc_cold gacc_hot age1 cent rtrunc' 
        print parentID, id, cent, sub['snap'], np.log10(ramp), np.log10(frest), np.log10(massr), rsize, np.log10(mgas_c),np.log10(mgas_h),Macc_c*1e-9,Macc_h*1e-9, age1,cent, rtrunc       
        main_id.append(id)
        sub_ID.append(sub['id'])
        snap.append(sub['snap'])
        ram_p.append(np.log10(ramp * 1e10))
        f_rest.append(np.log10(frest * 1e10))
        gmassr.append(massr)
        r_size.append(rsize)
        dens_r.append(np.log10(densr))
        densg.append(np.log10(dens))
        velwrtgas.append(np.sqrt(vel_wrt_gas2))
        gmass_cold.append(mgas_c)
        gmass_hot.append(mgas_h)
        gacc_cold.append(Macc_c*1e-9)# M_sun/year
        gacc_hot.append(Macc_h*1e-9)# M_sun/year
        hc_dist_center.append(gal_dist_center)
        hc_r200.append(hcr200)
        central.append(cent)
        redshift.append(redshift1)
        age.append(age1)
        r_trunc.append(rtrunc)
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#for other trace list:
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
ff.close()
# while sub['prog_sfid'] >= -1:
while sub['snap'] >= 50 :    
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
       
     #cutout_request = {'dm':'Coordinates' ,'gas':'Coordinates,Masses,Velocities,InternalEnergy,ElectronAbundance,GFM_Metallicity,Density','stars':'Coordinates,Masses,Velocities,GFM_StellarFormationTime'}
     #cutout =  get(sub['meta']['url'] + "cutout.hdf5" , cutout_request)
     cutout = "./cutout_ramp/cutout_"+str(sub['id'])+".hdf5"

    # make 2d histogram visualization of gas distribution
     with h5py.File(cutout,'r') as f:
        dx_DM          = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
        dy_DM          = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
        dz_DM          = f['PartType1']['Coordinates'][:,2] - sub['pos_z']
 
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
        
        temperature = (gama-1)* (u/KB)* mu* 1e10 
        #AGNRadiation        = f['PartType0']['GFM_AGNRadiation'][:]
        #CoolingRate         = f['PartType0']['GFM_CoolingRate'][:]
        Metallicity         = f['PartType0']['GFM_Metallicity'][:] / 0.0127
        #WindDMVelDisp       = f['PartType0']['GFM_WindDMVelDisp'][:]
        #StarFormationRate   = f['PartType0']['StarFormationRate'][:]
        Density             = f['PartType0']['Density'][:] # (10^10M/h)/(ckpc/h)3
        #NeutralHydrogenAbundance = f['PartType0']['NeutralHydrogenAbundance'][:]
        rsize = sub['halfmassrad_stars']
        #gas_profile = temperature
        radius_s      = np.sqrt(np.power(dx_star,2) +  np.power(dy_star,2) +  np.power(dz_star,2))
        radius      = np.sqrt(np.power(dx_gas,2) +  np.power(dy_gas,2) +  np.power(dz_gas,2))

        Msun = 1.989e+33
        kpc_in_cm=3.08568e+21
#----------------------------------------------------------------------------------------------              
        #star formation history test
        SFHp = f['PartType4']['GFM_StellarFormationTime'][:]
        SFH = np.mean(SFHp[(radius_s <= rsize)])
        age1 = cosmo.age((1/SFH)+1).value
#----------------------------------------------------------------------------------------------                
        #cold& hot gas
        mgas_c = np.sum(mgas[((radius <= rsize)& (temperature < temp_lim))])
        mgas_h = np.sum(mgas[((radius <= rsize)& (temperature > temp_lim))]) 
#----------------------------------------------------------------------------------------------        

        #ram pressure
        # massr = sub['massinhalfrad_gas'] * 1e10 # M/h
        massr = np.sum(mgas[(radius <= rsize)])
        densr = np.nanmedian(Density[(radius <= rsize)]* 1e10 * Msun/kpc_in_cm**3)   # (M/h)/(ckpc/h)3
        max_r       = round(max(radius))
        
        Grav = 4.3 * 1e-9 * 1000 # Km2 kpc Msun-1 s-2
        Pi = 3.14

        frest = (1e10 * (Pi/2.)*Grav) * (massr) * densr / (rsize)  # [g s-2 cm]  [g s-2 kpc-3]
        
        mask = ((radius <= 2.* rsize)&(radius > rsize))
        gasVelx = np.sum(mgas[mask]*dvx_gas[mask])/np.sum(mgas[mask])
        gasVely = np.sum(mgas[mask]*dvy_gas[mask])/np.sum(mgas[mask])
        gasVelz = np.sum(mgas[mask]*dvz_gas[mask])/np.sum(mgas[mask])
        vel_wrt_gas2 = (sub['vel_x'] - gasVelx)**2 + (sub['vel_y'] - gasVely)**2 + (sub['vel_z'] - gasVelz)**2
        dens = np.nanmedian(Density[mask]* 1e10 * Msun/kpc_in_cm**3)
        ramp = (1e10 * dens * vel_wrt_gas2)

#----------------------------------------------------------------------------------------------
#save as array in file FOR TRUNCATION RADII 
        bin_radii = -0.3
        rx = dx_gas
        ry = dy_gas
        rz = dz_gas
        
        rx_dm = dx_DM
        ry_dm = dy_DM
        rz_dm = dz_DM
        Grav = 4.3 * 1e-9 * 1000 # Km2 kpc Msun-1 s-2
        Pi = 3.14
        rtrunc_all = [];gap = [];ramii = [];gid = [];m200 = [];rsize1 = []
        for i in np.arange(2.0*(rsize),1,bin_radii):

            if i < rsize:
                maski = ( (np.sqrt((rx)**2 + (ry)**2 + (rz)**2)<= i))
                # maski_s = ( (np.sqrt((rx_s)**2 + (ry_s)**2 + (rz_s)**2)<= i))
                maski_dm = ( (np.sqrt((rx_dm)**2 + (ry_dm)**2 + (rz_dm)**2)<= i))
            else:
                maski = ( (np.sqrt((rx)**2 + (ry)**2 + (rz)**2)<= rsize))
                maski_dm = ( (np.sqrt((rx_dm)**2 + (ry_dm)**2 + (rz_dm)**2)<= rsize))
            densi =  np.nanmedian(Density[maski]* 1e10 * Msun/kpc_in_cm**3)
            #massdmi = np.sum(mass_DM[maski_dm])
            massdmi = len(rx_dm[maski_dm]) * TNG_dm_mass
            massgasi = np.sum(mgas[maski])
            massi = massgasi + massdmi
            rami = ((Pi/2.)*Grav) * (massi) * densi / (i)
            rtrunc_all.append(i)
            gap.append(abs(np.log10(ramp*1e10) - np.log10(rami*1e10)))
            ramii.append(rami)
            gid.append(id)
            m200.append(hcr200)
            rsize1.append(rsize)

        zz = np.array([gid,rtrunc_all,gap,np.log10(ramii)+10,rsize1,np.log10(m200)]).T
            
        np.savetxt(ff, zz, fmt='%i %1.8f %1.8f %1.8f %1.8f %1.8f',delimiter='\t')
        if len(zz[:,2]) > 0:
             w2 = np.where((zz[:,2] == np.nanmin(zz[:,2])))[0]
        if len(w2)>0:
             rtrunc =   zz[w2[0],1]
             frestor_trunc = zz[w2[0],3]
             gap1 = zz[w2[0],2]
        else:
             rtrunc = 2*int(rsize)
             frestor_trunc = np.log10(ram*1e10)        
#----------------------------------------------------------------------------------------------
        #Cold gas accretion
        mask_c = ((radius <= 2.* rsize)&(radius > rsize)&(temperature < temp_lim))
        if (mgas_c) > 0: 
            vx_gas = dvx_gas[mask_c] - sub['vel_x']
            vy_gas = dvy_gas[mask_c] - sub['vel_y']
            vz_gas = dvz_gas[mask_c] - sub['vel_z']
            rx = dx_gas[mask_c]
            ry = dy_gas[mask_c]
            rz = dz_gas[mask_c]
            Macc_c = np.sum((mgas[mask_c]/rsize) * (vx_gas*rx +vy_gas*ry + vz_gas*rz)/np.sqrt((rx)**2 + (ry)**2 + (rz)**2))
        else: 
            Macc_c = 0.0
        #Hot gas accretion
        
        mask_h = ((radius <= 2.* rsize)&(radius > rsize)&(temperature > temp_lim))
        if (mgas_h) > 0: 
            vx_gas = dvx_gas[mask_h] - sub['vel_x']
            vy_gas = dvy_gas[mask_h] - sub['vel_y']
            vz_gas = dvz_gas[mask_h] - sub['vel_z']
            rx = dx_gas[mask_h]
            ry = dy_gas[mask_h]
            rz = dz_gas[mask_h]
            Macc_h = np.sum((mgas[mask_h]/rsize) * (vx_gas*rx +vy_gas*ry + vz_gas*rz)/np.sqrt((rx)**2 + (ry)**2 + (rz)**2))
        else: 
            Macc_h = 0.0
#----------------------------------------------------------------------------------------------

        print '#parentid id central sub_id snap ramP  frest  massr  rsize  densr  dens velwrtgas gmass_cold gmass_hot gacc_cold gacc_hot age1 cent rtrunc ' 
        print parentID, id, cent, sub['snap'], np.log10(ramp), np.log10(frest), np.log10(massr), rsize, np.log10(mgas_c),np.log10(mgas_h),Macc_c*1e-9,Macc_h*1e-9, age1,cent,rtrunc         
        main_id.append(id)
        sub_ID.append(sub['id'])
        snap.append(sub['snap'])
        ram_p.append(np.log10(ramp))
        f_rest.append(np.log10(frest))
        gmassr.append(massr)
        r_size.append(rsize)
        dens_r.append(np.log10(densr))
        densg.append(np.log10(dens))
        velwrtgas.append(np.sqrt(vel_wrt_gas2))    
        gmass_cold.append(mgas_c)
        gmass_hot.append(mgas_h)
        gacc_cold.append(Macc_c*1e-9)# M_sun/year
        gacc_hot.append(Macc_h*1e-9)# M_sun/year
        hc_dist_center.append(gal_dist_center)
        hc_r200.append(hcr200)
        central.append(cent)
        redshift.append(redshift1)
        age.append(age1)
        r_trunc.append(rtrunc)

ff.close()
cc = np.array([main_id,sub_ID,snap,ram_p,f_rest,gmassr,r_size,dens_r,densg,velwrtgas,gmass_cold,gmass_hot,gacc_cold,gacc_hot,hc_dist_center,hc_r200,central,redshift,age,r_trunc]).T
np.savetxt(output, cc, fmt='%i %i %i %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %i %1.8f %1.8f %1.8f', delimiter='\t',header='ID\tSubID\tsnap\tramp\tfrest\tgmassr\tr_size\tdens_r\tdensg\tvelwrtgas\tgmass_cold\tgmass_hot\tgacc_cold\tgacc_hot\thc_dist_center\thc_r200\tcentral\tredshift\tage\tr_trunc')    
     
plt.plot(redshift, ram_p, lw= 3, color = 'blue',label = 'Real')
plt.plot(redshift, f_rest, lw= 3, color = 'red',label = 'Theory')
plt.xlabel('z')
plt.ylabel(r'$Ram\ Pressure\ [cgs]$')
plt.xlim(0.01, 0.5) 
plt.ylim(-18.5, -9)
plt.xscale('log')
plt.savefig('TNG_ramp_trace'+str(id) +'.pdf')
plt.close()



