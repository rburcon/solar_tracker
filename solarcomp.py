#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 10:14:36 2017

@author: rburcon
"""

from tkinter import *
from tkinter import StringVar, ttk
import matplotlib.image as mpimg
import numpy as np

#######################################

import matplotlib.pyplot as plt
try:
    import seaborn as sns
    sns.set(rc={"figure.figsize": (12, 6)})
except ImportError:
    print('We suggest you install seaborn using conda or pip and rerun this cell')
    
# built in python modules
import datetime

# python add-ons
import numpy as np
import pandas as pd

import pvlib
from pvlib.tools import cosd, sind
from pvlib.location import Location

#######################################


root = Tk()
root.wm_title("Solar Radiation - Tracker vs Fixed Mount")


def tracker():
    data_ini=(datainicio_e.get())
    data_fim=(datafim_e.get())
    #entalpia=IAPWS97(P=press,x=titulo)
    #print ("\nCom a pressão de %.2f MPA e o título de %.2f entalpia é de -->> %.2f kj/kg" %(press, titulo, entalpia.h)  )
    
    lat=float(latitude_e.get())
    long=float(longitude_e.get())
    alt=float(altitude_e.get())
    fuso=box.get()
    cidade=cidade_e.get()
    angle=float(tracker_angle.get())
    
    local1 = Location(lat, long, fuso, alt, cidade)
    print(local1)
    ano_i=int(data_ini[6:10])
    mes_i=int(data_ini[3:5])
    dia_i=int(data_ini[0:2])

    ano_f=int(data_fim[6:10])
    mes_f=int(data_fim[3:5])
    dia_f=int(data_fim[0:2])

    #dados do tracker##############################################################
    axis_tilt = float(eixo_e.get())
    axis_azimuth = float(inclina_azimute_e.get())
    latitude = local1.latitude
    max_angle = local1.longitude
    backtrack = True
    #razao entre o comprimento da placa e a distancia entre linhas.
    #gcr = 2.0/6.0
    gcr=float(altura_e.get())/float(entre_e.get())
    backtrack_ON=True
    
    
    print('  Axis tilt %.2f, Axis Azimuth %.2f, GCR %.2f' %(axis_tilt, axis_azimuth, gcr))
    
    
    
    ###############################################################################
    times = pd.date_range(start=datetime.datetime(ano_i,mes_i,dia_i), end=datetime.datetime(ano_f,mes_f,dia_f), freq='5Min')

    ephem_local1 = pvlib.solarposition.get_solarposition(times.tz_localize(local1.tz), local1.latitude, local1.longitude)
    ephemout = ephem_local1

    azimuth = ephemout['azimuth']
    apparent_azimuth = ephemout['azimuth']
    apparent_zenith = ephemout['apparent_zenith']

    times = azimuth.index

    az = apparent_azimuth - 180
    apparent_elevation = 90 - apparent_zenith
    x = cosd(apparent_elevation) * sind(az)
    y = cosd(apparent_elevation) * cosd(az)
    z = sind(apparent_elevation)

    axis_azimuth_south = axis_azimuth - 180
    xp = x*cosd(axis_azimuth_south) - y*sind(axis_azimuth_south);
    yp = (x*cosd(axis_tilt)*sind(axis_azimuth_south) + y*cosd(axis_tilt)*cosd(axis_azimuth_south) - z*sind(axis_tilt))
    zp = (x*sind(axis_tilt)*sind(axis_azimuth_south) + y*sind(axis_tilt)*cosd(axis_azimuth_south) + z*cosd(axis_tilt))

    # Calculate angle from x-y plane to projection of sun vector onto x-z plane
    # and then obtain wid by translating tmp to convention for rotation angles.
    wid = pd.Series(90 - np.degrees(np.arctan2(zp, xp)), index=times)

    # filter for sun above panel horizon
    wid[zp <= 0] = np.nan

    

    if backtrack:
        axes_distance = 1/gcr
        temp = np.minimum(axes_distance*cosd(wid), 1)

        # backtrack angle
        # (always positive b/c acosd returns values between 0 and 180)
        wc = np.degrees(np.arccos(temp))

        v = wid < 0
        widc = pd.Series(index=times)
        widc[~v] = wid[~v] - wc[~v]; # Eq 4 applied when wid in QI
        widc[v] = wid[v] + wc[v];    # Eq 4 applied when wid in QIV
    else:
        widc = wid


    times = pd.date_range(start=datetime.datetime(ano_i,mes_i,dia_i), end=datetime.datetime(ano_f,mes_f,dia_f), freq='5Min')

    ephem_local1 = pvlib.solarposition.get_solarposition(times.tz_localize(local1.tz), local1.latitude, local1.longitude)

    tracker_data = pvlib.tracking.singleaxis(ephem_local1['apparent_zenith'], ephem_local1['azimuth'],
                                         axis_tilt=axis_tilt, axis_azimuth=axis_azimuth, max_angle=angle,
                                         backtrack=backtrack_ON, gcr=gcr)
    #tracker_data.plot()
    #plt.ylim(-100,100)

    irrad_data = local1.get_clearsky(times.tz_localize(local1.tz))
    #print(irrad_data)
    dni_et = pvlib.irradiance.get_extra_radiation(irrad_data.index, method='asce')
    #plt.figure()
    
    #Plota as componentes da radiaçao
    #irrad_data.plot()
    #dni_et.plot(label='DNI ET')

    ground_irrad = pvlib.irradiance.get_ground_diffuse(tracker_data['surface_tilt'], irrad_data['ghi'], albedo=.25)
    
    
    #print(radiafixa)
    
    #ground_irrad.plot()
    ephem_data = ephem_local1

    haydavies_diffuse = pvlib.irradiance.haydavies(tracker_data['surface_tilt'], tracker_data['surface_azimuth'], 
                                                irrad_data['dhi'], irrad_data['dni'], dni_et,
                                                ephem_data['apparent_zenith'], ephem_data['azimuth'])
    #haydavies_diffuse.plot(label='haydavies diffuse')

    global_in_plane = cosd(tracker_data['aoi'])*irrad_data['dni'] + haydavies_diffuse + ground_irrad
    global_in_plane.plot(label='global in plane'+'(tracker)')
    
    radiatotal=global_in_plane.sum()
    
    
    #################################################
    dni_et = pvlib.irradiance.get_extra_radiation(times.dayofyear)
    sun_zen = ephem_data['apparent_zenith']
    AM = pvlib.atmosphere.get_relative_airmass(sun_zen)
    #models = ['isotropic', 'klucher', 'haydavies', 'reindl', 'king', 'perez']
    models = ['isotropic']
    totals = {}
    
    for model in models:
        total = pvlib.irradiance.get_total_irradiance(lat, 180,ephem_data['apparent_zenith'], ephem_data['azimuth'], dni=irrad_data['dni'], ghi=irrad_data['ghi'], dhi=irrad_data['dhi'], dni_extra=dni_et, airmass=AM, model=model, surface_type='urban')
        totals[model] = total
        #total.plot()
        #plt.title(model)
        #plt.ylabel('Irradiance (W/m^2)')
    
    #plt.figure()
    for model, total in totals.items():
       total['poa_global'].plot(lw=.5, label=model+'(fixed mount)')
    
    plt.legend()
    plt.ylabel('Irradiance (W/m^2)')
    
    
    radiafixa=total.poa_global.sum()
    #print(total)
    
    ################################################
    
    
    rendimento=((radiatotal/radiafixa)-1)*100
    
              
    print('Maximum theoretical yield-->%.1f'%rendimento,'%')
    
    #print(radiatotal)
    #print(ground_irrad)
    #plt.legend()
    #plt.ylabel('Irradiance $W/m^2$')
    plt.title(cidade+' Tracker Irradiance')
    plt.savefig(r'tracker.png')
 
      
    #f = open('data.csv', 'w')
    #f.write(str(global_in_plane))
    #f.close()
    plt.show()
    



def auxiliar():
    #img=mpimg.imread('auxiliar.png')
    #imgplot = plt.imshow(img)
    #plt.show()
    from PIL import Image        
    img = Image.open('auxiliar.png')
    img.show() 
    
def mapa():
    #img=mpimg.imread('auxiliar.png')
    #imgplot = plt.imshow(img)
    #plt.show()
    from PIL import Image        
    img = Image.open('mapasolar.png')
    img.show()    

def sobre():
    from tkinter import messagebox
    var = messagebox.showinfo("Info" , "Comparative solar radiation received - Tracker vs Fixed Mount - ROMAGNOLE PRODUTOS ELETRICOS S/A")


def posi():
    data_ini=(datainicio_e.get())
    data_fim=(datafim_e.get())
    #entalpia=IAPWS97(P=press,x=titulo)
    #print ("\nCom a pressão de %.2f MPA e o título de %.2f entalpia é de -->> %.2f kj/kg" %(press, titulo, entalpia.h)  )
    
    lat=float(latitude_e.get())
    long=float(longitude_e.get())
    alt=float(altitude_e.get())
    fuso=box.get()
    cidade=cidade_e.get()
    angle=float(tracker_angle.get())
    
    local1 = Location(lat, long, fuso, alt, cidade)
    print(local1)
    ano_i=int(data_ini[6:10])
    mes_i=int(data_ini[3:5])
    dia_i=int(data_ini[0:2])

    ano_f=int(data_fim[6:10])
    mes_f=int(data_fim[3:5])
    dia_f=int(data_fim[0:2])

    #dados do tracker##############################################################
    axis_tilt = float(eixo_e.get())
    axis_azimuth = float(inclina_azimute_e.get())
    latitude = local1.latitude
    max_angle = local1.longitude
    backtrack = True
    #razao entre o comprimento da placa e a distancia entre linhas.
    #gcr = 2.0/6.0
    gcr=float(altura_e.get())/float(entre_e.get())
    backtrack_ON=True
    
    
    print('  Axis tilt %.2f, Axis Azimuth %.2f, GCR %.2f' %(axis_tilt, axis_azimuth, gcr))
    
    
    
    ###############################################################################
    times = pd.date_range(start=datetime.datetime(ano_i,mes_i,dia_i), end=datetime.datetime(ano_f,mes_f,dia_f), freq='5Min')

    ephem_local1 = pvlib.solarposition.get_solarposition(times.tz_localize(local1.tz), local1.latitude, local1.longitude)
    ephemout = ephem_local1

    azimuth = ephemout['azimuth']
    apparent_azimuth = ephemout['azimuth']
    apparent_zenith = ephemout['apparent_zenith']

    times = azimuth.index

    az = apparent_azimuth - 180
    apparent_elevation = 90 - apparent_zenith
    x = cosd(apparent_elevation) * sind(az)
    y = cosd(apparent_elevation) * cosd(az)
    z = sind(apparent_elevation)

    axis_azimuth_south = axis_azimuth - 180
    xp = x*cosd(axis_azimuth_south) - y*sind(axis_azimuth_south);
    yp = (x*cosd(axis_tilt)*sind(axis_azimuth_south) + y*cosd(axis_tilt)*cosd(axis_azimuth_south) - z*sind(axis_tilt))
    zp = (x*sind(axis_tilt)*sind(axis_azimuth_south) + y*sind(axis_tilt)*cosd(axis_azimuth_south) + z*cosd(axis_tilt))

    # Calculate angle from x-y plane to projection of sun vector onto x-z plane
    # and then obtain wid by translating tmp to convention for rotation angles.
    wid = pd.Series(90 - np.degrees(np.arctan2(zp, xp)), index=times)

    # filter for sun above panel horizon
    wid[zp <= 0] = np.nan

    wid.plot(label='tracking angle')
    ephemout['apparent_elevation'].plot(label='apparent elevation')
    plt.legend()
    plt.title('Ideal tracking angle without backtracking')


    if backtrack:
        axes_distance = 1/gcr
        temp = np.minimum(axes_distance*cosd(wid), 1)

        # backtrack angle
        # (always positive b/c acosd returns values between 0 and 180)
        wc = np.degrees(np.arccos(temp))

        v = wid < 0
        widc = pd.Series(index=times)
        widc[~v] = wid[~v] - wc[~v]; # Eq 4 applied when wid in QI
        widc[v] = wid[v] + wc[v];    # Eq 4 applied when wid in QIV
    else:
        widc = wid

    widc.plot(label='tracking angle')
    #pyephemout['apparent_elevation'].plot(label='apparent elevation')
    plt.legend(loc=2)
    plt.title('Ideal tracking angle with backtracking')

    max_angle=angle
    tracking_angles = pd.DataFrame({'with backtracking':widc,'without backtracking':wid})
    tracker_theta = widc.copy()
    tracker_theta[tracker_theta > max_angle] = max_angle
    tracker_theta[tracker_theta < -max_angle] = -max_angle

    tracking_angles['with restriction'] = tracker_theta
    tracking_angles.plot()
    plt.title('Ideal tracking angle with backtracking and restriction')
    plt.show()


def anual():
    data_ini='01-01-2017'
    data_fim='31-12-2017'
    #entalpia=IAPWS97(P=press,x=titulo)
    #print ("\nCom a pressão de %.2f MPA e o título de %.2f entalpia é de -->> %.2f kj/kg" %(press, titulo, entalpia.h)  )
    
    lat=float(latitude_e.get())
    long=float(longitude_e.get())
    alt=float(altitude_e.get())
    fuso=box.get()
    cidade=cidade_e.get()
    angle=float(tracker_angle.get())
    
    local1 = Location(lat, long, fuso, alt, cidade)
    print(local1)
    ano_i=int(data_ini[6:10])
    mes_i=int(data_ini[3:5])
    dia_i=int(data_ini[0:2])

    ano_f=int(data_fim[6:10])
    mes_f=int(data_fim[3:5])
    dia_f=int(data_fim[0:2])

    #dados do tracker##############################################################
    axis_tilt = float(eixo_e.get())
    axis_azimuth = float(inclina_azimute_e.get())
    latitude = local1.latitude
    max_angle = local1.longitude
    backtrack = True
    #razao entre o comprimento da placa e a distancia entre linhas.
    #gcr = 2.0/6.0
    gcr=float(altura_e.get())/float(entre_e.get())
    backtrack_ON=True
    
    print('  Axis tilt %.2f, Axis Azimuth %.2f, GCR %.2f' %(axis_tilt, axis_azimuth, gcr))
    
    
    ###############################################################################
    times = pd.date_range(start=datetime.datetime(ano_i,mes_i,dia_i), end=datetime.datetime(ano_f,mes_f,dia_f), freq='1h')

    ephem_local1 = pvlib.solarposition.get_solarposition(times.tz_localize(local1.tz), local1.latitude, local1.longitude)
    ephemout = ephem_local1

    azimuth = ephemout['azimuth']
    apparent_azimuth = ephemout['azimuth']
    apparent_zenith = ephemout['apparent_zenith']

    times = azimuth.index

    az = apparent_azimuth - 180
    apparent_elevation = 90 - apparent_zenith
    x = cosd(apparent_elevation) * sind(az)
    y = cosd(apparent_elevation) * cosd(az)
    z = sind(apparent_elevation)

    axis_azimuth_south = axis_azimuth - 180
    xp = x*cosd(axis_azimuth_south) - y*sind(axis_azimuth_south);
    yp = (x*cosd(axis_tilt)*sind(axis_azimuth_south) + y*cosd(axis_tilt)*cosd(axis_azimuth_south) - z*sind(axis_tilt))
    zp = (x*sind(axis_tilt)*sind(axis_azimuth_south) + y*sind(axis_tilt)*cosd(axis_azimuth_south) + z*cosd(axis_tilt))

    # Calculate angle from x-y plane to projection of sun vector onto x-z plane
    # and then obtain wid by translating tmp to convention for rotation angles.
    wid = pd.Series(90 - np.degrees(np.arctan2(zp, xp)), index=times)

    # filter for sun above panel horizon
    wid[zp <= 0] = np.nan

    

    if backtrack:
        axes_distance = 1/gcr
        temp = np.minimum(axes_distance*cosd(wid), 1)

        # backtrack angle
        # (always positive b/c acosd returns values between 0 and 180)
        wc = np.degrees(np.arccos(temp))

        v = wid < 0
        widc = pd.Series(index=times)
        widc[~v] = wid[~v] - wc[~v]; # Eq 4 applied when wid in QI
        widc[v] = wid[v] + wc[v];    # Eq 4 applied when wid in QIV
    else:
        widc = wid


    times = pd.date_range(start=datetime.datetime(ano_i,mes_i,dia_i), end=datetime.datetime(ano_f,mes_f,dia_f), freq='5Min')

    ephem_local1 = pvlib.solarposition.get_solarposition(times.tz_localize(local1.tz), local1.latitude, local1.longitude)

    tracker_data = pvlib.tracking.singleaxis(ephem_local1['apparent_zenith'], ephem_local1['azimuth'],
                                         axis_tilt=axis_tilt, axis_azimuth=axis_azimuth, max_angle=angle,
                                         backtrack=backtrack_ON, gcr=gcr)
    #tracker_data.plot()
    #plt.ylim(-100,100)

    irrad_data = local1.get_clearsky(times.tz_localize(local1.tz))
    #print(irrad_data)
    dni_et = pvlib.irradiance.get_extra_radiation(irrad_data.index, method='asce')
    #plt.figure()
    #irrad_data.plot()
    #dni_et.plot(label='DNI ET')

    ground_irrad = pvlib.irradiance.get_ground_diffuse(tracker_data['surface_tilt'], irrad_data['ghi'], albedo=.25)
    
        
    #ground_irrad.plot()
    ephem_data = ephem_local1

    haydavies_diffuse = pvlib.irradiance.haydavies(tracker_data['surface_tilt'], tracker_data['surface_azimuth'], 
                                                irrad_data['dhi'], irrad_data['dni'], dni_et,
                                                ephem_data['apparent_zenith'], ephem_data['azimuth'])
    #haydavies_diffuse.plot(label='haydavies diffuse')

    global_in_plane = cosd(tracker_data['aoi'])*irrad_data['dni'] + haydavies_diffuse + ground_irrad
    #global_in_plane.plot(label='global in plane')
   
    radiatotal=global_in_plane.sum()
    
    
    #print(radiatotal)
    #################################################
    dni_et = pvlib.irradiance.get_extra_radiation(times.dayofyear)
    sun_zen = ephem_data['apparent_zenith']
    AM = pvlib.atmosphere.get_relative_airmass(sun_zen)
    #models = ['isotropic', 'klucher', 'haydavies', 'reindl', 'king', 'perez']
    models = ['isotropic']
    totals = {}
    
    for model in models:
        total = pvlib.irradiance.get_total_irradiance(lat, 180,ephem_data['apparent_zenith'], ephem_data['azimuth'], dni=irrad_data['dni'], ghi=irrad_data['ghi'], dhi=irrad_data['dhi'], dni_extra=dni_et, airmass=AM, model=model, surface_type='urban')
        totals[model] = total
        #total.plot()
        #plt.title(model)
        #plt.ylabel('Irradiance (W/m^2)')
    
    #plt.figure()
    #for model, total in totals.items():
       #total['poa_global'].plot(lw=.5, label=model)
    
    #plt.legend()
    #plt.ylabel('Irradiance (W/m^2)')
    
    
    radiafixa=total.poa_global.sum()
   
    #print(radiafixa)
    
    ################################################
    
    
    rendimento=((radiatotal/radiafixa)-1)*100
    #radiafixa=radiafixa/1000
    #radiatotal=radiatotal/1000
    #print('Maximum Solar Radiation for Fixed Mount-->%.1f'%radiafixa, 'kW.h/m2.year')
    #print('Maximum Solar Radiation for Tracer-->%.1f'%radiatotal, 'kW.h/m2')
    print('Maximum theoretical yearly yield-->%.1f'%rendimento,'%')
    
    #print(radiatotal)
    #print(ground_irrad)
    #plt.legend()
    #plt.ylabel('Irradiance $W/m^2$')
    #plt.title(cidade+' Tracker Irradiance')
   
 
      
    #f = open('data.csv', 'w')
    #f.write(str(global_in_plane))
    #f.close()
        
    
    #plt.show()


##################################################################

Label(root, text="Start Date [dd-mm-yyyy]").grid(row=0, sticky=W)
Label(root, text="End Date [dd-mm-yyyy]").grid(row=1, sticky=W)
Label(root, text="Latitude [degrees]").grid(row=2, sticky=W)
Label(root, text="Longitude [degrees]").grid(row=3, sticky=W)
Label(root, text="Elevation [m]").grid(row=4, sticky=W)
Label(root, text="Timezone").grid(row=5, sticky=W)
Label(root, text="City").grid(row=6, sticky=W)
Label(root, text="Shaft Tilt [degrees]").grid(row=7, sticky=W)
Label(root, text="Azimuth Slope [degrees]").grid(row=8, sticky=W)
Label(root, text="Plate Height [m]").grid(row=9, sticky=W)
Label(root, text="Length Between Axis [m]").grid(row=10, sticky=W)
Label(root, text="Sweep angle [degrees]").grid(row=11, sticky=W)




datainicio_e=Entry(root)
datainicio_e.grid(row=0, column=1, sticky=E)
datainicio_e.insert(0,"01-01-2018")

datafim_e=Entry(root)
datafim_e.grid(row=1, column=1, sticky=E)
datafim_e.insert(0,"02-01-2018")

latitude_e=Entry(root)
latitude_e.grid(row=2, column=1, sticky=E)
latitude_e.insert(0,"-23.54")

longitude_e=Entry(root)
longitude_e.grid(row=3, column=1, sticky=E)
longitude_e.insert(0,"-51.68")

altitude_e=Entry(root)
altitude_e.grid(row=4, column=1, sticky=E)
altitude_e.insert(0,"670")

value = StringVar()
box = ttk.Combobox(root, textvariable=value, state='readonly', width=18)
#box['values'] = ('America/Noronha', 'America/Belem','America/Fortaleza', 'America/Recife', 'America/Araguaina', 'America/Maceio', 'America/Bahia', 'America/Sao_Paulo','America/Campo_Grande', 'America/Cuiaba', 'America/Santarem', 'America/Porto_Velho', 'America/Boa_Vista', 'America/Manaus', 'America/Eirunepe', 'America/Rio_Branco')
box['values'] = ('America/Noronha', 'America/Bahia', 'America/Manaus', 'America/Rio_Branco')
box.current(1)
box.grid(column=1, row=5)

cidade_e=Entry(root)
cidade_e.grid(row=6, column=1, sticky=E)
cidade_e.insert(0,"Mandaguari")

eixo_e=Entry(root)
eixo_e.grid(row=7, column=1, sticky=E)
eixo_e.insert(0,"0")

inclina_azimute_e=Entry(root)
inclina_azimute_e.grid(row=8, column=1, sticky=E)
inclina_azimute_e.insert(0,"0")

altura_e=Entry(root)
altura_e.grid(row=9, column=1, sticky=E)
altura_e.insert(0,"2.0")

entre_e=Entry(root)
entre_e.grid(row=10, column=1, sticky=E)
entre_e.insert(0,"6.0")

tracker_angle=Entry(root)
tracker_angle.grid(row=11, column=1, sticky=E)
tracker_angle.insert(0,"60.0")

Button(root, text="Tracker ", command=tracker).grid(row=6, column=3, sticky=E)
Button(root, text="Position", command=posi).grid(row=7, column=3, sticky=E)
Button(root, text="A. yield ", command=anual).grid(row=8, column=3, sticky=E)

Button(root, text="Help      ", command=auxiliar).grid(row=9, column=3, sticky=E)
Button(root, text="Map       ", command=mapa).grid(row=10, column=3, sticky=E)
Button(root, text="About    ", command=sobre).grid(row=11, column=3, sticky=E)



root.mainloop()

