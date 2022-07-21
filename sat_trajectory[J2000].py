from skyfield.api import load
from datetime import timedelta, datetime
import numpy as np
from astropy.coordinates import SkyCoord, GCRS, FK5
import pyvista as pv
import math
import os
import time

def limit(val):
    if val > 1:
        val = 1
    elif val < -1:
        val = -1
    else:
        val = val
    return val

def cal_rortation(new_date):
    ts = load.timescale()
    julian = ts.utc(new_date.year, new_date.month, new_date.day, new_date.hour, new_date.minute, new_date.second)
    reference_solar_jul = ts.utc(new_date.year, 1,1,0,0,0)
    reference_earth_jul = ts.utc(new_date.year, new_date.month,new_date.day,0,0,0)
    rotation_solar = (julian.tdb - reference_solar_jul.tdb)/365*360
    rotation_earth = (julian.tdb - reference_earth_jul.tdb)/365*360
    
    return (rotation_solar - 90) - (rotation_earth - 180)

def z_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[c, s, 0],
                    [-s, c, 0],
                    [0, 0, 1]])
    
    return (DCM @ vector.T).T

def y_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[c, 0, s],
                    [0, 1, 0],
                    [-s, 0, c]])
    
    return (DCM @ vector.T).T

def x_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[1, 0, 0],
                    [0, c, s],
                    [0, -s, c]])
    
    return (DCM @ vector.T).T

def quaternion_to_dcm(q): # Transform quaternion to DCM
    q1, q2, q3, q4 = np.transpose(q)
   
    # passive roataion (rotatoin of frame)
    dcm = np.array([[q4**2 + q1**2 -q2**2-q3**2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4)],
                     [2*(q1*q2-q3*q4), q4**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q1*q4)],
                     [2*(q1*q3+q2*q4), 2*(q2*q3 - q1*q4), q4**2 - q1**2 - q2**2 + q3**2]])
    return dcm

def visualize_sat_j2k(stations_url, ground_station, vectors, quaternion):
    start = time.time()
    satellites = load.tle_file(stations_url, reload=True)
    print('Loaded', len(satellites), 'satellites')
    print('satellites', satellites)

    target = input("satellite name: ")
    if target not in {sat.name: sat for sat in satellites}:
        print("WRONG NAME")
    satellite = {sat.name: sat for sat in satellites}[target]
    print(satellite)
    
    sat_DCM = quaternion_to_dcm(quaternion)
    basis = np.array([[100,0,0],
                      [0,100,0],
                      [0,0,100]])
    sat_attitude = (sat_DCM @ basis).T
    
    earth = pv.Sphere(radius=earth_r, theta_resolution=360, phi_resolution=360, start_theta=270.001, end_theta=270)
    earth.active_t_coords = np.zeros((earth.points.shape[0], 2))

    for i in range(earth.points.shape[0]):
        x,y,z = earth.points[i,0]/earth_r, earth.points[i,1]/earth_r, earth.points[i,2]/earth_r
        x,y,z = limit(x), limit(y), limit(z)
        earth.active_t_coords[i] = [0.5 + math.atan2(-x, y)/(2 * math.pi), 0.5 + math.asin(z)/math.pi]

    new_date = datetime.now()
    print('Current Time: ',new_date)
    rot_angle = cal_rortation(new_date)
    earth.rotate_z(-90 + rot_angle)
    ts = load.timescale()
    tn = ts.now()
    geocentric = satellite.at(tn)
    sat_location = geocentric.position.km
       
    pl = pv.Plotter(window_size=(1600,1200))
    pl.add_background_image('starry-night-sky-fit.jpg', scale=1.001)
    pl.add_mesh(earth, texture = pv.read_texture("8k_earth_daymap.jpg"), smooth_shading=False)

    path = "/Users/jaewon/Desktop/J_one/WorkSpace"
    os.chdir(path)
    sat_mesh = pv.read("satellite.stl")
    sat_scale = 4
    gs_mesh = pv.read("gs.stl")
    gs_scale = 80
    colors = ['green', 'red' , 'blue', 'cyan', 'orange', 'pink', 'white','lightgreen']
      
    # 위성의 궤도
    steps =6000
    orbit = []

    for i in range(steps):
        new_date = new_date + timedelta(seconds = 1)
        yr, mon, day = new_date.year, new_date.month, new_date.day
        hour, min, sec = new_date.hour, new_date.minute, new_date.second
        
        t_now = ts.utc(yr,mon,day,hour,min,sec)
        location = satellite.at(t_now)
        orbit.append(np.array(location.position.km))
        
    xx, yy, zz = zip(*orbit)
    ################################################################### J2K frame ###################################################################
    gs_j2k = SkyCoord(x= ground_station[0], y= ground_station[1], z= ground_station[2], frame = GCRS(), unit= 'kpc', representation_type= 'cartesian').transform_to(FK5(equinox = 'j2000'))
    gs_j2k.representation_type = 'cartesian'

    loc_j2k = SkyCoord(x= xx, y= yy, z= zz, frame = GCRS(), unit= 'kpc', representation_type= 'cartesian').transform_to(FK5(equinox = 'j2000'))
    loc_j2k.representation_type = 'cartesian'
    orbit_j2k = np.array(list(zip(loc_j2k.x.value, loc_j2k.y.value, loc_j2k.z.value)))
    for i in range(len(orbit_j2k)):
        pl.add_mesh(orbit_j2k[i], point_size=2, color=colors[-1])
        distance = math.sqrt((loc_j2k.x.value[i]- ground_station[0])**2 + (loc_j2k.y.value[i]- ground_station[1])**2 +(loc_j2k.z.value[i]- ground_station[2])**2)
        if distance <= 3000:
            pl.add_mesh(orbit_j2k[i], point_size=7, color= colors[-3])
    # print("Current location[J2K]: ",loc_j2k.x.value, loc_j2k.y.value, loc_j2k.z.value)

    sat_j2k = SkyCoord(x= sat_location[0], y=sat_location[1], z= sat_location[2], frame = GCRS(), unit= 'kpc', representation_type= 'cartesian').transform_to(FK5(equinox = 'j2000'))
    sat_j2k.representation_type = 'cartesian'
    sat_location = np.array([sat_j2k.x.value, sat_j2k.y.value, sat_j2k.z.value])
    print("Current location[J2K]: ",sat_location[0], sat_location[1], sat_location[2])
    for i in range(gs_mesh.points.shape[0]):
        x_gs, y_gs, z_gs = gs_mesh.points[i,0]*gs_scale + ground_station[0], gs_mesh.points[i,1]*gs_scale + ground_station[1], gs_mesh.points[i,2]*gs_scale + ground_station[2] 
        gs_mesh.points[i,0], gs_mesh.points[i,1], gs_mesh.points[i,2] = x_gs, y_gs, z_gs

    for i in range(sat_mesh.points.shape[0]):
        sat_mesh.points[i] = (sat_DCM @ (sat_mesh.points[i].T)).T
        x_sat, y_sat, z_sat = sat_mesh.points[i,0]*sat_scale + sat_location[0], sat_mesh.points[i,1]*sat_scale + sat_location[1], sat_mesh.points[i,2]*sat_scale + sat_location[2] 
        sat_mesh.points[i,0], sat_mesh.points[i,1], sat_mesh.points[i,2] = x_sat, y_sat, z_sat

    pl.add_mesh(gs_mesh, color= colors[1])
    pl.add_mesh(dj, point_size= 10, color= colors[3], render_points_as_spheres= True)
    pl.add_mesh(sat_mesh, color= colors[-2])
    # 위성의 자세
    for i in range(3):
        pl.add_lines(np.array([sat_location + 3*sat_attitude[i], sat_location]), color = colors[i], width = 3)
    # 임의의 벡터
    for i in range(len(vectors)):
        loc = sat_DCM @ vectors[i]
        loc = loc / (math.sqrt(loc[0]**2 + loc[1]**2 + loc[2]**2)) * 200
        pl.add_lines(np.array([sat_location, sat_location + loc]), color= colors[i+3], width=4)
    # 좌표계
    for i in range(3):
        pl.add_lines(np.array([[0,0,0], basis[i]/100*earth_r*1.2]), color = colors[i], width = 2)
    pl.add_axes(line_width= 5, labels_off=False)
    
    print("Operation time:",round(time.time()-start, 3), "sec") # 실행시간 측정
    pl.show()
   
################################################ 지정 값 #################################################################
## tle 정보를 가져오는 url
stations_url = 'https://celestrak.org/NORAD/elements/gp.php?INTDES=2022-072'
# stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle'
## 지상국 위치
earth_r = 6378
ground_station = np.array([0,earth_r,0])
lon, lat = 127.3845475, 36.3504119
rot_angle = cal_rortation(datetime.now())
dj = np.array([earth_r,0,0])
dj = y_rotation(dj, -lat)
dj = z_rotation(dj, lon + rot_angle - 90)

## 임의의 지정 벡터
vectors = np.array([[1,1,1],[-1,-1,1]])
## 초기 자세에 대한 쿼터니온
# q1, q2, q3, q4 = 0.254887, -0.04494346, 0.16773126, 0.95125124
q1, q2, q3, q4 = 0,(1/4)**(1/2),(1/4)**(1/2),(1/2)**(1/2)
quaternion = np.array([q1, q2, q3, q4])

## main

# rot_angle = cal_rortation(datetime.now())
# print(-90 + rot_angle - lon)
visualize_sat_j2k(stations_url, ground_station, vectors, quaternion)
