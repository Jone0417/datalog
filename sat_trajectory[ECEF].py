from skyfield.api import load, wgs84
from datetime import datetime
import numpy as np
import pyvista as pv
import math
import os
import time

def z_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[c, s, 0],
                    [-s, c, 0],
                    [0, 0, 1]])
    
    return (DCM @ vector.T).T

def x_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[1, 0, 0],
                    [0, c, -s],
                    [0, s, c]])
    
    return (DCM @ vector.T).T

def y_rotation(vector, angle):
    
    s = math.sin(angle * math.pi / 180)
    c = math.cos(angle * math.pi / 180)
    
    DCM = np.array([[c, 0, s],
                    [0, 1, 0],
                    [-s, 0, c]])
    
    return (DCM @ vector.T).T

def limit(val):
    if val > 1:
        val = 1
    elif val < -1:
        val = -1
    else:
        val = val
    return val

def quaternion_to_dcm(q): # Transform quaternion to DCM
    q1, q2, q3, q4 = np.transpose(q)
   
    # passive roataion (rotatoin of frame)
    dcm = np.array([[q4**2 + q1**2 -q2**2-q3**2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4)],
                     [2*(q1*q2-q3*q4), q4**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q1*q4)],
                     [2*(q1*q3+q2*q4), 2*(q2*q3 - q1*q4), q4**2 - q1**2 - q2**2 + q3**2]])
    return dcm

def lat_lon_at(satellite, ti, earth_r): # load latittude, longitude, and height information of satellite
    geo = satellite.at(ti)
    lat, lon = wgs84.latlon_of(geo)
    lat, lon = lat.arcminutes()/60, lon.arcminutes()/60
    hei = wgs84.height_of(geo).km
    
    return lat_lon_rotation(lat, lon, np.array([earth_r + hei,0,0]))

def lat_lon_rotation(lat, lon, vector): # calculate rotation in y,z axis with each vectors
    s= math.sin(lat * math.pi/180)
    c= math.cos(lat * math.pi/180)
    dcm_y = np.array([[c,0,s],
                      [0,1,0],
                      [-s,0,c]])
    s= math.sin(lon * math.pi/180)
    c= math.cos(lon * math.pi/180)
    dcm_z = np.array([[c,-s,0],
                      [s,c,0],
                      [0,0,1]])
    mat = (dcm_z @ (dcm_y @ vector))
    
    return mat.T
    
def visualize_sat(stations_url, ground_station, vectors, quaternion):               # main function
    start = time.time()
    satellites = load.tle_file(stations_url, reload= True)                          # load tle data from stations_url(in network)
    print('Loaded', len(satellites), 'satellites')
    print('satellites', satellites)
    
    target = input("satellite name: ")                                              # choose target satellite with its name in satellite lists
    if target not in {sat.name: sat for sat in satellites}:
        print("WRONG NAME")
    satellite = {sat.name: sat for sat in satellites}[target]
    print(satellite)
    
    sat_DCM = quaternion_to_dcm(quaternion)                                         # convert input quaternions into DCM
    basis = np.array([[100,0,0],
                      [0,100,0],
                      [0,0,100]])
    sat_attitude = (sat_DCM @ basis).T                                              # decide satellite attitude with DCM
    
    earth = pv.Sphere(radius=earth_r, theta_resolution=360, phi_resolution=360, start_theta=270.001, end_theta=270) 
    earth.active_t_coords = np.zeros((earth.points.shape[0], 2))                    # load 3D earth

    for i in range(earth.points.shape[0]):
        x,y,z = earth.points[i,0]/earth_r, earth.points[i,1]/earth_r, earth.points[i,2]/earth_r
        x,y,z = limit(x), limit(y), limit(z)
        earth.active_t_coords[i] = [0.5 + math.atan2(-x, y)/(2 * math.pi), 0.5 + math.asin(z)/math.pi]

    earth.rotate_z(-90)                                                             # match earth with zeroa_latitude line
    
    pl = pv.Plotter(window_size=(1600,1200))
    pl.add_background_image('starry-night-sky-fit.jpg', scale=1.001)                
    pl.add_mesh(earth, texture = pv.read_texture("8k_earth_daymap.jpg"), smooth_shading=False) 

    path = "/Users/jaewon/Desktop/J_one/WorkSpace"
    os.chdir(path)
    sat_mesh = pv.read("satellite.stl")
    sat_scale = 4
    gs_mesh = pv.read("gs.stl")
    gs_scale = 80
    colors = ['green', 'red' , 'blue', 'cyan', 'magenta','orange', 'pink', 'white']
    
    for i in range(gs_mesh.points.shape[0]):                                        # load ground station
        x_gs, y_gs, z_gs = gs_mesh.points[i,0]*gs_scale + ground_station[0], gs_mesh.points[i,1]*gs_scale + ground_station[1], gs_mesh.points[i,2]*gs_scale + ground_station[2] 
        gs_mesh.points[i,0], gs_mesh.points[i,1], gs_mesh.points[i,2] = x_gs, y_gs, z_gs
    pl.add_mesh(gs_mesh, color= colors[1])
    print('Ground Station[Cartesian]:',ground_station[0], ground_station[1], ground_station[2])

    trac_li = []
    available_li = [[],[]]
    ts = load.timescale()
    now = datetime.now()
    time_now = ts.utc(now.year, now.month, now. day, now.hour, now.minute, range(now.second, now.second + 6000, 1))
    print('Currrnt time: ', now)
    
    for ti in time_now:                                                             # load 6000 sec of satellite orbit start from now 
        x,y,z = lat_lon_at(satellite, ti, earth_r)
        trac = np.array([x,y,z])
        trac_li.append(trac)
        distance = math.sqrt((x-ground_station[0])**2 + (y-ground_station[1])**2 + (z-ground_station[2])**2)
        if distance <= 3000:                                                        
            available_li[0].append(trac)                                            # communication availale if distance <= 3000
            available_li[1].append(ti.utc_strftime('%Y %b %d %H:%M:%S'))            # communication available period
        
    for ele in trac_li:
        pl.add_lines(ele, width = 2, color= colors[0])
    for ele in available_li[0]:
        pl.add_lines(ele, width = 5, color = colors[-2])
        
    if available_li[1] == []:
        print('[ Satellite is out of communication range ]')
    else:
        print('Satellite communication-available period: from {0} to {1}'.format(available_li[1][0], available_li[1][-1]))

    tn = time_now[0]                                                                # current location of satellite
    x_sat, y_sat, z_sat = lat_lon_at(satellite,tn, earth_r)
    print('Current location[Cartesian]: [', x_sat, y_sat, z_sat,']')
    sat_location = np.array([x_sat, y_sat, z_sat])
    for i in range(sat_mesh.points.shape[0]):
        sat_mesh.points[i] = (sat_DCM @ (sat_mesh.points[i].T)).T
        x_sat, y_sat, z_sat = sat_mesh.points[i,0]*sat_scale + sat_location[0], sat_mesh.points[i,1]*sat_scale + sat_location[1], sat_mesh.points[i,2]*sat_scale + sat_location[2] 
        sat_mesh.points[i,0], sat_mesh.points[i,1], sat_mesh.points[i,2] = x_sat, y_sat, z_sat
    pl.add_mesh(sat_mesh, color= colors[-3])
    
    for i in range(3):                                                              # orthonormal satellite body frame 
        pl.add_lines(np.array([sat_location + 3*sat_attitude[i], sat_location]), color = colors[i], width = 3)

    for i in range(len(vectors)):                                                   # some vector directions from satellite(ex, sun vector, LOS, aos, main camera lens ) 
        loc = sat_DCM @ vectors[i]
        loc = loc / (math.sqrt(loc[0]**2 + loc[1]**2 + loc[2]**2)) * 200
        pl.add_lines(np.array([sat_location, sat_location + loc]), color= colors[i+3], width=4)

    for i in range(3):                                                              # earth orthonormal frame (z: rotation axis, x: equinox, y: orthonormal to x,z)
        pl.add_lines(np.array([[0,0,0], basis[i]/100*earth_r*1.2]), color = colors[i], width = 2)
    pl.add_axes(line_width= 5, labels_off=False)
    print('Operation time',round(time.time()-start,3), 'sec')
    pl.show()    

########################################################################################################################

## tle 정보를 가져오는 url
stations_url = 'https://celestrak.org/NORAD/elements/gp.php?INTDES=2022-072'
# stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle'

## 지상국 위치
earth_r = 6378 # 지구 반지름
ground_station = np.array([0,earth_r,0])    
# gs_2 (Daejeon)
lon, lat = 127.3845475, 36.3504119
dj = np.array([earth_r,0,0])
dj = y_rotation(dj, -lat)
dj = z_rotation(dj, -lon)
# ground_station = dj

## 임의의 지정 벡터
vectors = np.array([[1,1,1],[-1,-1,1]])
## 임의의 자세 지정 quaternion     
q1, q2, q3, q4 = 0,(1/4)**(1/2),(1/4)**(1/2),(1/2)**(1/2)
quaternion = np.array([q1, q2, q3, q4])


#### main ####
visualize_sat(stations_url, ground_station, vectors, quaternion)


