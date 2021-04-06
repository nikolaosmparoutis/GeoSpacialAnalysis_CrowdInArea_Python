# to install the packages, on terminal run:   sh -e dependencies

import pandas as pd
import yaml
import geopandas as gpd
import os.path  # Ref lib line *os*
from os import path  # >>
from shapely.geometry import Point  # Ref lib line *shape*
from shapely.ops import nearest_points  # Ref lib line *near*
import math  # Ref lib line *math*
import geopy  # Ref lib line *Geo*
import geopy.distance as geod  # >>
import matplotlib.pyplot as plt  # Ref lib line *mat*
import folium  # Ref lib line *fol,web*
import webbrowser  # >>
from datetime import timedelta  # Ref lib line *timesdelta*

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', 10)


class ClassifyDevicesInArea:

    def __init__(self, absolute_path_devices_and_area):
        self.absolute_path = absolute_path_devices_and_area
        # self.devices = self.read_data()[0]
        # self.polygon = self.read_data()[1]

    def read_data(self):
        try:
            with open(self.absolute_path, 'r') as yml_file:
                cfg = yaml.load(yml_file, Loader=yaml.FullLoader)

                devices_df = pd.read_csv(cfg["path_csv"]["devices"],
                                         sep=",", encoding='latin-1', error_bad_lines=False, warn_bad_lines=False,
                                         low_memory=False, memory_map=True)

                polygon_geodf = gpd.GeoDataFrame(gpd.read_file(cfg["path_csv"]["polygon"]))

            return devices_df, polygon_geodf
        except FileNotFoundError:
            raise FileNotFoundError

    # optional for future use
    # def _swap_log_lat_polygon(geodf):
    #     import shapely.ops
    #     trans_coords = lambda Polygon: shapely.ops.transform(lambda x, y: (y, x), geom=Polygon)
    #     polygon_df = geodf.geometry.map(trans_coords)
    #     return polygon_df

    """This method converts the coordinates to Point keeping the same dataframe
    input: the overall devices
    output: the Points in a DataFrame"""

    @staticmethod
    def devices_location_to_points(overall_devices):
        # Ref lib line *shape*
        overall_devices.astype(float)
        P = []
        for long, lat in overall_devices[["latitude", "longitude"]].values:
            P.append(Point(lat, long))  # Point class by default  reverses its coordinates

        overall_devices.rename(columns={"latitude": "Point"}, inplace=True)
        overall_devices.drop(columns=['longitude'], inplace=True)
        overall_devices["Point"] = P
        return overall_devices

    """This method calculates the haversine formula
    input: the lat long of two points 
    returns: the distance in meters 
    """

    @staticmethod
    def distance_calc(lat1, lon1, lat2, lon2):
        # Ref lib line *math*
        R = 6378.137  # Radius of earth in KM

        dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180
        dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180

        a = math.sin(dLat / 2) * math.sin(dLat / 2) + \
            math.cos(lat1 * math.pi / 180) * \
            math.cos(lat2 * math.pi / 180) * \
            math.sin(dLon / 2) * math.sin(dLon / 2)
        theta = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d = R * theta
        return d * 1000  # meters

    """This method for each device considers the uncertainty factor in meters
    and for each device adds the uncertainty to the min distance from its nearest neighbor. 
    input:  the uncertainty in meters the min distance  
    output: the sum 
    """

    @staticmethod
    def add_uncertainty_to_min_distance(overall_devices, min_distances):
        new_distance = []
        for i in range(len(overall_devices)):
            sum_distance = overall_devices.loc[i, "uncertainty"] + min_distances[i]
            new_distance.append(sum_distance)

        overall_devices["new_distance"] = new_distance
        return overall_devices

    @staticmethod
    def find_nearest_point_and_distance(polyg, devices_ptns):
        # Ref lib line *near*
        nearest_coords = []
        min_distances = []
        for p in devices_ptns['Point']:
            p1, p2 = nearest_points(polyg, p)  # p1 is poly p2 is visitor's position approx for longitude
            nearest_coords.append(p1)
            # from geopy import distance # 0.10 millimetre in difference from our build in method <_distance_calc>
            # from geopy import Point # but converts the shapely.geometry.Point to geopy.Point and needs more work
            d = ClassifyDevicesInArea.distance_calc(p1.x, p1.y, p2.x, p2.y)
            min_distances.append(d)
        devices_ptns["nearest_coord"] = nearest_coords
        return min_distances, devices_ptns

    """This method interpolates the old points of devices to new positions.
        Using the  new distance of the points made by <add_uncertainty_to_min_distance>
        input: lat1, lon1, lat2, lon2, new_distance_i per point
        output: the new interpolated points to insert into <find_devices_in_polygon>
        """

    @staticmethod
    def interpolate_coordinates(devices_to_interpolate):
        new_point = []
        for i in range(0, len(devices_to_interpolate)):
            lat1_i = devices_to_interpolate["Point"].loc[i].y
            lon1_i = devices_to_interpolate["Point"].loc[i].x
            new_distance_i = devices_to_interpolate["new_distance"].loc[i]

            lat2_i = devices_to_interpolate["nearest_coord"].loc[i].y
            lon2_i = devices_to_interpolate["nearest_coord"].loc[i].x

            new_point.append(
                ClassifyDevicesInArea._helper_interpolate_coordinates(lat1_i, lon1_i, lat2_i, lon2_i, new_distance_i))
        devices_to_interpolate["new_point"] = new_point
        return devices_to_interpolate

    @staticmethod
    def _helper_interpolate_coordinates(lat1, lon1, lat2, lon2, new_distance_i):
        # Ref lib line *Geo*
        # given: lat1, lon1, default, bearing in degrees, default, distance in kilometres
        bearing = ClassifyDevicesInArea._get_bearing(lat1, lon1, lat2, lon2)
        origin = geopy.Point(lat1, lon1)
        destination = geod.distance(kilometers=new_distance_i / 1000).destination(origin, bearing)
        lon2, lat2 = destination.longitude, destination.latitude
        return Point(lon2, lat2)

    @staticmethod
    def _get_bearing(lat1, lon1, lat2, lon2):

        bearing = math.atan2(math.sin(lon2 - lon1) * math.cos(lat2), math.cos(lat1) * math.sin(lat2) - math.sin(lat1)
                             * math.cos(lat2) * math.cos(lon2 - lon1))
        bearing = math.degrees(bearing)
        bearing = (bearing + 360) % 360
        return bearing

    """This method finds the visitors (through their mobile device)
     who appeared inside the area and stores them in txt file.                                                       
     input : devices inside an area,                                                                                 
     output: the visitors as list"""

    @staticmethod
    def find_devices_in_polygon(overall_devices, poly):
        devices_in = pd.DataFrame()
        users_counter = 0
        for i in range(0, len(overall_devices)):
            if overall_devices["new_point"].loc[i].within(poly) is True:  # series object does not have within
                devices_in = devices_in.append(
                    overall_devices.loc[i])  # we use Point class The Point it self is class append to df is wrong
                users_counter += 1
        print("Traces inside the building:", users_counter)
        return devices_in

    @staticmethod
    def write_to_txt(df, txt_name, has_header):
        df.reset_index().to_csv(txt_name, sep='\t', header=has_header, index=False, mode="a")

    @staticmethod
    def plot_building_area(polygon):
        # Ref lib line *mat*
        polygon.plot(color='white', edgecolor='black')
        plt.title("building area", pad=20)
        plt.ion()  # keep the program after the show
        plt.draw()
        plt.pause(.001)  # for the GUI to take time to do the show

    @staticmethod
    def plot_points_map(polygon, chosen_devices_points, mapName):
        # Ref lib line *fol,web*
        chosen_devices_points.index = range(len(chosen_devices_points.index))
        fol_map = folium.Map(location=[59.12381366829599, 14.608602084061102],
                             tiles="OpenStreetMap", zoom_start=40)
        folium.GeoJson(polygon).add_to(fol_map)
        folium.LatLngPopup().add_to(fol_map)
        coords = [()]
        for i in range(1, len(chosen_devices_points)):  # position 0 is null, error for the Marker
            lat_long = (chosen_devices_points.loc[i].y,
                        chosen_devices_points.loc[i].x)
            coords.append(lat_long)
            folium.CircleMarker(location=coords[i],
                                radius=1, color='red').add_to(fol_map)
        fol_map.save(mapName)
        webbrowser.open(mapName)

    """This method groups the Points who found inside the area for each visitor
     input : the points inside the area 
     returns: the overall visits."""

    @staticmethod
    def group_devices_in_poly(devices_inside_poly):
        # Ref lib line os
        grouped_df = devices_inside_poly.groupby("hash_id")
        enterFirstTime = True

        if path.exists("../points_in_area.txt") is True:
            os.remove("../points_in_area.txt")

        for key, item in grouped_df:
            if enterFirstTime is True:
                ClassifyDevicesInArea.write_to_txt(grouped_df.get_group(key), 'points_in_area.txt', has_header=True)
            else:
                ClassifyDevicesInArea.write_to_txt(grouped_df.get_group(key), 'points_in_area.txt', has_header=False)
            enterFirstTime = False

        print("Num of visitors inside the area:", len(grouped_df))
        return grouped_df

    """ This method calculates the time inside the area for each visitor.
    The sampling time is per second so we do calculations on the groupedDF by HashId.
    Arguments: the classified points who found inside the area
    Returns: The time stay in format HH:MM:SS for each visitor inside the building. """

    @staticmethod
    def calculate_time_stay(visits):
        # Ref lib line *timesdelta*
        for key, item in visits["hash_id"]:
            sum_sampled_time = sum(item)
            print("Visitor {} stayed {:0>8}".format(key, str(timedelta(seconds=sum_sampled_time))))
