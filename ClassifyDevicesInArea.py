# to install the packages, on terminal run:   sh -e dependencies

import pandas as pd
import yaml
import geopandas as gpd
from shapely.geometry import Point

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

    """Convert the coordinates to Point keeping the same dataframe
    input: the overall devices
    output: the Points in a DataFrame"""

    @staticmethod
    def devices_location_to_points(overall_devices):
        overall_devices.astype(float)
        P = []
        for long, lat in overall_devices[["latitude", "longitude"]].values:
            P.append(Point(lat, long))  # Point class by default  reverses its coordinates

        overall_devices.rename(columns={"latitude": "Point"}, inplace=True)
        overall_devices.drop(columns=['longitude'], inplace=True)
        overall_devices["Point"] = P
        return overall_devices

    """Haversine formula
    input: the lat long of two points 
    returns: the distance in meters 
    """

    @staticmethod
    def distance_calc(lat1, lon1, lat2, lon2):
        import math
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

    """ For each device considers the uncertainty factor in meters
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
        from shapely.ops import nearest_points
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

    """Interpolates the old points of devices to new positions.
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
        import geopy
        import geopy.distance as geod
        # given: lat1, lon1, default, bearing in degrees, default, distance in kilometres
        bearing = ClassifyDevicesInArea._get_bearing(lat1, lon1, lat2, lon2)
        origin = geopy.Point(lat1, lon1)
        destination = geod.distance(kilometers=new_distance_i / 1000).destination(origin, bearing)
        lon2, lat2 = destination.longitude, destination.latitude
        return Point(lon2, lat2)

    @staticmethod
    def _get_bearing(lat1, lon1, lat2, lon2):
        import math as m
        bearing = m.atan2(m.sin(lon2 - lon1) * m.cos(lat2),
                          m.cos(lat1) * m.sin(lat2) - m.sin(lat1) * m.cos(lat2) * m.cos(lon2 - lon1))
        bearing = m.degrees(bearing)
        bearing = (bearing + 360) % 360
        return bearing

    """ finds the Points ( the visits of a visitor ) inside the area.
     input : all the devices 
     output: the overall visits."""

    @staticmethod
    def find_devices_in_polygon(overall_devices, poly):
        devices_in = pd.DataFrame()
        users_counter = 0
        for i in range(0, len(overall_devices)):
            if overall_devices["new_point"].loc[i].within(poly) is True:  # series object does not have within
                devices_in = devices_in.append(
                    overall_devices.loc[i])  # we use Point class The Point it self is class append to df is wrong
                users_counter += 1
        print("Points in polygon:", users_counter)
        return devices_in

    @staticmethod
    def write_to_txt(df, txt_name, has_header):
        df.reset_index().to_csv(txt_name, sep='\t', header=has_header, index=False, mode="a")

    @staticmethod
    def plot_points_map(polygon, chosen_devices_points, mapName):
        chosen_devices_points.index = range(len(chosen_devices_points.index))
        import matplotlib.pyplot as plt
        import folium

        polygon.plot(color='white', edgecolor='black')
        plt.title("building area", pad=20)
        plt.show()
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
        import webbrowser
        webbrowser.open(mapName)

    """finds the visitors (through their mobile device)
     who appeared inside the area and stores them in txt file.                                                       
     input : devices inside an area,                                                                                 
     output: the visitors as list"""

    @staticmethod
    def group_devices_in_poly(devices_inside_poly):
        devices_inside_poly.index = range(len(devices_inside_poly.index))
        visitors = []
        grouped_df = devices_inside_poly.groupby("hash_id")
        enterFirstTime = True
        for key, item in grouped_df:
            visitors.append(grouped_df.get_group(key))
            if enterFirstTime is True:
                ClassifyDevicesInArea.write_to_txt(grouped_df.get_group(key), 'points_in_area.txt', has_header=True)
            else:
                ClassifyDevicesInArea.write_to_txt(grouped_df.get_group(key), 'points_in_area.txt', has_header=False)
            enterFirstTime = False
        print("Num of visitors inside the area:", len(visitors))
        colm_visitor = pd.Series()
        colm_visitor["visitors"] = visitors
        return colm_visitor

    # @staticmethod
    # def calculate_time_stay(devices_inside_poly):
    #     import math
    #     for time_per_sec in devices_inside_poly["timestamp"]:
    #         math.fsum(time_per_sec)


if __name__ == '__main__':
    dp = ClassifyDevicesInArea(
        absolute_path_devices_and_area="/home/nikoscf/PycharmProjects/VisitorsInArea/paths")
    devices, polygon = dp.read_data()
    devices_points = dp.devices_location_to_points(devices)

    min_distances, devices_points = dp.find_nearest_point_and_distance(polygon.geometry[0], devices_points)
    devices_to_interpolate = dp.add_uncertainty_to_min_distance(devices_points, min_distances)
    new_interpolated_points = dp.interpolate_coordinates(devices_to_interpolate)
    # dp.plot_points_map(polygon, new_interpolated_points["new_point"], "devices_outside.html")
    devices_inside = dp.find_devices_in_polygon(new_interpolated_points, polygon.geometry[0])

    # dp.plot_points_map(polygon, devices_inside["new_point"], "devices_in_area.html")

    visitors = dp.group_devices_in_poly(devices_inside)
    # print(visitors)
    # dp.calculate_time_stay(visitors)
