# to install the packages, on terminal run:   sh -e dependencies

import pandas as pd
import yaml
import geopandas as gpd
from shapely.geometry import Point


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

        overall_devices["new_point_distance"] = new_distance
        return overall_devices


    @staticmethod
    def find_nearest_point_and_distance(polyg, devices_ptns):
        from shapely.ops import nearest_points
        nearest_coordinates = []
        min_distances = []
        for p in devices_ptns['Point']:
            p1, p2 = nearest_points(polyg, p)  # p1 is poly p2 is visitor's position approx for longitude
            nearest_coordinates.append(p1)
            # from geopy import distance # 0.10 millimetre in difference from our build in method <_distance_calc>
            # from geopy import Point # but converts the shapely.geometry.Point to geopy.Point and needs more work
            d = ClassifyDevicesInArea.distance_calc(p1.x, p1.y, p2.x, p2.y)
            min_distances.append(d)
        devices_ptns["nearest_coordinate"] = nearest_coordinates
        return min_distances, devices_ptns


    @staticmethod
    def helper_interpolate_coordinates(devices_to_interpolate):
        new_point = []
        for i in range(0, len(devices_to_interpolate)):
            lat1_i = devices_to_interpolate["Point"].loc[i].y
            lon1_i = devices_to_interpolate["Point"].loc[i].x
            new_distance_i = devices_to_interpolate["new_point_distance"].loc[i]

            lat2_i = devices_to_interpolate["nearest_coordinate"].loc[i].y
            lon2_i = devices_to_interpolate["nearest_coordinate"].loc[i].x

            new_point.append(
                ClassifyDevicesInArea._interpolate_coordinates(lat1_i, lon1_i, lat2_i, lon2_i, new_distance_i))
        devices_to_interpolate["new_interpolated_point"] = new_point
        return devices_to_interpolate


    """Interpolates the old points of devices to new positions.
       Using the  new distance of the points made by <add_uncertainty_to_min_distance>
       input: the new dataframe with the distance
       output: the new interpolated points to insert into <find_devices_in_polygon>
       """

    @staticmethod
    def _interpolate_coordinates(lat1, lon1, lat2, lon2, new_distance_i):

        # print("lat1")
        # print(lat1)
        # print("lon1")
        # print(lon1)
        # print("lat2")
        # print(lat2)
        # print("lon2")
        # print(lon2)

        import math
        bearing = ClassifyDevicesInArea._get_bearing(lat1, lon1, lat2, lon2)
        R = 6378.1  # Radius of the Earth
        brng = math.radians(bearing) # 1.57  Bearing is 90 degrees converted to radians.
        d = new_distance_i #15  Distance in km <!!!!!!!!!!!!

        lat1 = math.radians(lat1)  # Current lat point converted to radians
        lon1 = math.radians(lon1)  # Current long point converted to radians

        lat2 = math.asin(math.sin(lat1) * math.cos(d / R) +
                         math.cos(lat1) * math.sin(d / R) * math.cos(brng))

        lon2 = lon1 + math.atan2(math.sin(brng) * math.sin(d / R) * math.cos(lat1),
                                 math.cos(d / R) - math.sin(lat1) * math.sin(lat2))

        lat2 = math.degrees(lat2)
        lon2 = math.degrees(lon2)

        # print("--lat2--")
        # print("--lon2--")
        #
        # print(lat2)
        # print(lon2)

        return Point(lat2, lon2)

    @staticmethod
    def _get_bearing(lat1, lon1, lat2, lon2):
        import math as m
        bearing = m.atan2(m.sin(lon2 - lon1) * m.cos(lat2),
                          m.cos(lat1) * m.sin(lat2) - m.sin(lat1) * m.cos(lat2) * m.cos(lon2 - lon1))
        bearing = m.degrees(bearing)
        bearing = (bearing + 360) % 360
        return bearing

        # alternative way, but geopy is a heavy lib
        # import geopy
        # from geopy.distance as geod
        #
        # # given: lat1, lon1, b = bearing in degrees, d = distance in kilometers
        #
        # origin = geopy.Point(lat1, lon1)
        # destination = geod(kilometers=new_distance).destination(origin, b)
        # lat2, lon2 = destination.latitude, destination.longitude


    """ finds the Points ( the visits of a visitor ) inside the area.
     input : all the devices 
     output: the overall visits."""

    @staticmethod
    def find_devices_in_polygon(overall_devices, poly):
        devices_in = pd.DataFrame()
        for i in range(0, len(overall_devices)):
            if overall_devices["Point"].loc[i].within(poly) is True:  # series object does not have within
                devices_in = devices_in.append(
                    overall_devices.loc[i])  # # we use Point class The Point it self is class so cannot be ...
        return devices_in  # appended as is in a dataframe. Append all the dataframe


    """finds the visitors (through their mobile device) 
        who were inside the area at least one time.
        input : devices inside an area, 
        output: the visitors """

    @staticmethod
    def group_devices_in_poly(devices_inside_poly):
        global visitors
        grouped_df = devices_inside_poly.groupby("hash_id")
        print("Devices inside the polygon")
        for key, item in grouped_df:
            visitors = grouped_df.get_group(key)
        print("Num of visitors inside the area:", len(visitors))
        return visitors

    @staticmethod
    def write_to_txt(df):
        df.reset_index().to_csv('visitors_in_area.txt', sep='\t', header=True, index=False)

    @staticmethod
    def plot_polygon(gdf):
        import matplotlib.pyplot as plt
        import folium
        gdf.plot(color='white', edgecolor='black')
        plt.show()
        fol_map = folium.Map(location=[59.12381366829599, 14.608602084061102],
                             tiles="OpenStreetMap", zoom_start=20
                             )
        folium.GeoJson(gdf).add_to(fol_map)
        folium.LatLngPopup().add_to(fol_map)
        fol_map.save("polygon.html")
        import webbrowser
        webbrowser.open("polygon.html")


if __name__ == '__main__':
    dp = ClassifyDevicesInArea(
        absolute_path_devices_and_area="/home/nikoscf/PycharmProjects/VisitorsInArea/paths")

    devices, polygon = dp.read_data()
    # plot_polygon(polygon)
    devices_points = dp.devices_location_to_points(devices)

    min_distances, devices_points = dp.find_nearest_point_and_distance(polygon.geometry[0], devices_points)

    # print("min_distances=", min_distances)
    # print("devices_points=", devices_points)
    devices_to_interpolate = dp.add_uncertainty_to_min_distance(devices_points, min_distances)

    # print("devices_to_interpolate = ", devices_to_interpolate)
    new_interpolated_points = dp.helper_interpolate_coordinates(devices_to_interpolate)

    # print("new_interpolated_points=", new_interpolated_points["new_interpolated_point"])

    devices_inside = dp.find_devices_in_polygon(new_interpolated_points, polygon.geometry[0])

    visitors = dp.group_devices_in_poly(devices_inside)
    dp.write_to_txt(visitors)
