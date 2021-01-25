# to install the packages, on terminal run:   sh -e dependencies

import pandas as pd
import yaml
import geopandas as gpd
from shapely.geometry import Point


# noinspection SpellCheckingInspection
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

    # @staticmethod
    # def devices_location_to_points(overall_devices):
    #     devices_loc = overall_devices[["hash_id", "latitude", "longitude", "timestamp", "uncertainty"]] \
    #         .astype(float)
    #     P = []
    #     for lat, long in devices_loc[["latitude", "longitude"]].values:
    #         P.append(Point(lat, long))  # Point class by default  reverses its coordinates
    #
        # devices_loc["latitude"] = P
        # devices_loc.rename(columns={"latitude": "Point"}, inplace=True)
        # devices_loc.drop(columns=['longitude'], inplace=True)
    #     return devices_loc

    # @staticmethod
    # def devices_location_to_points(overall_devices):
    #     devices_loc = overall_devices[["hash_id", "latitude", "longitude", "timestamp", "uncertainty"]] \
    #         .astype(float)
    #     P = []
    #
    #     for lat, long in devices_loc[["latitude", "longitude"]].values:
    #         P.append(Point(lat, long))  # Point class by default  reverses its coordinates
    #
    #     devices_loc["latitude"] = P
    #     devices_loc.rename(columns={"latitude": "Point"}, inplace=True)
    #     devices_loc.drop(columns=['longitude'], inplace=True)
    #     return devices_loc

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

    @staticmethod
    def find_nearest_point_and_distance(polyg, devices_ptns):
        from shapely.ops import nearest_points
        min_distances = []
        for p in devices_ptns['Point']:
            p1, p2 = nearest_points(polyg, p)  # p1 is poly p2 is visitor's position approx for longitude
            # print(p1,p2)
            # from geopy import distance # 0.10 millimetre in difference from our build in method <_distance_calc>
            # from geopy import Point # but converts the shapely.geometry.Point to geopy.Point and needs more work
            d = ClassifyDevicesInArea.distance_calc(p1.x, p1.y, p2.x, p2.y)
            min_distances.append(d)
        return min_distances

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

    """Interpolates the old points of devices to new positions.
    Using the  new distance of the points made by <add_uncertainty_to_min_distance>
    input: the new dataframe with the distance
    output: the new interpolated points to insert into <find_devices_in_polygon>
    """

    @staticmethod
    def interpolate_new_points(points):
        pass
        # return interpolated_points

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
    min_distances = dp.find_nearest_point_and_distance(polygon.geometry[0], devices_points)
    print("min_distances", min_distances)
    devices_to_interpolate = dp.add_uncertainty_to_min_distance(devices_points, min_distances)
    print("devices_to_interpolate = ", devices_to_interpolate)
    new_interpolated_positions = dp.interpolate_new_points(devices_to_interpolate)

    devices_inside = dp.find_devices_in_polygon(devices_points, polygon.geometry[0])
    print("devices_inside = ", devices_inside.columns)

    visitors = dp.group_devices_in_poly(devices_inside)
    # dp.write_to_txt(visitors)
