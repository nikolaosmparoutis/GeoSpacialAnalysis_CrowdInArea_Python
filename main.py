# to install the packages, on terminal run:   sh -e dependencies

import pandas as pd
import yaml
import geopandas as gpd
from shapely.geometry import Point

absolute_path = "/home/nikoscf/PycharmProjects/VisitorsInArea/paths"


def read_data():
    try:
        with open(absolute_path, 'r') as yml_file:
            cfg = yaml.load(yml_file, Loader=yaml.FullLoader)

            devices_df = pd.read_csv(cfg["path_csv"]["devices"],
                                     sep=",", encoding='latin-1', error_bad_lines=False, warn_bad_lines=False,
                                     low_memory=False, memory_map=True)

            polygon_geodf = gpd.GeoDataFrame(gpd.read_file(cfg["path_csv"]["polygon"]))

        return devices_df, polygon_geodf
    except:
        raise FileNotFoundError


# optional for future use
# def _swap_log_lat_polygon(geodf):
#     import shapely.ops
#     trans_coords = lambda Polygon: shapely.ops.transform(lambda x, y: (y, x), geom=Polygon)
#     polygon_df = geodf.geometry.map(trans_coords)
#     return polygon_df

"""Convert the coordinates to Point
input: the overall devices
output: the Points in a DataFrame"""


def _devices_location_to_points(overall_devices):
    devices_loc = overall_devices[["hash_id", "latitude", "longitude"]].astype(float)

    points_of_devices = pd.DataFrame(columns=["hash_id"])
    hash_ids = []
    P = []

    for hash_id, long, lat in devices_loc.values:
        hash_id = int(hash_id)
        hash_ids.append(hash_id)
        P.append(Point(lat, long))  # Point class by default  reverses its coordinates

    points_of_devices["hash_id"] = hash_ids
    points_of_devices["Point"] = P
    return points_of_devices


""" finds if the Points ( the visits of a visitor ) inside the area.
 input : all the devices 
 output: the overall visits."""


def _find_nearest_point_and_distance(polygon, devices_inside):
    global p1, p2
    from shapely.ops import nearest_points
    distances = []
    for p in devices_inside['Point']:
        p1, p2 = nearest_points(polygon, p)  # p1 is poly p2 is visitor's position approx for longitude

        # from geopy import distance # 0.10 millimetre in difference from our build in method <_distance_calc>
        # from geopy import Point # but converts the shapely.geometry.Point to geopy.Point and needs more work
        distances = _distance_calc(p1.x, p1.y, p2.x, p2.y)
    return distances


def _distance_calc(lat1, lon1, lat2, lon2):
    import math
    R = 6378.137  # Radius of earth in KM

    dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180
    dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180

    a = math.sin(dLat / 2) * math.sin(dLat / 2) + \
        math.cos(lat1 * math.pi / 180) * math.cos(lat2 * math.pi / 180) * \
        math.sin(dLon / 2) * math.sin(dLon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c
    return d * 1000  # meters


def find_devices_in_polygon(overall_devices, poly):
    devices_in = pd.DataFrame()
    for i in range(0, len(overall_devices)):
        if overall_devices["Point"].loc[i].within(poly) is True:  # series object does not have within
            devices_in = devices_in.append(  # we use Point class
                overall_devices.loc[i])  # The Point it self is class so cannot be ...
    return devices_in  # appended as is in a dataframe. Append all the dataframe


"""finds the visitors who was inside the area at least one time.
input : devices inside an area, 
output: the visitors """


def group_devices_in_poly(devices_inside):
    global visitors
    grouped_df = devices_inside.groupby("hash_id")
    print("Devices inside the polygon")
    for key, item in grouped_df:
        visitors = grouped_df.get_group(key)
    print("Num of visitors inside the area:", len(visitors))
    return visitors


def write_to_txt(df):
    df.reset_index().to_csv('visitors_in_area.txt', sep='\t', header=True, index=False)


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
    devices, polygon = read_data()
    # plot_polygon(polygon)
    devices_points = _devices_location_to_points(devices)
    distance = _find_nearest_point_and_distance(polygon.geometry[0], devices_points)
    # distance = _distance_point_from_visitor(p1, p2)
    # devices_inside = find_devices_in_polygon(devices_points, polygon.geometry[0])
    # visitors = group_devices_in_poly(devices_inside)
    # write_to_txt(visitors)
