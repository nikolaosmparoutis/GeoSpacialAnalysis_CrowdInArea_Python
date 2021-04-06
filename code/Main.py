from FindDevicesInArea import ClassifyDevicesInArea
import matplotlib.pyplot as plt

if __name__ == '__main__':

    dp = ClassifyDevicesInArea(absolute_path_devices_and_area="/home/nikoscf/PycharmProjects/VisitorsInArea/code/paths")
    devices, polygon = dp.read_data()


    def _ploting(polygon,new_interpolated_points,devices_inside):
        try:
            dp.plot_building_area(polygon)
            dp.plot_points_map(polygon, new_interpolated_points["new_point"], "devices_outside.html")
            dp.plot_points_map(polygon, devices_inside["new_point"], "devices_in_area.html")
        except BaseException:  # to handle unknown chrome browser exception
            pass

    def main():
        devices_points = dp.devices_location_to_points(devices)
        min_distances, devices_points = dp.find_nearest_point_and_distance(polygon.geometry[0], devices_points)
        devices_to_interpolate = dp.add_uncertainty_to_min_distance(devices_points, min_distances)
        new_interpolated_points = dp.interpolate_coordinates(devices_to_interpolate)
        devices_inside = dp.find_devices_in_polygon(new_interpolated_points, polygon.geometry[0])
        visits = dp.group_devices_in_poly(devices_inside)
        dp.calculate_time_stay(visits)
        _ploting(polygon, new_interpolated_points, devices_inside)
        plt.ioff()
        plt.show()


main()