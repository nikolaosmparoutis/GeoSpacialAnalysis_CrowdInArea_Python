# IoT: GeoSpacial Analysis (Python)

Tracking the location per second of visitors in a small building area.
1. Considering the uncertainty of the sensors the program finds the potential position of the visitor.
2. The program classifies the geolocation (in lat,long) of each visitor regarding if he/she found to be in the area (a store). 
3. Writes the new data and the previous data with the header to an output file.
4. Finds the number of visitors inside the area and their unique number.
5. Calculate the time stay for each visitor inside the building. 
6. Plots the map with the points and the area, where every visitor was inside and outside the area.
7. Plots the map with the points and the area, where the visitors was inside the area.




=============== Output ===============

Number of visitors inside the building: 8

Visitor 1.0 stayed 00:00:23\
Visitor 2.0 stayed 00:31:44\
Visitor 3.0 stayed 00:00:57\
Visitor 5.0 stayed 00:10:40\
Visitor 6.0 stayed 00:23:30\
Visitor 7.0 stayed 00:02:55\
Visitor 8.0 stayed 00:00:16\
Visitor 9.0 stayed 01:23:15\

---------------------------
Visitor 4 did not entered the building.
