'''''
import csv
import googlemaps as gmaps
import numpy as np
import earth_points

API_key = 'AIzaSyDgprA-7ugcgDioCFG6CPO9_Uk5kgZETaE'  # enter the key you got from Google. I removed mine here
gmaps = gmaps.Client(key=API_key)

addresses = earth_points.get_addresses('Data/random100.csv')
markers = ["color:blue|size:mid|label:" + chr(65 + i) + "|"
           + r for i, r in enumerate(addresses)]

# path = "color:0x0000ff|weight:2"
# for k in enumerate(data):
#     path += "|" + str(k[1])

result_map = gmaps.static_map(
    center="Athens, Greece",
    scale=2,
    zoom=11,
    size=[1280, 1280],
    format="jpg",
    maptype="roadmap",
    markers=markers,
    path="color:0x0000ff|weight:2|" + "|".join(addresses))

with open("driving_route_map.jpg", "wb") as img:
    for chunk in result_map:
        img.write(chunk)

# max 15 points, undone
'''''
import matplotlib.pyplot as plt
import random

# Stations in use


def plot_2d(N, solution, stations, coordinates):
    stations_in_use = []
    stations_points_in_use = []
    stations_ids = [(i + N - stations) for i in range(stations)]
    for i in stations_ids:
        if i in solution:
            stations_in_use.append(i)
            stations_points_in_use.append(coordinates[i])
    nodes_in_use = list(range(N - stations + 1))
    nodes_in_use.extend(stations_in_use)

    sol_coord = []
    sol_x = []
    sol_y = []

    st_x = []
    st_y = []

    for i in solution:
        sol_coord.append(coordinates[i])
        sol_x.append(coordinates[i][0])
        sol_y.append(coordinates[i][1])

    for i in stations_in_use:
        st_x.append(coordinates[i][0])
        st_y.append(coordinates[i][1])

    fig, ax = plt.subplots(sharex=True, sharey=True)

    plt.scatter(st_x, st_y, marker='*', s=350, c='red')

    plt.plot(sol_x, sol_y, color='green', linestyle='dashed', linewidth=3,
             marker='o', markerfacecolor='blue', markersize=9)
    for i in solution:
        plt.annotate('(%s)' % i, xy=coordinates[i])

    # setting x and y axis range
    plt.ylim(min(sol_y) - 10, max(sol_y) + 10)
    plt.xlim(min(sol_x) - 10, max(sol_x) + 10)

    # naming the x-axis
    plt.xlabel('x - axis')
    # naming the y-axis
    plt.ylabel('y - axis')

    # giving a title to my graph
    plt.title('Routes')

    # function to show the plot
    plt.tight_layout()
    plt.show()


def plot_2d_diff_col(N, solution, stations, coordinates, vk):
    stations_in_use = []
    stations_points_in_use = []
    stations_ids = [(i + N - stations) for i in range(stations)]
    for i in stations_ids:
        if i in solution:
            stations_in_use.append(i)
            stations_points_in_use.append(coordinates[i])
    nodes_in_use = list(range(N - stations + 1))
    nodes_in_use.extend(stations_in_use)

    sol_coord = []
    sol_x = []
    sol_y = []

    st_x = []
    st_y = []

    vk_x = []
    vk_xx = []
    vk_y = []
    vk_yy = []

    for i in solution:
        sol_coord.append(coordinates[i])
        sol_x.append(coordinates[i][0])
        sol_y.append(coordinates[i][1])

    for i in stations_in_use:
        st_x.append(coordinates[i][0])
        st_y.append(coordinates[i][1])

    for i in vk:
        for k in i:
            vk_x.append(coordinates[k][0])
            vk_y.append(coordinates[k][1])
        vk_xx.append(vk_x)
        vk_yy.append(vk_y)
        vk_x = []
        vk_y = []

    fig, ax = plt.subplots(sharex=True, sharey=True)

    plt.scatter(st_x, st_y, marker='*', s=350, c='black')

    no_of_colours = len(vk)
    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
              for j in range(no_of_colours)]

    colour_list = ['blue', 'orange', 'purple', 'green', 'cyan', 'black', 'brown', 'pink', 'yellow', 'olive']

    for i in range(len(vk_xx)):
        plt.plot(vk_xx[i], vk_yy[i], color=colors[i], linestyle='dashed', linewidth=3,
                 marker='o', markerfacecolor='red', markersize=9)

    for i in solution:
        plt.annotate('(%s)' % i, xy=coordinates[i])

    # setting x and y axis range
    plt.ylim(min(sol_y) - 10, max(sol_y) + 10)
    plt.xlim(min(sol_x) - 10, max(sol_x) + 10)

    # naming the x-axis
    plt.xlabel('x - axis')
    # naming the y-axis
    plt.ylabel('y - axis')

    # giving a title to my graph
    plt.title('Routes')

    # function to show the plot
    plt.tight_layout()
    plt.show()