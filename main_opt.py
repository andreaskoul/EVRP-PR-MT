import time
from copy import deepcopy

from VeRyPy.cvrp_io import read_TSPLIB_CVRP
from local_search import *
from local_search.inter_route_operators import *
from local_search.intra_route_operators import *
from parallel_savings import parallel_savings_init, _refine_solution
from plot_2d import *

start_time = time.time()

Instance_path = r"Data/Instances/E-n29-k4-s7.evrp"   # insert instance here
problem = read_TSPLIB_CVRP(Instance_path)
coordinates = problem.coordinate_points
stations = problem.stations
distance_matrix = problem.distance_matrix
cust_demands = problem.customer_demands
N = len(distance_matrix)

# Time constrains
Time_all_min = problem.distance_matrix / 12 / 60
Service_time = 10
Reload_time = 30

# Problem Constrains
Max_L = problem.energy_capacity  # m   / Max Vehicle Mileage
Max_C = problem.capacity_constraint  # items / Max Vehicles Capacity
Max_T = 180  # min  | Max Time per route

# Initial Solution
solution = parallel_savings_init(
    D=distance_matrix,
    d=cust_demands,
    C=Max_C,
    stations=stations,
    L=Max_L,
    T=Time_all_min,
    # Route_time=Max_T,
    s=Service_time,
    minimize_K=True)

# Stations in use
stations_in_use = []
stations_ids = [(i + N - stations) for i in range(stations)]
for i in stations_ids:
    if i in solution:
        stations_in_use.append(i)

# Recharging Properties
p = 22  # kW    \ Charging Power (Average Charger in Athens)
e = 0.2  # kWh/km \ Energy Consumption
c = p / e  # km/hour \ Charging Speed

# Route Optimisation
ls_ops = [do_1point_move, do_2opt_move, do_2point_move, do_3opt_move, do_redistribute_move,
          do_2optstar_move, do_relocate_move, do_chain_move, do_insert_move, do_exchange_move]

# ls_inter = [do_2optstar_move,do_1point_move,do_2point_move,do_insert_move,do_redistribute_move,do_chain_move]
# ls_intra = [do_2opt_move,do_3opt_move,do_relocate_move,do_exchange_move]

rs_sol = do_local_search(ls_ops=[do_worst_rs_remove], sol=solution, D=distance_matrix, d=cust_demands, C=Max_C, L=Max_L,
                         stations_ids=stations_ids, max_iterations=10)

new_sol = do_local_search(ls_ops=ls_ops, sol=rs_sol, D=distance_matrix, d=cust_demands, C=Max_C, L=Max_L,
                          stations_ids=stations_ids, max_iterations=1000)

normal_sol = do_rs_insert(new_sol, distance_matrix, stations_ids, Max_L)   # Add Closest Recharging Stations

refined_sol = _refine_solution(normal_sol, D=distance_matrix, d=cust_demands, C=Max_C, L=Max_L,
                               stations_ids=stations_ids, minimize_K=True)

s1 = objf(solution, distance_matrix)
s2 = objf(rs_sol, distance_matrix)
s3 = objf(new_sol, distance_matrix)
s4 = objf(normal_sol, distance_matrix)
s5 = objf(refined_sol, distance_matrix)

# Vehicle Routing
no_vehicles = problem.vehicles
routes = []
routes_total_t = []
routes_total_m = []
routes_total_d = []
r = 0
route_prev = 0
route_list = []

for route_idx, route in enumerate(sol2routes(refined_sol)):
    routes.append(route)
    demand = totald(route, problem.customer_demands)
    routes_m = objf(route, problem.distance_matrix)
    routes_t = objf(route, Time_all_min) + Service_time * (len(route) - 2) + Reload_time
    r += 1
    if r <= no_vehicles:
        routes_total_t.append(routes_t)
    else:
        charging_time = 60 * route_prev / c
        # Charging time = (km in previous route)/ (charging speed) in minutes
        if charging_time > Reload_time:
            routes_t += charging_time
        else:
            routes_t += Reload_time
        routes_total_t.append(routes_t)
    route_prev = routes_m
    demand = totald(route, problem.customer_demands)
    routes_total_d.append(demand)
    routes_total_m.append(routes_m)
    # print("Route #%d : %s" % (route_idx + 1, route), ", Distance : ", routes_m, "km",
    #       ", Duration : ", round(routes_t / 60, 2), "hours", " Route Demands : ", int(demand))
#
# print("The total distance is : ", sum(routes_total_m), "km")
# print("The total duration is : ", sum(routes_total_t) / 60, "hours")

print("--- %s seconds ---" % (time.time() - start_time))

# Vehicle Route assignment
vk = []
while (len(routes) - no_vehicles) > 0:
    # Find two shortest routes (requiring min time of charging)
    routes2merge = []
    m_sort = deepcopy(routes_total_m)
    m_sort.sort()
    dist_1min = m_sort[0]
    dist_1min_id = routes_total_m.index(dist_1min)
    dist_2min = m_sort[1]
    dist_2min_id = routes_total_m.index(dist_2min)

    # Merge Shortest Routes
    routes2merge.append(routes[dist_1min_id])
    routes2merge.append(routes[dist_2min_id])
    routes[dist_1min_id] = routes2sol(routes2merge)

    # Calculate Charging time (Partial Charging)
    charging_time = 60 * routes_total_m[dist_2min_id] / c

    # Check if there is Capacity Constrain violation & add reload in merged
    routes_total_d[dist_1min_id] += routes_total_d[dist_2min_id]
    if routes_total_d[dist_1min_id] > Max_C:
        routes_total_t[dist_1min_id] += Reload_time

    # Merge Times & Distances
    routes_total_t[dist_1min_id] += routes_total_t[dist_2min_id] + charging_time
    routes_total_m[dist_1min_id] += routes_total_m[dist_2min_id]

    # Delete Merged Route & Update tables
    del (routes[dist_2min_id])
    del (routes_total_t[dist_2min_id])
    del (routes_total_m[dist_2min_id])
    del (routes_total_d[dist_2min_id])

for k in range(no_vehicles):
    vk.append(routes[k])
    print("Vehicle #%d : %s" % (k + 1, routes[k]), ", Distance : ", routes_total_m[k], "km",
          ", Duration : ", round(routes_total_t[k] / 60, 2), "hours", " Route Demands : ", routes_total_d[k])

total_time = sum(routes_total_t) / 60  # in hours
total_distance = sum(routes_total_m)  # in km

print("The total distance is : ", total_distance, "km")
print("The total duration is : ", round(total_time, 2), "hours")

# plotting the points
plot_2d_diff_col(N, refined_sol, stations, coordinates, vk)
