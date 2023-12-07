#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementation of the Clarke and Wright (1964)
savings functions for parallel (as in multiple route) savings heuristic. The
provided algorithm is built to be extendable and many other savings variants
implemented in VeRyPy use the generic parallel savings heuristic implemented 
here.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has minimal dependencies: only numpy and scipy
are needed for reading and preparing the problem instance."""
###############################################################################

from __future__ import division
from __future__ import print_function

import logging
from builtins import range
from logging import log, DEBUG

from VeRyPy.config import CAPACITY_EPSILON as C_EPS
from VeRyPy.config import COST_EPSILON as S_EPS
from VeRyPy.config import TIME_EPSILON as T_EPS
from VeRyPy.util import routes2sol, objf, without_empty_routes
from local_search import LSOPT, do_local_search
from local_search.inter_route_operators import *
from local_search.intra_route_operators import *
from VeRyPy.routedata import RouteData


def clarke_wright_savings_function(D, stations=0):
    N = len(D) + 1 - stations
    n = N - 1
    savings = [None] * int((n * n - n) / 2)
    idx = 0
    for i in range(1, N):
        for j in range(i + 1, N):
            s = D[i, 0] + D[0, j] - D[i, j]
            savings[idx] = (s, -D[i, j], i, j)
            idx += 1
    savings.sort(reverse=True)
    return savings


def parallel_savings_init(D, d, C, stations=None, s=None, T=None, L=None, Route_time=None, minimize_K=False,
                          savings_callback=clarke_wright_savings_function):
    """
    Implementation of the basic savings algorithm / construction heuristic for
    capaciated vehicle routing problems with symmetric distances (see, e.g.
    Clarke-Wright (1964)). This is the parallel route version, aka. best
    feasible merge version, that builds all of the routes in the solution in
    parallel making always the best possible merge (according to varied savings
    criteria, see below).
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.
    
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring a merge that would increase the total distance.
       WARNING: This only works when the solution from the savings algorithm is
       final. With postoptimimization this non-improving move might have still
       led to improved solution.
   
    * optional savings_callback is a function of the signature:
        sorted([(s_11,x_11,i_1,j_1)...(s_ij,x_ij,i,j)...(s_nn,x_nn,n,n) ]) =
            savings_callback(D)
        where the returned (sorted!) list contains savings (that is, how much
        solution cost approximately improves if a route merge with an edge
        (i,j) is made). This should be calculated for each i in {1..n},
        j in {i+1..n}, where n is the number of customers. The x is a secondary
        sorting criterion but otherwise ignored by the savings heuristic.
        The default is to use the Clarke Wright savings criterion.
        
    See clarke_wright_savings.py, gaskell_savings.py, yellow_savings.py etc.
    to find specific savings variants. They all use this implementation to do 
    the basic savings procedure, and they differ only by the savings
    calculation. There is also the sequential_savings.py, which builds the
    routes one by one.


    Clarke, G. and Wright, J. (1964). Scheduling of vehicles from a central
    depot to a number of delivery points. Operations Research, 12, 568-81.
    """

    N = len(D) + 1 - stations
    ignore_negative_savings = not minimize_K
    # 1. make route for each customer
    routes = [[i] for i in range(1, N)]
    route_demands = d[1:] if C else [0] * N
    if L:
        route_costs = [D[0, i] + D[i, 0] for i in range(1, N)]
    if Route_time:
        route_times = [T[0, i] + T[i, 0] for i in range(1, N)]
    try:
        # 2. compute initial savings
        savings = savings_callback(D, stations)

        # zero based node indexing!
        endnode_to_route = [0] + list(range(0, N - 1))

        # 3. merge
        # Get potential merges the best savings first (second element is secondary
        # sorting criterion, and it is ignored)
        for best_saving, _, i, j in savings:
            # if __debug__:
            #     log(DEBUG - 1, "Popped savings s_{%d,%d}=%.2f" % (i, j, best_saving))

            if ignore_negative_savings:
                cw_saving = D[i, 0] + D[0, j] - D[i, j]
                if cw_saving < 0.0:
                    break
            left_route = endnode_to_route[i]
            right_route = endnode_to_route[j]

            # the node is already an internal part of a longer segment
            if ((left_route is None) or
                    (right_route is None) or
                    (left_route == right_route)):
                continue

            # if __debug__:
            #     log(DEBUG - 1, "Route #%d : %s" %
            #         (left_route, str(routes[left_route])))
            #     log(DEBUG - 1, "Route #%d : %s" %
            #         (right_route, str(routes[right_route])))

            # check capacity constraint validity
            if C:
                merged_demand = route_demands[left_route] + route_demands[right_route]
                if merged_demand - C_EPS > C:
                    # if __debug__:
                    #     log(DEBUG - 1, "Reject merge due to " +
                    #         "capacity constraint violation")
                    continue

            # if there is time constraint, check its validity
            if Route_time:
                merged_time = route_times[left_route] - T[0, i] + \
                              route_times[right_route] - T[0, j] + \
                              T[i, j] + s * (len(routes[left_route]) + len(routes[right_route]))
                if merged_time - T_EPS > Route_time:
                    # if __debug__:
                    #     log(DEBUG - 1, "Reject merge due to " +
                    #         "maximum time constraint violation")
                    continue

            # check range constraint validity
            if L:
                merged_cost = route_costs[left_route] - D[0, i] + \
                              route_costs[right_route] - D[0, j] + \
                              D[i, j]
                if merged_cost - S_EPS > L:
                    # Check if there is a charging station nearby
                    dist_no_rs = D[routes[left_route][-1], 0] + D[0, routes[right_route][0]]
                    dist_cl_st = 10000
                    closest_station = 0
                    for k in range(1,stations):
                        st_id = N + k - 1
                        d_st = D[routes[left_route][-1], st_id] + D[st_id, routes[right_route][0]]
                        if d_st < dist_cl_st:
                            dist_cl_st = d_st
                            closest_station = N + k - 1
                    if dist_cl_st < dist_no_rs:
                        # Check capacity constrain
                        if C:
                            merged_demand = route_demands[left_route] + route_demands[right_route]
                            if merged_demand - C_EPS > C:
                                # if __debug__:
                                #     log(DEBUG - 1, "Reject merge due to " +
                                #         "capacity constraint violation")
                                continue

                        # Check time constrain and add charging for km done in left route
                        if Route_time:
                            merged_time = route_times[left_route] - T[0, i] + \
                                          route_times[right_route] - T[0, j] + \
                                          T[i, j] + s * (len(routes[left_route]) + len(routes[right_route])) + \
                                          objf(routes[left_route], T)
                            if merged_time - T_EPS > Route_time:
                                # if __debug__:
                                #     log(DEBUG - 1, "Reject merge due to " +
                                #         "maximum time constraint violation")
                                continue
                        if Route_time: route_times[left_route] = merged_time

                        # update bookkeeping only on the receiving (left) route
                        route_costs[left_route] = merged_cost
                        route_demands[left_route] = merged_demand

                        # reverse the merged routes as necessary
                        if routes[left_route][0] == i:
                            routes[left_route].reverse()
                        if routes[right_route][-1] == j:
                            routes[right_route].reverse()

                        # the nodes that become midroute points cannot be merged
                        if len(routes[left_route]) > 1:
                            endnode_to_route[routes[left_route][-1]] = None
                        if len(routes[right_route]) > 1:
                            endnode_to_route[routes[right_route][0]] = None

                        # all future references to right_route are to merged route
                        endnode_to_route[routes[right_route][-1]] = left_route

                        # merge with list concatenation
                        routes[left_route].extend([closest_station])
                        routes[left_route].extend(routes[right_route])
                        route_costs[left_route] = 0
                        routes[right_route] = None
                        continue
                    else:
                        # if __debug__:
                        #     log(DEBUG - 1, "Reject merge due to " +
                        #         "maximum route length constraint violation")
                        continue

            # update bookkeeping only on the receiving (left) route
            if C: route_demands[left_route] = merged_demand
            if L: route_costs[left_route] = merged_cost
            if Route_time: route_times[left_route] = merged_time

            # merging is done based on the joined endpoints, reverse the 
            # merged routes as necessary
            if routes[left_route][0] == i:
                routes[left_route].reverse()
            if routes[right_route][-1] == j:
                routes[right_route].reverse()

            # the nodes that become midroute points cannot be merged
            if len(routes[left_route]) > 1:
                endnode_to_route[routes[left_route][-1]] = None
            if len(routes[right_route]) > 1:
                endnode_to_route[routes[right_route][0]] = None

            # all future references to right_route are to merged route
            endnode_to_route[routes[right_route][-1]] = left_route

            # merge with list concatenation
            routes[left_route].extend(routes[right_route])
            routes[right_route] = None
            # if __debug__:
            #     dbg_sol = routes2sol(routes)
            #     log(DEBUG - 1, "Merged, resulting solution is %s (%.2f)" %
            #         (str(dbg_sol), objf(dbg_sol, D)))
    except KeyboardInterrupt:  # or SIGINT
        interrupted_sol = routes2sol(routes)
        raise KeyboardInterrupt(interrupted_sol)
    return routes2sol(routes)


def _refine_solution(sol, D, d, C, L, stations_ids, minimize_K):
    # refine until stuck at a local optima
    local_optima_reached = False
    while not local_optima_reached:
        sol = without_empty_routes(sol)
        if not minimize_K:
            sol.append(0)  # make sure there is an empty route to move the pt to

        # improve with relocation and keep 2-optimal
        sol = do_local_search((do_1point_move, do_2opt_move, do_2point_move, do_relocate_move), sol, D, d, C, L,
                              stations_ids, max_iterations=1000, iteration_strategy=LSOPT.BEST_ACCEPT)

        # try to redistribute the route with the smallest demand
        sol = without_empty_routes(sol)
        routes = RouteData.from_solution(sol, D, d)
        min_rd = min(routes, key=lambda rd: rd.demand)
        routes.remove(min_rd)

        if not minimize_K:
            routes.append(RouteData())

        if __debug__:
            log(DEBUG, "Applying do_redistribute_move on %s (%.2f)" % (str(sol), objf(sol, D)))

        redistribute_result = do_redistribute_move(min_rd, routes,
                                                   D, stations_ids, L, d, C, strategy=LSOPT.FIRST_ACCEPT,
                                                   # Note: Mole and Jameson do not specify exactly
                                                   # how the redistribution is done (how many
                                                   # different combinations are tried).
                                                   # Increase the recombination_level if the
                                                   # for more agressive and time consuming search
                                                   # for redistributing the customers on other
                                                   # routes.
                                                   recombination_level=0)
        redistribute_delta = redistribute_result[-1]

        if (redistribute_delta is not None) and \
                (minimize_K or redistribute_delta < 0.0):

            updated_sol = RouteData.to_solution(redistribute_result[:-1])
            if __debug__:
                log(DEBUG - 1, ("Improved from %s (%.2f) to %s (%.2f)" %
                                (sol, objf(sol, D), updated_sol, objf(updated_sol, D)))
                    + "using inter route heuristic do_redistribute_move\n")
            sol = updated_sol
        else:
            local_optima_reached = True
            if __debug__:
                log(DEBUG - 1, "No move with do_redistribute_move\n")
    return sol


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ps_algorithm():
    def call_init(D, d, C, L, minimize_K):
        return parallel_savings_init(D, d, C, L, minimize_K)

    return call_init


if __name__ == "__main__":
    from VeRyPy.shared_cli import cli

    cli(*get_ps_algorithm())
