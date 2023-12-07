# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and implements intra route (i.e. within one route) local 
search improvement heuristics such as 2-opt, one-point-move etc. All operators 
assume we start from a feasible solution. Also, all functions implementing the
operations have the following signature:

do_X_move(route, D, strategy, best_delta), where

* route is a list of nodes forming a route.The route must start and end to the
    depot (index 0), e.g. [0,1,2,3,0].
* D is the symmetric full distance matrix, given as numpy-like 2darray.
    Thus, D[i,j] gives the distance between nodes i and j. 
* strategy is either FIRST_ACCEPT (default) or BEST_ACCEPT. 
    First accept returns a modified route as soon as the first improvement is  
    encoutered. Best accept tries all possible combinations and returns the
    best one.
* best_delta is the required level of change in the in route cost. It can be
   used to  set the upper bound (requirement) for the improvement. Usually it 
   is None, which sets the level to 0.0. In that case only improving deltas are
   accepted. If set to (large) positive value, best worsening move can also be
   returned.
   
All intra route improvement operators return the new improved route and the
improvement (delta) as 2-tuple or (None,None) if no improvement was found."""
###############################################################################

from __future__ import division
# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function

import copy
from builtins import range

import VeRyPy.util
from VeRyPy.util import *
from VeRyPy.config import COST_EPSILON as S_EPS
from VeRyPy.routedata import RouteData
from local_search import LSOPT


def do_2opt_move(route, D, stations_ids=None, L=None, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """ 2-opt local search operation for the symmetric distances D
    Remove 2 edges from the route and check if swapping then endpoints
    would yield an improvement.
    """
    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(0, rN - 1):
        for j in range(i + 1, rN - 1):
            if route[j] in stations_ids:
                continue
            a = route[i]
            b = route[i + 1]
            c = route[j]
            d = route[j + 1]

            # a->c b->d( 2-opt move )
            #        _____
            #       /     \
            # >--a  b-<-c  -->d
            #     \____/ 

            # For edges a-b and c-d, try if savings could be made by traveling 
            #  from a-c and b-d (effectively reversing the chain from
            #  b to c).
            delta = D[a, c] + D[b, d] \
                    - D[a, b] - D[c, d]
            if delta + S_EPS < best_delta:
                best_move = (i, j)
                best_delta = delta
                if strategy == LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break  # j loop
        if accept_move:
            break  # i loop

    if best_move:
        i, j = best_move
        # print("REMOVEME:","best_move", i,j, route, best_delta)
        return route[:i + 1] + route[j:i:-1] + route[j + 1:], best_delta
    return None, None


def do_3opt_move(route, D, stations_ids=None, L=None, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """ 3-opt local search operation for the symmetric distances D """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False

    for i in range(0, rN - 1):
        for j in range(i + 1, rN - 1):
            for k in range(j + 1, rN - 1):
                # the edge endpoints
                a = route[i]
                b = route[i + 1]
                c = route[j]
                d = route[j + 1]
                e = route[k]
                f = route[k + 1]
                # print("search abcdef", a,b,c,d,e,f)

                # After removing edges a-b, c-d, and e-f, try all of the 7
                # combinations in which the segments can be reconnected.
                removed_weights = D[a, b] + D[c, d] + D[e, f]

                # The combinations could be iterated, but for simplicity
                #  (and speed), the loop has been unrolled below.

                ## 2-opt moves

                # #a->b# c->e d->f ( 2-opt move )
                #               _____
                #              /     \
                # >-a--b->-c  d-<-e  f-->
                #          \_____/  
                delta = (D[a, b] + D[c, e] + D[d, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (i + 1, j + 1, 1),
                                 (k, j, -1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                # a->c b->d #e->f#  ( 2-opt move )
                #        _____
                #       /     \
                # >--a  b-<-c  d->-e--f-->
                #     \____/ 
                delta == (D[a, c] + D[b, d] + D[e, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (j, i, -1),
                                 (j + 1, k + 1, 1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                # a->e #d->c# b->f  (2-opt move)
                #        ____________
                #       /            \
                # >--a  b-<-c--d-<-e  f-->
                #     \___________/ 
                delta = (D[a, e] + D[d, c] + D[b, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (k, j, -1),
                                 (j, i, -1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                ## 3-opt moves

                # a->c b->e d->f  ( 3-opt move )
                #         _________
                #        /         \
                # >--a  b-<-c  d-<-e  f-->
                #     \____/    \____/
                delta = (D[a, c] + D[b, e] + D[d, f]) - removed_weights;
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (j, i, -1),
                                 (k, j, -1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                # a->d e->b c->f  (3-opt move)
                #         __________
                #        /          \
                # >--a  b->-c   d->-e  f-->
                #     \______\_/      /
                #             \______/
                delta = (D[a, d] + D[e, b] + D[c, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (j + 1, k + 1, 1),
                                 (i + 1, j + 1, 1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                # a->d e->c b->f  (3-opt move)
                #          __________
                #         /   _____  \
                #        /   /     \  \  
                # >--a  b-<-c  d->-e  f-->
                #     \_______/ 
                delta = (D[a, d] + D[e, c] + D[b, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (j + 1, k + 1, 1),
                                 (j, i, -1), (k + 1, None, 1))
                    best_delta = delta

                    if strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop

                # a->e b->d c->f  (3-opt move)
                #             _______
                #            /       \
                # >--a  b->-c  d-<-e  f-->
                #     \  \____/   / 
                #      \_________/  
                delta = (D[a, e] + D[d, b] + D[c, f]) - removed_weights
                if delta + S_EPS < best_delta:
                    best_move = ((None, i + 1, 1), (k, j, -1),
                                 (i + 1, j + 1, 1), (k + 1, None, 1))
                    best_delta = delta

                    if best_move and strategy == LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break  # k loop
            if accept_move:
                break  # j loop
        if accept_move:
            break  # i loop

    if best_move:
        sgmt1, sgmt2, sgmt3, sgmt4 = best_move

        # the head containing leaving the depot +
        #  the first segment, reversed or not +
        #  the second segment, reversed or not +
        #  the tail containing return to the depot
        return route[sgmt1[0]:sgmt1[1]:sgmt1[2]] + \
               route[sgmt2[0]:sgmt2[1]:sgmt2[2]] + \
               route[sgmt3[0]:sgmt3[1]:sgmt3[2]] + \
               route[sgmt4[0]:sgmt4[1]:sgmt4[2]], best_delta

    return None, None


def do_relocate_move(route, D, stations_ids=None, L=None, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """Relocate local search operation for the symmetric distances D.
    Check if a node on the route can be moved to another position on the same
    route.
    
    Please note that relocate search space is a subset of 3-opt.
    However, the search space of 3-opt is larger.
    """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(1, rN - 1):
        for j in range(1, rN):
            if i == j or j == i - 1:
                continue

            a = route[i - 1]
            b = route[i]
            c = route[i + 1]

            d = route[j - 1]
            e = route[j]

            # check the no op position
            if d == b or e == b:
                continue

            # move b from between a and c to between d and e
            delta = -D[a, b] - D[b, c] + D[a, c] \
                    - D[d, e] + D[d, b] + D[b, e]

            if delta + S_EPS < best_delta:
                best_move = (i, j)
                best_delta = delta
                if strategy == LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break  # j loop
        if accept_move:
            break  # i loop

    if best_move:
        i, j = best_move
        if i < j:
            return route[:i] + route[i + 1:j] + [route[i]] + route[j:], best_delta
        else:
            return route[:j] + [route[i]] + route[j:i] + route[i + 1:], best_delta

    return None, None


def do_exchange_move(route, D, stations_ids=None, L=None, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """(Node) exchange local search operation for the symmetric distances D.
    Checks if two nodes on the route can be swapped. 
    
    Please note that exchange search space is a subset of 4-opt.
    However, the search space of 4-opt is significantly larger.
    """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(1, rN - 1):
        for j in range(i + 1, rN - 1):
            if i == j:
                continue
            a = route[i - 1]
            b = route[i]
            c = route[i + 1]

            d = route[j - 1]
            e = route[j]
            f = route[j + 1]

            if c == e:
                delta = -D[a, b] - D[b, e] - D[e, f] \
                        + D[a, e] + D[e, b] + D[b, f]
            else:
                # swap b and e from between a and c to between d and f
                delta = -D[a, b] - D[b, c] + D[a, e] + D[e, c] \
                        - D[d, e] - D[e, f] + D[d, b] + D[b, f]

            # print("REMOVEME:", i,j, "(", delta, ")", "best =", best_delta)
            if delta + S_EPS < best_delta:
                best_move = (i, j)
                best_delta = delta

                if strategy == LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break  # j loop
        if accept_move:
            break  # i loop

    if best_move:
        i, j = best_move
        return route[:i] + [route[j]] + route[i + 1:j] + [route[i]] + route[j + 1:], best_delta
    return None, None


def do_worst_rs_remove(route, D, stations_ids=None, L=None, strategy=LSOPT.BEST_ACCEPT, best_delta=None):
    st_inroute = []
    vehicle_charge = []

    initial_delta = VeRyPy.util.objf(route, D)
    if not best_delta:
        best_delta = 0

    for k in stations_ids:
        if k in route:
            tempk = route.index(k)
            st_inroute.append(tempk)
            # print("route.index(k)", route.index(k))
            vehicle_charge.append(VeRyPy.util.objf(route[0:tempk], D))
    if not st_inroute:
        return None, None
    else:
        min_ch = min(vehicle_charge)
        min_ch_id = vehicle_charge.index(min_ch)
        min_ch_st_id = st_inroute[min_ch_id]
        del (route[min_ch_st_id])
        best_delta = VeRyPy.util.objf(route, D) - initial_delta
        return route, best_delta


# def do_rs_insert_move(route, D, stations_ids, L, d=None, C=None,
#                       best_delta=None):
#     """ Try to insert (unrouted) recharging station(s) on the existing receiving route.
#     #The
#    # operation succeeds only if all customers can be inserted on the recieving
#    # route. The difference to the one point move is that there is no delta for
#    # removing the unrouted node from its RouteData / list.
#
#     Returns a 3-tuple: An empty RouteData, A copy of the updated RouteData and
#      the route cost (usually route length) delta; or (None,None,None) if the
#      insertion fails.
#
#     * unrouted
#        Can be a single station (int), a list of those, or a route (RouteData).
#        If RouteData is given its .demand field must be set.
#     * recieving_route_data
#        is the route (RouteData) that the customers to insert are inserted to.
#        Its .demand field must be set.
#     * D,d,C,L
#        the problem definition and constraints
#     * strategy
#        is the applied loca search strategy to use. Please note that if there is
#        more than one customer to insert, the maximum route cost constraint L
#        is set, and strategy is set to first accept, the best insertion position
#        for the customers are still searched for to leave the most amount of
#        room for subsequent insertions.
#     """
#     # test
#     initial_route_cost = VeRyPy.util.objf(route, D)
#     total_demand = VeRyPy.util.totald(route, d)
#     if not best_delta:
#         best_delta = 0
#
#     best_node = None
#     route_td = VeRyPy.util.objf(route, D)
#     station_dist = []
#     for i in range(1, len(route) - 1):
#         route_disti = VeRyPy.util.objf(route[0:i + 1], D)
#         insert_after = route[i - 1]
#         insert_before = route[i]
#         for k in stations_ids:
#             delta = D[insert_after, k - 1] + D[k - 1, insert_before] - D[insert_after, insert_before]
#             station_dist.append(delta)
#             if k in route[1:i]:
#                 # rs_index = route[1:i].index(k)
#                 route_disti = 0  # -= VeRyPy.util.objf(route[0:rs_index], D)
#         if route_disti > L:  # if distance is violated add rs
#             best_insert_pos = i
#             # find closest recharging station
#             closest_rs_dist = min(station_dist)
#             closest_rs = station_dist.index(closest_rs_dist)
#
#             insert_delta = closest_rs_dist
#             remove_delta = route_disti + D[insert_after, k - 1]
#
#             # print("remove_delta",remove_delta)
#             route_td += insert_delta
#             if route_td - remove_delta <= L:  # Check if full charge can complete the following route
#                 if initial_route_cost > route_td:
#                     best_insert_pos, best_node = i, closest_rs
#                     best_insert_delta = insert_delta
#                     best_route_dist = route_td
#         if best_node is None:
#             print("none")
#             continue
#
#         elif best_node is not None:
#             initial_route_cost = best_route_dist
#             # we found a valid insertion location accept it
#             route = route[:best_insert_pos] \
#                     + [best_node] \
#                     + route[best_insert_pos:]
#             best_delta += best_insert_delta
#             return route, best_delta
#     return None, None


def do_rs_insert(solution, D, stations_ids, L):
    routes = sol2routes(solution)
    totalR = []
    new_routes = []
    for R in routes:
        totalR.append(objf(R, D))
        dist_i = 0
        for i in range(1, len(R)):
            dist_i += objf(R[i - 1:i + 1], D)
            insert_after = R[i - 1]
            insert_before = R[i]

            # Check if there are stations
            for s in stations_ids:
                if s in R:
                    s_id = R.index(s)
                    dist_i -= objf(R[0:s_id + 1], D)

            # Calculate Distance off hypothetically added station
            station_dist = []
            for s in stations_ids:
                s_delta = D[insert_after, s] \
                          + D[s, insert_before] \
                          - D[insert_after, insert_before]
                station_dist.append(s_delta)

            if dist_i > L:
                # Check the node at which there is no charge
                s_index = station_dist.index(min(station_dist))
                closest_st = stations_ids[s_index]
                new_tolald = objf(R,D) + int(min(station_dist))

                # Check the previous node as well
                prev_node = i - 1
                insert_after = R[prev_node - 1]
                insert_before = R[prev_node]
                station_dist = []
                for s in stations_ids:
                    s_delta = D[insert_after, s] + D[s, insert_before] - D[insert_after, insert_before]
                    station_dist.append(s_delta)
                station_dist_prev = min(station_dist)
                s_index_prev = station_dist.index(station_dist_prev)
                closest_st_prev = stations_ids[s_index_prev]
                new_tolald_prev = objf(R,D) + station_dist_prev

                if new_tolald_prev > new_tolald:
                    R = R[:i] + [closest_st] + R[i:]

                elif new_tolald_prev < new_tolald:
                    R = R[:i - 1] + [closest_st_prev] + R[i - 1:]
        new_routes.append(R)
    new_sol = routes2sol(new_routes)
    return new_sol


# def do_C_constrain_fix(solution, D, C, d):
#     # ###
#     routes = sol2routes(solution)
#     totald_R = []
#     best_remove = []
#     for R in routes:
#         totald_R.append(totald(R, d))
#         cap_i = 0  # Check capacity
#         for i in range(1, len(R)-1):
#             cap_i += d[R[i]]
#             print(cap_i)
#             print([R[i]])
#             remove_i = R[i]
#         # Check if there is visit to depot
#             if cap_i > C:  # Check if we have reached load capacity
#                 print("MPHKAME")
#                 dist_else = []
#                 total_dist_fix = []
#                 other_route = []
#                 other_route_option = []
#                 for r in routes:  # Check routes
#                     new_routes = copy.deepcopy(routes)
#                     if totald(r,d)+d[remove_i] < C:  # Check routes that can take the extra node
#                         print(True)
#                         for k in r:
#                             dist_else.append(D[remove_i, k])
#                         # dist_else.append(D[remove_i, k] for k in r)
#                         print("dist_else",dist_else)
#                         min_else = min(dist_else)  # Find the best point to enter the node
#                         print("min_else",min_else)
#                         templist =D[remove_i, :].tolist()
#                         other_node = templist.index(min_else)
#                         print("other_node",other_node)
#
#                         other_route_option.append(r)
#                         updated_route = r[:other_node] + [remove_i] + r[other_node:]
#                         other_route.append(updated_route)
#
#                         new_routes[routes.index(r)] = updated_route
#                         total_dist_fix.append(objf(routes2sol(new_routes), D))
#                     else:
#                         print(False)
#                         continue
#                 min_fix = min(total_dist_fix)
#                 min_fix_id = total_dist_fix.index(min_fix)
#                 min_fix_route = other_route[min_fix_id]
#                 routes[routes.index(other_route_option[min_fix_id])]=min_fix_route
#                 routes[routes.index(R)]=R[:remove_i]+R[remove_i+1:]
#                 print("routes[routes.index(R)]",routes[routes.index(R)])
#                 # Check the closest node that does not belong to the route
#             if  R[i] == 0:
#                 # s_id = R.index(s)
#                 cap_i = 0
#     new_sol = routes2sol(routes)
#     return new_sol
