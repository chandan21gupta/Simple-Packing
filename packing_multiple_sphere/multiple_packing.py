import numpy as np
import math
from random import randint

def shift_centers(multi_sphere, optimal_result):
    pdb_id = optimal_result['pdb_id']
    shifted_centers_x = optimal_result['x']
    shifted_centers_y = optimal_result['y']
    shifted_centers_z = optimal_result['z']

    x = []
    y = []
    z = []

    radius_list = []  

    for ii in range(len(pdb_id)):
        for keys in multi_sphere[pdb_id[ii]]:
            np.array(multi_sphere[pdb_id[ii]][keys]['sphere_center'])
            shifted_centers = np.array([shifted_centers_x[ii], shifted_centers_y[ii], shifted_centers_z[ii]])
            #print("protein "+pdb_id[ii]+" : "+str(shifted_centers[0])+" "+str(shifted_centers[1])+" "+str(shifted_centers[2]))
            #print("cluster "+keys+" : "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][0])+" "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][1])+" "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][2]))
            multi_sphere[pdb_id[ii]][keys]['sphere_center'] = multi_sphere[pdb_id[ii]][keys]['sphere_center'] - shifted_centers
            #print("cluster "+keys+" : "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][0])+" "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][1])+" "+str(multi_sphere[pdb_id[ii]][keys]['sphere_center'][2]))

# def shift_centers(multi_sphere, optimal_result):
#   # shifted_centers_x = optimal_result['x']
#  #    shifted_centers_y = optimal_result['y']
#  #    shifted_centers_z = optimal_result['z']
#   ''' Shifting the centers of the multiple spheres with respect to the original single sphere '''
#   pdb_id = optimal_result['pdb_id']
#     # for i in range(len(shifted_centers_x)):
#     #     print('protein '+i+' center - '+str(shifted_centers_x[i])+' '+str(shifted_centers_y[i])+' '+str(shifted_centers_z[i]))
#     shifted_centers_x = optimal_result['x']
#     shifted_centers_y = optimal_result['y']
#     shifted_centers_z = optimal_result['z']

#     x = []
#     y = []
#     z = []                    

#     radius_list = []

#     for ii in range(len(pdb_id)):
#         for keys in multi_sphere[pdb_id[ii]]:
#             np.array(multi_sphere[pdb_id[ii]][keys]['sphere_center'])
#             shifted_centers = np.array([shifted_centers_x[ii], shifted_centers_y[ii], shifted_centers_z[ii]])
#             #print(str(shifted_centers[0])+' '+str(shifted_centers[1])+' '+str(shifted_centers[2]))
#             print(shifted_centers)
#             #print(dict[keys]['sphere_center'])
#             multi_sphere[pdb_id[ii]]['sphere_center'] = multi_sphere[pdb_id[ii]]['sphere_center'] - shifted_centers
#     #         print(dict[keys]['sphere_center'])
#     #         x.append(dict[keys]['sphere_center'][0])
#     #         y.append(dict[keys]['sphere_center'][1])
#     #         z.append(dict[keys]['sphere_center'][2])
#     #         radius_list.append(dict[keys]['sphere_radius'])

#     # x_loc = np.array(x)
#     # y_loc = np.array(y)
#     # z_loc = np.array(z)

#     # center_list = [x_loc, y_loc, z_loc]

def overlap_detection_multiple(pdb_id, dict):
    current_clusters = dict[pdb_id]
    for keys in dict:
        if keys != pdb_id:
            clusters = dict[keys]
            for i in range(len(current_clusters)):
                tempR1 = current_clusters[i][1]
                x1 = current_clusters[i][0][0]
                y1 = current_clusters[i][0][1]
                z1 = current_clusters[i][0][2]
                for j in range(len(clusters)):
                    tempR2 = clusters[j][1]
                    x2 = clusters[j][0][0]
                    y2 = clusters[j][0][1]
                    z2 = clusters[j][0][2]
                    tempDis = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)
                    Distance = math.sqrt(tempDis)
                    if (tempR1 + tempR2) > Distance:
                        return 1
    return 0

def getSum(key, dict):
    sumx = 0
    sumy = 0
    sumz = 0
    for keys in dict:
        if keys != key:
            for center in dict[keys]:
                sumx = sumx + center[0][0]
                sumy = sumy + center[0][1]
                sumz = sumz + center[0][2]
    return sumx,sumy,sumz

def getDistance(pdb_id, cluster, dict):
    current_centers = dict[pdb_id][cluster][0]
    distance = 0
    for keys in dict:
        if keys != pdb_id:
            for i in range(len(dict[keys])):
                centers = dict[keys][i][0]
                distance = distance+math.sqrt((current_centers[0]-centers[0])**2 + (current_centers[1]-centers[1])**2 + (current_centers[2]-centers[2])**2)
    return distance

def do_multiple_packing(multi_sphere_data, single_sphere_data, iteration = 10001, step = 1):
    # Shift the centers
    # Make dictionary of the form {'pdb_id' : [[centers,radius], [centers,radius], ... ,[centers,radius]], 'pdb_id' : ...}


    shift_centers(multi_sphere_data,single_sphere_data)

    dict = {}

    N = 0

    learning_rate = step

    sum_list = []

    for key in multi_sphere_data:
        dict[key] = []
        for keys in multi_sphere_data[key]:
            dict[key].append([multi_sphere_data[key][keys]['sphere_center'],multi_sphere_data[key][keys]['sphere_radius']])

    for key in dict:
        N = N+len(dict[key])

    for i in range(1,iteration):
        for key in dict:
            current_protein = dict[key]  #[[centers,radius],[centers,radius]...]
            random_cluster = randint(0,len(current_protein)-1) # integer
            choosen_cluster = current_protein[random_cluster]  #[[center,radius]]
            print("Intial distance")
            print(getDistance(key,random_cluster,dict))
            n = N - len(current_protein)
            sumx = 0
            sumy = 0
            sumz = 0
            sumx,sumy,sumz = getSum(key,dict)
            GradX = n*choosen_cluster[0][0] - sumx
            GradY = n*choosen_cluster[0][1] - sumy
            GradZ = n*choosen_cluster[0][2] - sumz
            tempsum = GradX * GradX + GradY * GradY + GradZ * GradZ
            tempsum = math.sqrt(tempsum)

            #print(tempsum)

            #print("Initial distance: "+str(tempsum))

            # if show_log != 0:
            #     if (ii == iteration - 1 or ii == math.ceil(iteration * 4 / 5) or ii == math.ceil(
            #             iteration * 3 / 5) or ii == math.ceil(iteration * 2 / 5) or ii == math.ceil(
            #             iteration * 1 / 5)) and jj == 0:
            #         print('setp', ii, 'tempsum:', tempsum, 'GradX:', GradX)

            GradX2 = GradX / tempsum
            GradY2 = GradY / tempsum
            GradZ2 = GradZ / tempsum

            # choosen_cluster[0][0] = choosen_cluster[0][0] - GradX2*learning_rate
            # choosen_cluster[0][1] = choosen_cluster[0][1] - GradY2*learning_rate
            # choosen_cluster[0][2] = choosen_cluster[0][2] - GradZ2*learning_rate
            
            for j in range(len(current_protein)):
                current_protein[j][0][0] = current_protein[j][0][0] - GradX2*learning_rate
                current_protein[j][0][1] = current_protein[j][0][1] - GradY2*learning_rate
                current_protein[j][0][2] = current_protein[j][0][2] - GradZ2*learning_rate

            mark = overlap_detection_multiple(key,dict)
            if mark == 1:
                for j in range(len(current_protein)):
                    current_protein[j][0][0] = current_protein[j][0][0] + GradX2*learning_rate
                    current_protein[j][0][1] = current_protein[j][0][1] + GradY2*learning_rate
                    current_protein[j][0][2] = current_protein[j][0][2] + GradZ2*learning_rate
            dict[key] = current_protein
            sum_list.append(tempsum)
            print("Final distance")
            print(getDistance(key,random_cluster,dict))

    x = []
    y = []
    z = []

    radius_list = []

    for key in dict:
        protein = dict[key]
        for i in range(len(protein)):
            x.append(protein[i][0][0])
            y.append(protein[i][0][1])
            z.append(protein[i][0][2])
            radius_list.append(protein[i][1])

    x_loc = np.array(x)
    y_loc = np.array(y)
    z_loc = np.array(z)

    center_list = [x_loc,y_loc,z_loc]

    show_dict = {}
    show_dict['sum'] = tempsum
    show_dict['grad'] = GradX
    show_dict['radius_list'] = radius_list
    show_dict['centers'] = center_list

    return show_dict




