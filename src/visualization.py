from mpl_toolkits.mplot3d import Axes3D, axes3d
import matplotlib.pyplot as plt
import numpy as np
import json
import math

def show_center_img(x, y, z):
    print("Show distribution img of protein's center")
    ax = plt.subplot(111, projection='3d')  # create project
    ax.scatter(x, y, z, c='r')  # draw

    ax.set_zlabel('Z')  # aix
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()


def show_sum_img(sumlist, protein_num, pdb_id):
    ##process sumlist
    stride = int(len(sumlist) / protein_num)
    newlist = process_sumlist(sumlist, protein_num, stride)
    print('Show sum img')
    ax = plt.subplot()
    for line in range(protein_num):
        ax.plot(newlist[stride * line:stride * (line + 1)], label=pdb_id[line].upper())
    plt.legend()
    plt.title('Loss of protein')
    plt.xlabel('Iteration')
    plt.ylabel('Loss')


def process_sumlist(sumlist, protein_num, stride):
    newlist = [0 for _ in range(len(sumlist))]
    for ii in range(len(sumlist)):
        temp_remainder = ii % protein_num
        temp_quotient = int((ii - temp_remainder) / protein_num)
        newlist[temp_remainder * stride + temp_quotient] = sumlist[ii]
    return newlist


def get_packing_and_plot_ball(optimal_result, boundary_shpere):
    # optimal_result is data read from json file by instruction "json.load(file)"
    # boundary_shpere is a list of the radius of each protein
    pdb_id = optimal_result['pdb_id']
    x_loc = optimal_result['x']
    y_loc = optimal_result['y']
    z_loc = optimal_result['z']

    radius_list = []
    center_list = [x_loc, y_loc, z_loc]
    for ii in range(len(pdb_id)):
        radius_list.extend([boundary_shpere[pdb_id[ii]]['radius'][0]])
    print(radius_list)
    print(center_list)

    drawing_center_with_ball(radius_list, center_list)

def checkOverlap(multi_sphere):
    # multi_sphere_info = {  '1f1b': {   0: {  'sphere_center': array([36.968, 47.965, -9.969], dtype=float32),
    #                                                  'sphere_id': 0,
    #                                                  'sphere_radius': 39.7004,
    #                                                  'vmd_radius': 26.820240783691407},

    #                                            1: {  'sphere_center': ***,
    #                                                  'sphere_id': ***,
    #                                                  'sphere_radius': ***,
    #                                                  'vmd_radius': ****},

    #                                            2: {  sphere_info }

    #                                            ...
    #                                      },

    #                              'pdb_id': { 0:{}, 1:{}, 2:{}, ... n:{} },
    #                              'pdb_id': { 0:{}, 1:{}, 2:{}, ... n:{} },
    #                              ...
    #                            }

    flag = False

    dict = {}
    for key in multi_sphere:
        dict[key] = []
        for keys in multi_sphere[key]:
            dict[key].append([multi_sphere[key][keys]['sphere_center'], multi_sphere[key][keys]['sphere_radius']])
    
    for key in dict:
        print(key+" "+str(len(dict[key])))

    for key in dict:
        for i in range(len(dict[key])):
            x1 = dict[key][i][0][0]
            y1 = dict[key][i][0][1]
            z1 = dict[key][i][0][2]
            tempR1 = dict[key][i][1]
            for keys in dict:
                if(keys != key):
                    for j in range(len(dict[keys])):
                        x2 = dict[keys][j][0][0]
                        y2 = dict[keys][j][0][1]
                        z2 = dict[keys][j][0][2]
                        tempR2 = dict[keys][j][1]
                        Distance = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                        if tempR1 + tempR2 > Distance:
                            print("Overlap!!")
                            print(key+" "+str(i))
                            print(keys+" "+str(j))
                            flag = True
    return flag



def get_multiple_packing_and_plot_ball(multi_sphere, optimal_result):
    pdb_id = optimal_result['pdb_id']

    multi = {}

    for ids in pdb_id:
        multi[ids] = multi_sphere[ids]

    shifted_centers_x = optimal_result['x']
    shifted_centers_y = optimal_result['y']
    shifted_centers_z = optimal_result['z']

    # for i in range(len(shifted_centers_x)):
    #     print('protein '+i+' center - '+str(shifted_centers_x[i])+' '+str(shifted_centers_y[i])+' '+str(shifted_centers_z[i]))

    x = []
    y = []
    z = []

    radius_list = []

    for ii in range(len(pdb_id)):
        for keys in multi[pdb_id[ii]]:
            np.array(multi[pdb_id[ii]][keys]['sphere_center'])
            shifted_centers = np.array([shifted_centers_x[ii], shifted_centers_y[ii], shifted_centers_z[ii]])
            #print(str(shifted_centers[0])+' '+str(shifted_centers[1])+' '+str(shifted_centers[2]))
            # print(shifted_centers)
            # print(multi[pdb_id[ii]][keys]['sphere_center'])
            multi[pdb_id[ii]][keys]['sphere_center'] = multi[pdb_id[ii]][keys]['sphere_center'] - shifted_centers
            # print(multi[pdb_id[ii]][keys]['sphere_center'])
            x.append(multi[pdb_id[ii]][keys]['sphere_center'][0])
            y.append(multi[pdb_id[ii]][keys]['sphere_center'][1])
            z.append(multi[pdb_id[ii]][keys]['sphere_center'][2])
            radius_list.append(multi[pdb_id[ii]][keys]['sphere_radius'])

    x_loc = np.array(x)
    y_loc = np.array(y)
    z_loc = np.array(z)

    center_list = [x_loc, y_loc, z_loc]

    if checkOverlap(multi):
        for key in multi:
            print(key)

    drawing_center_with_ball(radius_list, center_list)

    print(len(x))

def drawing_center_with_ball(radius_list, loc_list):
    # for loc_list row 0,1,2 of the matrix represent x,y,z, colum 0,1,... represent location of each protein

    # draw one ball
    center = [loc_list[0][0], loc_list[1][0], loc_list[2][0]]
    radius = radius_list[0]
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x_list = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y_list = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z_list = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

    # draw all ball
    for k in range(1, len(radius_list)):
        center = [loc_list[0][k], loc_list[1][k], loc_list[2][k]]
        radius = radius_list[k]
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

        x_list = np.concatenate([x, x_list], axis=1)
        y_list = np.concatenate([y, y_list], axis=1)
        z_list = np.concatenate([z, z_list], axis=1)

    X = x_list
    Y = y_list
    Z = z_list

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #ax.set_aspect('equal')
    ax.plot_surface(x_list, y_list, z_list, rstride=4, cstride=4, color='b')
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
# Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    print(len(radius_list))
    plt.show()

def get_json_and_plot_ball(file):
    # input the path of the json file, then the function will read information in
    # json file and call the drawing function
    with open(file) as f:
        packing_result = json.load(f)
        pdb_id = packing_result['optimal_result']['pdb_id']  # list of the protein id
        x_loc = packing_result['optimal_result']['x']  # location of each protein
        y_loc = packing_result['optimal_result']['y']
        z_loc = packing_result['optimal_result']['z']

    radius_list = []
    center_list = [x_loc, y_loc, z_loc]  # combine x,y,z lists into one list
    for ii in range(len(pdb_id)):
        radius_list.extend(
            [packing_result['boundary_shpere'][pdb_id[ii]]['radius'][0]])  # list of radius of each protein
    print(radius_list)
    print(center_list)
    drawing_center_with_ball(radius_list, center_list)

# ######################### test part ####################
# file='../IOfile/packing.json'
# get_json_and_plot_ball(file)
#
