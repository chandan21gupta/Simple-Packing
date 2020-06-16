import sys
sys.path.insert(1,'../src/')
sys.path.insert(2,'../packing_single_sphere/')
import pdb2ball_multiple as P2BM
import pdb2ball_single as P2B
import visualization as DR
import json_saver as js
import simulate as simu
import json
import multiple_packing as mp

#sys.path.append("..")
import pprint


op_p2mb = { 'target_protein': '1bxn',
            'PDB_ori_path': '../IOfile/pdbfile/',
            'savepath' :'../IOfile/pdb_multi_sphere/',
            'k_use' : 1,
            'k' : 3,
            'saveORnot' : 1,
            'show_info': 1}

packing_op = simu.packing_op

packing_op['iteration'] = 5001

def packing_with_target_mtsp( op_p2mb ):

    # convert pdb file into single ball and get the center and radius of this ball.
    multi_sphere_info = P2BM.pdb2ball_multiple(op_p2mb)
    dic_print = pprint.PrettyPrinter(indent=4)
    dic_print.pprint(multi_sphere_info)

    # set target protein
    print('target protein is', op_p2mb['target_protein'],'\n\n')
    protein_name = []
    protein_name.append(op_p2mb['target_protein'])
    sphere_list = []
    sphere_list.append(multi_sphere_info[protein_name[0]])

    packing_result = simu.packing_with_target(packing_op)

    filenames = ['single_sphere','multiple_sphere']

    js.save_to_json(packing_result['optimal_result'],filenames[0])
    js.save_to_json(multi_sphere_info,filenames[1])

    savepath = '../IOfile/json/'


    with open(savepath+filenames[0]+'.json') as jsonfile:
        single_sphere_data = json.load(jsonfile)

    with open(savepath+filenames[1]+'.json') as jsonfile:
        multi_sphere_data = json.load(jsonfile)

    # #print(data)
    # DR.get_multiple_packing_and_plot_ball(data,min_dict)

    DR.get_multiple_packing_and_plot_ball(multi_sphere_data,single_sphere_data)

    result = mp.do_multiple_packing(multi_sphere_data, single_sphere_data)

    DR.drawing_center_with_ball(result['radius_list'], result['centers'])

    # print(sphere_list)

    # # select random proteins
    # random_protein = RS.get_random_protein(boundary_shpere,protein_number = random_protein_number)
    #
    # # get important info
    # info = RS.get_radius_and_id(random_protein, radii_list = radii_list, protein_name = protein_name, show_log = show_log)
    # radius_list = info['radius_list']
    # protein = info['protein_key']

    # set box

    # initialization

    # packing


if __name__ == '__main__':
    packing_op['boundary_shpere'] =  P2B.pdb2ball_single(PDB_ori_path='../IOfile/pdbfile/', show_log=0)
    try:
        packing_op = sys.argv[1]
        packing_with_target_mtsp(sys.argv[2])
    except:
        packing_with_target_mtsp(op_p2mb)




