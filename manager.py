#! /usr/bin/python3

import random
import sys
from datetime import datetime
import pickle
import os

from turbine import *



gene_params = ["N", "tt_1o", "tt_2o", "tt_3o", "tt_4o", "tt_5o", \
               "a_3_1", "a_3_2", "a_3_3", "a_3_4", "a_3_5", \
               "a_2_1", "a_2_2", "a_2_3", "a_2_4", "a_2_5", "a_2_6", \
               "c_x3_1", "c_x3_2", "c_x3_3", "c_x3_4", "c_x3_5", "c_x3_6", \
               "d_theta_i_1", "d_theta_i_2", "d_theta_i_3", "d_theta_i_4", "d_theta_i_5", "d_theta_i_6"]



def cross(gene_1, gene_2):
    new_gene = gene_1.copy()
    for param in gene_params:
        r = random.random()
        if (r < 0.3):
            new_gene[param] = gene_2[param]
        if (0.3 <= r) and (r < 0.7):
            new_gene[param] = (gene_1[param] + gene_2[param]) / 2

    return new_gene


def mutate(gene):
    num = random.randint(0, len(gene_params)-1)
    gene_param = gene_params[num]

    new_gene = gene.copy()

    if random.random() > 0.7:
        return new_gene

    if gene_param == "N":
        new_gene[gene_param] = random.randint(8000, 15000)
    if gene_param == "tt_1o":
        new_gene[gene_param] = get_tt_1i() - random.randint(int(get_dt_stg_ideal() * 0.3), int(get_dt_stg_ideal() * 1.7))
    if gene_param == "tt_2o":
        new_gene[gene_param] = get_tt_1i() - 2 * random.randint(int(get_dt_stg_ideal() * 0.4), int(get_dt_stg_ideal() * 1.6))
    if gene_param == "tt_3o":
        new_gene[gene_param] = get_tt_1i() - 3 * random.randint(int(get_dt_stg_ideal() * 0.6), int(get_dt_stg_ideal() * 1.4))
    if gene_param == "tt_4o":
        new_gene[gene_param] = get_tt_1i() - 4 * random.randint(int(get_dt_stg_ideal() * 0.8), int(get_dt_stg_ideal() * 1.2))
    if gene_param == "tt_5o":
        new_gene[gene_param] = get_tt_1i() - 5 * random.randint(int(get_dt_stg_ideal() * 0.9), int(get_dt_stg_ideal() * 1.1))
    if gene_param == "a_3_1" or gene_param == "a_3_2" or gene_param == "a_3_3" or gene_param == "a_3_4" or gene_param == "a_3_5":
        new_gene[gene_param] = m.radians(random.randint(0, 30))
    if gene_params == "a_2_1" or gene_param == "a_2_2" or gene_param == "a_2_3" or gene_param == "a_2_4" or gene_param == "a_2_5" or gene_param == "a_2_6":
        new_gene[gene_param] = m.radians(random.randint(30, 90))
    if gene_param == "d_theta_i_1" or gene_param == "d_theta_i_2" or gene_param == "d_theta_i_3" or gene_param == "d_theta_i_4" or gene_param == "d_theta_i_5" or gene_param == "d_theta_i_6":
        new_gene[gene_param] = m.radians(random.randint(0, 15))
    if gene_param == "c_x3_1":
        new_gene[gene_param] = random.randint(int(138 * 0.95), int(139 * 1.05))
    if gene_param == "c_x3_2":
        new_gene[gene_param] = random.randint(int(gene["c_x3_1"] * 0.9), int(gene["c_x3_1"] * 1.1))
    if gene_param == "c_x3_3":
        new_gene[gene_param] = random.randint(int(gene["c_x3_2"] * 0.9), int(gene["c_x3_2"] * 1.1))
    if gene_param == "c_x3_4":
        new_gene[gene_param] = random.randint(int(gene["c_x3_3"] * 0.9), int(gene["c_x3_3"] * 1.1))
    if gene_param == "c_x3_5":
        new_gene[gene_param] = random.randint(int(gene["c_x3_4"] * 0.9), int(gene["c_x3_4"] * 1.1))
    if gene_param == "c_x3_6":
        new_gene[gene_param] = random.randint(int(gene["c_x3_5"] * 0.9), int(gene["c_x3_5"] * 1.1))
    
    return new_gene


individuals_num = 3000
shuffle_times = 3000
max_loop_times = 10000
g_max_point = 0
same_count = 0


best_turbine = None

def main_loop(acm, generations):
    global best_turbine
    result_list = []
    for (i, turbine) in enumerate(generations):
        val = turbine.evaluate()
        if val > (total_max_point() - 10):
            print("")
            print("[FOUND!]")
            print("")
            turbine.print()
        result_list.append(val)
    result_sum = sum(result_list)

    max_point = max(result_list)
    if max_point < 1:
        print("max point {}".format(max_point))
        raise ZeroDivisionError()
    for (i, val) in enumerate(result_list):
        if val == max_point:
            best_turbine = generations[i]
            break


    if (acm % 100) == 50:
        print("\n[STEP{}({})] max point:{}\n".format(acm, datetime.now().strftime("%Y/%m/%d %H:%M:%S"), max_point))
        try:
            os.remove("./data.pickle")
        except:
            pass
        with open("./data.pickle", 'wb') as f:
            pickle.dump(best_turbine.get_gene(), f)

        
    if acm > max_loop_times:
        return

    selected = []
    for i in range(int(individuals_num / 10)):
        target_val = random.random() * result_sum
        target_id = 0
        for (j, n) in enumerate(result_list):
            target_val = target_val - n
            if target_val < 0:
                target_id = j
                break
        selected.append((result_list[target_id], generations[target_id]))
    selected_sum = sum([i[0] for i in selected])
    new_generations = []
    for target in selected:
        num = 0.5 * individuals_num * (target[0] / selected_sum)
        for i in range(int(num)):
            new_generations.append(Turbine(target[1].get_gene()))
    for i in range(individuals_num - len(new_generations)):
        target_val_1 = random.random() * selected_sum
        target_val_2 = random.random() * selected_sum
        target_id_1 = 0
        target_id_2 = 0
        for (j, n) in enumerate(selected):
            target_val_1 = target_val_1 - n[0]
            if target_val_1 < 0:
                target_id_1 = j
                break
        for (j, n) in enumerate(selected):
            target_val_2 = target_val_2 - n[0]
            if target_val_2 < 0:
                target_id_2 = j
                break
        new_gene = cross(selected[target_id_1][1].get_gene(), selected[target_id_2][1].get_gene())
        new_generations.append(Turbine(new_gene))

    for i in range(100):
        for j in range(len(new_generations)):
            if random.random() > 0.8:
                target_id = int(random.random() * 99)
                new_generations[target_id] = Turbine(mutate(new_generations[target_id].get_gene()))

    global g_max_point
    global same_count
    if max_point == g_max_point:
        same_count = same_count + 1
        if same_count > 100:
            best_turbine.evaluate(is_debug=True)
            best_turbine.print()
            same_count = 0
            for i in range(shuffle_times):
                for j in range(len(new_generations)):
                    if random.random() > 0.8:
                        target_id = int(random.random() * 99)
                        new_generations[target_id] = Turbine(mutate(new_generations[target_id].get_gene()))
    else:
        g_max_point = max_point
        same_count = 0


    main_loop(acm+1, new_generations)


if __name__ == '__main__':
    sys.setrecursionlimit(max_loop_times + 10)

    a_3 = m.radians(15)
    a_2 = m.radians(70)
    d_theta = m.radians(8)

    data_exist = os.path.exists("./data.pickle")

    while 1:
        generations = []
        for i in range(individuals_num):
            gene = {"N":11788, "tt_1o":1209, "tt_2o":1150, "tt_3o":1094, "tt_4o":1036, "tt_5o":933, \
                    "a_3_1":0.02, "a_3_2":0.04, "a_3_3":0.014, "a_3_4":0.04, "a_3_5":a_3, \
                    "a_2_1":1.1, "a_2_2":1.1, "a_2_3":1.1, "a_2_4":1.07, "a_2_5":a_2, "a_2_6":a_2, \
                    "c_x3_1":143.6, "c_x3_2":140.7, "c_x3_3":142.2, "c_x3_4":151, "c_x3_5":139, "c_x3_6":139, \
                    "d_theta_i_1":0.03, "d_theta_i_2":0.08, "d_theta_i_3":0.121, "d_theta_i_4":0.35, "d_theta_i_5":d_theta, "d_theta_i_6":d_theta}
            if data_exist:
                with open("./data.pickle", 'rb') as f:
                    gene = pickle.load(f)
            for j in range(shuffle_times):
                gene = mutate(gene)
            generations.append(Turbine(gene))
        try:
            main_loop(0, generations)
            best_turbine.print()
            exit()
        except ZeroDivisionError as err:
            print("Zero Division")
            continue
        except KeyboardInterrupt as err:
            best_turbine.print()
            with open("./data.pickle", 'wb') as f:
                pickle.dump(best_turbine.get_gene(), f)
            exit()
