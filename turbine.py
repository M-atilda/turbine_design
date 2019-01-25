#! /usr/bin/python3

import math as m
import random

# consts
gamma = 1.33
R = 287
Cp = 1250
ata = 0.91

# util
def calc_MFP(M):
    return M * m.sqrt(gamma / R) / (1 + ((gamma - 1) / 2) * (M**2))**((gamma + 1) / 2 / (gamma - 1))

def calc_next_Pt(pt_i, tt_i, tt_o, m_0):
    it_3 = 635767.3
    it_i_d = ((Cp * tt_o * (m_0 + m_cool)) - ((m_cool * 2) * it_3)) / m_0
    tt_i_d = it_i_d / Cp
    return pt_i * (1 + ((tt_i_d / tt_i) - 1) / ata)**(gamma / (gamma - 1))




# sharing params
Ns = 6
AR = 4

M_1i = 0.2

alpha_1i = 0
alpha_6o = 0

# N = 15000

m_0_1i = 3.13
m_0_6o = 3.29
m_cool = (m_0_6o - m_0_1i) / (2 * Ns)
m_0_2i = m_0_1i + (2 * m_cool)
m_0_3i = m_0_2i + (2 * m_cool)
m_0_4i = m_0_3i + (2 * m_cool)
m_0_5i = m_0_4i + (2 * m_cool)
m_0_6i = m_0_5i + (2 * m_cool)


tt_1i = 1263
def get_tt_1i():
    return tt_1i
tt_6o = 867.3
dt_stg_ideal = (tt_1i - tt_6o) / Ns
def get_dt_stg_ideal():
    return dt_stg_ideal
d_i_ideal = Cp * dt_stg_ideal
pt_1i = 195718.5
pt_6o = 33815.76

t_start = tt_1i / (1 + (((gamma - 1) / 2) * (M_1i**2)))
c_x_start = int(M_1i * m.sqrt(gamma * R * t_start))

x_start = 0.12
d_0_start = 0.34
A_1_start = m_0_1i * m.sqrt(tt_1i) / pt_1i / calc_MFP(M_1i) / m.cos(0)
d_1_start = m.sqrt(d_0_start**2 + A_1_start * 4 / m.pi)



# limits
d_a_max = m.radians(120)
u_t_max = 300
psi_max = 1.5


# GA params
point_1_clear = 1
point_2_clear = 4
point_3_clear = 16
point_4_clear = 64
point_5_clear = 256
point_6_clear = 512

point_theta_clear = 2

point_psi_clear = 2
point_stg_clear = 5


def total_max_point():
    return (point_1_clear + point_2_clear + point_3_clear + point_4_clear + point_5_clear + point_6_clear) * (1 + point_theta_clear + point_stg_clear)


def fix_gene(stg, gene, result):
    new_gene = gene.copy()
    if result["kind"] == "wide":
        new_gene["d_theta_i_{}".format(stg)] = gene["d_theta_i_{}".format(stg)] - m.radians(random.randint(0, 15))
    if result["kind"] == "large":
        if random.random() < 0.3:
            new_gene["N"] = gene["N"] * 0.95
        else:
            new_gene["d_theta_i_{}".format(stg)] = gene["d_theta_i_{}".format(stg)] - m.radians(random.randint(0, 15))
    if result["kind"] == "small":
        if random.random() < 0.3:
            new_gene["N"] = gene["N"] * 1.05
        else:
            new_gene["d_theta_i_{}".format(stg)] = gene["d_theta_i_{}".format(stg)] + m.radians(random.randint(0, 15))
    if result["kind"] == "bend_2" or result["kind"] == "slow":
        new_gene["a_2_{}".format(stg)] = result["a_2"]
    if result["kind"] == "fast" or result["kind"] == "bend_3":
        if stg == 6:
            pass
        else:
            new_gene["a_3_{}".format(stg)] = result["a_3"]
    if result["kind"] == "match":
        r = random.random()
        if r < 0.3:
            new_gene["a_2_{}".format(stg)] = result["a_2"]
        if 0.3 <= r and r < 0.6:
            if stg == 6:
                pass
            else:
                new_gene["a_3_{}".format(stg)] = result["a_3"]
        if 0.6 <= r and r < 0.9:
            if stg == 6:
                pass
            else:
                new_gene["tt_{}o".format(stg)] = result["tt_3"]
    if stg == 1:
        if new_gene["d_theta_i_{}".format(stg)] > m.radians(15):
            new_gene["d_theta_i_{}".format(stg)] =  m.radians(15)
        if new_gene["d_theta_i_{}".format(stg)] < m.radians(0):
            new_gene["d_theta_i_{}".format(stg)] = m.radians(0)
    else:
        if new_gene["d_theta_i_{}".format(stg)] > (gene["d_theta_i_{}".format(stg-1)] + m.radians(15)):
            new_gene["d_theta_i_{}".format(stg)] = gene["d_theta_i_{}".format(stg-1)] + m.radians(15)
        if new_gene["d_theta_i_{}".format(stg)] < (gene["d_theta_i_{}".format(stg-1)] - m.radians(15)):
            new_gene["d_theta_i_{}".format(stg)] = gene["d_theta_i_{}".format(stg-1)] - m.radians(15)
        
    return new_gene


def calc_stage(N, d_0_1, d_0_3, tt_1, tt_3, pt_1, pt_3, m_1, a_1, a_2, a_3, c_x1, c_x3, is_debug):
    c_x2 = (c_x1 + c_x3) / 2
    c_1 = c_x1 / m.cos(a_1)
    c_2 = c_x2 / m.cos(a_2)
    c_3 = c_x3 / m.cos(a_3)

    t_1 = tt_1 - (c_1**2 / 2 / (gamma / (gamma - 1)) / R)
    t_2 = tt_1 - (c_2**2 / 2 / (gamma / (gamma - 1)) / R)
    t_3 = tt_3 - (c_3**2 / 2 / (gamma / (gamma - 1)) / R)
    M_1 = c_1 / m.sqrt(gamma * R * t_1)
    M_2 = c_2 / m.sqrt(gamma * R * t_2)
    M_3 = c_3 / m.sqrt(gamma * R * t_3)
        
    m_2 = m_1 + m_cool
    m_3 = m_2 + m_cool

    A_1 = m_1 * m.sqrt(tt_1) / pt_1 / calc_MFP(M_1) / m.cos(a_1)
    A_2 = m_2 * m.sqrt(tt_1) / pt_1 / calc_MFP(M_2) / m.cos(a_2)
    A_3 = m_3 * m.sqrt(tt_3) / pt_3 / calc_MFP(M_3) / m.cos(a_3)

    d_m_1 = m.sqrt(d_0_1**2 + (A_1 / 2) * 4 / m.pi)
    d_1_3 = m.sqrt(d_0_3**2 + A_3 * 4 / m.pi)
    d_m_3 = m.sqrt(d_0_3**2 + (A_3 / 2) * 4 / m.pi)
    d_m_mean = d_m_1 + 3 / 4 * (d_m_3 - d_m_1)

    u_0_3 = d_0_3 * m.pi * N / 60
    u_1_3 = d_1_3 * m.pi * N / 60
    u_mean = d_m_mean * m.pi * N / 60
    if u_1_3 > u_t_max:
        if is_debug:
            print("u_1_3:{}".format(u_1_3))
        return {"result":False, "point":0, "kind":"large"}

    dt_stg = tt_1 - tt_3
    psi = Cp * dt_stg / (u_0_3**2)
    if psi > psi_max:
        if is_debug:
            print("psi:{}".format(psi))
        return {"result":False, "point":(psi_max - psi), "kind":"small"}

    d_a_s = a_2 + a_1
    if d_a_s > d_a_max:
        if is_debug:
            print("das:{}".format(d_a_s))
        a_2 = 0.9 * (d_a_max - a_1)
        return {"result":False, "point":point_psi_clear, "kind":"bend_2", "a_2":a_2}
    c_u2 = c_x2 * m.tan(a_2)

    w_u2 = c_u2 - u_mean
    if w_u2 < 0:
        if is_debug:
            print("w_u2:{}".format(w_u2))
        a_2 = 1.1 * m.atan(u_mean / c_x2)
        return {"result":False, "point":point_psi_clear + 0.01, "kind":"slow", "a_2":a_2}
    b_2 = m.atan(w_u2 / c_x2)

    c_u3 = c_x3 * m.tan(a_3)
    w_u3 = c_u3 + u_mean

    w_3 = m.sqrt(w_u3**2 + c_x3**2)
    M_w3 = w_3 / m.sqrt(gamma * R * t_3)
    if M_w3 > 0.9:
        if is_debug:
            print("M_w3:{}".format(M_w3))
        w_3 = 0.85 * m.sqrt(gamma * R * t_3)
        w_u3 = m.sqrt(w_3**2 - c_x3**2)
        c_u3 = max(0, w_u3 - u_mean)
        a_3 = m.atan(c_u3 / c_x3)
        return {"result":False, "point":point_psi_clear + 0.02, "kind":"fast", "a_3":a_3}
        
    b_3 = m.atan(w_u3 / c_x3)
    d_a_r = b_3 + b_2
    if d_a_r > d_a_max:
        if is_debug:
            print("dar:{}".format(d_a_r))
        b_3 = d_a_max - b_2
        w_u3 = c_x3 * m.tan(b_3)
        c_u3 = max(w_u3 - u_mean, 10)
        a_3 = m.atan(c_u3 / c_x3)
        return {"result":False, "point":point_psi_clear + 0.03, "kind":"bend_3", "a_3":a_3}

    d_c_u = c_u2 + c_u3
    d_i = u_mean * d_c_u
    d_t_calc = d_i / Cp
    if (abs(d_t_calc / (tt_1 - tt_3) - 1) > 0.05):
        if is_debug:
            print("match:{}".format(d_t_calc / (tt_1 - tt_3)))
        c_u2_ideal = max(0, ((tt_1 - tt_3) * Cp / u_mean) - c_u3)
        a_2 = m.atan(c_u2_ideal / c_x2)
        c_u3_ideal = max(0, ((tt_1 - tt_3) * Cp / u_mean) - c_u2)
        a_3 = m.atan(c_u3_ideal / c_x3)
        tt_3 = tt_1 - d_t_calc
        return {"result":False, "point":point_psi_clear + 0.04, "kind":"match", "a_2":a_2, "a_3":a_3, "tt_3":tt_3}

    # d_h = 0.5 * (w_u3**2 - w_u2**2)
    # r = d_h / d_i

    return {"result":True, "point":point_stg_clear, "d_1_3":d_1_3}



class Turbine:
    def __init__(self, new_gene):
        # object "gene" {"N", "tt_1o", "tt_2o", "tt_3o", "tt_4o",
        #                "a_3_1", "a_3_2", "a_3_3", "a_3_4",
        #                "a_2_1", "a_2_2", "a_2_3", "a_2_4", a_2_5,
        #                "c_x3_1", "c_x3_2", "c_x3_3", "c_x3_4", "c_x3_5",
        #                "d_theta_i_1", "d_theta_i_2", "d_theta_i_3", "d_theta_i_4", "d_theta_i_5"}

        self.gene = new_gene.copy()
        self.point = 0


    def calc(self, is_debug):
        # initialize
        self.point = 0
        N = self.gene["N"]
        
        global tt_1i
        global pt_1i
        global m_0_1i
        a_2_1 = self.gene["a_2_1"]
        a_3_1 = self.gene["a_3_1"]
        c_x3_1 = self.gene["c_x3_1"]
        d_theta_i_1 = self.gene["d_theta_i_1"]
        
        tt_1o = self.gene["tt_1o"]
        pt_1o = calc_next_Pt(pt_1i, tt_1i, tt_1o, m_0_1i)
        global point_1_clear
        
        tt_2i = tt_1o
        pt_2i = pt_1o
        global m_0_2i
        a_2_2 = self.gene["a_2_2"]
        a_3_2 = self.gene["a_3_2"]
        
        tt_2o = self.gene["tt_2o"]
        pt_2o = calc_next_Pt(pt_2i, tt_2i, tt_2o, m_0_2i)
        global point_2_clear
        c_x3_2 = self.gene["c_x3_2"]
        d_theta_i_2 = self.gene["d_theta_i_2"]

        pt_3i = pt_2o
        tt_3i = tt_2o
        global m_0_3i
        a_2_3 = self.gene["a_2_3"]
        a_3_3 = self.gene["a_3_3"]

        tt_3o = self.gene["tt_3o"]
        pt_3o = calc_next_Pt(pt_3i, tt_3i, tt_3o, m_0_3i)
        global point_3_clear
        c_x3_3 = self.gene["c_x3_3"]
        d_theta_i_3 = self.gene["d_theta_i_3"]
        
        pt_4i = pt_3o
        tt_4i = tt_3o
        global m_0_4i
        a_2_4 = self.gene["a_2_4"]
        a_3_4 = self.gene["a_3_4"]

        tt_4o = self.gene["tt_4o"]
        pt_4o = calc_next_Pt(pt_4i, tt_4i, tt_4o, m_0_4i)
        global point_4_clear
        c_x3_4 = self.gene["c_x3_4"]
        d_theta_i_4 = self.gene["d_theta_i_4"]

        pt_5i = pt_4o
        tt_5i = tt_4o
        global m_0_5i
        a_2_5 = self.gene["a_2_5"]
        a_3_5 = 0

        tt_5o = self.gene["tt_5o"]
        pt_5o = calc_next_Pt(pt_5i, tt_5i, tt_5o, m_0_5i)
        global point_5_clear
        c_x3_5 = self.gene["c_x3_5"]
        d_theta_i_5 = self.gene["d_theta_i_5"]

        pt_6i = pt_5o
        tt_6i = tt_5o
        global m_0_6i
        a_2_6 = self.gene["a_2_6"]
        a_3_6 = 0

        global tt_6o
        global pt_6o
        global point_6_clear
        c_x3_6 = self.gene["c_x3_6"]
        d_theta_i_6 = self.gene["d_theta_i_6"]
        
        acm_i = [{"x":0, "y":(0.334 / 2)}, {"x":0.06, "y":(0.32 / 2)}, {"x":x_start, "y":(d_0_start / 2)}]
        acm_o = [{"x":0, "y":(0.384 / 2)}, {"x":0.06, "y":(0.396 / 2)}, {"x":x_start, "y":(d_1_start / 2)}]
        dy_dx_i_0 = (acm_i[-1]["y"] - acm_i[-2]["y"]) / (acm_i[-1]["x"] - acm_i[-2]["x"]) # inner gradient @ HPT
        dy_dx_o_0 = (acm_o[-1]["y"] - acm_o[-2]["y"]) / (acm_o[-1]["x"] - acm_o[-2]["x"]) # outer gradient @ HPT
        d_theta_i_0 = m.atan(dy_dx_i_0)



        dy_dx_i_1 = m.tan(d_theta_i_0 + m.radians(d_theta_i_1))
        chord_1 = 0.015 * 2
        d_0_1_1 = acm_i[-1]["y"] * 2
        d_0_3_1 = (acm_i[-1]["y"] + (chord_1 * 2) * dy_dx_i_1) * 2
        result_1 = calc_stage(N, d_0_1_1, d_0_3_1, tt_1i, tt_1o, pt_1i, pt_1o, m_0_1i, 0, a_2_1, a_3_1, c_x_start, c_x3_1, is_debug)
        self.point = self.point + (result_1["point"] * point_1_clear)
        if not result_1["result"]:
            self.gene = fix_gene(1, self.gene, result_1)
            return
        if is_debug:
            print("stg 1 clear")
        self.point = self.point + point_1_clear
        d_1_3_1 = result_1["d_1_3"]
        dy_dx_o_1 = (d_1_3_1 / 2 - acm_o[-1]["y"]) / chord_1
        d_theta_o_0 = m.atan(dy_dx_o_0)
        d_theta_o_1 = m.atan(dy_dx_o_1)
        if abs(d_theta_o_1 - d_theta_o_0) > m.radians(15):
            if is_debug:
                print(d_1_3_1 / 2)
                print(d_0_3_1 / 2)
                print(acm_i[-1])
                print(acm_o[-1])
                print(dy_dx_o_0)
                print(dy_dx_o_1)
                print(m.degrees(d_theta_o_0))
                print(m.degrees(d_theta_o_1))
                print("d_theta_o_1:{}".format(m.degrees(d_theta_o_1 - d_theta_o_0)))
            self.gene = fix_gene(1, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 1 clear")
        self.point = self.point + (point_theta_clear * point_1_clear)
        acm_i.append({"x":(acm_i[-1]["x"] + chord_1), "y":(d_0_3_1 / 2)})
        acm_o.append({"x":(acm_o[-1]["x"] + chord_1), "y":(d_1_3_1 / 2)})



        dy_dx_i_2 = m.tan(d_theta_i_1 + m.radians(d_theta_i_2))
        chord_2 = (acm_o[-1]["y"] - acm_i[-1]["y"]) / AR * 2
        d_0_1_2 = d_0_3_1
        d_0_3_2 = (acm_i[-1]["y"] + (chord_1 * 2) * dy_dx_i_2) * 2
        result_2 = calc_stage(N, d_0_1_2, d_0_3_2, tt_2i, tt_2o, pt_2i, pt_2o, m_0_2i, a_3_1, a_2_2, a_3_2, c_x3_1, c_x3_2, is_debug)
        self.point = self.point + (result_2["point"] * point_2_clear)
        if not result_2["result"]:
            self.gene = fix_gene(2, self.gene, result_2)
            return
        if is_debug:
            print("stg 2 clear")
        self.point = self.point + point_2_clear
        d_1_3_2 = result_2["d_1_3"]
        dy_dx_o_2 = (d_1_3_2 / 2 - acm_o[-1]["y"]) / chord_2
        d_theta_o_2 = m.atan(dy_dx_o_2)
        if abs(d_theta_o_2 - d_theta_o_1) > m.radians(15):
            if is_debug:
                print(d_1_3_2 / 2)
                print(d_0_3_2 / 2)
                print(acm_i)
                print(acm_o)
                print(dy_dx_o_1)
                print(dy_dx_o_2)
                print(m.degrees(d_theta_o_1))
                print(m.degrees(d_theta_o_2))
                print("d_theta_o_2:{}".format(m.degrees(d_theta_o_2 - d_theta_o_1)))
            self.gene = fix_gene(2, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 2 clear")
        self.point = self.point + (point_theta_clear * point_2_clear)
        acm_i.append({"x":(acm_i[-1]["x"] + chord_2), "y":(d_0_3_2 / 2)})
        acm_o.append({"x":(acm_o[-1]["x"] + chord_2), "y":(d_1_3_2 / 2)})



        dy_dx_i_3 = m.tan(d_theta_i_2 + m.radians(d_theta_i_3))
        chord_3 = (acm_o[-1]["y"] - acm_i[-1]["y"]) / AR * 2
        d_0_1_3 = d_0_3_2
        d_0_3_3 = (acm_i[-1]["y"] + (chord_3 * 2) * dy_dx_i_3) * 2
        result_3 = calc_stage(N, d_0_1_3, d_0_3_3, tt_3i, tt_3o, pt_3i, pt_3o, m_0_3i, a_3_2, a_2_3, a_3_3, c_x3_2, c_x3_3, is_debug)
        self.point = self.point + (result_3["point"] * point_3_clear)
        if not result_3["result"]:
            self.gene = fix_gene(3, self.gene, result_3)
            return
        if is_debug:
            print("stg 3 clear")
        self.point = self.point + point_3_clear
        d_1_3_3 = result_3["d_1_3"]
        dy_dx_o_3 = (d_1_3_3 / 2 - acm_o[-1]["y"]) / chord_3
        d_theta_o_3 = m.atan(dy_dx_o_3)
        if abs(d_theta_o_3 - d_theta_o_2) > m.radians(15):
            if is_debug:
                print(d_1_3_3 / 2)
                print(d_0_3_3 / 2)
                print(acm_i)
                print(acm_o)
                print(dy_dx_o_2)
                print(dy_dx_o_3)
                print(m.degrees(d_theta_o_2))
                print(m.degrees(d_theta_o_3))
                print("d_theta_o_3:{}".format(m.degrees(d_theta_o_3 - d_theta_o_2)))
            self.gene = fix_gene(3, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 3 clear")
        self.point = self.point + (point_theta_clear * point_3_clear)
        acm_i.append({"x":(acm_i[-1]["x"] + chord_3), "y":(d_0_3_3 / 2)})
        acm_o.append({"x":(acm_o[-1]["x"] + chord_3), "y":(d_1_3_2 / 3)})



        dy_dx_i_4 = m.tan(d_theta_i_3 + m.radians(d_theta_i_4))
        chord_4 = (acm_o[-1]["y"] - acm_i[-1]["y"]) / AR * 2
        d_0_1_4 = d_0_3_3
        d_0_3_4 = (acm_i[-1]["y"] + (chord_4 * 2) * dy_dx_i_4) * 2
        result_4 = calc_stage(N, d_0_1_4, d_0_3_4, tt_4i, tt_4o, pt_4i, pt_4o, m_0_4i, a_3_3, a_2_4, a_3_4, c_x3_3, c_x3_4, is_debug)
        self.point = self.point + (result_4["point"] * point_4_clear)
        if not result_4["result"]:
            self.gene = fix_gene(4, self.gene, result_4)
            return
        if is_debug:
            print("stg 4 clear")
        self.point = self.point + point_4_clear
        d_1_3_4 = result_4["d_1_3"]
        dy_dx_o_4 = (d_1_3_4 / 2 - acm_o[-1]["y"]) / chord_4
        d_theta_o_4 = m.atan(dy_dx_o_4)
        if abs(d_theta_o_4 - d_theta_o_3) > m.radians(15):
            self.gene = fix_gene(4, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 4 clear")
        self.point = self.point + (point_theta_clear * point_4_clear)
        acm_i.append({"x":(acm_i[-1]["x"] + chord_4), "y":(d_0_3_4 / 2)})
        acm_o.append({"x":(acm_o[-1]["x"] + chord_4), "y":(d_1_3_4 / 2)})



        dy_dx_i_5 = m.tan(d_theta_i_4 + m.radians(d_theta_i_5))
        chord_5 = (acm_o[-1]["y"] - acm_i[-1]["y"]) / AR * 2
        d_0_1_5 = d_0_3_4
        d_0_3_5 = (acm_i[-1]["y"] + (chord_5 * 2) * dy_dx_i_5) * 2
        result_5 = calc_stage(N, d_0_1_5, d_0_3_5, tt_5i, tt_5o, pt_5i, pt_5o, m_0_5i, a_3_4, a_2_5, a_3_5, c_x3_4, c_x3_5, is_debug)
        self.point = self.point + (result_5["point"] * point_5_clear)
        if not result_5["result"]:
            self.gene = fix_gene(5, self.gene, result_5)
            return
        if is_debug:
            print("stg 5 clear")
        self.point = self.point + point_5_clear
        d_1_3_5 = result_5["d_1_3"]
        dy_dx_o_5 = (d_1_3_5 / 2 - acm_o[-1]["y"]) / chord_5
        d_theta_o_5 = m.atan(dy_dx_o_5)
        if abs(d_theta_o_5 - d_theta_o_4) > m.radians(15):
            self.gene = fix_gene(5, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 5 clear")
        self.point = self.point + (point_theta_clear * point_5_clear)
        acm_i.append({"x":(acm_i[-1]["x"] + chord_5), "y":(d_0_3_5 / 2)})
        acm_o.append({"x":(acm_o[-1]["x"] + chord_5), "y":(d_1_3_5 / 2)})



        dy_dx_i_6 = m.tan(d_theta_i_5 + m.radians(d_theta_i_6))
        chord_6 = (acm_o[-1]["y"] - acm_i[-1]["y"]) / AR * 2
        d_0_1_6 = d_0_3_5
        d_0_3_6 = (acm_i[-1]["y"] + (chord_6 * 2) * dy_dx_i_6) * 2
        result_6 = calc_stage(N, d_0_1_6, d_0_3_6, tt_6i, tt_6o, pt_6i, pt_6o, m_0_6i, a_3_5, a_2_6, a_3_6, c_x3_5, c_x3_6, is_debug)
        self.point = self.point + (result_6["point"] * point_6_clear)
        if not result_6["result"]:
            self.gene = fix_gene(6, self.gene, result_6)
            return
        if is_debug:
            print("stg 6 clear")
        self.point = self.point + point_6_clear
        d_1_3_6 = result_6["d_1_3"]
        dy_dx_o_6 = (d_1_3_6 / 2 - acm_o[-1]["y"]) / chord_6
        d_theta_o_6 = m.atan(dy_dx_o_6)
        if abs(d_theta_o_6 - d_theta_o_5) > m.radians(15):
            self.gene = fix_gene(6, self.gene, {"kind":"wide"})
            return
        if is_debug:
            print("theta 6 clear")
        self.point = self.point + (point_theta_clear * point_6_clear)
        # acm_i.append({"x":(acm_i[-1]["x"] + chord_4), "y":(d_0_3_4 / 2)})
        # acm_o.append({"x":(acm_o[-1]["x"] + chord_4), "y":(d_1_3_4 / 2)})

        
        return

    def get_gene(self):
        return self.gene.copy()

    def print(self):
        print(self.gene)

    def evaluate(self, is_debug=False):
        try:
            self.calc(is_debug)
        except ValueError as err:
            self.point = 0
        except ZeroDivisionError as err:
            self.point = 0
        except TypeError:
            self.point = 0
        return self.point


if __name__ == '__main__':
    a_3 = m.radians(15)
    a_2 = m.radians(70)
    d_theta = m.radians(8)

    gene = {"N":10340, \
            "tt_1o":1200, "tt_2o":1133, "tt_3o":1067, "tt_4o":1000, "tt_5o":933, \
            "a_3_1":0.0, "a_3_2":0.0, "a_3_3":0.28, "a_3_4":a_3, "a_3_5":a_3, \
            "a_2_1":0.99, "a_2_2":0.89, "a_2_3":0.82, "a_2_4":a_2, "a_2_5":a_2, "a_2_6":a_2, \
            "c_x3_1":140, "c_x3_2":147, "c_x3_3":147, "c_x3_4":139, "c_x3_5":139, "c_x3_6":139, \
            "d_theta_i_1":0.244, "d_theta_i_2":0.2, "d_theta_i_3":2, "d_theta_i_4":d_theta, "d_theta_i_5":d_theta, "d_theta_i_6":d_theta}
    t = Turbine(gene)
    print(t.evaluate(is_debug=True))
    print(t.get_gene())
