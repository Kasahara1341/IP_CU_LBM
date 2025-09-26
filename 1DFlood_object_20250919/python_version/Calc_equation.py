from Element import Element
from Node import Node

import numpy as np

def calc_mainloop(elements_list, nodes_list, dt, H, Q, stages):

    # Iwasaki 修正 ルンゲクッタ 各ステージの流出する流量(q)のリスト作成
    for stage in range(stages):
        for target_elements in elements_list:
            for target_element in target_elements:
                target_element.solve_mass_equation(dt)# 質量保存側
        for target_nodes in nodes_list:
            for target_node in target_nodes:
                target_node.solve_momentum_equation()
    # 出力用
    for target_elements in elements_list:
        for target_element in target_elements:
            H.append(target_element.get_variable_depth())
    for target_nodes in nodes_list:
        for target_node in target_nodes:
            Q.append(target_node.get_variable_q())

