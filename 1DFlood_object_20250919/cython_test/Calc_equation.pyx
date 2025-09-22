# mainloop.pyx
from Element cimport Element
from Node cimport Node

cpdef void calc_mainloop(list elements_list, list nodes_list,
                  double dt, list H, list Q, int stages):
    cdef int stage
    cdef Element target_element
    cdef Node target_node
    cdef list target_elements
    cdef list target_nodes

    # 二重リストを平坦化
    cdef list flat_elements = [e for sublist in elements_list for e in sublist]
    cdef list flat_nodes    = [e for sublist in nodes_list for e in sublist]


    # Iwasaki 修正ルンゲクッタ 各ステージ
    for stage in range(stages):
        for target_element in flat_elements:
            target_element.solve_mass_equation(dt)   # 質量保存側
        for target_node in flat_nodes:
            target_node.solve_momentum_equation()

    # 出力用
    for target_element in flat_elements:
        H.append(target_element.get_variable_depth())
    for target_node in flat_nodes:
        Q.append(target_node.get_variable_q())
