#ifndef ITEMS_INDEX_H
#define ITEMS_INDEX_H

#define IDX_dz 0
#define IDX_dt 1
#define IDX_c 2
#define IDX_nx 3
#define IDX_ny 4
#define IDX_nz 5
#define IDX_Q 6
#define IDX_ratiox 7
#define IDX_ratioy 8
#define IDX_PFthick 9
#define IDX_num_calc 10
#define IDX_num_wall 11
#define IDX_num_IBPoints 12
#define IDX_dIBM 13
#define IDX_nu 14
#define IDX_alpha 15
#define IDX_sigma 16
#define IDX_tau 17
#define IDX_taus 18
#define IDX_wall1_size 19
#define IDX_wall2_size 20
#define IDX_wall3_size 21
#define IDX_wall4_size 22
#define IDX_wall5_size 23
#define IDX_wall6_size 24
#define IDX_vector_start 25
#define IDX_w(k)   (IDX_vector_start + (int)items[IDX_Q]*0 + (k))
#define IDX_cx(k)  (IDX_vector_start + (int)items[IDX_Q]*1 + (k))
#define IDX_cy(k)  (IDX_vector_start + (int)items[IDX_Q]*2 + (k))
#define IDX_cz(k)  (IDX_vector_start + (int)items[IDX_Q]*3 + (k))

#endif // ITEMS_INDEX_H
