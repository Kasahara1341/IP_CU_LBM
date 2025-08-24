index_list = [
    "IDX_dz",   
    "IDX_dt",   
    "IDX_c",        
    "IDX_nx",
    "IDX_ny",
    "IDX_nz",
    "IDX_Q",
    "IDX_ratiox",
    "IDX_ratioy",
    "IDX_PFthick",
    "IDX_num_calc",    
    "IDX_num_wall",    
    "IDX_num_IBPoints",
    "IDX_dIBM",
    "IDX_nu",    
    "IDX_alpha",    
    "IDX_sigma", 
    "IDX_tau",
    "IDX_taus",
    "IDX_wall1_size",
    "IDX_wall2_size",
    "IDX_wall3_size",
    "IDX_wall4_size",
    "IDX_wall5_size",
    "IDX_wall6_size",
    "IDX_vector_start",
    # ...追加可能
]

vector_list = [
    "IDX_w(k)   (IDX_vector_start + (int)items[IDX_Q]*0 + (k))",
    "IDX_cx(k)  (IDX_vector_start + (int)items[IDX_Q]*1 + (k))",
    "IDX_cy(k)  (IDX_vector_start + (int)items[IDX_Q]*2 + (k))",
    "IDX_cz(k)  (IDX_vector_start + (int)items[IDX_Q]*3 + (k))",
]

with open("itemIndex.hpp", "w") as f:
    f.write("#ifndef ITEMS_INDEX_H\n#define ITEMS_INDEX_H\n\n")
    for i, name in enumerate(index_list):
        f.write(f"#define {name} {i}\n")
        print(name,i)
    for i, name in enumerate(vector_list):
        f.write(f"#define {name}\n")
    f.write("\n#endif // ITEMS_INDEX_H\n")
