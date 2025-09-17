arch = 89
# オプション
NVCC_FLAGS = -dc -arch=sm_$(arch)

# ディレクトリ
SRC_DIR = ../function
OBJ_DIR = ../library

# ソースファイルとオブジェクトファイル
SRC_FILES := $(wildcard $(SRC_DIR)/*.cu)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(SRC_FILES))

low:
	-@rm lib.a a.out
	ar r lib.a ../library/*
	nvcc Boegman2005run18.cu lib.a -arch=sm_$(arch) -rdc=true
	./a.out
	# python3 make_figure.py
	python3 make_VTK.py
precomp:
	-@rm ../library/*.o
comp: $(OBJ_FILES)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	nvcc $(NVCC_FLAGS) -o $@ $<
rm:
	-@rm data/*
	-@rm vtu_0/*
	-@rm figure/*
rmR:
	-@rm make_result/1_*.csv
	-@rm make_result/interface*.png
fig:
	python3 make_figure.py
vtk:
	python3 make_VTK.py
csv:
	python3 make_result/csv_height_A.py
	python3 make_result/csv_height_B.py
png:
	python3 make_result/png_height.py
result:
	python3 make_result/csv_height_B.py
	python3 make_result/png_height.py
suikou2025:
	python3 make_result/suikou2025/png_height.py
san:
	compute-sanitizer --tool memcheck --log-file error_memcheck.txt ./a.out
san2:
	-@rm lib.a
	ar r lib.a ../library/*
	nvcc internal_wave.cu lib.a 
	compute-sanitizer ./a.out
