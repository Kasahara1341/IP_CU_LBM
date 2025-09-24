run:
	g++ -O3 MAIN.cpp src/ElemNod.cpp src/Solver.cpp src/ReadWrite.cpp && ./a.out 
gif:
	python3 make_Qgif.py
	python3 make_Hgif.py

bind:
g++ -O3 -shared -std=c++17 -fPIC \
    $(python3 -m pybind11 --includes) \
    src/*.cpp -o example$(python3-config --extension-suffix)