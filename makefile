
all:
	g++ -O3 COW2.cpp -o COW2
	./COW2
	python3 ./viz2.py ./results/COW2_chain.txt