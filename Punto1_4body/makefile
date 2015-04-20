all: orbitas_1year.txt orbitas_1000year.txt 4body.x orbitas-1000yr.png orbitas_1yr.png

orbitas-1000yr.png orbitas_1yr.png: trayectorias.py
	python trayectorias.py

orbitas_1000year.txt: ./4body.x ic.txt
	./4body.x ic.txt 0.1 1000

orbitas_1year.txt: ./4body.x ic.txt
	./4body.x ic.txt 0.0001 1

4body.x: 4body.c
	cc 4body.c -o 4body.x -lm

clean:
	rm orbitas_1year.txt rm orbitas_1000year.txt