all: orbitas_1year.txt orbitas_1000year.txt 4body.x orbitas-1000yr.png orbitas_1yr.png

orbitas-1000yr.png: trayectorias.py orbitas_1000year.txt
	python trayectorias.py

orbitas_1yr.png: trayectorias.py orbitas_1year.txt	

orbitas_1000year.txt: ./4body.x ic.txt
	./4body.x ic.txt 0.01 1000

orbitas_1year.txt: ./4body.x ic.txt
	./4body.x ic.txt 0.00001 1

4body.x: 4body.c
	cc 4body.c -o 4body.x -lm

clean:
	rm -f orbitas_1year.txt orbitas_1000year.txt