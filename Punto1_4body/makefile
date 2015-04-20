all: orbitas.txt 4body.x

orbitas.txt: ./4body.x
	./4body.x ic.txt 0.001

4body.x: 4body.c
	cc 4body.c -o 4body.x -lm

clean:
	rm orbitas.txt