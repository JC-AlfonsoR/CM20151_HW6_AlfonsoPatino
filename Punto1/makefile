all: orbitas.txt valores.txt 4body.x

orbitas.txt valores.txt: 
	./4body.x ic.txt 100 > valores.txt

4body.x: 4body.c
	cc 4body.c -o 4body.x -lm

clean:
	rm orbitas.txt valores.txt