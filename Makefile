all: 	
	clear
	gcc -o MDFOpenMP mdf_openmp.c -fopenmp -lm -lpapi
	gcc -o MDF mdf.c -lm -lpapi
	gcc -o teste_pthread mdf_pthread.c -lpthread -lm -lpapi

	echo comandos: ./MDF ./teste_pthread ./MDFOpenMP

sequencial:
	clear
	gcc -o MDF mdf.c -lm -lpapi

	echo "comando: ./MDF"

openmp:
	clear
	gcc -o MDFOpenMP mdf_openmp.c -fopenmp -lm -lpapi

	echo comando: ./MDFOpenMP

pthread:
	clear
	gcc -o teste_pthread mdf_pthread.c -lpthread -lm -lpapi

	echo "comando: ./teste_pthread"