all: 	
	clear
	gcc -o MDFOpenMP mdf_openmp.c -fopenmp -lm -lpapi -fsanitize=thread
	gcc -o MDF mdf.c -lm -lpapi -fsanitize=thread
	gcc -o teste_pthread mdf_pthread.c -lpthread -lm -lpapi -fsanitize=thread

	echo comandos: ./MDF ./teste_pthread ./MDFOpenMP

sequencial:
	clear
	gcc -o MDF mdf.c -lm -Ofast

	echo "comando: ./MDF"

openmp:
	clear
	gcc -o MDFOpenMP mdf_openmp.c -fopenmp -lm -O1

	echo comando: ./MDFOpenMP

pthread:
	clear
	gcc -o teste_pthread mdf_pthread.c -lpthread -lm -lpapi -fsanitize=thread

	echo "comando: ./teste_pthread"