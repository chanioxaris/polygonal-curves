OBJS = structs.o main.o metric_functions.o hashtable.o functions.o preprocessing.o clustering.o clustering_init.o clustering_assignment.o clustering_update.o output_functions.o binary_tree.o
EXEC = clustering
GCCC = gcc -c
GCCO = gcc -o
LIBS = -lm

$(EXEC): $(OBJS)
	$(GCCO) $(EXEC) $(OBJS) $(LIBS)

structs.o: structs.c
	$(GCCC) structs.c

main.o: main.c
	$(GCCC) main.c

metric_functions.o: metric_functions.c
	$(GCCC) metric_functions.c

hashtable.o: hashtable.c
	$(GCCC) hashtable.c

functions.o: functions.c
	$(GCCC) functions.c

preprocessing.o: preprocessing.c
	$(GCCC) preprocessing.c
	
clustering.o: clustering.c
	$(GCCC) clustering.c
	
clustering_init.o: clustering_init.c
	$(GCCC) clustering_init.c
	
clustering_assignment.o: clustering_assignment.c
	$(GCCC) clustering_assignment.c
	
clustering_update.o: clustering_update.c
	$(GCCC) clustering_update.c

output_functions.o: output_functions.c
	$(GCCC) output_functions.c
	
binary_tree.o: binary_tree.c
	$(GCCC) binary_tree.c
	
clean:
	rm -rf $(OBJS) $(EXEC)