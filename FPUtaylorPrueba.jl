N = 32  #Número de partículas
n = 30000 #numero de pasos
alpha = 0.25 #Parámetro de interacción (0.25 en el trabajo FPU)

x = generarics(N)
campo(vector) = campoFPU(vector, alpha)
runFPU(x,campo,n,N)