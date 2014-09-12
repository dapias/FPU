include("FPUjulia.jl")

#Se recomienda incrementar el tamaño de paso conforme se aumenta el número de partículas. Probado para N = 32. h = 1., y n= 10000.

N = 2  #Número de partículas
h = 0.1 #Stepsize
n = 1000 #numero de pasos
alpha = 0.5 #Parámetro de interacción (0.25 en el trabajo FPU)

x = generarics(N)
campo(vector) = campoFPU(vector, alpha)
runFPU(x,campo, h,n,N)

