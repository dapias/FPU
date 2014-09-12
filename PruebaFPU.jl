include("/FPUjulia.jl")

N = 2
h = 0.1
n = 10000 #numero de pasos

x = generarics(N)
campo = campoFPU(x,N)
runFPU(x,campo, h,n,N)

