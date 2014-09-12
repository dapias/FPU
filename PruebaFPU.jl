include("FPUjulia.jl")

N = 2
h = 0.1
n = 1000 #numero de pasos

x = generarics(N)
campo = campoFPU
runFPU(x,campo, h,n,N)

