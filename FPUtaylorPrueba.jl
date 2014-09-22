include("FPUtaylor.jl")
using PyCall
#ordenTaylor = 28  #Valor por defecto: 28
#epsAbs = 1.0e-20  Valor por defecto: 1.0e-20
N = 2  #Número de partículas (FPU:32)
n = 300 #numero de pasos, número de pasos recomendado para 32 partículas: 30000
alpha = 0.25 #Parámetro de interacción (0.25 en el trabajo FPU)

x = generarics(N)
campo(vector) = campoFPU(vector, alpha)
runFPU(x,campo,n,N)
