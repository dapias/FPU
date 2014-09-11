# # Integrador que utiliza el tiempo inicial, el arreglo inicial, el campo vectorial, el paso y el número de pasos para integrar la ecuación ordinario dada por \dot{x} = F(x), considerando que es una ecuación diferencial autónoma.
using PyPlot
pygui(true)


function rungeK(x_init,campo,h)
	state = x_init
		k1 = campo(state)
		k2 = campo(state + .5*h*k1)
		k3 = campo(state + .5*h*k2)
		k4 = campo(state + h*k3)
		phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
		state = state+h*phi
	return state
end

# Campo vectorial para FPU

function campoFPU(vector, alpha=0.25, N =2)
  y = vector
  D = zeros(2N)
    D[1] = vector[N+1]
    D[N] = vector[2N]
    D[N+1] = y[2]-2y[1]+alpha*((y[2]-y[1])^2 -(y[1])^2)
    D[2N] = y[N-1]-2*y[N]+alpha*((y[N])^2-(y[N]-y[N-1])^2)
    for i in 2:N-1
        D[N+i] = y[i+1]+y[i-1]-2*y[i]+alpha*((y[i+1]-y[i])^2-(y[i]-y[i-1])^2)
        D[i] = y[N+i]
  end
    return D
end


#Generador de condiciones iniciales en el primer modo normal

function generarics(N)
  ics = zeros(2N)
  for i in 1:N
    ics[i] = (2/(N+1))^(1/2)*sin(i*pi/(N+1))
    ics[N+i] = 0
  end
  ics
end



function runFPU(x_init, scampo, h, n,N)
  q = x_init
  for i in 1:n
    q = hcat(q,rungeK(q[:,end],campo,h))
  end
  coordq =  q[1:N,:]
  veloq =  q[N+1:2N,:]
  coordQ = transformada(coordq[:,1],N)
  for i in 2:n+1
    w =transformada(coordq[:,i],N)
    coordQ = hcat(coordQ,w)
  end
  veloQ = transformada(veloq[:,1],N)
  for i in 2:n+1
    w =transformada(veloq[:,i],N)
    veloQ = hcat(veloQ,w)
  end

  tiempo = [h*i for i in 0:n]
  omegak = [2*(sin(pi*k/(2*(N+1)))) for k in 1:N] # Frecuencias de modo
  energia = [((coordQ[k,1]).^2*omegak[k]^2 + (veloQ[k,1]).^2)/2. for k in 1:N]
    for i in 2:n+1
    energia = hcat(energia, [((coordQ[k,i]).^2*omegak[k]^2 + (veloQ[k,i]).^2)/2. for k in 1:N])
  end
  plotEnergia(energia,tiempo,N)
end

function plotEnergia(energia,tiempo, N)
for i in 1:N
    plot(tiempo, [x::Float64 for x in energia[i,:]])
end
  end


#Transformada discreta con seno.

function transformada(vector,N)
  l = length(vector)
  matriz = eye(l)
  for j in 1:l
      for i in 1:l
            matriz[j,i] = sqrt(2/(N+1))*sin(i*j*pi/(N+1))
        end
    end
    matriz*vector
end




x = generarics(2)
campo = campoFPU
runFPU(x,campo, 0.001, 1000,2)
