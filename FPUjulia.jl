# # Integrador que utiliza el tiempo inicial, el arreglo inicial, el campo vectorial, el paso y el número de pasos para integrar la ecuación ordinario dada por \dot{x} = F(x), considerando que es una ecuación diferencial autónoma.
using PyPlot

function rungeK(t_init,x_init,campo,h,n)
	state = x_init
	t = t_init
	for i in 1:n
		k1 = campo(state)
		k2 = campo(state + .5*h*k1)
		k3 = campo(state + .5*h*k2)
		k4 = campo(state + h*k3)
		phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
		state = state+h*phi
		t=t+h
  end
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




function generarics(N)
  ics = zeros(2N)
  for i in 1:N
    ics[i] = (2/(N+1))^(1/2)*sin(i*pi/(N+1))
    ics[N+i] = 0
  end
  ics
end
    #omegak[I] = 4*(sin(pi*(I+1)/2/N))**2 # Frecuencias de modo

function runFPU(x_init,t_init, campo, h, n,N)
  q = x_init
  t = t_init
  for i in 1:n
    q = hcat(q,rungeK(t,q[:,end],campo,h,1))
  end
  Q = transformada(q[:,1],N)
  for i in 2:n+1
    w =transformada(q[:,i],N)
    Q = hcat(Q,w)
  end
  tiempo = [h*i for i in 0:n]
  coord =  Q[1:N,:]
  velo =  Q[N+1:2N,:]
  omegak = [2*(sin(pi*k/(2*(N+1)))) for k in 1:N] # Frecuencias de modo
  energia = [((coord[k,1]).^2*omegak[k]^2 + (velo[k,1]).^2)/2. for k in 1:N]
    for i in 2:n+1
    energia = hcat(energia, [((coord[k,i]).^2*omegak[k]^2 + (velo[k,i]).^2)/2. for k in 1:N])
  end
  plotEnergia(energia,tiempo,N)
end

function plotEnergia(energia,tiempo, N)
for i in 1:N
    plot(tiempo, [x::Float64 for x in energia[i,:]])
end
  end



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




