# # Integrador que utiliza el tiempo inicial, el arreglo inicial, el campo vectorial, el paso y el número de pasos para integrar la ecuación ordinario dada por \dot{x} = F(x), considerando que es una ecuación diferencial autónoma.
using PyCall
pygui(:qt)
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

function campoFPU(vector, alpha)
  y = vector
  len = length(y)
  N = iceil(len/2)
  D = zeros(2N)
  D[1] = y[N+1]
  D[N] = y[2N]
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
    #ics[i] = (2/(N+1))^(1/2)*sin(i*pi/(N+1)) Dar toda la energía al primer modo
    ics[i] = sin(pi*(i)/N) #Energía repartida pero concentrada en el primer modo
    ics[N+i] = 0
  end
  ics
end



function runFPU(x_init, campo, h, n,N)
  q = x_init
  omegak = [2*(sin(pi*k/(2*(N+1)))) for k in 1:N] # Frecuencias de modo
  tiempo = h*[i for i in 0:n]
  # tiempo = omegak[1]/(2*pi)*h*[i for i in 0:n] En términos de ciclos de la primera frecuencia

  matriz = crearMatriz(N)
  coordQ = transformada(q[1:N],matriz)
  veloQ = transformada(q[N+1:end],matriz)
  energia = [((coordQ[k]).^2*omegak[k]^2 + (veloQ[k]).^2)/2. for k in 1:N]

  for i in 1:n
    q = rungeK(q,campo,h)
    coordQ = transformada(q[1:N],matriz)
    veloQ = transformada(q[N+1:end],matriz)
    energia = hcat(energia,[((coordQ[k]).^2*omegak[k]^2 + (veloQ[k]).^2)/2. for k in 1:N])
  end

  plotEnergia(energia,tiempo,N )
end


#Función para dibujar la energía de cada modo.
function plotEnergia(energia,tiempo, N)
  for i in 1:N
    plot(tiempo, [x::Float64 for x in energia[i,:]])
  end
  xlabel("Tiempo")
  ylabel("Energia")
  title(L"FPU-$\alpha$")
end


#Transformada discreta con seno.

function crearMatriz(N)
  matriz = eye(N)
  for j in 1:N
    for i in 1:N
      matriz[j,i] = sqrt(2/(N+1))*sin(i*j*pi/(N+1))
      #matriz[j,i] = sin(i*j*pi/(N+1))
    end
  end
  matriz
end



function transformada(vector,matriz)
  matriz*vector
end



