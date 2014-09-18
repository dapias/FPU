# Solución del problema FPU usando un integrador de Taylor.

using TaylorSeries
using PyPlot
pygui(true)

const ordenTaylor = 28
const epsAbs = 1.0e-20

#Condiciones iniciales
function generarics(N)
  ics = zeros(2N)
  for i in 1:N
    #ics[i] = (2/(N+1))^(1/2)*sin(i*pi/(N+1)) Dar toda la energía al primer modo
    ics[i] = sin(pi*(i)/N) #Energía repartida pero concentrada en el primer modo
    ics[N+i] = 0
  end
  ics
end



function campoFPU{T<:Number}(x ::Array{T,1},alpha)
  #Esta función recibe un arreglo, lo convierte en uno tipo Taylor, resuelve las ecuaciones de movimiento y regresa los coeficientes de la solución.

  order = ordenTaylor
  len = length(x)
  vec0T = [ Taylor([x[i]], order) for i in 1:len]
  N = iceil(len/2)
  y = [ Taylor(0., order) for i=1:len ]
  D = [ Taylor(0., order) for i=1:len ]

  for k = 0:order-1
    knext = k+1

    for i in 1:len
      y[i] = Taylor( vec0T[i].coeffs[1:k+1], k)
    end

    #Ecuaciones de movimiento
    D[1] = y[N+1]
    D[N] = y[2N]
    D[N+1] = y[2]-2y[1]+alpha*((y[2]-y[1])^2 -(y[1])^2)
    D[2N] = y[N-1]-2*y[N]+alpha*((y[N])^2-(y[N]-y[N-1])^2)
    for i in 2:N-1
      D[N+i] = y[i+1]+y[i-1]-2*y[i]+alpha*((y[i+1]-y[i])^2-(y[i]-y[i-1])^2)
      D[i] = y[N+i]
    end

    #Actualización de los coeficientes
    for i in 1:len
      vec0T[i].coeffs[knext+1]  = D[i].coeffs[knext] / knext
    end

  end
  return vec0T
end

# Integrador
function taylorStepper{T<:Number}( jetEqs::Function, vec0::Array{T,1}, order::Int64, epsilon::T)
  #Esta función recibe un arreglo correspondiente a los datos iniciales, resuelve las ecuaciones de movimiento y evalúa la solución en el tiempo óptimo (determinado mediantela función stepsize).

  n = length( vec0 )
  vec1T = jetEqs( vec0 )

  # Step-size
  hh = Inf
  for i=1:n
    h1 = stepsize( vec1T[i], epsilon )
    hh = min( hh, h1 )
  end

  # Values at t0+h
  for i=1:n
    vec0[i] = evalTaylor( vec1T[i], hh )
  end

  return hh, vec0
end


# Returns the maximum step size from epsilon and the last two coefficients of the x-Taylor series
function stepsize{T<:Real}(x::Taylor{T}, epsilon::Float64)
  ord = x.order
  h = Inf
  for k in [ord-1, ord]
    kinv = 1.0/k
    aux = abs( x.coeffs[k+1] )
    h = min(h, (epsilon/aux)^kinv)
  end
  return h
end



function runFPU(x_init, campo, n,N)
  q = x_init
  omegak = [2*(sin(pi*k/(2*(N+1)))) for k in 1:N] # Frecuencias de modo
  matriz = crearMatriz(N)
  coordQ = transformada(q[1:N],matriz)
  veloQ = transformada(q[N+1:end],matriz)
  energia = [((coordQ[k]).^2*omegak[k]^2 + (veloQ[k]).^2)/2. for k in 1:N]
  tiempo = zeros(n+1)


  for i in 1:n
    t, q = taylorStepper(campo, x, ordenTaylor, epsAbs )
    if i > 1
      tiempo[i] = t + tiempo[i-1]
    end

    coordQ = transformada(q[1:N],matriz)
    veloQ = transformada(q[N+1:end],matriz)
    energia = hcat(energia,[((coordQ[k]).^2*omegak[k]^2 + (veloQ[k]).^2)/2. for k in 1:N])
  end


  plotEnergia(energia,tiempo,N )
end


#Función para dibujar la energía de cada modo.
function plotEnergia(energia,tiempo, N)
  for i in 1:N
    plot(tiempo[1:end-1], [x::Float64 for x in energia[i,:]][1:end-1])
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
