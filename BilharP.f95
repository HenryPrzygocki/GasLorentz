program BilharP
  use Dados
  implicit none
  integer(sd)              :: i!, kmax = 101, kc
  integer(id)              :: m, n, mmax
  integer(kind=4)          :: s1=132, s2=44534, s3=234231
  real(d)                  :: t1, v0, tm, Fx, Dlt, J0, dt, vx0, vy0, c1, vf(2) = 0.0_d
  real(d), allocatable     :: Visco(:), Dif(:), J(:), C(:)
  character(len=3)         :: IDENTIDADE

  open(5,file=path//"u00Finais.dat",Status="Unknown")


  r = (real(raio,d)/1000.0_d)*ld
  tm = (ld**2.0_d - pi*r**2.0_d)/(2.0_d*r)
  tmax = 50.0_d*tm
  mmax = 150*int(tmax/tm)
  dt = tmax/(real(mmax,d))
  allocate(J(0:mmax),Visco(0:mmax),Dif(0:mmax),C(0:mmax))
  write(IDENTIDADE,"(I3.3)") raio
  open(0,file=path//IDENTIDADE//"u00Info.dat",Status="Unknown")
  write(0,*) "Informacoes sobre a simulacao"
  write(0,*) "Simulação com raio constante e variacao de velocidade"
  write(0,*) "Numero de particulas: ", nmax
  write(0,*) "Velocidade: ", u
  write(0,*) "Tempo medio entre colisoes: ", tm

  C = 0.0_d
  J = 0.0_d
  Visco = 0.0_d
  Dif = 0.0_d


  open(1,file=path//IDENTIDADE//"u00Velocidade.dat",Status="Unknown")
  open(2,file=path//IDENTIDADE//"u00Momentos.dat",Status="Unknown")
  open(3,file=path//IDENTIDADE//"u00Difusao.dat",Status="Unknown")
  open(4,file=path//IDENTIDADE//"u00Viscosidade.dat",Status="Unknown")


  call CPU_TIME(c1)
  Particula: do i = 1, nmax
     call CONDini(s1, s2, s3)
     t = 0.0_d; l = 0; disco = .false.; k = 0.0_d; v = 1.0_d
     vx = v*cos(theta); vy = v*sin(theta); Dlt = 0.0_d!; kc = 0
     n = 0; vf = 0.0_d
     J0 = vx*vy
     J(0) = J(0) + J0**2.0_d
     C(0)= C(0) + 2.0_d*vx**2.0_d
     vx0 = vx; vy0 = vy
     Colisoes: do
        x0 = x; y0 = y; theta0 = theta; v0 = v; vxc = vx
        call CollisionsCisal(v0)
        t1 = (abs(x - x0) + abs(y - y0))/abs(v0)
        Dlt = Dlt + t1                                                     ! Tempo entre duas colisões consecutivas com discos em unidades do tempo médio
        t = t + t1
        if(t > tmax) t = tmax
        m = int(t/dt)
        if(disco(1) .or. disco(2)) then
           Fx = (vx - vxc)/Dlt
           J(n:m) = J(n:m) + (vx*vy + y*Fx)*J0
           C(n:m) = C(n:m) + 2.0_d*vxc*vx0
           Dlt = 0.0_d
           n = m
        else
           if(abs(x) >= z0) then
             ! dr(1) = dr(1) + 2.0_d*x
              if(x > 0.0_d) then
                 x = -z0
              else
                 x = z0
              endif
           else
             ! dr(2) = dr(2) + 2.0_d*y
             if(y >= z0) then
               y = -z0
               k = k + 1
             else
               y = z0
               k = k - 1
             endif
           endif
        endif
        if(t == tmax) exit
     enddo Colisoes
     if(mod(real(i,d),real(nmax)*0.001_d) == 0.0) call Status(i, c1)

     vf(1) = vf(1) + v**2.0_d
     enddo Particula
     vf(1) = vf(1)/real(nmax,d)
     J(0:mmax) = J(0:mmax)/real(nmax,d)
     C(0:mmax) = C(0:mmax)/real(nmax,d)

    do i = 1, mmax
      Visco(i) = Visco(i) + (J(i-1)+J(i))*dt/2.0_d
      Dif(i) = Dif(i) + (C(i-1)+C(i))*dt/2.0_d
      If(i+1<=mmax) Visco(i + 1) = Visco(i)
      If(i+1<=mmax) Dif(i + 1) = Dif(i)
      write(1,*) i*dt/tm, abs(C(i))
      write(2,*) i*dt/tm, abs(J(i))
      write(3,*) i*dt/tm, Dif(i)
      write(4,*) i*dt/tm, Visco(i)
   enddo
   write(5,*) r, Visco(mmax), vf(1)
  call Status(0, c1)

  close(1)
  close(2)
  close(3)
  close(4)
  deallocate(J, Visco, C, Dif)
enddo
close(0)


end program BilharP
