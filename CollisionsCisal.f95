! Colisoes em um bilhar de Sinai com cisalhamento e com condicoes de contorno periodicas
Subroutine CollisionsCisal(v0)
  Use Dados
  Implicit none
  Real(d)                :: aa, bb, den, alpha, cosine, xc(2), yc(2), phi, csine(4), thetaq, xf
  Real(d)                :: rx(2), delta(2), x1, x2, va, v0
  Logical                :: parcial(2)

  x0 = x; y0 = y; theta0 = theta                              ! guarda coordenadas antes da colisão
  parcial(1) = .false.                                              ! teste de colisao com espalhador
  parcial(2) = .false.
  If(abs(vy) >= 1.0E-10) vx = vx - k*u
  va = sqrt(vx**2.0_d + vy**2.0_d)                                        ! velocidade relativa do partícula com o espalhador
  thetaq = asin(vy/va)                                            ! angulo da velocidade relativa
  if(vx < 0.0_d) thetaq = pi - thetaq
  if(thetaq < 0.0_d) thetaq = 2.0_d*pi + thetaq

  rx(1) = mod(k*u*t, z0)                                   ! posicao 1 do centro do espalhador
  if(rx(1) >= 0.0_d) rx(2) = rx(1) - ld                        ! posicao 2 do centro do espalhador
  if(rx(1) < 0.0_d) then
    rx(2) = rx(1)
    rx(1) = rx(1) + ld
  endif

  if( k*u == 0.0_d ) then
    rx(1) = 0.0_d
    rx(2) = ld + r
    thetaq = theta
  endif



  aa = tan(thetaq);  bb = y0 - x0*aa

  den = 1.0_d + aa**2.0_d
  delta(1) = r**2.0_d*den - (bb + aa*rx(1))**2.0_d
  delta(2) = r**2.0_d*den - (bb + aa*rx(2))**2.0_d

  If(delta(1) >= 0.0_d .and. .not.disco(1)) then                      ! Verifica o disco(1)

    If(abs(x-rx(1))>abs(y)) then
      alpha = acos((x-rx(1))/sqrt((x-rx(1))**2.0_d + y**2.0_d))
      If(y<0.0_d) alpha = 2.0_d*pi - alpha                        ! alpha e a posicao angular da particula em relacao ao centro do espalhador 1
    else
      alpha = asin(y/sqrt(y**2.0_d + (x-rx(1))**2.0_d))
      If((x-rx(1))<0.0_d) alpha = pi - alpha
    endif

    If(alpha<0.0_d) alpha = 2.0_d*pi + alpha
    cosine = cos(alpha - theta)

    If(cosine < 0.0_d) then       ! verifica se a trajetoria da particula esta orientada em direcao ao espalhador 1
      xc(1) = (rx(1)-aa*bb + sqrt(delta(1)))/den
      xc(2) = (rx(1)-aa*bb - sqrt(delta(1)))/den
      if(abs(x-xc(1))<=abs(x-xc(2))) then
        x1 = xc(1)
      else
        x1 = xc(2)
      Endif
      if(abs(x1)<=z0) then
        parcial(1) = .true.
      else
        parcial(1) = .false.
      endif
    Endif
  Endif

  If(delta(2) >= 0.0_d .and..not.disco(2)) then            ! repete os calculos para o espalhador 2

    If(abs(x-rx(2))>abs(y)) then
      alpha = acos((x-rx(2))/sqrt((x-rx(2))**2.0_d + y**2.0_d))
      If(y<0.0_d) alpha = 2.0_d*pi - alpha
    else
      alpha = asin(y/sqrt(y**2.0_d + (x-rx(2))**2.0_d))
      If((x-rx(2))<0.0_d) alpha = pi - alpha
    endif

    If(alpha<0.0_d) alpha = 2.0_d*pi + alpha
    cosine = cos(alpha - theta)

    If(cosine < 0.0_d) then
      xc(1) = (rx(2)-aa*bb + sqrt(delta(2)))/den
      xc(2) = (rx(2)-aa*bb - sqrt(delta(2)))/den
      if(abs(x-xc(1))<abs(x-xc(2))) then
        x2 = xc(1)
      else
        x2 = xc(2)
      Endif
      if(abs(x2)<=z0) then
        parcial(2) = .true.
      else
        parcial(2) = .false.
      endif
    Endif
  Endif


  If(parcial(1) .and. .not.parcial(2)) then
    xr = rx(1)
    xf = x1
    disco(1) = .true.
    disco(2) = .false.
  Elseif(.not.parcial(1) .and. parcial(2)) then
    xr = rx(2)
    xf = x2
    disco(1) = .false.
    disco(2) = .true.
  ElseIf(parcial(1) .and. parcial(2)) then
    If(abs(x-x1) < abs(x-x2)) then                                ! caso ainda nao tenha decidido em qual espalhador houve a colisão, escolhe-se a mais proxima
      xr = rx(1)
      xf = x1
      disco(1) = .true.
      disco(2) = .false.
    Else
      xr = rx(2)
      xf = x2
      disco(1) = .false.
      disco(2) = .true.
    Endif
  Elseif(.not.parcial(1).and..not.parcial(2)) then
    disco(1)= .false.
    disco(2)= .false.
  Endif

  If(disco(1) .or. disco(2)) then
    y = aa*xf + bb
    if(abs((xf-xr)/r) < 1.0_d) then
      phi = acos((xf-xr)/r)
      If(y < 0.0_d) phi = 2.0_d*pi - phi
    else
      phi = asin(y/r)
      if(xf-xr < 0.0_d) phi = pi - phi
    endif
    x = (y-y0)/tan(theta) + x0 + xf - r*cos(phi) - xr
    xr = (y-y0)/tan(theta) + x0 - r*cos(phi)
    if(abs(xr) < 1.0E-10) xr = 0.0_d
  Endif

  If(abs(y) > z0 .or. abs(x) > z0) then
    disco(1) = .false.
    disco(2) = .false.
  Endif

  If(disco(1) .or. disco(2)) then                                       ! aceitou uma colisao com o disco
    l=l+1
    vx = -v0*cos(2.0_d*phi-theta) + 2.0_d*k*u*cos(phi)**2.0_d          ! troca de velocidade na colisao com o espalhador
    vy = -v0*sin(2.0_d*phi-theta) + k*u*sin(2.0_d*phi)
    v = sqrt(vx**2.0_d + vy**2.0_d)

    theta = atan(vy/vx)
    If(vx < 0.0_d) theta = pi + theta                               ! angulo de saida apos a colisao

    theta = dmod(theta,2.0_d*pi)
    If(theta<0.0_d) theta = 2.0_d*pi + theta
    If(xr > 0.0_d) then
      disco(1) = .true.
      disco(2) = .false.
    else
      disco(1) = .false.
      disco(2) = .true.
    endif
    ! if(track) then
    !   write(9,*) xr + dr(1), dr(2)
    !   write(8,*) x + dr(1), y + dr(2)
    ! endif

    return

  Else                                                            ! Paredes
    aa = tan(theta0);  bb = y - x*aa
    yc(1) = aa*z0 + bb
    yc(2) = -aa*z0 + bb
    xc(1) = (z0 - bb)/aa
    xc(2) = -(z0 + bb)/aa
    csine(1) = cos(theta0)
    csine(2) = cos(pi/2.0_d - theta0)
    csine(3) = cos(pi - theta0)
    csine(4) = cos(3.0_d*pi/2.0_d - theta0)
    If(abs(yc(1))<z0.and.csine(1)>0.0_d) then
      x = z0
      y = yc(1)
    elseif(abs(xc(1))<z0.and.csine(2)>0.0_d) then
      x = xc(1)
      y = z0
    elseif(abs(yc(2))<z0.and.csine(3)>0.0_d) then
      x = -z0
      y = yc(2)
    elseif(abs(xc(2))<z0.and.csine(4)>0.0_d) then
      x = xc(2)
      y = -z0
    else
      !Erro. Sem colisão
      write(6,*) "Erro: Nenhuma colisao detectada"
      write(6,*) " Nenhuma colisao detectada "
      write(6,*) " Colisao = ", l
      write(6,*) "x = ",x, "  y = ",y,"   Theta = ", (180.0_d/pi)*theta
      stop
    Endif
    disco(1) = .false. ; disco(2) = .false.

  Endif

End Subroutine CollisionsCisal
