subroutine CONDini(s1, s2, s3)
use Dados
implicit none
  Integer(kind=4), intent(inout)  :: s1, s2, s3
  Real(d)                         :: r4_random, aux
    Do
    x = (r4_random(s1, s2, s3) - 0.5_d)*ld
    y = (r4_random(s1, s2, s3) - 0.5_d)*ld
    aux = sqrt(x**2.0_d + y**2.0_d)
    if(aux > r) exit
    Enddo
theta = r4_Random(s1, s2, s3)*2.0_d*pi

return
End subroutine CONDini
