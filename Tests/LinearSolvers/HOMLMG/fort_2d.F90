
subroutine fort_set_coef (lo, hi, exact, elo, ehi, alpha, alo, ahi, beta, blo, bhi, &
     rhs, rlo, rhi, dx, prob_lo, prob_hi, a, b, sigma, w, bct) bind(c)
  use amrex_fort_module, only : amrex_real
  use iso_c_binding, only : c_char
  implicit none
  integer, dimension(3), intent(in) :: lo, hi, elo, ehi, alo, ahi, blo, bhi, rlo, rhi
  real(amrex_real), intent(in) :: dx(2), prob_lo(2), prob_hi(2), a, b, sigma, w
  character(kind=c_char), intent(in) :: bct
  real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2))
  real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2))
  real(amrex_real), intent(inout) :: beta (blo(1):bhi(1),blo(2):bhi(2))
  real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2))

  integer :: i,j
  double precision x, y, xc, yc
  double precision r, theta, dbdrfac
  double precision pi, fpi, tpi, fac

  pi = 4.d0 * atan(1.d0)
  tpi = 2.0d0 * pi
  fpi = 4.0d0 * pi
  fac = 12.d0 * pi**2

  xc = (prob_hi(1) + prob_lo(1))/2.d0
  yc = (prob_hi(2) + prob_lo(2))/2.d0

  theta = 0.5d0*log(3.d0) / (w + 1.d-50)
      
  do j = lo(2)-1, hi(2)+1
     y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
     do i = lo(1)-1, hi(1)+1
        x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
         
        r = sqrt((x-xc)**2 + (y-yc)**2)
             
        beta(i,j) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
     end do
   end do
  
   do j = lo(2), hi(2)
      y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
      do i = lo(1), hi(1)
         x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
         
         r = sqrt((x-xc)**2 + (y-yc)**2)
         
         dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
         dbdrfac = dbdrfac * b
         
         alpha(i,j) = 1.d0

         if (bct .eq. 'p' .or. bct .eq. 'n') then
            exact(i,j) = 1.d0 * cos(tpi*x) * cos(tpi*y)   &
                 &      + .25d0 * cos(fpi*x) * cos(fpi*y)

            rhs(i,j) = beta(i,j)*b*fac*(cos(tpi*x) * cos(tpi*y)    &
                 &                        + cos(fpi*x) * cos(fpi*y))   &
                 &   + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y)    &
                 &                     + pi*sin(fpi*x) * cos(fpi*y))   &
                 &            + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y)    &
                 &                     + pi*cos(fpi*x) * sin(fpi*y)))  &
                 &                   + a * (cos(tpi*x) * cos(tpi*y)    &
                 &               + 0.25d0 * cos(fpi*x) * cos(fpi*y) )
         else
            exact(i,j) = 1.d0 * sin(tpi*x) * sin(tpi*y)  &
                 &      + .25d0 * sin(fpi*x) * sin(fpi*y) 

            rhs(i,j) = beta(i,j)*b*fac*(sin(tpi*x) * sin(tpi*y)   &
                 &                        + sin(fpi*x) * sin(fpi*y))  &
                 &  + dbdrfac*((x-xc)*(-tpi*cos(tpi*x) * sin(tpi*y)   &
                 &                     - pi*cos(fpi*x) * sin(fpi*y))  &
                 &           + (y-yc)*(-tpi*sin(tpi*x) * cos(tpi*y)   &
                 &                     - pi*sin(fpi*x) * cos(fpi*y))) &
                 &                   + a * (sin(tpi*x) * sin(tpi*y)   &
                 &               + 0.25d0 * sin(fpi*x) * sin(fpi*y))
         end if
      end do
   end do

end subroutine fort_set_coef
