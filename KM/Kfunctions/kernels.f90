subroutine double_Bernoulli(x, degre,y)
  implicit none
  real ( kind = 8 ), intent(in)    ::  x
  integer ( kind = 4 ), intent(in) ::  degre
  real ( kind = 8 ), intent(out)   ::  y
  if      (degre == 0) then 
    y = 1.
  else if (degre == 1) then
    y = x - 0.5
  else if (degre == 2) then
    y = x * x - x + 1./6.
  else if (degre == 3) then
    y = x * (x * (x - 1.5) + 0.5);
  else if (degre == 4) then
    y = x * (x * (x * (x - 2.0) + 1.0)) - 1./30.
  else if (degre == 5) then
    y = x * (x * (x * (x * (x - 2.5) + 5./3.)) - 1./6.)
  else if (degre == 6) then
    y = x * (x * (x * (x * (x * (x - 3.0) + 2.5)) - 0.5)) + 1./42.
  end if
  return
end
