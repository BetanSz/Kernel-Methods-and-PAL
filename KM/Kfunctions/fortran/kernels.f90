subroutine bernoulli(x, degre,s) 
  implicit none
  real ( kind = 8 ), intent(in)    ::  x
  integer ( kind = 4 ), intent(in) ::  degre
  real ( kind = 8 ) , intent(out) ::s
  if      (degre == 0) then 
    s = 1.
  else if (degre == 1) then
    s = x - 0.5
  else if (degre == 2) then
    s = x * x - x + 1./6.
  else if (degre == 3) then
    s = x * (x * (x - 1.5) + 0.5);
  else if (degre == 4) then
    s = x * (x * (x * (x - 2.0) + 1.0)) - 1./30.
  else if (degre == 5) then
    s = x * (x * (x * (x * (x - 2.5) + 5./3.)) - 1./6.)
  else if (degre == 6) then
    s = x * (x * (x * (x * (x * (x - 3.0) + 2.5)) - 0.5)) + 1./42.
  end if
  return 
end

subroutine kernel_spline(x, y, degre,z)
    real ( kind = 8 ), intent(in)    ::  x
    real ( kind = 8 ), intent(in)    ::  y
    integer ( kind = 4 ), intent(in) :: degre
    real ( kind = 8 ), intent(out)   :: z
    real ( kind = 8 ) u
    real ( kind = 8 ) v
    if (x<y) then
        u=x
        v=y
    else 
        u=y
        v=x
    end if

    if (degre == 0) then
        z = 1.
    else if (degre == 1) then
        z = 1. + u
    else if (degre == 2) then
        z = 1. + u * v + u * u * (v - u / 3. ) / 2.
    else 
        z = 0.
    end if
    return
end

subroutine kernel_bernoulli(x, y, degre,s)
  real ( kind = 8 ), intent(in)    ::  x
  integer ( kind = 4 ), intent(in) ::  degre
  real ( kind = 8 ), intent(in)   ::  y
  real ( kind = 8 ), intent(out)   ::   s
  real ( kind = 8 )   aux1
  real ( kind = 8 )   aux2
  integer ( kind = 4 ) i
    s=0.0
    aux1=1000
    aux2=0
    do i = 0, degre, 1
        call bernoulli(x,i,aux1)
        call bernoulli(y,i,aux2)
        s=s+( aux1*aux2 )/ (gamma(real(i+1,8))**2)
        !write(*,*) 'i,s,coso',i,s,gamma(real(i+1,8))**2
    end do
    call bernoulli(abs(x-y),2 * degre,aux1)
    !write(*,*) 'coso F',gamma(real(2*degre+1,8))
    s=s + ( ((-1.)**(degre+1.)) / gamma(real(2*degre+1,8)))* aux1
    return
end



subroutine kprod_bernoulli(X, Y, K, n ,p)
  real ( kind = 8 ), intent(in)    ::  X(n)
  real ( kind = 8 ), intent(in)    ::  Y(n)
  integer ( kind = 4 ), intent(in)    ::  K(n)
  integer ( kind = 4 ), intent(in) ::  n
  real ( kind = 8 ), intent(out)   ::  p
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  p=1.0
  do i = 1, n, 1
   call kernel_bernoulli(X(i),Y(i),K(i),aux1)
   p = p*aux1
   !write(*,*) 'i,X(i),Y(i),K(i),p',i,X(i),Y(i),K(i),p
  end do
  return
end

subroutine kprod_spline(X, Y, K, n ,p)
  real ( kind = 8 ), intent(in)    ::  X(n)
  real ( kind = 8 ), intent(in)    ::  Y(n)
  integer ( kind = 4 ), intent(in)    ::  K(n)
  integer ( kind = 4 ), intent(in) ::  n
  real ( kind = 8 ), intent(out)   ::  p
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  p=1.0
  do i = 1, n, 1
   call kernel_spline(X(i),Y(i),K(i),aux1)
   p = p*aux1
   !write(*,*) 'i,X(i),Y(i),K(i),p',i,X(i),Y(i),K(i),p
  end do
  return
end

subroutine predict(kernel_flag,nflag,XX,alpha,Y,K,n,nXX,ret_val)
  implicit none
  integer ( kind = 4 ), intent(in) ::  n
  integer ( kind = 4 ), intent(in) ::  nxx
  integer ( kind = 4 ), intent(in) ::  K(n)
  real ( kind = 8 ), intent(in) ::  Y(n)
  real ( kind = 8 ), intent(in) ::  alpha(nxx)
  real ( kind = 8 )   , intent(in) ::  XX(nXX,n)
  integer ( kind = 4 ), intent(in) ::  nflag
  character ( len = *), intent(in) ::  kernel_flag
  real ( kind = 8 ), intent(out)   ::  ret_val
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  ret_val=0.0
  if (kernel_flag=='KBN') then
        do i = 1, nXX, 1
            call kprod_bernoulli(XX(i,:),Y(:),K(:),n,aux1)
            ret_val=ret_val+alpha(i)*aux1
        end do
  end if
  if (kernel_flag=='KSP') then
        do i = 1, nXX, 1
            call kprod_spline(XX(i,:),Y(:),K(:),n,aux1)
            ret_val=ret_val+alpha(i)*aux1
        end do
  end if
  !K=K +regul*np.eye(len(K))
  return
end

subroutine predict2(Y,K,n,ret_val)
  implicit none
  integer ( kind = 4 ), intent(in) ::  n
!  integer ( kind = 4 ), intent(in) ::  nxx
!  real ( kind = 8 ), intent(in) ::  alpha(nxx)
  integer ( kind = 4 ), intent(in) ::  K(n)
  real ( kind = 8 ), intent(in) ::  Y(n)
  real ( kind = 8 ), intent(out)   ::  ret_val
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  ret_val=n
  !K=K +regul*np.eye(len(K))
  return
end

subroutine predict3(alpha,nXX,ret_val)
  implicit none
!  integer ( kind = 4 ), intent(in) ::  ny
  integer ( kind = 4 ), intent(in) ::  nxx
  real ( kind = 8 ), intent(in) ::  alpha(nxx)
!  integer ( kind = 4 ), intent(in) ::  K(n)
!  real ( kind = 8 ), intent(in) ::  Y(ny)
  real ( kind = 8 ), intent(out)   ::  ret_val
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  ret_val=nxx
  !K=K +regul*np.eye(len(K))
  write (*,*) 'nxx, funziona',nxx
  return
end

subroutine predict4(alpha,Y,K,n,nXX,ret_val)
  implicit none
  integer ( kind = 4 ), intent(in) ::  n
  integer ( kind = 4 ), intent(in) ::  nxx
  integer ( kind = 4 ), intent(in) ::  K(n)
  real ( kind = 8 ), intent(in) ::  Y(n)
  real ( kind = 8 ), intent(in) ::  alpha(nxx)
!  real ( kind = 8 )   , intent(in) ::  XX(nXX,n)
!  integer ( kind = 4 ), intent(in) ::  nflag
!  character ( len = *), intent(in) ::  kernel_flag
  real ( kind = 8 ), intent(out)   ::  ret_val
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  ret_val=0.0
  write (*,*) 'n,Y funziona',n,Y
  return
end

subroutine predict5(alpha,Y,K,n,nXX,ret_val)
  implicit none
  integer , intent(in) ::  n
  integer , intent(in) ::  nxx
  integer , intent(in) ::  K(n)
  real ( kind = 8 ), intent(in) ::  Y(n)
  real ( kind = 8 ), intent(in) ::  alpha(nxx)
!  real ( kind = 8 )   , intent(in) ::  XX(nXX,n)
!  integer ( kind = 4 ), intent(in) ::  nflag
!  character ( len = *), intent(in) ::  kernel_flag
  real ( kind = 8 ), intent(out)   ::  ret_val
  integer  i
  real ( kind = 8 )   aux1
  ret_val=0.0
  write (*,*) 'n,Y',n,Y
  return
end

subroutine vect_kprod(kernelF, XX, Y, K, nXX , nY, nK, P)
  implicit none
  real ( kind = 8 ), intent(in)    ::  XX(nXX,nK)
  real ( kind = 8 ), intent(in)    ::  Y(nY,nK)
  integer ( kind = 4 ), intent(in) ::  K(nK)
  integer ( kind = 4 ), intent(in) ::  nXX
  integer ( kind = 4 ), intent(in) ::  nY
  integer ( kind = 4 ), intent(in) ::  nK
  character ( len = *), intent(in) ::  kernelF
  real ( kind = 8 ), intent(out)   ::  P(nY)
  integer ( kind = 4 ) i
  real ( kind = 8 )   aux1
  aux1=0
  if (kernelF=='KBN') then
    do i = 1, nY, 1
       call kprod_bernoulli(XX(i,:),Y(i,:),K,nK,aux1)
       P(i) = aux1
    !write(*,*) 'F',i
    end do
  end if

  if (kernelF=='KSP') then
    do i = 1, nY, 1
       call kprod_spline(XX(i,:),Y(i,:),K,nK,aux1)
       P(i) = aux1
    !write(*,*) 'F',i
    end do
  end if
  return
end

subroutine inner_product(w, v, n, ret)
  implicit none
  real ( kind = 8 ), intent(in)    ::  w(n)
  real ( kind = 8 ), intent(in)    ::  v(n)
  integer ( kind = 4 ), intent(in) ::  n
  real ( kind = 8 ), intent(out)   ::  ret
  integer ( kind = 4 ) i
  ret=0
  do i = 1, n, 1
      ret= ret+w(i)*v(i)
    end do
  return
end

subroutine calcul_k(kernel_flag,nflag,XX,K,n,nXX,ret_k)
! Note that no profiting whatosever from the fact is symetric!
  implicit none
  integer ( kind = 4 ), intent(in) ::  n
  integer ( kind = 4 ), intent(in) ::  nxx
  integer ( kind = 4 ), intent(in) ::  K(n)
  real ( kind = 8 )   , intent(in) ::  XX(nXX,n)
  integer ( kind = 4 ), intent(in) ::  nflag
  character ( len = *), intent(in) ::  kernel_flag
  real ( kind = 8 ), intent(out)   ::  ret_k(nXX,nXX)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 )   aux1

  if (kernel_flag=='KBN') then
        do i = 1, nXX, 1
            do j = 1, nXX, 1
                call kprod_bernoulli(XX(i,:),XX(j,:),K,n,aux1)
                ret_k(i,j)=aux1
            end do
        end do
  end if
  if (kernel_flag=='KSP') then
        do i = 1, nXX, 1
            do j = 1, nXX, 1
                call kprod_spline(XX(i,:),XX(j,:),K,n,aux1)
                ret_k(i,j)=aux1
            end do
        end do
  end if
  return
end

subroutine calcul_k_generic(kernel_flag,nflag,K,XX,YY,n,nXX,nYY,ret_k)
! Note that no profiting whatosever from the fact is symetric!
  implicit none
  integer ( kind = 4 ), intent(in) ::  n
  integer ( kind = 4 ), intent(in) ::  nXX
  integer ( kind = 4 ), intent(in) ::  nYY
  integer ( kind = 4 ), intent(in) ::  K(n)
  real ( kind = 8 )   , intent(in) ::  XX(nXX,n)
  real ( kind = 8 )   , intent(in) ::  YY(nYY,n)
  integer ( kind = 4 ), intent(in) ::  nflag
  character ( len = *), intent(in) ::  kernel_flag
  real ( kind = 8 ), intent(out)   ::  ret_k(nXX,nYY)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 )   aux1
  if (kernel_flag=='KBN') then
        do i = 1, nXX, 1
            do j = 1, nYY, 1
                call kprod_bernoulli(XX(i,:),YY(j,:),K,n,aux1)
                ret_k(i,j)=aux1
            end do
        end do
  end if
  if (kernel_flag=='KSP') then
        do i = 1, nXX, 1
            do j = 1, nYY, 1
                call kprod_spline(XX(i,:),YY(j,:),K,n,aux1)
                ret_k(i,j)=aux1
            end do
        end do
  end if
  return
end


