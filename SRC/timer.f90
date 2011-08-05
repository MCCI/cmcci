subroutine timer(tw,tu,ts)

  ! We use normal reals here, because that's the only kind the compiler seems
  ! to accept for this intrinsic.
  real, intent(out)  :: tw,tu,ts

  real  :: tarray(2),etime

  tw = etime(tarray)
  tu = tarray(1)
  ts = tarray(2)

  return
end subroutine timer
