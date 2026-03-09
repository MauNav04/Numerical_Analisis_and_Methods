function matrix = pentadiagonal(m, a, b, c, d, e)
  if m < 5
    error("pentadiagonal: m debe ser mayor o igual a 5")

  endif

  if length(a) != m
    error("pentadiagonal: el vector a debe ser de largo m")
  endif

  if length(b) != (m-1)
    error("pentadiagonal: el vector b debe ser de largo m-1")
  endif

  if length(c) != (m-1)
    error("pentadiagonal: el vector c debe ser de largo m-1")
  endif

  if length(d) != (m-2)
    error("pentadiagonal: el vector d debe ser de largo m-2")
  endif

  if length(e) != (m-2)
    error("pentadiagonal: el vector e debe ser de largo m-2")
  endif

  matrix = diag(a) ...
  + diag(b, -1) ...
  + diag(c, 1) ...
  + diag(d, -2) ...
  + diag(e, 2)

 endfunction
