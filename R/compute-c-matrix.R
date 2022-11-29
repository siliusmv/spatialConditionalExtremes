
#' Compute the matrix C that solves the equation
#' C H^-1 C^T = H^-1 J H^-1
#' @export
get_C = function(H, J) {
  if (length(H) == 1 && length(J) == 1) return(get_C_univ(H, J))
  tmp = svd(solve(H))
  M_1 = tmp$u %*% diag(sqrt(tmp$d)) %*% t(tmp$v)
  tmp = svd(solve(H, J) %*% solve(H))
  M_2 = tmp$u %*% diag(sqrt(tmp$d)) %*% t(tmp$v)
  C = t(solve(M_1, M_2))
  C
}

# Compute C when H and J are scalars
get_C_univ = function(H, J) {
  M_1 = sqrt(1 / H)
  M_2 = sqrt(J / H^2)
  C = M_2 / M_1
  C
}
