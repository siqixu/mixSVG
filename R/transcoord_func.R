transcoord_func = function (si, transfunc = "gaussian", q = 0.2, c = 0)
{
  si <- scale(si)
  l <- quantile(abs(si), probs = q)

  if (transfunc == "gaussian") {
    out <- exp(-(si + c)^2/(2 * l^2))
  }
  if (transfunc == "cosine") {
    out <- cos(2 * pi * si/l + pi*c)
  }
  return(out)
}
