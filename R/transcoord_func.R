transcoord_func = function (si, transfunc = "gaussian", l = 0.5, c = 0)
{
  si <- scale(si)
  if (transfunc == "gaussian") {
    out <- exp(-(si + c)^2/(2*l))
  }
  if (transfunc == "cosine") {
    out <- cos(pi*si/l + pi*c/2)
  }
  return(out)
}
