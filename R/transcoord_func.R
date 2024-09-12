transcoord_func = function (si, transfunc = "gaussian", l = 0.2, c = 0)
{
  si <- scale(si)
  if (transfunc == "gaussian") {
    out <- exp(-(si + c)^2/(2*l^2))
  }
  if (transfunc == "cosine") {
    out <- cos(pi*(si + l*c/2)/l)
  }
  return(out)
}
