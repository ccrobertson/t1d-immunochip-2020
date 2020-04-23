prettyColumn = function(x) {
  if (is.numeric(x)) {
    if (max(abs(x),na.rm=TRUE)>1) {
      x_out = format(x, digits=1)
    } else {
      x_out = format(x, digits=2)
    }
  } else {
    x_out = x
  }
  return(x_out)
}
prettyTable = function(df) {
  data.frame(lapply(df, prettyColumn))
}
