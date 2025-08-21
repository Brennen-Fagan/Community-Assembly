# Need to use stackoverflow.com/a/47012149 to convert our arguments to a list.
callMeMaybe2 <- function(listofcharargs) {
  Charargs = unlist(listofcharargs)
  if(is.null(Charargs)) return(alist())
  eval(parse(
    text = paste0("alist(",
                  paste(parse(text = Charargs),
                        collapse = ","),")")
  ))
}
