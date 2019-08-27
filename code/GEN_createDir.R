
## ---- mkdir

mkdir <- function(file.name){
  if (!file.exists(file.name)){
    dir.create(file.name)
  } 
}




mkdir("code")
mkdir("figures")
mkdir("data")

rm(mkdir)





