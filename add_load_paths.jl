
## Add paths of local modules and scripts to current load paths list

add_to_load_path!(path; load_path=LOAD_PATH) = path in load_path ? () : push!(load_path,path) 

add_to_load_path!(".\\src\\")

