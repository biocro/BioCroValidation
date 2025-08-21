# Helping function for overwriting elements of a base list with elements of a
# second list; here we assume that both lists have names, and that the names of
# new_list is a subset of the names of base_list
combine_lists <- function(base_list, new_list) {
    for (element_name in names(new_list)) {
        base_list[[element_name]] <- new_list[[element_name]]
    }
    base_list
}
