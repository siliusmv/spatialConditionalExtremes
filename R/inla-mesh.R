#' @export
create_mesh = function(args, timeout = 5) {
  p = parallel:::mcfork()
  if (inherits(p, "masterProcess")) {
    mesh = do.call(INLA::inla.mesh.2d, args)
    mesh_file = tempfile()
    saveRDS(mesh, mesh_file)
    parallel:::mcexit(0L, mesh_file)
  }
  stdout = parallel:::readChildren(timeout = timeout)
  if (is.raw(stdout)) {
    child_pid = attr(stdout, "pid")
    mesh_file = unserialize(stdout)
    # for some reason, mcexit() does not always seem to close the process,
    # so we do it once more to be on the safe side
    tryCatch(parallel:::rmChild(child_pid), error = function(e) NULL)
    tryCatch(parallel:::rmChild(child_pid), error = function(e) NULL)
    mesh = readRDS(mesh_file)
    file.remove(mesh_file)
    return(mesh)
  }
  for (child in parallel:::children()) {
    tryCatch(parallel:::rmChild(child), error = function(e) NULL)
    tryCatch(parallel:::rmChild(child), error = function(e) NULL)
  }
  NULL
}
