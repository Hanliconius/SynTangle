# Convenience loader for a source checkout. Installed users: library(SynTangleR).
local({
    here <- dirname(normalizePath(sys.frame(1)$ofile))
    source(file.path(here,'..','..','r','SynTangleR','R','flip_assist.R'),local=parent.env(environment()))
})
