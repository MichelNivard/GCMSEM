# R CMD check runs this top-level script in tests/.
source(file.path('r', 'test-core.R'), chdir=FALSE)
source(file.path('r', 'test-presentation.R'), chdir=FALSE)
source(file.path('r', 'test-transform.R'), chdir=FALSE)
source(file.path('r', 'test-estimation.R'), chdir=FALSE)
