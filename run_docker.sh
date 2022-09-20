#!/bin/bash
## version
tag=0.2.9

# 		--entrypoint=/bin/bash \
docker run --rm -it \
		--entrypoint=/bin/bash \
        -v "$(pwd)"/input:/flywheel/v0/input \
        -v "$(pwd)"/output:/flywheel/v0/output \
        -v "$(pwd)"/work:/flywheel/v0/work \
        pereanez/mrtrix3preproc2sh:$tag
