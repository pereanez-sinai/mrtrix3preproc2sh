#!/bin/bash
## version
tag=0.0.4

# 		--entrypoint=/bin/bash \
docker run --rm -it \
		--entrypoint=/bin/bash \
        -v "$(pwd)"/input:/flywheel/v0/input \
        -v "$(pwd)"/output:/flywheel/v0/output \
        pereanez/mrtrix3preproc2sh:$tag
