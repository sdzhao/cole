#!/bin/bash
# how to use the script:
#
# to build image, run:
# ./scripts/cole.sh build
#
# to run rstudio in Docker, run:
# ./scripts/cole.sh r

ACTION=$1
IMAGE="cole/base"

if [[ -z $ACTION ]]; then
    echo "Usage: cole.sh [build/r]"
    exit 1
fi

if [[ $ACTION = 'build' ]];
then
    echo "Building Docker image: $IMAGE"
    docker build . -t $IMAGE -f Dockerfiles/Dockerfile
elif [[ $ACTION = 'r' ]];
then
    echo "Running R"
    docker run -v $(pwd)/local:/home --rm -p 8787:8787 -e PASSWORD=cole $IMAGE
else
    echo "Usage: cole.sh [build/r]"
fi
