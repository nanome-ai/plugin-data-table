#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

if [ -z "$(docker network ls -qf name=data-table-network)" ]; then
    echo "creating network"
    docker network create --driver bridge data-table-network
fi

existing=$(docker ps -aqf name=data-table)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name data-table \
--restart unless-stopped \
--network data-table-network \
--env no_proxy=data-table-server \
--env NO_PROXY=data-table-server \
-e ARGS="$*" \
data-table

docker run -d \
--name data-table-server \
--restart unless-stopped \
--network data-table-network \
-p 80:80 \
data-table-server
