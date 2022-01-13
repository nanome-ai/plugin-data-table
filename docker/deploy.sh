#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

if [ -z "$(docker network ls -qf name=web-table-network)" ]; then
    echo "creating network"
    docker network create --driver bridge web-table-network
fi

existing=$(docker ps -aqf name=web-table)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name web-table \
--restart unless-stopped \
--network web-table-network \
--env no_proxy=web-table-server \
--env NO_PROXY=web-table-server \
-e ARGS="$*" \
web-table

docker run -d \
--name web-table-server \
--restart unless-stopped \
--network web-table-network \
-p 80:80 \
web-table-server
